#include <linux/ratelimit.h>
#include <linux/spinlock.h>
#include <linux/sched.h>
#include <linux/printk.h>
#include <linux/random.h>
#include <linux/kernel.h>

#include "core_rl.h"
#include "sched.h"

spinlock_t mLock = __SPIN_LOCK_UNLOCKED(mLock);
long QVals[NR_STATE];
long NR_visit[NR_STATE];

int prev_state;

// dummy variable for lock acquire and release
unsigned long flags;

//Function prototype declaration
int calculate_Qstate(int cpu, int pid);
int calculate_state(void);
void update_weights(int hit);

int max_index(int states[]);
long calculate_reward(void);

void printfp(char *s, long a);
long mul(long x, long y);
long div(long x, long y);

int best_action(int pid, int next_state[]);
int least_action(int states[]);

long R;
# define precision 1000 // 0.001
# define prscale 100 // precision/10

// All in frac * precision; eg 0.2 --> 200
#define	alpha 100
#define	gamma 950
#define epsilon 100

int cpu_lock = -1;
int timer_tick = 0;

#define time_quanta 1

int select_task_rq_rl(struct task_struct *p, int cpu_given, int sd_flags, int wake_flags){

	int next_state[NR_CPU];
	int cpu_chosen;

	int is_acquire = 0;

	unsigned int coin;

	// Find a better way to do this
	spin_lock_irqsave(&mLock,flags);
	if(cpu_lock == -1){
		cpu_lock = -2;
		is_acquire = 1;
		update_weights(p->pid);
	}
	spin_unlock_irqrestore(&mLock,flags);
	
	cpu_chosen = best_action(p->pid,next_state);
	
	if(is_acquire == 1){

		// With prob of epsilon, make random action
		coin = get_random_int() % precision;
		if(coin <= epsilon){
			cpu_chosen =  least_action(next_state);
		}

		// Store prevState [ To be used in next weight update ] 
		prev_state = next_state[cpu_chosen];

		NR_visit[prev_state]++;				
		cpu_lock = cpu_chosen;

	}
		
	return cpu_chosen;
}

int best_action(int pid, int next_state[]){
	
	int cpu;

	for(cpu=0;cpu<NR_CPU;cpu++){
		// Get the current state f(s,a)		
		next_state[cpu] = calculate_Qstate(cpu,pid);		
	}

	// Max over all Q(s,a)
	return max_index(next_state);
}

void get_nr_process(long proc[]){
	int cpu;
	struct rq *rq;
	
	for(cpu=0; cpu<NR_CPU; cpu++){
		rq = cpu_rq(cpu);
		proc[cpu] = rq->nr_running;
	}
}

int get_proc_variance(long proc[]){
	int cpu,min,max,num;
		
	min = proc[0];
	max = min;
	
	for(cpu=1; cpu<NR_CPU; cpu++){
		num = proc[cpu];
		if(num > max) max = num;
		if(num < min) min = num;	
	}
	
	num = (max - min);

	return (num >= 4);
}

int calculate_Qstate(int cpu, int pid){
	long proc[NR_CPU];
	int is_hit, is_ghost;
	
	is_hit = (pid % NR_CPU == cpu);

	get_nr_process(proc);
	
	// Going to add a new process to this
	proc[cpu]++;

	is_ghost = get_proc_variance(proc);

	proc[cpu]--;
	
	return 2*is_hit + is_ghost;
}

long calculate_reward(){
	int cpu,min,max,num,R;
	long is_hit;

	struct rq *rq;

	rq = cpu_rq(0);
	min = rq->nr_running;
	max = min;

	for(cpu=1;cpu<NR_CPU;cpu++){
		rq = cpu_rq(cpu);
		num = rq->nr_running;
		if(num > max) max = num;
		if(num < min) min = num;
	}
	
	/*
		CPU hit reward for previous action 
		1 -> if hit correctly 
		0 -> wrong hit
	*/
	is_hit = (prev_state == 2) || (prev_state == 3);
	R = is_hit * precision;
	
	/*
		A ghost reward of -5 if cpus highly imbalanced
	*/

	num = max - min;
	if(num > 3){
		R = -5*precision;
	}

	return R;
}

int least_action(int states[]){
	int index = 0;
	long min = NR_visit[states[0]];
	long visitN;	
	int cpu;

	for(cpu=1;cpu<NR_CPU;cpu++){
		visitN = NR_visit[states[cpu]];
		if(visitN < min){
			min = visitN;
			index = cpu;
		}
	}

	return index;

}

int max_index(int states[]){
	int index = 0;
	long maxVal = QVals[states[0]];
	long Qv;	
	int cpu;


	for(cpu=1;cpu<NR_CPU;cpu++){
		Qv = QVals[states[cpu]]; 
		if(Qv > maxVal){
			maxVal = Qv;
			index = cpu;
		}
	}

	return index;
}

// Weight update functions

void update_weights(int hit){
	long delta, Vval;
	int action;
	int next_state[NR_CPU];
	int state;

	R = calculate_reward();	

	action = best_action(hit, next_state);
	Vval = QVals[next_state[action]];
	
	delta = (R + mul(gamma, Vval)) - QVals[prev_state];

	QVals[prev_state] += mul(alpha, delta);
}

void scheduler_tick_rl(int cpu){
	if(cpu != cpu_lock)
		return;

	timer_tick++;

	if(timer_tick == time_quanta){
		R = calculate_reward();
		cpu_lock = -1;
		timer_tick = 0;
	}
}

// Fixed point simulation functions

long mul(long x, long y){
	long res;
	res = (x * y)/ precision;
	return res;
}

long div(long x, long y){
	long res;
	res = (x * precision)/y;
	return res;
}

// Warning - change format specifier with precision

void printfp(char *s, long a){
	long dec = a/1000;
	long fra = a % 1000;
	
	if(a < 0)
		trace_printk("%s -%ld.%03ld\n",s,-1*dec,-1*fra);
	else
		trace_printk("%s %ld.%03ld \n",s,dec,fra);
}

