#include <linux/ratelimit.h>
#include <linux/spinlock.h>
#include <linux/sched.h>
#include <linux/printk.h>
#include <linux/random.h>

#include "core_rl.h"
#include "sched.h"

spinlock_t mLock = __SPIN_LOCK_UNLOCKED(mLock);
long QVals[NR_STATE];

int prev_state;

// dummy variable for lock acquire and release
unsigned long flags;

//Function prototype declaration
int calculate_Qstate(int cpu, struct task_struct *p);
void update_weights(struct task_struct *p);
int best_action(struct task_struct *p, int next_state[]);

int max_index(int Qval[]);
long calculate_reward(void);

void printfp(char *s, long a);
long mul(long x, long y);
long div(long x, long y);

# define precision 1000 // 0.001
# define prscale 100 // precision/10

// All in frac * precision; eg 0.02 --> 20
#define	alpha 200
#define	gamma 950
#define epsilon 10

int cpu_lock = -1;
int timer_tick = 0;

#define time_quanta 1

int select_task_rq_rl(struct task_struct *p, int cpu_given, int sd_flags, int wake_flags){
	
	int cpu_chosen;
	int next_state[NR_CPU];
	int is_acquire = 0;

	unsigned int coin;
	
	// Find a better way to do this
	spin_lock_irqsave(&mLock,flags);
	if(cpu_lock == -1){
		trace_printk("Lock Acquired\n");
		cpu_lock = -2;
		is_acquire = 1;
		update_weights(p);
	}
	spin_unlock_irqrestore(&mLock,flags);

	cpu_chosen = best_action(p,next_state);
	
	if(is_acquire == 1){
		// With prob of epsilon, make random action
		coin = get_random_int() % precision;
		if(coin <= epsilon){
			trace_printk("Random Action\n");
			cpu_chosen =  get_random_int() % NR_CPU;
		}else{
			trace_printk("Deterministic Action\n");
		}

		// Store prevState [ To be used in next weight update ] 
		prev_state = next_state[cpu_chosen];

		cpu_lock = cpu_chosen;
				
		trace_printk("CPU chosen: %d\n",cpu_chosen);			
	}

	return cpu_chosen;
}

// Weight update functions

void update_weights(struct task_struct *p){
	long R, delta, Vval;
	int action, next_state[NR_CPU], state;

	R = calculate_reward();	

	// Max Q(s',a')
	action = best_action(p, next_state);
	Vval = QVals[next_state[action]];
	
	delta = (R + mul(gamma, Vval)) - QVals[prev_state];

	printfp("Reward", R);
	printfp("Vval", Vval);
	printfp("prev_Vval", QVals[prev_state]);	
	printfp("Delta", delta);

	QVals[prev_state] += mul(alpha, delta);

	for(state = 0; state < NR_STATE;state++){
		printfp("Qval", QVals[state]);
	}	
}

int best_action(struct task_struct *p, int next_state[]){
	
	int cpu;

	for(cpu=0;cpu<NR_CPU;cpu++){
		// Get the current state f(s,a)		
		next_state[cpu] = calculate_Qstate(cpu,p);		
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

int calculate_Qstate(int cpu, struct task_struct *p){
		
	long proc[NR_CPU];
	int is_hit, is_ghost;
	
	is_hit = ((p->pid) % NR_CPU == cpu);

	get_nr_process(proc);
	
	// Going to add a new process to this
	proc[cpu]++;

	is_ghost = get_proc_variance(proc);

	proc[cpu]--;

	if(cpu_lock == -2){
		trace_printk("Nr_running %ld\n", proc[cpu]);		
		printfp("is_hit", is_hit);
		printfp("is_ghost", is_ghost);				
		trace_printk("\n");		
	}
	
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
		A ghost reward of -1 if cpus highly imbalanced
	*/

	num = max - min;
	if(num > 3){
		R = -5*precision;
	}

	return R;
}

int max_index(int state[]){
	int index = 0;
	long maxVal = QVals[state[0]];
	long Qv;	
	int cpu;

	for(cpu=1;cpu<NR_CPU;cpu++){
		Qv = QVals[state[cpu]]; 
		if(Qv > maxVal){
			maxVal = Qv;
			index = cpu;
		}
	}

	return index;
}


void scheduler_tick_rl(int cpu){
	if(cpu != cpu_lock)
		return;

	timer_tick++;

	if(timer_tick == time_quanta){
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


