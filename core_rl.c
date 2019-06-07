#include <linux/ratelimit.h>
#include <linux/spinlock.h>
#include <linux/sched.h>
#include <linux/printk.h>

#include "core_rl.h"
#include "sched.h"

spinlock_t mLock = __SPIN_LOCK_UNLOCKED(mLock);
int weights[NR_FEAT];

struct param_rl prev_param;
long prev_Vval;

// dummy variable for lock acquire and release
unsigned long flags;

//Function prototype declaration
void calculate_Qparam(struct param_rl *par, int cpu, struct task_struct *p);
void calculate_param(struct param_rl *par);
void update_weights(void);
long calculate_Qval(struct param_rl p);

int max_index(long Qval[]);
void store_param(struct param_rl p);
long calculate_reward(void);

void printfp(char *s, long a);
long mul(long x, long y);
long div(long x, long y);

# define precision 1000 // 0.001
# define prscale 100 // precision/10

// All in frac * precision; eg 0.02 --> 20
#define	alpha 200
#define	gamma 980
#define epsilon 10

int cpu_lock = -1;
int timer_tick = 0;

#define time_quanta 2

int select_task_rq_rl(struct task_struct *p, int cpu_given, int sd_flags, int wake_flags){
	struct param_rl parameters[NR_CPU];
	long Qval[NR_CPU];
	int cpu, cpu_chosen;

	int is_acquire = 0;

	// Find a better way to do this
	spin_lock_irqsave(&mLock,flags);
	if(cpu_lock == -1){
		trace_printk("Lock Acquired\n");
		cpu_lock = -2;
		is_acquire = 1;
	}
	spin_unlock_irqrestore(&mLock,flags);
	
	for(cpu=0;cpu<NR_CPU;cpu++){
		// Get the current state f(s,a)		
		calculate_Qparam(&parameters[cpu],cpu,p);
		// Get Q(s,a) 
		Qval[cpu] = calculate_Qval(parameters[cpu]);		
	}

	// Max over all Q(s,a)
	cpu_chosen = max_index(Qval);
	
	if(is_acquire == 1){
		// Store max{a : Q(s,a)}. [ To be used in next weight update ] 
		prev_Vval = Qval[cpu_chosen];
		
		// Store f(s,a)
		store_param(parameters[cpu_chosen]);
		
		cpu_lock = cpu_chosen;

		trace_printk("CPU chosen: %d\n",cpu_chosen);	
	}
		
	return cpu_chosen;
}

long calculate_Qval(struct param_rl p){
	
	long res1,res2,res3,qval;

	res1 = mul((p.bias), weights[0]);
	res2 = mul((p.nr_running), weights[1]);
	res3 = mul((p.is_hit), weights[2]);

	qval = res1+res2+res3;
	return qval;
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
	int cpu,min,max,num,var;
		
	min = proc[0];
	max = min;
	
	for(cpu=1; cpu<NR_CPU; cpu++){
		num = proc[cpu];
		if(num > max) max = num;
		if(num < min) min = num;	
	}
	
	num = (max - min);

	if(num != 0)
		var = ((num - 1) * (num - 1)) * 10;
	else
		var = -10;

	return var;
}

void calculate_Qparam(struct param_rl *par, int cpu, struct task_struct *p){
		
	long proc[NR_CPU];
	
	par->bias = prscale;
	par->is_hit = ((p->pid) % NR_CPU == cpu) * prscale;

	get_nr_process(proc);
	
	// Going to add a new process to this
	proc[cpu]++;

	par->nr_running = get_proc_variance(proc);

	proc[cpu]--;

	if(cpu_lock == -2){
		trace_printk("Nr_running %ld\n", proc[cpu]);		
		printfp("is_hit", par->is_hit);
		printfp("is_ghost", par->nr_running);				
		trace_printk("\n");		
	}	
}

void calculate_param(struct param_rl *par){
	
	long proc[NR_CPU];
	int cpu;

	par->bias = prscale;
	par->is_hit = prev_param.is_hit;

	get_nr_process(proc);

	for(cpu=0; cpu<NR_CPU; cpu++){
		printfp("Nr_pr",proc[cpu]*precision);
	}

	par->nr_running = get_proc_variance(proc);

	printfp("is_hit", par->is_hit);
	printfp("is_ghost", par->nr_running);
}

long calculate_reward(){
	int cpu,min,max,num,R;

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
		0.1 -> if hit correctly 
		0 -> wrong hit
	*/
	R = prev_param.is_hit;

	/*
		A ghost reward of -1 if cpus highly imbalanced
	*/

	num = max - min;
	if(num > 3){
		R = -1*precision;
	}

	return R;
}


void store_param(struct param_rl p){
	prev_param.bias = p.bias;
	prev_param.nr_running = p.nr_running;
	prev_param.is_hit = p.is_hit;
}

int max_index(long Qval[]){
	int index = 0;
	long maxVal = Qval[0];
	int cpu;

	for(cpu=1;cpu<NR_CPU;cpu++){
		if(Qval[cpu] > maxVal){
			maxVal = Qval[cpu];
			index = cpu;
		}
	}

	return index;
}

// Weight update functions

void update_weights(){
	long R, delta, alphaDelta, Vval;
	struct param_rl parameter;

	R = calculate_reward();	
	
	calculate_param(&parameter);

	Vval = calculate_Qval(parameter);
	
	delta = (R + mul(gamma, Vval)) - prev_Vval;

	printfp("Reward", R);
	printfp("Vval", Vval);
	printfp("prev_Vval", prev_Vval);	
	printfp("Delta", delta);
	
	alphaDelta = mul(alpha, delta);

	printfp("weights[0]", weights[0]);
	printfp("weights[1]", weights[1]);
	printfp("weights[2]", weights[2]);

	weights[0] += mul(alphaDelta, prev_param.bias);	
	weights[1] += mul(alphaDelta, prev_param.nr_running);
	weights[2] += mul(alphaDelta, prev_param.is_hit);

	printfp("weights[0]", weights[0]);
	printfp("weights[1]", weights[1]);
	printfp("weights[2]", weights[2]);

}

void scheduler_tick_rl(int cpu){
	if(cpu != cpu_lock)
		return;

	timer_tick++;

	if(timer_tick == time_quanta){
		update_weights();
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


