#include <linux/ratelimit.h>
#include <linux/spinlock.h>
#include <linux/sched.h>
#include <linux/printk.h>

#include "core_rl.h"
#include "sched.h"

spinlock_t mLock = __SPIN_LOCK_UNLOCKED(mLock);
int weights[NR_FEAT];

struct param_rl prev_param;
long prev_Qval;

// dummy variable for lock acquire and release
unsigned long flags;

//Function prototype declaration
void calculate_param(struct param_rl *p, int cpu);
void update_weights(long curr_Qval);
long calculate_Qval(struct param_rl p);

int max_index(long Qval[]);
void store_param(struct param_rl p);
long calculate_reward(void);

void printfp(char *s, long a);
long mul(long x, long y);
long div(long x, long y);


# define precision 1000 // 0.001

// All in frac * precision; eg 0.02 --> 20
#define	alpha 200
#define	gamma 980
#define epsilon 10

int select_task_rq_rl(struct task_struct *p, int cpu_given, int sd_flags, int wake_flags){
	struct param_rl parameters[NR_CPU];
	long Qval[NR_CPU]; long maxQval;
	int cpu, cpu_chosen;

	spin_lock_irqsave(&mLock,flags);
	
	trace_printk("\n");
	
	for(cpu=0;cpu<NR_CPU;cpu++){
		// Get the current state f(s',a')		
		calculate_param(&parameters[cpu],cpu);
		// Get Q(s',a') with current weights
		Qval[cpu] = calculate_Qval(parameters[cpu]);	
	}

	// Finding max Q(s',a')	
	cpu = max_index(Qval);
	maxQval = Qval[cpu];

	// Update weights
	update_weights(maxQval);

	// Current state is s from now onwards

	// Calculate Q(s,a) [ with updated weights ]
	for(cpu=0;cpu<NR_CPU;cpu++){
		Qval[cpu] = calculate_Qval(parameters[cpu]);	
	}

	// Max over all Q(s,a)
	cpu_chosen = max_index(Qval);
	
	// Store Q(s,a) to be used in next weight update 
	prev_Qval = Qval[cpu_chosen];
	
	// Store f(s,a)
	store_param(parameters[cpu_chosen]);
 

	spin_unlock_irqrestore(&mLock,flags);
	
	trace_printk("CPU chosen: %d\n",cpu_chosen);	
	return cpu_chosen;
}
 

void update_weights(long curr_Qval){
	long R, delta;

	R = calculate_reward();	
	
	delta = (R + mul(gamma, curr_Qval)) - prev_Qval;

	printfp("Reward", R);	
	printfp("Delta", delta);
	
	weights[0] = weights[0] + mul(mul(alpha, delta), prev_param.bias);	
	weights[1] = weights[1] + mul(mul(alpha, delta), prev_param.nr_running);

	printfp("weights[0]", weights[0]);
	printfp("weights[1]", weights[1]);
}

long calculate_Qval(struct param_rl p){
	
	long res1,res2,qval;

	res1 = mul((p.bias), weights[0]);
	res2 = mul((p.nr_running), weights[1]);

	qval = res1+res2;
	return qval;
}

void calculate_param(struct param_rl *par, int cpu){
		
	struct rq *rq;
	long proc[NR_CPU];
	int cpuN,mean,total,var;

	total = 0;

	for(cpuN=0; cpuN<NR_CPU; cpuN++){
		rq = cpu_rq(cpuN);
		proc[cpuN] = rq->nr_running;
		total += proc[cpuN];	
	}

	// Going to add a new process to this
	proc[cpu]++;
	total += 1;

	// bias is 0.1
	par->bias = precision / 10;
	
	mean = total/NR_CPU;

	var = 0;
	
	for(cpuN=0;cpuN<NR_CPU;cpuN++){
		rq = cpu_rq(cpuN);
		var += (proc[cpuN] - mean) * (proc[cpuN] - mean);	
	}

	// var/10
	par->nr_running = var*10;
	
	printfp("variance", par->nr_running);
}

long calculate_reward(){
	int cpu,min,max,num;
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

	num = min - max;

	return 10 * ( 3 - (num*num) );
}


void store_param(struct param_rl p){
	prev_param.bias = p.bias;
	prev_param.nr_running = p.nr_running;
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


