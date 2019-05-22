
struct param_rl{
	long bias;
	long nr_running;
};

#define NR_CPU 4
#define NR_FEAT 2

int select_task_rq_rl(struct task_struct *p, int cpu, int sd_flags, int wake_flags);
