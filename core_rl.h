#define NR_CPU 4
#define NR_STATE 4

int select_task_rq_rl(struct task_struct *p, int cpu, int sd_flags, int wake_flags);
void scheduler_tick_rl(int cpu);
