from desipipe import Queue, Environment, TaskManager, spawn, setup_logging

setup_logging()

queue = Queue('LRGsysnet')
queue.clear()

environ = Environment('nersc-cosmodesi')
tm = TaskManager(queue=queue, environ=environ)
# 24 tasks per node, 72 workers total = 3 nodes, for 3 hours each
tm_sysnet = tm.clone(scheduler=dict(max_workers=75), provider=dict(provider='nersc', time='03:00:00', mpiprocs_per_worker=1, nodes_per_worker=1/25, output='_sbatch/slurm-%j.out', error='_sbatch/slurm-%j.err'))


@tm_sysnet.bash_app
def sysnet(i, zrange):
    return ['./scripts/Y1ab2ndgen_sysnet_tracer_zbin.sh', i, 'LRG_ffa', f'{zrange[0]:.1f}_{zrange[1]:.1f}']


if __name__ == '__main__':

    setup_logging()
    for zrange in [(0.4, 0.6), (0.6, 0.8), (0.8, 1.1)]:
        for i in range(0, 25):
            sysnet(i, zrange)
    spawn(queue, spawn=True)  # that'll launch all jobs