#!/usr/bin/env python
""" 
    A high-level code for running the SYSNet software
"""
import sysnet

if __name__ == '__main__':
    config = sysnet.parse_cmd_arguments('sysnet_config.yaml')
    pipeline = sysnet.SYSNetSnapshotMPI(config)
    # T_0 sets the number frequency of snapshots (T_0=50 saves a snapshot every 50 epochs)
    pipeline.config.scheduler_kwargs.update(T_0=50, T_mult=1) 
    pipeline.run()
