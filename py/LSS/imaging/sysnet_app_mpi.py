#!/usr/bin/env python
""" 
    A high-level code for running the SYSNet software
"""
import sysnet

if __name__ == '__main__':
    config = sysnet.parse_cmd_arguments('sysnet_config.yaml')
    pipeline = sysnet.SYSNetMPI(config)
    pipeline.run()
