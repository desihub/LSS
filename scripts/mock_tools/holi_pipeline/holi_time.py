#!/usr/bin/env python3

import subprocess
import re
import shlex
import numpy as np
import os

import matplotlib.pyplot as plt


def get_time_real():
    # cmd = "grep real *.log| awk '{print $2}'"
    cmd = "grep real *.log"
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        print("Error running command:", res.stderr)
        return None
    # convert the output to a list of floats
    l_val = []
    for line in res.stdout.split("\n"):
        # print(line)
        if line == "":
            continue
        str_t = line.split("\t")[1]
        # print(str_t)
        lst = re.findall(r"(\d*\.*\d+)", str_t)
        # print(lst)
        lf = [float(t) for t in lst]
        vf = lf[0] + lf[1] / 60
        l_val.append(vf)
        # print(lst, vf)
    return l_val


def process_logs_directory():
    os.chdir("logs")
    t_logs = get_time_real()
    t_tasks = np.array(t_logs).reshape(-1, 10)
    # print(np.mean(t_tasks, axis=0))
    ts1 = np.sum(t_tasks[:, :3], axis=1)
    # print(ts1)
    ts4 = np.sum(t_tasks[:, 3:6], axis=1)
    # print("4",ts4)
    ts5 = np.sum(t_tasks[:, 6:8], axis=1)
    # print("5:",ts5)
    ts6 = t_tasks[:, 8].copy()
    # print("6:",ts6)
    ts7 = t_tasks[:, 9].copy()
    # print("7:",ts7)
    t_tasks = np.stack((ts1, ts4, ts5, ts6, ts7)).T
    # print("t_tasks shape:", t_tasks.shape)
    # print(np.mean(t_tasks, axis=0))
    return t_tasks


def process_logs_roots():
    t_logs = get_time_real()
    t_tasks = np.array(t_logs).reshape(-1, 3)
    # print(np.mean(t_tasks, axis=0))
    return t_tasks


def print_stat(a_tasks, a_stage, full=True):
    print("=====================================")
    print(np.mean(a_tasks, axis=0))
    print(np.max(a_tasks, axis=0))
    print(np.mean(a_stage, axis=0))
    print(np.max(a_stage, axis=0))
    if full:
        print("=====================================")
        np.set_printoptions(precision=1, suppress=True)
        print(a_tasks)
        print(a_stage)


def create_plot(a_tasks, a_stage, nb_cpu):
    fig, ax = plt.subplots(layout="constrained")
    print(a_tasks.shape)
    # plt.boxplot(a_tasks, tick_labels=["1", "4", "5", "6", "7"])
    ax.set_title("Time for step 1-7 of Holi pipeline\nFor 18 seeds on /cfs file system with 1 CPU per step")
    vs = a_tasks.T
    plt.boxplot(
        (vs[0], a_stage[:, 1]/nb_cpu, vs[1], vs[2], vs[3], vs[4]),
        #tick_labels=["1: simu cat", "3: BRICKMASK", "4: apply mask", "5: contaminant", "6: join cat", "7: init AltMTL"],
        tick_labels=["1", "3", "4", "5", "6", "7"],
        
    )
    # for label in ax.get_xticklabels():
    #     label.set_rotation(45)
    #     label.set_rotation_mode("anchor")
    ax.set_xlabel("Step number")
    ax.set_ylabel("Time (minutes)")
    ax.grid()


a_stage, a_tasks = process_logs_roots(), process_logs_directory()
print_stat(a_tasks, a_stage)
create_plot(a_tasks, a_stage, nb_cpu=10)
plt.show()

# print(process_logs_roots())
# print("=====================================")
# print(process_logs_directory())
