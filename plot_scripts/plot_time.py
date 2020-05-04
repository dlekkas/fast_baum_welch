import numpy as np
import pandas as pd
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

def plot(plot_time, input_file, save):
    times = {}
    with open(input_file) as fp:
        line = fp.readline()
        while line:
            tokens = line.split(",")
            if (not(tokens[0] in times)):
                times[tokens[0]] = {}
            key_pair = (int(tokens[2]), int(tokens[3]))
            times[tokens[0]][key_pair] = float(tokens[-1])
            line = fp.readline()

    for k in times:
        d = times[k]
        od = collections.OrderedDict(sorted(d.items()))
        times[k] = od

    num_points = len(od)
    x=range(num_points)
    x_l = [str(i) for i in od]
    fig, ax = plt.subplots(figsize=(11,9))

    labs=[]
    lns=[]
    for i in times:
        print(i)
        d = times[i]
        y = [d[j] for j in d]
        print(y)
        line=plt.plot(x, y, marker='o', linestyle='-', label=i)
        lns += line

    ax.set_xticks(x)
    ax.set_xticklabels(x_l)

    for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(11)

    for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(10)

    if (plot_time):
        ax.set_title('Time benchmark', fontsize=14)
        ax.set_xlabel('HMM size(M, N)', fontsize=14)
        ax.set_ylabel('Time(msec)', fontsize=14)

    else:
        ax.set_title('Cycles benchmark', fontsize=14)
        ax.set_xlabel('HMM size(M, N)', fontsize=14)
        ax.set_ylabel('Number of Cycles', fontsize=14)

    labs = [l.get_label() for l in lns]
    lgd = ax.legend(lns, labs, ncol=1, loc=2, fontsize=12)
    fig.savefig(save)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print ('Usage: python preprocess.py <time/cycles>')
        exit()
    if (sys.argv[1] == "time"):
        plot(1,"../build/results_time.txt", "times.png")
    elif (sys.argv[1] == "cycles"):
        plot(0,"../build/results_cycles.txt", "cycles.png")
    else:
        print ('Usage: python preprocess.py <time/cycles>')
        exit()

