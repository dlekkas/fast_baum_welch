import numpy as np
import pandas as pd
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys


def compute_flops(M, N, T, iterations):
    M_sq = pow(M,2)
    #flops_per_it = 13*M_sq*(T-1) + M_sq*(T+1) + 10*M*T + T - M + M*N*(T+2)
    flops_per_it = M*T*(M + N)
    flops = flops_per_it*iterations
    return flops

def plot(input_file):
    times = {}
    with open(input_file) as fp:
        line = fp.readline()
        while line:
            tokens = line.split(",")
            if (not(tokens[0] in times)):
                times[tokens[0]] = {}
            key_pair = (int(tokens[2]), int(tokens[3]))
            times[tokens[0]][key_pair] = float(tokens[5])
            T = int(tokens[4])
            iterations = int(tokens[6])
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
        y = []
        for j in d:
            m = j[0]
            n = j[1]
            flops = compute_flops(m, n, T, iterations)
            y.append(flops/d[j])
        print(y)
        line=plt.plot(x, y, marker='o', linestyle='-', label=i)
        lns += line

    ax.set_xticks(x)
    ax.set_xticklabels(x_l)

    for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(11)

    for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(10)

    ax.set_title('Baum-Welch Performance', fontsize=14)
    ax.set_xlabel('HMM size(M, N)', fontsize=14)
    ax.set_ylabel('Performance(flops/cycle)', fontsize=14)

    labs = [l.get_label() for l in lns]
    lgd = ax.legend(lns, labs, ncol=1, loc=2, fontsize=12)
    fig.savefig("perf.png")

if __name__ == '__main__':
    if len(sys.argv) != 1:
        print ('Usage: python preprocess.py')
        exit()
    plot("../build/results_cycles.txt")
