import numpy as np
import pandas as pd
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys


# def compute_flops(M, N, T, iterations):
#     M_sq = pow(M,2)
#     #flops_per_it = 13*M_sq*(T-1) + M_sq*(T+1) + 10*M*T + T - M + M*N*(T+2)
#     flops_per_it = M*T*(M + N)
#     flops = flops_per_it*iterations
#     return flops

def plot(input_file):
    times = {}
    flops = {}
    with open(input_file) as fp:
        line = fp.readline()
        tokens = line.split(",")
        if 'var_M' in input_file:
            xlabel = 'M'
            key_idx = 2
            title = 'N = ' + tokens[3] + ', T = ' + tokens[4]
        elif 'var_N' in input_file:
            xlabel = 'N'
            key_idx = 3
            title = 'M = ' + tokens[2] + ', T = ' + tokens[4]
        elif 'var_T' in input_file:
            xlabel = 'T'
            key_idx = 4
            title = 'M = ' + tokens[2] + ', N = ' + tokens[3]
        else:
            print ('Error: Expected file whose name contains \'var_M\' or \'var_N\' or \'var_T\'')
            exit()
        while line:
            tokens = line.split(",")
            if (not(tokens[0] in times)):
                times[tokens[0]] = {}
            if not(tokens[0] in flops):
                flops[tokens[0]] = {}
            key = int(tokens[key_idx])
            times[tokens[0]][key] = float(tokens[5])
            flops[tokens[0]][key] = float(tokens[9])
            # T = int(tokens[4])
            iterations = int(tokens[6])
            line = fp.readline()

    for k in times:
        d = times[k]
        od = collections.OrderedDict(sorted(d.items()))
        times[k] = od
        d = flops[k]
        od = collections.OrderedDict(sorted(d.items()))
        flops[k] = od

    num_points = len(od)
    x=range(num_points)
    x_l = [str(i) for i in od]
    fig, ax = plt.subplots(figsize=(11,9))

    labs=[]
    lns=[]
    for i in times:
        print(i)
        d = times[i]
        d_flops = flops[i]
        y = []
        for j in d:
            y.append(d_flops[j]/d[j])
        print(y)
        line=plt.plot(x, y, marker='o', linestyle='-', label=i)
        lns += line

    ax.set_xticks(x)
    ax.set_xticklabels(x_l)

    for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(11)

    for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(10)

    processor = 'Intel(R) Core(TM) i7-8750H CPU @ 2.20GHz'
    compiler = 'icc 19.1'
    flags = '-mavx2 -mfma'
    subtitle = processor + '\n' + '\n' + compiler + ' (' + flags + ')'

    ax.set_title('Baum-Welch Performance (' + title + ')' + '\n' + subtitle, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=14)
    plt.ylabel('Performance [flop/cycle]', rotation=0, fontsize=14)
    ax.yaxis.set_label_coords(0.06,1.02)
    ax.set_ylim(0, 6)

    labs = [l.get_label() for l in lns]
    lgd = ax.legend(lns, labs, ncol=3, loc=2, fontsize=12)
    fig.savefig("perf_" + xlabel + ".png")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print ('Usage: python plot_performance.py <cycles_file>')
        exit()
    plot(sys.argv[1])
