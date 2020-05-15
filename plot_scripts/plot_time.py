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
            key = int(tokens[key_idx])
            times[tokens[0]][key] = float(tokens[5])
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
        ax.set_title('Time benchmark (' + title + ')', fontsize=14)
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel('Time (msec)', fontsize=14)

    else:
        ax.set_title('Cycles benchmark (' + title + ')', fontsize=14)
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel('Number of Cycles', fontsize=14)

    labs = [l.get_label() for l in lns]
    lgd = ax.legend(lns, labs, ncol=1, loc=2, fontsize=12)
    fig.savefig(save + '_' + xlabel + '.png')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print ('Usage: python plot_time.py <file> <time/cycles>')
        exit()
    if (sys.argv[2] == "time"):
        plot(1, sys.argv[1], "times")
    elif (sys.argv[2] == "cycles"):
        plot(0, sys.argv[1], "cycles")
    else:
        print ('Usage: python preprocess.py <file> <time/cycles>')
        exit()
