import numpy as np
import pandas as pd
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

def plot(plot_time, input_file, save_name, show_ci=False):

    times = {}
    if show_ci:
        cis = {}
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
            if not(tokens[0] in times):
                times[tokens[0]] = {}
            if show_ci:
                if not(tokens[0] in cis):
                    cis[tokens[0]] = {}

            key = int(tokens[key_idx])
            times[tokens[0]][key] = float(tokens[5])
            if show_ci:
                cis[tokens[0]][key] = (float(tokens[7]), float(tokens[8]))
            line = fp.readline()

    for k in times:
        d = times[k]
        od = collections.OrderedDict(sorted(d.items()))
        times[k] = od
        if show_ci:
            d = cis[k]
            od = collections.OrderedDict(sorted(d.items()))
            cis[k] = od

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
        print(save_name,":")
        print(y)
        if not show_ci:
            line=plt.plot(x, y, marker='o', linestyle='-', label=i)
        else:
            d = cis[i]
            ci_lows = [d[j][0] for j in d]
            ci_highs = [d[j][1] for j in d]
            print("conf. intervals :")
            print(list(zip(ci_lows, ci_highs)))
            ci_lows_err = list(np.array(y) - np.array(ci_lows))
            ci_highs_err = list(np.array(ci_highs) - np.array(y))
            print("ci_lows_err:")
            print(ci_lows_err)
            print("ci_highs_err:")
            print(ci_highs_err)

            conf_intervals = [ci_lows_err, ci_highs_err]
            line=plt.errorbar(x, y, xerr=0.0, yerr=conf_intervals, marker='o', capsize=6)
            labs += [i]
        lns += line

    if not show_ci:
        labs = [l.get_label() for l in lns]

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

    if show_ci:
        ln = []
        cnt = 0
        for i in lns:
            if cnt % 3 == 0:
                ln.append(i)
            cnt += 1
        lgd = ax.legend(ln, labs, ncol=1, loc=2, fontsize=12)
    else:
        lgd = ax.legend(lns, labs, ncol=1, loc=2, fontsize=12)
    fig.savefig(save_name + '_' + xlabel + '.png')

if __name__ == '__main__':
    if (len(sys.argv) != 2) and (len(sys.argv) != 3):
        print ('Usage: python plot_time.py <file> [OPTIONAL: <show_ci>]')
        exit()
    show_ci = False
    if len(sys.argv) == 3:
        show_ci = True

    if ("time" in sys.argv[1]):
        plot(1, sys.argv[1], "times", show_ci)
    elif ("cycles" in sys.argv[1]):
        plot(0, sys.argv[1], "cycles", show_ci)
    else:
        print ('Usage: python plot_time.py <file> [OPTIONAL: <show_ci>]')
        exit()
