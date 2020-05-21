import sys
import matplotlib.pyplot as plt
import collections
import numpy as np
import math

def parse_threads_and_time(file_name):
	res_d = {}
	with open(file_name, 'r') as f:
		for line in f:
			line = line.replace(' ', '')
			fields = line.split(',')

			threads_no = int(fields[1].split(':')[1])
			time = float(fields[2].split(':')[1])

			if threads_no in res_d:
				res_d[threads_no].append(time)
			else:
				res_d[threads_no] = [time]

	return res_d


# plot a bar plot showing the median execution time
# for different number of cores used


def plot_bar_diagram():
    fig_name = 'cycles_compilers.png'
    implemenation = 'VectLoopUnrolling'

    cycles = [1.1, 2.3, 3.5, 4.4]
    comps  = ['gcc [-O3]', 'gcc [-O3 -march=native]', 'icc [-O3]', 'clang [-O3]']

    yticks = []

    fig, ax = plt.subplots(figsize=(11,9))

    plt.xticks(range(len(comps)), comps)
    plt.xlabel("\ncompiler [flags]", fontsize=14)
    ax.yaxis.set_label_coords(0.029,1.02)
    plt.ylabel('cycles',rotation=0, fontsize=14)
    plt.title("Compilers comparison (M = , N = , T = )\n\n" + implemenation, fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(cycles + yticks, fontsize=12)
    plt.bar(range(len(cycles)), cycles, color='darkred', zorder=2)
    for y in cycles:
        plt.axhline(y=y, linewidth=0.8, color='black', linestyle='--', zorder=1)
    plt.savefig(fig_name)

plot_bar_diagram()
