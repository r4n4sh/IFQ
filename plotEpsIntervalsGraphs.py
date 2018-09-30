#!/usr/bin/python
import subprocess
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import ticker
import numpy as np
import re
import cPickle
import sys

def plot_speeds(speeds, dataset):
    average_speeds = speeds
    MS = 12
    LW = 2
    intervals = [1,5,10,15,30,50]
    print intervals
    print average_speeds["acc1"]
    print average_speeds["acc2"]
    print average_speeds["acc4"]
    print average_speeds["acc8"]
    print average_speeds["hit"]
    print average_speeds["cm"]

    plt.plot(intervals, average_speeds["hit"],":*" ,label="$HIT$", markersize=MS, markeredgecolor="k",linewidth=LW, c="purple")
    plt.plot(intervals, average_speeds["acc1"],"->",label="$ACC_1$", markersize=MS, markeredgecolor="k",linewidth=LW, c="blue")
    plt.plot(intervals, average_speeds["acc2"],"-.s",label="$ACC_2$", markersize=MS, markeredgecolor="k",linewidth=LW, c="red")
    plt.plot(intervals, average_speeds["acc4"],"--^",label="$ACC_4$", markersize=MS, markeredgecolor="k",linewidth=LW, c="yellow")
    plt.plot(intervals, average_speeds["acc8"],":<",label="$ACC_8$", markersize=MS, markeredgecolor="k",linewidth=LW, c="orange")
    plt.plot(intervals, average_speeds["cm"],"-D",label="$ECM$", markersize=MS, markeredgecolor="k",linewidth=LW, c="darkcyan")

    #plt.xscale("log",basex=2)
    
    #plt.yscale("log",basex=2)
    #plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=2))

    plt.xlabel("Interval size (% of window size)", fontsize=20)
    ylabel_str = "Queries/seconds [Millions]"
    plt.ylabel(ylabel_str, fontsize=20)
    plt.tick_params(labelsize=10)
    plt.xlim(1, 50)
    plt.ylim(0, 2) #For Queries
    plt.legend(loc="best", numpoints=2) # keys of the graphs
    plt.tight_layout()
    plt.savefig('vary_intervals_sizes.png')
    plt.clf()


csvFiles = dict()
#csvFiles['Chicago16'] = sys.argv[1] + '/MyChicago_' + sys.argv[1]+ '_epsilon.txt'
#csvFiles['SanJose14'] = sys.argv[1] + '/MySanJose_' + sys.argv[1]+ '_epsilon.txt'
#csvFiles['Univ1'] = sys.argv[1] + '/MyUniv1_' + sys.argv[1]+ '_epsilon.txt'
#csvFiles['Univ2'] = sys.argv[1] + '/MyUniv2_' + sys.argv[1]+ '_epsilon.txt'
#csvFiles['UCLA_TCP'] = sys.argv[1] + '/MyUCLATCP_' + sys.argv[1]+ '_epsilon.txt'
csvFiles['UCLA_UDP'] = 'test_intervals.txt'

minEps = 4
maxEps = 14
speeds = dict()

for dataset in csvFiles:
    lines = open(csvFiles[dataset]).readlines()
    speeds[dataset] = dict()
    for line in lines:
        alg = line.split(' ')[0].replace('./','') # algo = hhh2RSS
        if (alg not in speeds[dataset]):
            speeds[dataset][alg]=[]
        ops = re.search(r'(\d+) pairs', line).groups()[0]
        timing = re.search("\s(\d+\.\d+)s\s", line).groups()[0] # all timing values
        speeds[dataset][alg].append((float(ops) / float(timing)) / 1e6) #For Queries in Millions
    plot_speeds(speeds[dataset], dataset)
