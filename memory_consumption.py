#!/usr/bin/python
import subprocess
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import ticker
import numpy as np
import re
import cPickle
import numpy
from math import log
from math import sqrt


FIRST_PHINDEX = 4
LAST_PHINDEX = 14
window_size = 65536
epsilonsrange = range(FIRST_PHINDEX,LAST_PHINDEX)
epsilons = [2**(-x) for x in epsilonsrange] #actually its counters parameters (counters = 1/epsilon)



def plot_memory(memory):
    average_memory= memory
    phis = [2**(-x) for x in epsilonsrange]
    MS = 10
    LW = 4

    plt.plot(phis, average_memory["RAW"],"-.D"	,label="$RAW$", markersize=MS, linewidth=LW, c="black")
    plt.plot(phis, average_memory["hhh2RSS"],"-.D"	,label="HIT", markersize=MS, linewidth=LW, c="red")
    plt.plot(phis, average_memory["acc"],"-^",label="ACC_1", markersize=MS, linewidth=LW, c="cyan")
    plt.plot(phis, average_memory["acc1"],"-o",label="ACC_2", markersize=MS, linewidth=LW, c="blue")
    plt.plot(phis, average_memory["baseWRSS"],"-v",label="WCSS", markersize=MS, linewidth=LW, c="green")
    #plt.plot(phis, average_memory["k"],"-o",label="K", markersize=MS, linewidth=LW, c="yellow")


    ##for algorithm in speed:
    ##    plt.plot(phis, speed[algorithm],label=algorithm)
    plt.xscale("log",basex=2)
    plt.yscale("log",basex=2)

    #ticks,labels = plt.xticks()
    #plt.xticks(ticks[::2],labels[::2])
    plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=2))

    plt.xlabel("Accuracy Guarantee ($\epsilon$)", fontsize=36)
    ylabel_str = "Memory [MB]"
    plt.ylabel(ylabel_str, fontsize=24)
    plt.tick_params(labelsize=20)
    plt.xlim(left = 2**-(LAST_PHINDEX-1), right=2**-FIRST_PHINDEX)
    plt.ylim(0, 80000000)
    plt.legend(loc="best") # keys of the graphs
    plt.tight_layout()
    plt.savefig('memory.png')
    plt.clf()

def calc_memory(eps):
    RSS_WCSS = 108 * (4/eps)
    RSS_HIT = 108 * (5/eps)
    RSS_ACC1 = 108 * (5/eps)

    B_WCSS = (4/eps) * (4 + ((log((4/eps), 2))/8))

    b_WCSS = 4 * (4/eps)

    ACC1 = (5/eps) * (5/eps) * (4 + log(5/eps, 2)/8)
    ACC2 = (5/eps) * 2 * sqrt((5/eps)) * (4 + log(5/eps, 2)/8)
    acc = []
    for i in range(3, 32):
        acc.append((5/eps)**(1 + (1.0/ i)) * i  * (4 + log(5/eps, 2)/8) + RSS_ACC1)
    HIT = (5/eps) * log(5/eps, 2) * (4 + log(5/eps, 2)/8)

    WCSS = RSS_WCSS + B_WCSS + b_WCSS
    ACC1 = RSS_ACC1 + ACC1
    ACC2 = RSS_ACC1 + ACC2
    HIT = RSS_ACC1  + HIT
    RAW = (4/eps) * WCSS

    memory["RAW"].append((RAW/ 1e6))
    memory["hhh2RSS"].append((HIT/ 1e6))
    memory["baseWRSS"].append((WCSS/ 1e6))
    memory["acc"].append((ACC1 / 1e6))
    memory["acc1"].append((ACC2/ 1e6))
    memory["k"].append((acc[5])/ 1e6)

memory = dict()
memory["RAW"]=[]
memory["hhh2RSS"]=[]
memory["baseWRSS"]=[]
memory["acc"]=[]
memory["acc1"]=[]
memory["k"]=[]

for epsilon in epsilons:
	out = calc_memory(epsilon)

print "memory of y HIT", memory["hhh2RSS"]
print "memory of y ACC", memory["acc"]
print "memory of y ACC1", memory["acc1"]
print "memory of y baseWRSS", memory["baseWRSS"]
print "memory of y RAW", memory["RAW"]

plot_memory(memory)
