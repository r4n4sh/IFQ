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


FIRST_PHINDEX = 3
LAST_PHINDEX = 23
window_size = 65536
epsilonsrange = range(FIRST_PHINDEX,LAST_PHINDEX)
epsilons = [2**(-x) for x in epsilonsrange] #actually its counters parameters (counters = 1/epsilon)



def plot_memory(memory):
    average_memory= memory
    phis = [2**(-x) for x in epsilonsrange]
    MS = 12
    LW = 4

    plt.plot(phis, average_memory["RAW"],"-.8"  ,label="$RAW$",markersize=MS, linewidth=LW, c="black")
    plt.plot(phis, average_memory["hhh2RSS"],":*"   ,label="$HIT$" ,markersize=MS, linewidth=LW, c="purple")
    plt.plot(phis, average_memory["acc"],"->",label="$ACC_1$", markersize=MS, linewidth=LW, c="blue")
    plt.plot(phis, average_memory["acc1"],"-.s",label="$ACC_2$", markersize=MS, linewidth=LW, c="red")
    plt.plot(phis, average_memory["baseWRSS"],"-v",label="$WCSS$", markersize=MS, linewidth=LW, c="sienna")
    plt.plot(phis, average_memory["acc4"],"--^",label="$ACC_4$", markersize=MS, linewidth=LW, c="yellow")
    plt.plot(phis, average_memory["acc8"],":<",label="$ACC_8$",markersize=MS, linewidth=LW, c="orange")
    plt.plot(phis, average_memory["ECM"],"-D",label="$ECM$", markersize=MS, linewidth=LW, c="darkcyan")


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
    plt.savefig('test_memory.png')
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
    for i in range(1, 32):
        acc.append((5/eps)**(1 + (1.0/ i)) * i  * (4 + log(5/eps, 2)/8) + RSS_ACC1)
    HIT = (5/eps) * log(5/eps, 2) * (4 + log(5/eps, 2)/8)

    WCSS = RSS_WCSS + B_WCSS + b_WCSS
    ACC1 = RSS_ACC1 + ACC1
    ACC2 = RSS_ACC1 + ACC2
    HIT = RSS_ACC1  + HIT
    RAW = (4/eps) * WCSS

    memory["RAW"].append((RAW/ 1e6)) #10^6 is for Mega
    memory["hhh2RSS"].append((HIT/ 1e6))
    memory["baseWRSS"].append((WCSS/ 1e6))
    memory["acc"].append((ACC1 / 1e6))
    memory["acc1"].append((ACC2/ 1e6))
    memory["acc4"].append((acc[4])/ 1e6)
    memory["acc8"].append((acc[8])/ 1e6)

    logOneOverDelta = 10 #corresponds to delta ~ 0.1%
    memory["ECM"].append((logOneOverDelta * eps**(-2) * 32)/ 1e6)


def myprint():
    i = 0
    print "##########   HIT   #####################"
    for x in epsilonsrange:
        print 2**(-x), "**", -x,"**", memory["hhh2RSS"][i]
        i = i + 1
    i = 0
    print "##########   ACC_1   #####################"
    for x in epsilonsrange:
        print 2**(-x), "**", -x,"**", memory["acc"][i]
        i = i + 1
    i = 0
    print "##########   ACC_2   #####################"
    for x in epsilonsrange:
        print 2**(-x), "**", -x,"**", memory["acc1"][i]
        i = i + 1
    i = 0
    print "##########   ACC_4   #####################"
    for x in epsilonsrange:
        print 2**(-x), "**", -x,"**", memory["acc4"][i]
        i = i + 1
    i = 0
    print "##########   ACC_8   #####################"
    for x in epsilonsrange:
        print 2**(-x), "**", -x,"**", memory["acc8"][i]
        i = i + 1
    i = 0
    print "##########   WCSS   #####################"
    for x in epsilonsrange:
        print 2**(-x), "**",-x,"**", memory["baseWRSS"][i]
        i = i + 1
    i = 0
    print "##########   RAW   #####################"
    for x in epsilonsrange:
        print 2**(-x), "**", -x,"**", memory["RAW"][i]
        i = i + 1

memory = dict()
memory["RAW"]=[]
memory["hhh2RSS"]=[]
memory["baseWRSS"]=[]
memory["acc"]=[]
memory["acc1"]=[]
memory["acc4"]=[]
memory["acc8"]=[]
memory["ECM"]=[]

for epsilon in epsilons:
    out = calc_memory(epsilon)
myprint()
#print "memory of HIT", memory["hhh2RSS"]
#print "memory of ACC", memory["acc"]
#print "memory of ACC1", memory["acc1"]
#print "memory of ACC4", memory["acc4"]
#print "memory of ACC8", memory["acc8"]
#print "memory of baseWRSS", memory["baseWRSS"]
#print "memory of RAW", memory["RAW"]
#print "memory of ECM", memory["ECM"]

plot_memory(memory)
