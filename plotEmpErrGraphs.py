#!/usr/bin/python
import subprocess
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import ticker
import numpy as np
import re
import cPickle
import sys
import numpy
from math import log
from math import sqrt


FIRST_PHINDEX = 3
LAST_PHINDEX = 14
window_size = 65536
epsilonsrange = range(FIRST_PHINDEX,LAST_PHINDEX)
epsilons = [2**(-x) for x in epsilonsrange] #actually its counters parameters (counters = 1/epsilon)


def calc_memory(eps, alg):
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

    if (alg == "RAW"):
        memory["RAW"].append((RAW/ 1e6)) #10^6 is for Mega
    if (alg == "hit"):
        memory["hit"].append((HIT/ 1e6))
    if (alg == "baseWRSS"):
        memory["baseWRSS"].append((WCSS/ 1e6))
    if (alg == "acc1"):
        memory["acc1"].append((ACC1 / 1e6))
    if (alg == "acc2"):
        memory["acc2"].append((ACC2/ 1e6))
    if (alg == "acc4"):
        memory["acc4"].append((acc[4])/ 1e6)
    if (alg == "acc8"):
        memory["acc8"].append((acc[8])/ 1e6)
    if (alg == "ECM"):
        logOneOverDelta = 10 #corresponds to delta ~ 0.1%
        memory["ECM"].append((logOneOverDelta * eps**(-2) * 32)/ 1e6)



def memory_consumption(alg):
    epsilons = [2**(-x) for x in eprange[alg]] #actually its counters parameters (counters = 1/epsilon)

    for epsilon in epsilons:
        out = calc_memory(epsilon, alg)


def plot_speeds(speeds, dataset):
    average_speeds = speeds
    MS = 12
    LW = 2


    memory_consumption("hit")
    memory_consumption("acc1")
    memory_consumption("acc2")
    memory_consumption("acc4")
    memory_consumption("acc8")
    #memory_consumption("ECM")

    print "values of x acc8", memory["acc8"]
    print "value of y acc8", average_speeds["acc8"]
    print "values of x acc4", memory["acc4"]
    print "value of y acc4", average_speeds["acc4"]
    print "values of x acc2", memory["acc2"]
    print "value of y acc2", average_speeds["acc2"]
    print "values of x acc1", memory["acc1"]
    print "value of y acc1", average_speeds["acc1"]
    print "values of x hit", memory["hit"]
    print "value of y hit", average_speeds["hit"]


    plt.plot(memory["hit"], average_speeds["hit"],":*",label="$HIT$", markersize=MS, linewidth=LW, markeredgecolor="k", c="purple")
    plt.plot(memory["acc1"], average_speeds["acc1"],"->",label="$ACC_1$", markersize=MS, linewidth=LW,markeredgecolor="k",c="blue")
    plt.plot(memory["acc2"], average_speeds["acc2"],"-.s",label="$ACC_2$", markersize=MS, linewidth=LW, markeredgecolor="k",c="red")
    plt.plot(memory["acc4"], average_speeds["acc4"],"--^",label="$ACC_4$", markersize=MS, linewidth=LW,  markeredgecolor="k",c="yellow")
    plt.plot(memory["acc8"], average_speeds["acc8"],":<",label="$ACC_8$", markersize=MS, linewidth=LW, markeredgecolor="k",c="orange")

    #plt.xscale("log",basex=2)
    #plt.yscale("log",basex=2)
    plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=2))

    plt.xlabel("Memory (Mbytes)", fontsize=36)
    ylabel_str = "Observed Error"
    plt.ylabel(ylabel_str, fontsize=24)
    plt.tick_params(labelsize=10)
    plt.xlim(1, 15)
    plt.ylim(0, 3000) #For Queries
    plt.legend(loc="best", numpoints=2) # keys of the graphs
    plt.tight_layout()
    plt.savefig('emp_error.png')
    plt.clf()



csvFiles = dict()
csvFiles['Chicago16'] ='merged_emp_error.txt'
#csvFiles['SanJose14'] = sys.argv[1] + '/MySanJose_epsilon.txt'
#csvFiles['Univ1'] = sys.argv[1] + '/MyUniv1_epsilon.txt'
#csvFiles['Univ2'] = sys.argv[1] + '/MyUniv2_epsilon.txt'
#csvFiles['UCLA_TCP'] = sys.argv[1] + '/MyUCLA_TCP_epsilon.txt'
#csvFiles['UCLA_UDP'] = sys.argv[1] + '/MyUCLA_UDP_epsilon.txt'

#csvFiles['Chicago16'] = 'Queries/MyChicago_epsilon_query.txt'
#csvFiles['Chicago16'] = 'Updates/MyChicago_epsilon_update.txt'

minEps = 3
maxEps = 14
speeds = dict()
eprange = dict()
eprange["hit"]= range(10,14)
eprange["acc1"]=range(6,10)
eprange["acc2"]=range(8,12)
eprange["acc4"]=range(9,14)
eprange["acc8"]=range(9,15)


for dataset in csvFiles:
    memory = dict()
    memory["RAW"]=[]
    memory["hit"]=[]
    memory["baseWRSS"]=[]
    memory["acc1"]=[]
    memory["acc2"]=[]
    memory["acc4"]=[]
    memory["acc8"]=[]
    memory["ECM"]=[]


    lines = open(csvFiles[dataset]).readlines()
    speeds[dataset] = dict()
    for line in lines:
        alg = line.split(' ')[0].replace('./','') # algo = hhh2RSS
        print alg
        if (alg not in speeds[dataset]):
            speeds[dataset][alg]=[]
        timing = re.search("\s(\d+\.\d+)\s", line).groups()[0]
        print timing
        speeds[dataset][alg].append(float(timing))
    plot_speeds(speeds[dataset], dataset)
