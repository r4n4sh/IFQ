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
    phindices = range(minEps,maxEps)
    acc_phindices = range(4,13)
    phis = [2**(-x) for x in phindices]
    acc_phis = [2**(-x) for x in acc_phindices]
    MS = 12
    LW = 4


    print "values of x", phis
    print "values of acc x", acc_phis
    print "value of y", average_speeds["acc4"]
    #print "values of y HIT", average_speeds["hhh2RSS"]
    #print "values of y ACC", average_speeds["acc"]
    #print "values of y ACC1", average_speeds["acc1"]
    #print "values of y baseWRSS", average_speeds["baseWRSS"]

    plt.plot(acc_phis, average_speeds["acc4"],"-o",label="$ACC_4$", markersize=MS, linewidth=LW, c="red")

    #plt.xscale("log",basex=2)
    #plt.yscale("log",basex=2)
    plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=2))

    plt.xlabel("Accuracy Guarantee ($\epsilon$)", fontsize=36)
    ylabel_str = "Error"
    plt.ylabel(ylabel_str, fontsize=24)
    plt.tick_params(labelsize=20)
    plt.xlim(left = 2**-(maxEps-1), right=2**-minEps)
    plt.ylim(0, 200) #For Queries
    plt.legend(loc="best") # keys of the graphs
    plt.tight_layout()
    plt.savefig('test.png')
    plt.clf()



csvFiles = dict()
csvFiles['Chicago16'] ='tmp.txt'
#csvFiles['SanJose14'] = sys.argv[1] + '/MySanJose_epsilon.txt'
#csvFiles['Univ1'] = sys.argv[1] + '/MyUniv1_epsilon.txt'
#csvFiles['Univ2'] = sys.argv[1] + '/MyUniv2_epsilon.txt'
#csvFiles['UCLA_TCP'] = sys.argv[1] + '/MyUCLA_TCP_epsilon.txt'
#csvFiles['UCLA_UDP'] = sys.argv[1] + '/MyUCLA_UDP_epsilon.txt'

#csvFiles['Chicago16'] = 'Queries/MyChicago_epsilon_query.txt'
#csvFiles['Chicago16'] = 'Updates/MyChicago_epsilon_update.txt'

minEps = 4
maxEps = 13
speeds = dict()


for dataset in csvFiles:
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
