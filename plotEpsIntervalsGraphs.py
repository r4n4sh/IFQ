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
    HIT_values = sorted(average_speeds["hhh2RSS"])
    yErr = []
    yVals = []
    t_n_1 = 2.262 # t statistic for 9 degrees of freedom @ 0.95 confidence two tail
    for eps in HIT_values:
        	print eps
    		mean = sum(average_speeds["hhh2RSS"][eps]) / len(average_speeds["hhh2RSS"][eps])
    		yVals.append(mean)
    		sos = sum([(x - mean)**2 for x in average_speeds["hhh2RSS"][eps]])
    		variance = (sos / (len(average_speeds["hhh2RSS"][eps]) - 1) )# ** 0.5
    		err = variance * t_n_1 / pow(len(average_speeds["hhh2RSS"][eps]), 0.5)
    		yErr.append(err)
    algLabel = 'HIT'
    print '\n',yVals,dataset
    plt.errorbar(HIT_values, yVals, yerr=yErr,fmt='-.D' ,label=algLabel, color="red", capsize=10, elinewidth=3, capthick=1.5)


    ACC_values = sorted(average_speeds["acc"])
    yErr = []
    yVals = []
    t_n_1 = 2.262 # t statistic for 9 degrees of freedom @ 0.95 confidence two tail
    for eps in ACC_values:
        	print eps
    		mean = sum(average_speeds["acc"][eps]) / len(average_speeds["acc"][eps])
    		yVals.append(mean)
    		sos = sum([(x - mean)**2 for x in average_speeds["acc"][eps]])
    		variance = (sos / (len(average_speeds["acc"][eps]) - 1) )# ** 0.5
    		err = variance * t_n_1 / pow(len(average_speeds["acc"][eps]), 0.5)
    		yErr.append(err)
    algLabel = 'ACC_1'
    print '\n',yVals,dataset
    plt.errorbar(ACC_values, yVals, yerr=yErr,fmt='-^' ,label=algLabel, color= "cyan", capsize=10, elinewidth=3, capthick=1.5)


    ACC2_values = sorted(average_speeds["acc1"])
    yErr = []
    yVals = []
    t_n_1 = 2.262 # t statistic for 9 degrees of freedom @ 0.95 confidence two tail
    for eps in ACC2_values:
        	print eps
    		mean = sum(average_speeds["acc1"][eps]) / len(average_speeds["acc1"][eps])
    		yVals.append(mean)
    		sos = sum([(x - mean)**2 for x in average_speeds["acc1"][eps]])
    		variance = (sos / (len(average_speeds["acc1"][eps]) - 1) )# ** 0.5
    		err = variance * t_n_1 / pow(len(average_speeds["acc1"][eps]), 0.5)
    		yErr.append(err)
    algLabel = 'ACC_2'
    print '\n',yVals,dataset
    plt.errorbar(ACC2_values, yVals, yerr=yErr,fmt='-o' ,label=algLabel, color = "blue", capsize=10, elinewidth=3, capthick=1.5)

    ACC2_values = sorted(average_speeds["acc4"])
    yErr = []
    yVals = []
    t_n_1 = 2.262 # t statistic for 9 degrees of freedom @ 0.95 confidence two tail
    for eps in ACC2_values:
        	print eps
    		mean = sum(average_speeds["acc4"][eps]) / len(average_speeds["acc4"][eps])
    		yVals.append(mean)
    		sos = sum([(x - mean)**2 for x in average_speeds["acc4"][eps]])
    		variance = (sos / (len(average_speeds["acc4"][eps]) - 1) )# ** 0.5
    		err = variance * t_n_1 / pow(len(average_speeds["acc4"][eps]), 0.5)
    		yErr.append(err)
    algLabel = 'ACC_4'
    print '\n',yVals,dataset
    plt.errorbar(ACC2_values, yVals, yerr=yErr,fmt='-o' ,label=algLabel, color = "blue", capsize=10, elinewidth=3, capthick=1.5)


    WCSS_values = sorted(average_speeds["baseWRSS"])
    yErr = []
    yVals = []
    t_n_1 = 2.262 # t statistic for 9 degrees of freedom @ 0.95 confidence two tail
    for eps in WCSS_values:
        	print eps
    		mean = sum(average_speeds["baseWRSS"][eps]) / len(average_speeds["baseWRSS"][eps])
    		yVals.append(mean)
    		sos = sum([(x - mean)**2 for x in average_speeds["baseWRSS"][eps]])
    		variance = (sos / (len(average_speeds["baseWRSS"][eps]) - 1) )# ** 0.5
    		err = variance * t_n_1 / pow(len(average_speeds["baseWRSS"][eps]), 0.5)
    		yErr.append(err)
    algLabel = 'WCSS'
    print '\n',yVals,dataset
    plt.errorbar(WCSS_values, yVals, yerr=yErr,fmt='-v' ,label=algLabel, color="green", capsize=10, elinewidth=3, capthick=1.5)


    plt.xscale("log",basex=2)
    plt.yscale("log",basex=2)
    plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=2))
    plt.xlabel("Accuracy Guarantee ($\epsilon$)", fontsize=36)
    ylabel_str = sys.argv[1] + "/second [Millions]" #For Queries in Millions
    plt.ylabel(ylabel_str, fontsize=24)
    plt.tick_params(labelsize=20)
    plt.xlim(left = 2**-(maxEps-1), right=2**-minEps)
    plt.ylim(0, 85) #For Queries
    #plt.legend(loc="best") # keys of the graphs
    #plt.legend(loc="best",prop={'size':16},ncol=20) # keys of the graphs

    plt.tight_layout()
    plt.savefig('graphs/' + sys.argv[1] + '/' + dataset + sys.argv[1] + '_' + str(minEps) + '-' + str(maxEps) +'_epsilon.png')
    plt.clf()



csvFiles = dict()
#csvFiles['Chicago16'] = sys.argv[1] + '/MyChicago_' + sys.argv[1]+ '_epsilon.txt'
#csvFiles['SanJose14'] = sys.argv[1] + '/MySanJose_' + sys.argv[1]+ '_epsilon.txt'
#csvFiles['Univ1'] = sys.argv[1] + '/MyUniv1_' + sys.argv[1]+ '_epsilon.txt'
#csvFiles['Univ2'] = sys.argv[1] + '/MyUniv2_' + sys.argv[1]+ '_epsilon.txt'
#csvFiles['UCLA_TCP'] = sys.argv[1] + '/MyUCLATCP_' + sys.argv[1]+ '_epsilon.txt'
csvFiles['UCLA_UDP'] = sys.argv[1] + '/MyUCLAUDP_' + sys.argv[1]+ '_epsilon.txt'

minEps = 4
maxEps = 14
speeds = dict()


for dataset in csvFiles:
	lines = open(csvFiles[dataset]).readlines()
	speeds[dataset] = dict()
	for line in lines:
		alg = line.split(' ')[0].replace('./','') # algo = hhh2RSS
		if (alg not in speeds[dataset]):
			speeds[dataset][alg]=dict()
        	counters = re.search(r'(\d+) counters', line).groups()[0]
            ops = re.search(r'(\d+) pairs', line).groups()[0]
            	counters = 1 / float(counters);
    		if (counters not in speeds[dataset][alg]):
    			speeds[dataset][alg][counters]=[]
		timing = re.search("\s(\d+\.\d+)s\s", line).groups()[0] # all timing values
        	#speeds[dataset][alg].append((1e7 / float(timing)) / 1e3)
                if (alg == "raw" and sys.argv[1]=="Queries"):
                    timing = float(timing) / 2
                    timing = str(timing)
            	if (alg == "raw"):
                	speeds[dataset][alg].append((ops / float(timing)) / 1e6) #For Queries in Millions
            	else:
        		speeds[dataset][alg][counters].append((ops / float(timing)) / 1e6) #For Queries in Millions
	print speeds[dataset][alg]
	plot_speeds(speeds[dataset], dataset)
