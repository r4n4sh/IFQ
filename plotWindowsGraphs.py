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
    #phindices = np.arange(21-len(average_speeds["ALS"]),21)
    windindices = range(minWind,maxWind)
    #phis = [2**(-x) for x in reversed(phindices)]
    phis = [2**(x) for x in windindices]
    #plt.plot(phis, average_speeds["ALS"], "-D", label="IM-SUM")
    #plt.plot(phis, average_speeds["LS"],"-d",label="DIM-SUM")
    MS = 12
    LW = 4
    # try:
        # del average_speeds["CMH"]
    # except:
        # pass
    #print "Plotting !", len(phis), len(average_speeds["full2w"])
    #plt.plot(phis, average_speeds["full2w"],"-->"	,label="Full Ancestry", markersize=MS, linewidth=LW, c="brown")
    #plt.plot(phis, average_speeds["ancestry2w"],"--s"	,label="Partial Ancestry", markersize=MS, linewidth=LW, c="black")
    #plt.plot(phis, average_speeds["hhh2w"],"-o"	,label="HSSH", markersize=MS, linewidth=LW, c="red")
    print "values of x"
    print phis
    print "values of y"
    print average_speeds["hhh2RSS"]
    plt.plot(phis, average_speeds["hhh2RSS"],"-.D"	,label="Interval_Queries", markersize=MS, linewidth=LW, c="blue")
    #plt.plot(phis, average_speeds["CMH"],"-^",label="CMH", markersize=MS, linewidth=LW, c="purple")
    #plt.plot(phis, average_speeds["CCFC"],"-v",label="CS", markersize=MS, linewidth=LW, c="sienna")


    ##for algorithm in speed:
    ##    plt.plot(phis, speed[algorithm],label=algorithm)
    plt.xscale("log",basex=2)
    #ticks,labels = plt.xticks()
    #plt.xticks(ticks[::2],labels[::2])
    plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=2))

    plt.xlabel("$W$", fontsize=36)
    ylabel_str = sys.argv[1] + "/second [Thousands]"
    plt.ylabel(ylabel_str, fontsize=24)
    plt.tick_params(labelsize=20)
    plt.xlim(right = 2**(maxWind-1), left=2**minWind)
    #plt.ylim(0.0, 2500)
    plt.ylim(500.0, 20000)
    #plt.ylim(0.0, 45.)
    #plt.ylim(down = 0)
    #topy = int((max(average_speeds["ALS"])*2)/50000+1)*50000
    #plt.ylim(top = topy)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig('graphs/' + dataset + '_' + str(minWind) + '-' + str(maxWind) + sys.argv[1] +'_windows.png')
    plt.clf()

    #fig = plt.gcf()
    #fig.canvas.set_window_title("{} {} Runs Gamma 4".format(dataset,RUNS))


csvFiles = dict()
csvFiles['mychicago'] = 'MyChicago_windows.txt'
#csvFiles['mychicago'] = 'mytestfile.txt'
#csvFiles['SJ14'] = 'SJ14_RSS_HHH_measurement.HeavyHitter.txt'
#csvFiles['SJ13'] = 'SJ13_RSS_HHH_measurement.HeavyHitter.txt'
#csvFiles['CH16'] = 'Chicago16_RSS_HHH_measurement.HeavyHitter.txt'
#csvFiles['CH15'] = 'Chicago15_RSS_HHH_measurement.HeavyHitter.txt'

minWind = 10
maxWind = 20
speeds = dict()


for dataset in csvFiles:
	lines = open(csvFiles[dataset]).readlines()
	speeds[dataset] = dict()
	for line in lines:
		alg = line.split(' ')[0].replace('./','') # algo = hhh2RSS
		if (alg not in speeds[dataset]):
			speeds[dataset][alg]=[]
		timing = re.search("\s(\d+\.\d+)s\s", line).groups()[0] # all timing values
        	speeds[dataset][alg].append((1e7 / float(timing)) / 1e3)
        	#speeds[dataset][alg].append((1e3 * float(timing)) / 1e7) #TODO: the correct calc i think!!
		print (1e7 / float(timing)) / 1e3
	plot_speeds(speeds[dataset], dataset)
