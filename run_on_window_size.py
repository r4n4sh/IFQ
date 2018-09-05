#!/usr/bin/python
import subprocess
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import ticker
import numpy as np
import re
import cPickle
import numpy

def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
            return None
    if len(lst) %2 == 1:
            return lst[((len(lst)+1)/2)-1]
    else:
            return float(sum(lst[(len(lst)/2)-1:(len(lst)/2)+1]))/2.0

# dataset = raw_input(
    # '''
    # Please pick dataset from the following options:
    # UCLA-UDP
    # UCLA-TCP
    # YouTube
    # Zipf1
    # CAIDA
    # SanJose
    # '''
    # )

#datasets = ['Zipf1',  'YouTube', 'CAIDA', 'Chicago16', 'SanJose14', 'Chicago15', 'SanJose13']
#datasets = [ 'SanJose13','Chicago16', 'SanJose14', 'Chicago15',  'CAIDA']
#datasets = ['Zipf1.3','Zipf0.7','Zipf1']
#datasets = ['Chicago15','SanJose14','Chicago16','SanJose13','Zipf1.3','Zipf0.7','Zipf1',  'YouTube']
#datasets = ['SanJose13','Zipf1.3','Zipf0.7','Zipf1',  'YouTube']
#datasets = ['YouTube']
#datasets = ['SanJose13', 'YouTube','Chicago16', 'SanJose14', 'Chicago15']
#datasets = ['Chicago16']
#datasets = ['SanJose14']
#datasets = ['Univ_1']
datasets = ['Univ_2']
datasets = ['Chicago16', 'SanJose14', 'Univ_1', 'Univ_2']
datasets = ['UCLA_TCP', 'UCLA_UDP']

for dataset in datasets:

	target_file = None
	M = 1
	f = open(dataset+"_window"+".txt", "w+")
	if dataset == "Chicago16":
		target_file = "/home/ranas/wrss/Chicago16Weighted.txt"
		M = 1
	elif dataset == "SanJose14":
		target_file = "/home/ranas/wrss/SJ14Weighted.txt"
		M = 1
	elif dataset == "Univ_1":
		target_file = "/home/ranas/wrss/univ1.ifqFormat"
		M = 1
	elif dataset == "Univ_2":
		target_file = "/home/ranas/wrss/univ2.ifqFormat"
		M = 1
    elif dataset == "UCLA_TCP":
    	target_file = "/home/ranas/wrss/UCLA_TCP_trace.ifqFormat"
    	M = 1
	elif dataset == "UCLA_UDP":
		target_file = "/home/ranas/wrss/UCLA_UDP_trace.ifqFormat"
		M = 1

	##answer = raw_input("run on real trace? (y/n)")
	##if answer == 'y':
	##    target_file = "../../data/UCLA_UDP_trace"
	##else:
	##    target_file = None

	RUNS = 11
	FIRST_WINDINDEX = 10
	LAST_WINDINDEX = 20
	#LAST_PHINDEX = 12 #for ACC
	#LAST_PHINDEX = 9 #for RAW
	gamma = 4
    	n = 10000000
    	threshold = 1
    	counters = 256
	windowsrange = range(FIRST_WINDINDEX,LAST_WINDINDEX)
	windows = [2**x for x in windowsrange]
	fn = "./hhh2RSS"
	speeds = []
    	range = 1

	results = {}
	first_iteration = True

	for window in windows:
        timing_list = []
        for run in xrange(RUNS):
		    command = [fn,"-w", str(window)]
            print command
            if target_file:
				command += ["-np", str(n), "-t", str(threshold), "-f", target_file, "-M", str(M), "-gamma", str(gamma), "-c", str(counters)]
            print " ".join(command)
            out = subprocess.check_output(command)
            if (range == 0):
                    timing = re.search("\s(\d+\.\d+)s\s", out).groups()[0]
                    timing_list.append(timing)
            else:
                    f.write(out)
                    print out
        if (range == 0):
                medianv = median(timing_list)
                result = out.split();
                result[4] = str(medianv) +"s"
                output = ' '.join(map(str, result))
                print output
                f.write(output)
                f.write('\n')
