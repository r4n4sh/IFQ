#!/usr/bin/python
import subprocess
import random
import math
import re
import cPickle
import sys


def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
            return None
    if len(lst) %2 == 1:
            return lst[((len(lst)+1)/2)-1]
    else:
            return float(sum(lst[(len(lst)/2)-1:(len(lst)/2)+1]))/2.0

#datasets = ['Zipf1',  'YouTube', 'CAIDA', 'Chicago16', 'SanJose14', 'Chicago15', 'SanJose13']
#datasets = [ 'SanJose13','Chicago16', 'SanJose14', 'Chicago15',  'CAIDA']
#datasets = ['Zipf1.3','Zipf0.7','Zipf1']
#datasets = ['Chicago15','SanJose14','Chicago16','SanJose13','Zipf1.3','Zipf0.7','Zipf1',  'YouTube']
#datasets = ['SanJose13','Zipf1.3','Zipf0.7','Zipf1',  'YouTube']
#datasets = ['YouTube']
#datasets = ['SanJose13', 'YouTube','Chicago16', 'SanJose14', 'Chicago15']
datasets = ['Chicago16']
#datasets = ['SanJose14']
#datasets = ['Univ_1']
#datasets = ['Univ_2']
#datasets = ['Chicago16', 'SanJose14', 'Univ_1', 'Univ_2']
#datasets = ['Chicago16', 'SanJose14', 'Univ_1', 'Univ_2', 'UCLA_TCP', 'UCLA_UDP']
#datasets = ['UCLA_TCP', 'UCLA_UDP']

for dataset in datasets:

        target_file = None
        M = 1
        f = open("result/" + dataset+ sys.argv[1] + "_emp_error"+".txt", "w+")
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

        FIRST_PHINDEX = 4
        LAST_PHINDEX = 12
        #LAST_PHINDEX = 14 #for RAW
        gamma = 1
        #n = 10000000
        n = 34095360
        threshold = 1
        window_size = 1048576
        eprange = dict()

        eprange["hit"]= range(10,15)
        eprange["acc1"]=range(6,9)
        eprange["acc2"]=range(8,12)
        eprange["acc4"]=range(9,14)
        eprange["acc8"]=range(9,15)

        epsilonsrange = eprange[sys.argv[1]];
        print epsilonsrange
        epsilons = [2**x for x in epsilonsrange] #actually its counters parameters (counters = 1/epsilon)
        fn = "./hhh2RSS"
        speeds = []
        range = 1
        k_algo = str(1)

        results = {}
        first_iteration = True
        for epsilon in epsilons:
                timing_list = []
                if (sys.argv[1].find("acc") != -1):
                    k_algo = re.sub(r'\D', "", sys.argv[1])

                command = [fn,"-c", str(epsilon)]
                if target_file:
                    command += ["-np", str(n), "-t", str(threshold), "-f", target_file, "-M", str(M), "-gamma", str(gamma), "-w", str(window_size), "-k", k_algo]
                print " ".join(command)
                out = subprocess.check_output(command)
                if (range == 0):
                    timing = re.search("\s(\d+\.\d+)\s", out).groups()[0]
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
