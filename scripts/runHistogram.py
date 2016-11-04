import matplotlib.pyplot as plt
import numpy as np
import os
from os import listdir, path
import sys


####### Usage ###########
#Cool fancy version 2.0
#python2.7 runHistogram.py folder1ctflogs folder2ctflogs folder3ctflogs + .... + foldernctflogs
#These folders MUST be in the CWD (can be symlinked)
#eg: python2.7 runHistogram.py test_liz_simple_driftcorr_5_ctf control_leginon_dfcorr_corrected_ctf
##########################

#### CONFIG #######
TOTAL_BINS = 50 #For Histogram
NORMALIZE_HIST = False #Also for Histogram
###############

def gatherCtfLogs(directory):
    resolutions = list()
    logFiles = [a for a in listdir(os.getcwd() + '/' + directory) if '.txt' in a and 'avrot' not in a]
    for logFile in logFiles:

        logFullText = open(directory + '/' + logFile, 'r').read()
        try:
            resolution = float(logFullText[-15:-6]) #Dont' forget negatives you idiot
            resolutions.append(resolution)
        except Exception, err:
            print('failed to get resolution on file: ' + str(logFile) + ', so just skippin it for now, cuz ' + str(err))
            continue
    return resolutions

if len(sys.argv) < 3:
    print('bad arguments you goof')
    sys.exit()

directories = list()
for i in range(1, len(sys.argv)):
    directories.append(sys.argv[i])

results = list()

for d in directories:
    results.append(gatherCtfLogs(d))

for i in range(0, len(sys.argv) - 1):
    plt.hist(results[i], bins=TOTAL_BINS, histtype='stepfilled', normed=NORMALIZE_HIST, label=directories[i])
#plt.hist(results2, bins=TOTAL_BINS, histtype='stepfilled', normed=NORMALIZE_HIST, color='r', label=ctfDirectory2)
plt.title("CTF Resolutions Comparison")
plt.xlabel("Resolutions (A)")
plt.ylabel("Number of Results")
plt.legend()
plt.show()