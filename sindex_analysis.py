#! /usr/bin/env python

import sys
sys.path.append("/Users/jaburt/Desktop/219134_Stuff/HD219134_SINDEX/SystPy")
from SystPy import *
import numpy as np

import matplotlib.pyplot as plt

keck_vels_file = 'HD219134_KECK.vels'
apf_vels_file  = 'HD219134_APF.vels'

keck_data = []
apf_data = []

# Parse the Keck vels file
with open(keck_vels_file,'r') as f:
    for line in f:
        ls = line.split()
        # Don't take lines with an invalid sindex value
        if float(ls[3]) == -1: continue
        row = [float(ls[0]), float(ls[3]), float(ls[2]), float(ls[1])]
        keck_data.append(row)

with open(apf_vels_file,'r') as f:
    for line in f:
        ls = line.split()
        # Don't take lines with an invalid sindex value
        if float(ls[3]) == -1: continue
        row = [float(ls[0]), float(ls[3]), float(ls[2]), float(ls[1])]
        apf_data.append(row)

# Convert to numpy arrays
keck_data = np.array(keck_data)
apf_data = np.array(apf_data)

# Try to normalize the data (might not want to?)
# The shorter baseline of APF data combined with the low errors sort of blow the 
# correlation away. Not sure if thats because it isn't there, or there just aren't
# enough APF points yet.
keck_data[:,1] = keck_data[:,1] - np.median(keck_data[:,1])
apf_data[:,1] = apf_data[:,1] - np.median(apf_data[:,1])


k = Kernel()

k.addDataArray(keck_data)
k.addDataArray(apf_data)

# Generare a periodogram for the kernel k, args = [kernel, num_samples, min_period, max_period, type, time_col, rv_col, err_col]
per = periodogram_ls(getCompiledMatrix(k), 1000000, 0.5, 10000, 0, 0,1,2)

# Plot the periodogram
plt.plot(per[:,0], per[:,1], c='black', lw=2)

peak = per[:,1].argsort()[-1]

print "Period of highest peak", per[peak,0], "(days)"
print "Power at highest peak", per[peak,1]

plt.text(per[peak,0], per[peak,1], "%6.2f" % per[peak,0], ha='right')

plt.xlim([0.5,10000])
plt.ylim([0,100])

plt.xscale('log')

plt.xlabel('Period', fontsize=18)
plt.ylabel("Power",fontsize=18)

plt.savefig("sindex_periodogram.pdf", bbox_inches='tight')

plt.clf()

# Combined and sorted Data
tot_time = np.concatenate(([keck_data[:,0], apf_data[:,0]]))
tot_s  = np.concatenate(([keck_data[:,1], apf_data[:,1]]))
tot_rv = np.concatenate(([keck_data[:,3], apf_data[:,3]]))
sort = tot_time.argsort()
tot_s = tot_s[sort]
tot_rv = tot_rv[sort]


# Correlation for only the APF Data
# x = np.correlate(apf_data[:,3], apf_data[:,1], mode='full')

# Correlation for only KECK Data
# x = np.correlate(keck_data[:,3], keck_data[:,1], mode='full')

# Correlation for combined data sets
x = np.correlate(tot_rv, tot_s, mode='full')


lag = np.arange(len(x)) - (len(x)-1)/2
plt.plot(lag, x)
plt.plot([0,0],[max(x), min(x)], 'r--')
plt.plot([min(lag),max(lag)], [0,0], 'g--')

plt.show()










