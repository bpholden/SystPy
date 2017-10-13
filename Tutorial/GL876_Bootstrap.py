# generate bootstrap members for GL876

import numpy as np

from SystPy import * 


# Load the previously created fit file
k = loadKernel('GL876.fit')
k.calculate()

# Perform the bootstrap analysis
kl = bootstrap(k, 2, 0)

KLSave(kl, 'GL876Bootstrap.kl')

print kl.getElementsStats(K_STAT_MEDIAN)

print kl.getElementsStats(K_STAT_STDDEV)
