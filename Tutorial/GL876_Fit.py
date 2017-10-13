# Using the published parameters, finds best fit for GL876

# Assumes access to SystPy - files are in PATH

import numpy as np
import matplotlib.pyplot as plt

from SystPy import *

minimize= False

k = Kernel()

# Load the sys file for 876
# Data from Keck, Harps, APF, PFS (In that order)
k.addDataFile("GL876.sys", directory="")


# Set the epoch to match that of the paper
k.setEpoch(2450602.093)

k.setPar(K_P_DATA1, 5.0)


# Params for 1.9 day planet
k.addPlanet([K_PER, 1.93778, K_MASS, 6.83/317.828, K_MA, 355, K_ECC, 0.207, K_LOP, 234, K_INC, 59.0, K_DONE])

# Params for 30 day planet
k.addPlanet([K_PER, 30.0881, K_MASS, 0.7142, K_MA, 294.59, K_ECC, 0.25591, K_LOP, 48.76, K_INC, 59.0, K_DONE])

# Params for 60 day planet
k.addPlanet([K_PER, 61.1166, K_MASS, 2.2756, K_MA, 325.7, K_ECC, 0.0324, K_LOP, 50.3, K_INC, 59.0, K_DONE])

# Params for 124 day planet
k.addPlanet([K_PER, 124.26, K_MASS, 14.6/317.828, K_MA, 335, K_ECC, 0.055, K_LOP, 239, K_INC, 59.0, K_DONE])

# Update kernel params for the new system we've added
k.calculate()


# initial minimize of the offsets
for i in range(k.getNsets()):
    k.minimize1d(-1, K_P_DATA1 + i)
k.calculate()


# Set up for integrated fit
# How do we tell whats avilible?
print k.setIntMethod.__doc__
k.setElementType(K_JACOBI)
k.setIntMethod(K_RK89)



# Turn off minimization for all the planets
for i in range(1, k.nPlanets()+1):
    for param in range(K_ELEMENTS_SIZE):
        k.setElementFlag(i, param, K_ACTIVE & K_MINIMIZE)

# Turn on minimization for data offset
for i in range(k.getNsets()):
    k.setParFlag(K_P_DATA1 + i, K_ACTIVE | K_MINIMIZE)

# Minimize the data offsets
if minimize:
    k.minimize()

# Turn on Minimization for all planets (basic elements + inclination)
for i in range(1, k.nPlanets() + 1):
    for param in range(K_PER,K_INC+1):
        k.setElementFlag(i, param, K_ACTIVE | K_MINIMIZE)

# Minimize all the system parameters
if minimize:
    k.minimize()

# Save the resulting fit 
k.save("GL876.fit")


# Generate a periodogram of residuals.
# Args: data matrix, samples, min per, max per, n/a, time col, val col, err col
# getCompiledMatrix(k)
per = periodogram_ls(getResidualMatrix(k), 20000, 0.5, 10000, 0, K_T_TIME, K_T_SVAL, K_T_ERR)

plt.plot(per[:,K_PS_TIME], per[:,K_PS_Z], c='black')

plt.xscale('log')

plt.xlabel("Period", fontsize=18)
plt.ylabel("Power", fontsize=18)

plt.xlim([0.5, 10000])

plt.savefig("ResidualPeriodogram.pdf", bbox_inches='tight')












