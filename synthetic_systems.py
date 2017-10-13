from random import shuffle

import numpy as np

from SystPy import *


stellar_mass = []
error_term = []
obs_time = []
period = []
mass = []

# Load Kepler Systems
# This is where we get our planets, stellar mass, and observation error
with open('kepSystem','r') as f:
    for line in f:
        ls = line.split()
        if ls[0] == 'System':
            stellar_mass.append(float(ls[1]))
            error_term.append(float(ls[2]))
            obs_time.append(float(ls[3]))
            period.append([])
            mass.append([])
        else:
            period[-1].append(float(ls[0]))
            mass[-1].append(float(ls[1]))

# Create the synth system
# I'm using the kepler system to determine the planets to add
epoch = 2452833.85
kernels = []
for i in range(len(stellar_mass)):
    # Create the kernel
    k = Kernel()
    # Set the stellar mass and epoch
    k.setMstar(stellar_mass[i])
    k.setEpoch(epoch)
    # Add some planets to the system
    for j in range(len(period[i])):
        k.addPlanet([K_PER, period[i][j], K_MASS, mass[i][j], K_DONE])
    k.calculate()
    # Get the predicted RV due to the planets we added
    # Restricts to 1 RV per day
    predRV = k.stellarRV(epoch, epoch+4*365.25, int(8*365.25))

    # Add the data in shuffled for some reason?
    idx = range(int(8*365.25))
    shuffle(idx)
    err = error_term[i] * np.random.randn(100)
    data = []
    for j in range(len(err)):
        data.append([predRV[idx[j],0], predRV[idx[j],1] + err[j], error_term[i]])

    # We just created the data set, so add it to the kernel
    k.addDataArray(np.array(data)) 
    k.calculate()
    # Remove the starting planets we added
    k.removePlanet(-1)
    k.calculate()
    kernels.append(k)
     
print "Finished Data Set Creation"


correct_planets = 0
false_planets = 0
    
for i in range(len(stellar_mass)):
    # This figures we know the period from a transit
    # if we didn't, would just use fitSystem(k[i])
    for j in range(len(period[i])):
        kernels[i].addPlanet([K_PER, period[i][j], K_DONE])
        for w in range(K_ELEMENTS_SIZE):
            kernels[i].setElementFlag(j, w, K_ACTIVE & K_MINIMIZE)
        kernels[i].setElementFlag(j, K_MASS, K_ACTIVE | K_MINIMIZE)
        kernels[i].setElementFlag(j, K_MA, K_ACTIVE | K_MINIMIZE)
        kernels[i].calculate()
        kernels[i].minimize1d(K_SIMPLEX, 5000, j, K_MA)
        kernels[i].minimize1d(K_SIMPLEX, 5000, j, K_MASS)

    kernels[i].setParFlag(K_DATA_NOISE1, K_ACTIVE | K_MINIMIZE)

    kernels[i].calculate()
    kernels[i].minimize()

    for j in range(len(period[i])):
        m1 = mass[i][j]
        m2 = kernels[i].getElement(j,K_MASS)
        up = m1 + m1*0.1
        down = m1 - m1*0.1
        if m2 < up and m2 > down: correct_planets += 1
        else: false_planets += 1

    kernels[i].removePlanet(-1)
    kernels[i].calculate()

            
print correct_planets
print false_planets            
        
        

