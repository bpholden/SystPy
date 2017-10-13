# Integrating the GL876 System

import numpy as np
import matplotlib.pyplot as plt

from SystPy import *



k = loadKernel("GL876.fit")

# Set up color array for the RV points
c = ['b','g','r','Gold','black','purple']
colors = []
times = np.empty(0)
for i in range(k.getNsets()):
    data = matrix_to_array(k.getData(i))[:,K_T_TIME]
    times = np.append(times, data)
    colors += [c[i]]*len(data)


colors = np.array(colors)
if len(colors) != len(times):
    print "Size miss-match"
else:
    idx = times.argsort()
    colors = colors[idx]
    
    


# Integrate the system for 1000 years

# Use the SWIFT package for integration
k.setIntMethod(K_SWIFTRMVS)

# Integrate from the initial epoch to epoch + 1000 years, with 2000 samples
# Integrate returns a chain of length(samples) with a snapshot of the system parameters at each time
times = np.linspace(k.getEpoch(), k.getEpoch() + 1000.0 * 365.25, 2000)
(system, length) = k.integrate(times)

# Convert the chain for system states to a radial velocity
rv = sys_to_els(system, len(times), 0)

# Plot the eccentricity of each planet of the integration time
for i in range(1,k.nPlanets() + 1):
    plt.plot((times-k.getEpoch())/365.25,rv[:,i*K_ELEMENTS_SIZE+K_ECC+1],c=c[i-1])


# Plot apperance

plt.tick_params(labelsize=16)

plt.xlabel("Time after JD " + repr(round(k.getEpoch(),2)) + " $[Years]$", fontsize=20)
plt.ylabel("Eccentricity",fontsize=20)

plt.savefig("IntegrationPlot.pdf", bbox_inches='tight')

plt.clf()




# Generate a model RV curve and overplot the data points and error bars

k.setIntMethod(K_RK89)

# Get the time range of the RV data
start, finish = k.getRange()

# Add a bit to the start and finish so we don't lay right on the edge of a data point
start -= (finish-start) * 0.1
finish += (finish-start) * 0.1

# Get a combined array with all the RV data sorted by JD
data = matrix_to_array(getCompiledMatrix(k))

# Add the error bars for the data points
plt.errorbar(data[:,K_T_TIME],data[:,K_T_SVAL],yerr= data[:,K_T_ERR],fmt=None,ecolor='gray',elinewidth=2.,zorder=25,capthick=2)

# Plot the actual datapoints colored by vel file
plt.scatter(data[:,K_T_TIME],data[:,K_T_SVAL], c=colors, s=100.0,edgecolor='none',zorder=20,alpha=0.75)

# Calculate the predicted RV from time start to time finish with 1000 points
rvCurve = k.stellarRV(start,finish,1000)

# Plot the RV curve (Below the data points)
plt.plot(rvCurve[:,0],rvCurve[:,1],color='black',zorder=10)

# Plot appearance
plt.xlim([start,finish])

plt.xlabel('Julian Days [JD]', fontsize=20)
plt.ylabel('Radial Velocity $[m/s]$', fontsize=20)

plt.locator_params(nbins=5)
plt.tick_params(labelsize=18)

plt.savefig('RVCurve.pdf', bbox_inches='tight')


