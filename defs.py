from ctypes import *
import numpy as np
import os
import sys

from constants import *

lib = None

libpaths = [ '/Users/holden/Applications/mac/Systemic/libsystemic.dylib', '/Users/jaburt/Dropbox/Research/Systemic/libsystemic.dylib', '/Users/bradholden/Dropbox/APF_TESS/SystPy/libsystemic.dylib']
libpath = ''
for l in libpaths:
    if os.path.exists(l):
        libpath = l

try:
    lib = CDLL(libpath)
#    lib = CDLL('/Users/holden/Applications/mac/Systemic/libsystemic.dylib')    
except OSError as e:
    print "Couldn't locate libsystemic.dylib in the current directory. " , e



if lib == None:
    print "Couldn't load the systemic library."
    sys.exit(1)

def vector_to_array(v_ptr):
    """Takes a C pointer to a GSL_VECTOR and returns a copy of the vector as a Numpy array."""
    lib.ok_vector_len.argtypes = [c_long]
    length = lib.ok_vector_len(v_ptr)
    v = np.empty(length)

    lib.ok_vector_block.argtypes = [c_long]
	lib.ok_vector_block.restype = c_long
    lib.ok_block_to_ptr(lib.ok_vector_block(v_ptr), v.ctypes.data_as(POINTER(c_double)))
    return v

def matrix_to_array_alt(m_ptr):
    """Takes a C pointer to a GSL_MATRIX and returns a copy of that matrix as a Numpy array."""
    lib.ok_matrix_rows.argtypes = [c_long]
    lib.ok_matrix_cols.argtypes = [c_long]
    w,h = lib.ok_matrix_rows(m_ptr) , lib.ok_matrix_cols(m_ptr)
    length = w * h

    v = np.empty(length)
    lib.ok_matrix_block.argtypes = [c_long]
	lib.ok_matrix_block.restype = c_long
    lib.ok_block_to_ptr(lib.ok_matrix_block(m_ptr), v.ctypes.data_as(POINTER(c_double)))

    arr = v.reshape(w,h).copy()
    return arr

def matrix_to_array(m_ptr):
    """Takes a C pointer to a GSL_MATRIX and returns a copy of that matrix as a Numpy array."""
    lib.ok_matrix_rows.argtypes = [c_long]
    lib.ok_matrix_cols.argtypes = [c_long]
    w,h = lib.ok_matrix_rows(m_ptr) , lib.ok_matrix_cols(m_ptr)

    v = np.empty((w,h))
    for i in range(w):
        for j in range(h):
            v[i,j] = mGet(m_ptr, i, j)

    return v


def vGet(vector, i):
    """Returns the i-th element of a C GSL_VECTOR."""
    lib.gsl_vector_get.restype = c_double
    lib.gsl_vector_get.argtypes = [c_long, c_int]
    return lib.gsl_vector_get(vector, c_int(i))

def vSet(vector, i, val):
    """Sets the i-th element of a C GSL_VECTOR."""
    lib.gsl_vector_set.argtypes = [c_long, c_int, c_double]
    lib.gsl_vector_set(vector, c_int(i), c_double(val))

def mGet(matrix, i, j):
    """Returns the (i, j)-th element of a C GSL_MATRIX."""
    lib.gsl_matrix_get.restype = c_double
    lib.gsl_matrix_get.argtypes = [c_long, c_int, c_int]
    return lib.gsl_matrix_get(matrix, c_int(i), c_int(j))

def mSet(matrix, i, j, val):
    """Sets the (i, j)-th element of a C GSL_MATRIX."""
    lib.gsl_matrix_set.argtypes = [c_long, c_int,c_int, c_double]
    lib.gsl_matrix_set(matrix, c_int(i), c_int(j), c_double(val))

def array_to_matrix(arr):
    """Converts a Numpy array to a pointer to a C GSL_MATRIX."""
    r_num = len(arr)
    c_num = len(arr[0])
    lib.gsl_matrix_alloc.restype = c_long
    m = lib.gsl_matrix_alloc(c_int(r_num), c_int(c_num))
    for i in range(r_num):
        for j in range(c_num):
            mSet(m, i, j, arr[i,j])
    return m

def array_to_vector(arr):
    """Converts a 1-D Numpy array to a pointer to a C GSL_VECTOR."""
    if len(arr) != len(arr.flatten()):
        print "Array Conversion Warning: array is not 1-D"
    array = arr.flatten()
    lib.gsl_vector_alloc.restype = c_long
    v = lib.gsl_vector_alloc(c_int(len(array)))
    for i,val in enumerate(array):
        vSet(v, i, val)
    return v

class Kernel(Structure):
    """A Class wrapping the functionality of the C Kernel struct """

    # Python overrides
    def __init__(self, alloc=True):
        """Initialize an empty Kernel.

        Keyword Arg:
        alloc = True : Determines whether to allocate memory for a kernel
        """
        
        if alloc:
            lib.K_alloc.restype = c_long
            self.ptr = lib.K_alloc()
        else:
            self.ptr = None

    def __del__(self):
        lib.K_free.argtypes = [c_long]
        lib.K_free(self.ptr)

    # Kernel manipulation

    def copy(self):
        """Returns a new Kernel object that is a copy of the caller."""
        lib.K_clone.restype = c_long
        lib.K_clone.argtypes = [c_long]
        new_k = lib.K_clone(self.ptr)
        k2 = Kernel(alloc=False)
        k2.ptr = new_k
        return k2

    def calculate(self):
        """Updates all the values of the Kernel (Model values, Chi^2, RMS, etc.)
        Not intended to be required by user."""
        lib.K_calculate.argtypes = [c_long]
        lib.K_calculate(self.ptr)

    def addDataFile(self, location, directory=''):
        r"""Add a .sys or .vels file to the current Kernel
        
        Keyword Arg:
        directory = '' : Specify the directory of the file (ex. "datafiles/")
        """  
        lib.K_addDataFile.argtypes = [c_long, c_char_p, c_int]   
        if "sys" in location:
            vels = []
            trans = []
            with open(directory+location,'r') as f:
                for line in f:
                    if "RV[]" in line:
                        ls = line.split()
                        vels.append(ls[1][1:-1])
                    if "TD[]" in line:
                        ls = line.split()
                        trans.append(ls[1][1:-1])
                    if "Mass" in line:
                        ls = line.split()
                        sm = float(ls[1])
            for v in vels:
                try:
                    f = open(directory + v, 'r')
                except IOError:
                    print "Warning! File %s could not be found. It will not be added to the Kernel." % (directory + v)
                else:
                    f.close()
                    lib.K_addDataFile(self.ptr, c_char_p(directory + v), c_int(K_T_RV))
            for t in trans:
                try:
                    f = open(directory + v, 'r')
                except IOError:
                    print "Warning! File %s could not be found. It will not be added to the Kernel." % (directory + v)
                else:
                    f.close()
                    lib.K_addDataFile(self.ptr, c_char_p(directory + t), c_int(K_T_TIMING))
            self.setMstar(sm)
        else:
            try:
                f = open(directory + location,'r')
            except IOError:
                print "Warning! File %s could not be found. It will not be added to the Kernel." % (directory + location)
            else:
                f.close()
                if '.tds' in location:
                    lib.K_addDataFile(self.ptr, c_char_p(directory+location), c_int(K_T_TIMING))
                else:
                    lib.K_addDataFile(self.ptr, c_char_p(directory+location), c_int(K_T_RV))
        self.calculate()

    def addDataArray(self, arr):
        """Add a Numpy array of the form [ [Time, RV, Error] ... ] to the Kernel."""
        r_num = len(arr)
        lib.gsl_matrix_alloc.restype = c_long
        m = lib.gsl_matrix_alloc(c_int(r_num), c_int(3))
        for i in range(r_num):
            mSet(m, i, 0, arr[i,0])
            mSet(m, i, 1, arr[i,1])
            mSet(m, i, 2, arr[i,2])
        lib.K_addDataTable.argtypes = [c_long, c_long, c_char_p, c_int]
        lib.K_addDataTable(self.ptr, m,c_char_p("Radial Set"), K_T_RV)
        self.calculate()

    def removeData(self, i=-1):
        """Remove all datasets from the Kernel

        Keyword Arg
        i = -1 : Remove the ith dataset rather then every dataset
        """
        lib.K_removeData.argtypes = [c_long, c_int]
        lib.K_removeData(self.ptr, c_int(i))
        self.calculate()

    def addPlanet(self, params):
        """Add a planet to the system.
        Input is a list of values in the form [Param, Value, ... , K_DONE]
        For example this will add a Jupiter type planet
        k.addPlanet([K_PER, 4000., K_MASS, 1.0, K_ECC, 0.01, K_DONE])
        """
        if len(params) % 2 != 1:
            print "Incorrect argument length. Did you forget K_DONE?"
            print "No planet was added."
            return
        p = np.array(params).astype(float)
        lib.K_addPlanet.argtypes = [c_long, POINTER(c_double)]
        lib.K_addPlanet(self.ptr, p.ctypes.data_as(POINTER(c_double)))
        self.calculate()

    def removePlanet(self, i=-1):
        """Removes all planets from the system.

        Keyword Arg:
        i = -1 : remove the i-th planet rather then all planets.
        """
        lib.K_removePlanet.argtypes = [c_long, c_int]
        lib.K_removePlanet(self.ptr, c_int(i))
        self.calculate()

    def setElement(self, i, j, val):
        """Sets the i-th planet's orbital parameter j to the specified value."""
        lib.K_setElement.argtypes = [c_long, c_int, c_int, c_double]
        lib.K_setElement(self.ptr, c_int(i), c_int(j), c_double(val))
        self.calculate()

    def setElementFlag(self, i, j, val):
        """Sets the i-th planet's orbital element j's minimization flag.
        Use K_MINIMIZE | K_ACTIVE to turn minimization on
        or  K_MINIMIZE & K_ACTIVE to turn minimization off.
        """
        lib.K_setElementFlag.argtypes = [c_long, c_int, c_int, c_int]
        lib.K_setElementFlag(self.ptr, c_int(i), c_int(j), c_int(val))

    def setElementType(self, i):
        """Sets the coordinate system.
        Options are K_ASTROCENTRIC or K_JACOBI.
        """
        lib.K_setElementType.argtypes = [c_long, c_int]
        lib.K_setElementType(self.ptr, c_int(i))
        self.calculate()

    def setMstar(self, val):
        """Sets the stellar mass (in solar mass units) to use in the Kernel."""
        lib.K_setMstar.argtypes = [c_long, c_double]
        lib.K_setMstar(self.ptr, c_double(val))
        self.calculate()

    def setEpoch(self, val):
        """Sets the initial epoch to use in the Kernel."""
        lib.K_setEpoch.argtypes = [c_long, c_double]
        lib.K_setEpoch(self.ptr, c_double(val))
        self.calculate()

    def setIntMethod(self, i):
        """Sets the integration method to use.
        Options are one of:
        K_KEPLER
        K_RK45
        K_RK89
        K_ADAMS
        K_BULIRSCHSTOER
        K_SWIFTRMVS (Used for integrating only)
        """
        lib.K_setIntMethod.argtypes = [c_long, c_int]
        lib.K_setIntMethod(self.ptr, c_int(i))
        self.calculate()

    def setPar(self, i, val):
        """Sets the i-th parameter to val. Usually a data offset."""
        lib.K_setPar.argtypes = [c_long, c_int, c_double]
        lib.K_setPar(self.ptr, c_int(i), c_double(val))
        self.calculate()

    def setParFlag(self, i, idx):
        """Sets the i-th parameter minimization flag.
        Use K_MINIMIZE | K_ACTIVE to turn minimization on
        or  K_MINIMIZE & K_ACTIVE to turn minimization off.
        """
        lib.K_setParFlag.argtypes = [c_long, c_int, c_int]
        lib.K_setParFlag(self.ptr, c_int(i), c_int(idx))
        self.calculate()

    def setTrend(self, val):
        self.setPar(20, val)

    def setTrendFlag(self, flag):
        self.setParFlag(20, flag)
    
    def minimize(self, algo=K_SIMPLEX, maxiter=5000, params=[]):
        """Minimize the Kernel using previously set element and par flags to minimize parameters.

        Keyword Args:
        algo = K_SIMPLEX : minimization algorithm (Options K_LM, K_DIFFEVOL, K_SA)
        maxiter = 5000 : maximum iterations before stopping
        params = [] : additional parameters to minimize
        """
        
        if params == []:
            lib.K_minimize.argtypes = [c_long, c_int, c_int, c_void_p]
            lib.K_minimize(self.ptr, c_int(algo), c_int(maxiter), None)
        else:
            lib.K_minimize.argtypes = [c_long, c_int, c_int, POINTER(c_double)]
            p = np.array(params)
            lib.K_minimize(self.ptr, c_int(algo), c_int(maxiter), p.ctypes.data_as(POINTER(c_double)))

    def minimize1d(self, row, col, algo=K_SIMPLEX, maxiter=5000, params=[]):
        """Minimizes planet i's orbital element j only.
            To minimize a Velocity offset use -1, K_P_DATA(# of set)

        Keyword Args:
        algo = K_SIMPLEX : minimization algorithm (Options K_LM, K_DIFFEVOL, K_SA)
        maxiter = 5000 : maximum iterations before stopping
        params = [] : additional parameters to minimize
        """
        
        if params == []:
            lib.K_1dminimize.argtypes = [c_long, c_int, c_int, c_int, c_int,c_void_p]
            lib.K_1dminimize(self.ptr, c_int(algo), c_int(maxiter), c_int(row), c_int(col), None)
        else:
            lib.K_1dminimize.argtypes = [c_long, c_int, c_int, c_int, c_int,POINTER(c_double)]
            p = np.array(params)
            lib.K_1dminimize(self.ptr, c_int(algo), c_int(maxiter), c_int(row), c_int(col), p.ctypes.data_as(POINTER(c_double)))

    def clone(self):
        """Returns a copy of this Kernel object."""
        knew = Kernel(alloc=False)
        lib.K_clone.restype = c_long
        lib.K_clone.argtypes = [c_long]
        knew.ptr = lib.K_clone(self.ptr)
        return knew
        


    # Get Kernel Information
    def getElement(self, i, j):
        """Returns the i-th planet's orbital element j."""
        lib.K_getElement.restype = c_double
        lib.K_getElement.argtypes = [c_long, c_int, c_int]
        return lib.K_getElement(self.ptr, c_int(i), c_int(j))

    def getElementFlag(self, i, j):
        """Returns the i-th planets orbital element j's flag."""
        lib.K_getElementFlag.argtypes = [c_long, c_int, c_int]
        return lib.K_getElementFlag(self.ptr, c_int(i), c_int(j))

    def getElements(self):
        """ Returns an array of system elements for this Kernel."""
        lib.K_getElements.argtypes = [c_long]
        lib.K_getElements.restype = c_long
        return matrix_to_array(lib.K_getElements(self.ptr))

    def getAllElements(self):
        """ Returns a more detailed array of elements than a simple getElements() call. Also includes Semi-major axis, RV half amplitude."""
        lib.K_getAllElements.argtypes = [c_long]
        lib.K_getAllElements.restype = c_long
        return matrix_to_array(lib.K_getAllElements(self.ptr))

    def getActivePars(self):
        lib.K_getActivePars.argtypes = [c_long]
        return lib.K_getActivePars(self.ptr)

    def getActiveElements(self):
        lib.K_getActiveElements.argtypes = [c_long]
        return lib.K_getActiveElements(self.ptr)

    def getTrend(self):
        return self.getPar(20)

    def getTrendFlag(self):
        return self.getParFlag(20)

    def getNrPars(self):
        lib.K_getNrPars.argtypes = [c_long]
        return lib.K_getNrPars(self.ptr)

    def getElementType(self):
        """ Returns the current internal element representation. Either K_JACOBI, or K_ASTROCENTRIC."""
        lib.K_getElementType.argtypes = [c_long]
        return lib.K_getElementType(self.ptr)

    def getXYZ(self):
        lib.K_getXYZ.argtypes = [c_long]
        return matrix_to_array(lib.K_getXYZ(self.ptr))

    def getMstar(self):
        """ Returns the current mass of the star for this Kernel."""
        lib.K_getMstar.restype = c_double
        lib.K_getMstar.argtypes = [c_long]
        return lib.K_getMstar(self.ptr)

    def getEpoch(self):
        """ Returns the current initial epoch of the system. """
        lib.K_getEpoch.restype = c_double
        lib.K_getEpoch.argtypes = [c_long]
        return lib.K_getEpoch(self.ptr)

    def getChi2(self):
        """ Returns the current Chi^2 for the current model fit to the data."""
        lib.K_getChi2.restype = c_double
        lib.K_getChi2.argtypes = [c_long]
        return lib.K_getChi2(self.ptr)

    def Chi2(self):
        """ Returns the current Chi^2 for the current model fit to the data."""
        self.getChi2()

    def getChi2_nr(self):
        lib.K_getChi2_nr.restype = c_double
        lib.K_getChi2_nr.argtypes = [c_long]
        return lib.K_getChi2_nr(self.ptr)
    
    def getChi2_ttv(self):
        """Return the chi squared for transit timing variations."""
        lib.K_getChi2_tts.restype = c_double
        lib.K_getChi2_tts.argtypes = [c_long]
        return lib.K_getChi2_tts(self.ptr)

    def getTTVPeriod(self):
        times = np.empty(0)
        for i in range(self.getNsets()):
            data = matrix_to_array(self.getData(i))
            if data[0,K_T_FLAG] == 2.:
                times = np.concatenate((times, data[:,K_T_VAL]))

        diff = times[1:] - times[:-1]
        if len(diff) == 0:
            print "No transit data found."
            return None
        return min(diff)

    def getLoglik(self):
        lib.K_getLoglik.restype = c_double
        lib.K_getLoglik.argtypes = [c_long]
        return lib.K_getLoglik(self.ptr)

    def getRMS(self):
        """ Returns the current RMS of the model fit to the data."""
        lib.K_getRms.restype = c_double
        lib.K_getRms.argtypes = [c_long]
        return lib.K_getRms(self.ptr)

    def RMS(self):
        """ Returns the current RMS of the model fit to the data."""
        return self.getRMS()

    def getRms(self):
        """ Returns the current RMS of the model fit to the data."""
        return self.getRMS()

    def getJitter(self):
        """ Returns the current jitter for the model fit to the data."""
        lib.K_getJitter.restype = c_double
        lib.K_getJitter.argtypes = [c_long]
        return lib.K_getJitter(self.ptr)

    def getNrvs(self):
        """ Returns the number of radial velocity data points associated with the current Kernel."""
        lib.K_getNrvs.argtypes = [c_long]
        return lib.K_getNrvs(self.ptr)

    def getNtts(self):
        """ Returns the number of transit timing data points associated with the current Kernel."""
        lib.K_getNtts.argtypes = [c_long]
        return lib.K_getNtts(self.ptr)

    def nPlanets(self):
        """ Returns the number of planets in the Kernel."""
        lib.K_getNplanets.argtypes = [c_long]
        return lib.K_getNplanets(self.ptr)

    def getNsets(self):
        """ Returns the current number of datasets in the Kernel."""
        lib.K_getNsets.argtypes = [c_long]
        return lib.K_getNsets(self.ptr)

    def nData(self):
        """Return the number of data points in the kernel."""
        lib.K_getNdata.argtypes = [c_long]
        return lib.K_getNdata(self.ptr)

    def getNdata(self):
        """Return the number of data points in the kernel."""
        return self.nData()

    def getIntMethod(self):
        """Return the current Kernel Integration Method."""
        lib.K_getIntMethod.argtypes = [c_long]
        return lib.K_getIntMethod(self.ptr)

    def getPars(self):
        """Return an array of current Parameters."""
        lib.K_getPars.argtypes = [c_long]
        return vector_to_array(lib.K_getPars(self.ptr))

    def getPar(self, i):
        """Return the parameter value of par i."""
        lib.K_getPar.restype = c_double
        lib.K_getPar.argtypes = [c_long, c_int]
        return lib.K_getPar(self.ptr, c_int(i))

    def getParFlag(self, i):
        """Return the minimization flag for par i."""
        lib.K_getParFlag.argtypes = [c_long, c_int]
        return lib.K_getParFlag(self.ptr, c_int(i))

    def getCompiled(self):
        """Return a pointer to a GSLMatrix holding combined datasets."""
        lib.K_getCompiled.argtypes = [c_long]
        lib.K_getCompiled.restype = c_long
        return lib.K_getCompiled(self.ptr)

    def getData(self, i):
        """Return a GSL Matrix pointer of dataset i."""
        lib.K_getData.argtypes = [c_long, c_int]
        lib.K_getData.restype = c_long   # <-- that's the line I'd try adding
        return lib.K_getData(self.ptr, c_int(i))

    def getRange(self):
        """Return the range over which the RV data spans as (start,finish)."""
        start = c_double()
        finish = c_double()
        lib.K_getRange.argtypes = [c_long, POINTER(c_double), POINTER(c_double)]
        lib.K_getRange(self.ptr, byref(start), byref(finish))

        return start.value, finish.value

    def getElementType(self):
        """Return the current coordiante system of the Kernel."""
        lib.K_getElementType.argtypes = [c_long]
        return lib.K_getElementType(self.ptr)

    def stellarRV(self, start=0, finish=0, samples=2000):
        """Calculate the predicted stellar RV for the current planetary system of the Kernel. Returns an array with times and RVs.""" 
        if start == 0 and finish==0:
            start = self.getEpoch()
            finish = start + 4*365.25
        lib.K_integrateStellarVelocity.argtypes = [c_long, c_double, c_double, c_int, c_void_p, c_void_p]
        lib.K_integrateStellarVelocity.restype = c_long
        return matrix_to_array(lib.K_integrateStellarVelocity(self.ptr, c_double(start), c_double(finish), c_int(samples), None, None))

    def periodogram(self, ret_type, algo=K_SIMPLEX, circ=False, sample=0, samples=20000, Pmin=0.5, Pmax=10000.0):
        """Calculate a periodogram of the data in the kernel."""
        matrix_to_array(lib.ok_periodogram_full(self.ptr, c_int(ret_type), c_int(algo), c_int(circ), c_int(sample), c_int(samples), c_double(Pmin), c_double(Pmax)))
        
    def integrate(self, times, intMethod=K_SWIFTRMVS):
        """Integrate the system in the Kernel to each time in times, returning system snapshots for each time. Returns (systems, length)."""
        old_int = self.getIntMethod()
        self.setIntMethod(intMethod)
        self.calculate()
        lib.gsl_vector_alloc.restype = c_long
        v = lib.gsl_vector_alloc(c_int(len(times)))
        for i,t in enumerate(times):
            vSet(v,i,t)

        lib.K_integrate.argtypes = [c_long, c_long, c_void_p, c_void_p]
        system = lib.K_integrate(self.ptr, v, None, None)

        self.setIntMethod(old_int)
        self.calculate()

        return (system,len(times))

    def integrateRange(self, start, finish, samples, intMethod=K_SWIFTRMVS):
        """Integrates the system in the Kernel from start to finish with samples. Returns tuple (pointer to ok_system**, length) """
        old_int = self.getIntMethod()
        self.setIntMethod(intMethod)
        self.calculate()
        args = [c_long, c_double, c_double, c_int, c_void_p, c_void_p]
        lib.K_integrateRange.argtypes = args
        lib.K_integrateRange.restype = c_long
        system = lib.K_integrateRange(self.ptr, c_double(start), c_double(finish), c_int(samples), None, None)
        self.setIntMethod(old_int)
        self.calculate()

        return (system,samples)

    def save(self, filename):
        """Save a snapshot of the current Kernel in filename."""
        lib.fopen.restype = c_long
        fp = lib.fopen(c_char_p(filename), c_char_p("w"))
        lib.K_save.argtypes = [c_long, c_long]
        lib.K_save(self.ptr, fp)



def loadKernel(filename, skip=0):
    """Load a Kernel snapshot in filename."""
    k = Kernel(alloc=False)
    lib.fopen.restype = c_long
    lib.K_load.restype = c_long
    fp = lib.fopen(c_char_p(filename),c_char_p("r"))
    lib.K_load.argtypes = [c_long, c_int]
    k.ptr = lib.K_load(fp, c_int(skip))
    k.calculate()
    return k


# Kernel related functions not suited to be methods

def getCompiledMatrix(k):
    lib.ok_buf_to_matrix.restype = c_long
    lib.ok_buf_to_matrix.argtypes = [c_long, c_int, c_int]
    return lib.ok_buf_to_matrix(k.getCompiled(), k.getNdata(), K_DATA_SIZE)

def getResidualMatrix(k):
    #return lib.K_getResidualMatrix(k.ptr)
    lib.ok_buf_to_matrix.restype = c_long
    lib.ok_buf_to_matrix.argtypes = [c_long, c_int, c_int]
    
    arr = matrix_to_array(lib.ok_buf_to_matrix(k.getCompiled(), k.getNdata(), K_DATA_SIZE))
    arr[:,K_T_SVAL] = arr[:,K_T_SVAL] - arr[:,K_T_PRED]
    arr[:,K_T_VAL] = arr[:,K_T_VAL] - arr[:,K_T_PRED]
    return array_to_matrix(arr)


def periodogram_ls(data, samples, Pmin, Pmax, method, timecol, valcol, sigcol):
    lib.ok_periodogram_ls.restype = c_long
    args = [c_long, c_int, c_double, c_double, c_int, c_int, c_int, c_int, c_void_p]
    lib.ok_periodogram_ls.argtypes = args
    return matrix_to_array(lib.ok_periodogram_ls(data, c_int(samples), c_double(Pmin), c_double(Pmax), c_int(method), c_int(timecol), c_int(valcol), c_int(sigcol), None))

def periodogram_boot(data, trials, samples, Pmin, Pmax, method, timecol, valcol, sigcol, seed=lib.time(None)):
    lib.ok_periodogram_boot.restype = c_long
    args = [c_long, c_int,c_int, c_double, c_double, c_int, c_int, c_int, c_int, c_int, c_void_p, c_void_p]
    lib.ok_periodogram_boot.argtypes = args
    mat = lib.ok_periodogram_boot(data, c_int(trials), c_int(samples), c_double(Pmin), c_double(Pmax), c_int(method), c_int(timecol), c_int(valcol), c_int(sigcol), c_int(seed), None, None)
    return matrix_to_array(mat)

def sys_to_rvs(system, length):
    lib.ok_get_rvs.restype = c_long
    lib.ok_get_rvs.argtypes = [c_long, c_int]
    return matrix_to_array(lib.ok_get_rvs(system, c_int(length)))

def sys_to_xyzs(system, length):
    lib.ok_get_xyzs.restype = c_long
    lib.ok_get_xyzs.argtypes = [c_long, c_int]
    return matrix_to_array(lib.ok_get_xyzs(system, c_int(length)))

def sys_to_els(system, length, internal):
    lib.ok_get_els.restype = c_long
    lib.ok_get_els.argtypes = [c_long, c_int, c_int]
    return matrix_to_array(lib.ok_get_els(system, c_int(length), c_int(internal)))


def findTransits(system, length, planetIDX, intMethod, eps, flags):
    """ Find transits nearest the times of the system snapshots for planet.
        Args:
        system (pointer) - pointer to a system struct as returned by k.integrate()[0]
        length (int) - Number of snapshots in system as returned by k.integrate()[1]
        planet (int) - Planet index 1-N
        intMeth - Integration method to use [One of the built-ins, i.e. K_RK89, K_SWIFTRMVS, etc. ]
        eps (float) - Desired precision to find transits with
    """
    lib.ok_find_transits.restype = c_long
    lib.ok_find_transits.argtypes = [c_long, c_int, c_int, c_int, c_double, c_void_p, c_void_p]
    v = lib.ok_find_transits(system, length, planetIDX, intMethod, eps, None, None)
    return vector_to_array(v)

def getPeakIndx(array):
    a = np.array(array)
    x = a[1:] - a[:-1]
    y = a[:-1] - a[1:]
    x = np.insert(x, 0, 0)
    y = np.insert(y, -1, 0)

    peaks = (x > 0) & (y > 0)

    return np.arange(0,len(a))[peaks]

def fitSystem(k, fapLim=0.001, minPer=1.5):
    if k.getNrvs() < 20:
        print "Too few data points"
        return
    #print "Beginning Fit"

    k.removePlanet(-1)
    k.calculate()

    for i in range(k.getNsets()):
        k.setParFlag(K_P_DATA1+i, K_MINIMIZE | K_ACTIVE)

    k.setTrendFlag(K_ACTIVE | K_MINIMIZE)

    while True:
        topP = -1
        k.calculate()
        p = periodogram_ls(getResidualMatrix(k), 100000, minPer, 10000, 0, K_T_TIME, K_T_SVAL, K_T_ERR)
        amp = p[:,K_PS_Z]
        peaks = getPeakIndx(amp)
        amp = amp[peaks]
        per = p[:,K_PS_TIME][peaks]
        fap = p[:,K_PS_FAP][peaks]
        sort = amp.argsort()[::-1]
        
        for idx in sort[:10]:
            if per[idx] > minPer and fap[idx] < fapLim:
                topP = idx
                break
        if topP == -1:
            break
        # Check if the planet we found has a similar period to an existing planet
        # This can suggest an interacting system, or perhaps low data points
        # We don't want the fitter stuck in a loop of adding numerous planets
        # with the same period
        pl_per = k.getElements()[1:,K_PER]
        cont = False
        for tp in pl_per:
            if np.abs(tp - per[topP]) / per[topP]  < 0.1:
                print tp, per[topP], "Found period similar to existing planet."
                print "Simple Periodogram analysis will not be sufficent"
                break
        else:
            cont = True
        if not cont: break
        #print "Found Planet with Period =", per[topP]
        k.addPlanet([K_PER, per[topP], K_DONE])
        lastPl = k.nPlanets()
        for i in range(2):
            k.minimize1d(lastPl,K_MA)
            k.minimize1d(lastPl,K_MASS)
            k.minimize1d(lastPl,K_PER)
            for j in range(K_ELEMENTS_SIZE):
                k.setElementFlag(lastPl,j,R_STEADY)
            k.setElementFlag(lastPl,K_PER,R_MINIMIZE)
            k.setElementFlag(lastPl,K_MASS,R_MINIMIZE)
            k.setElementFlag(lastPl,K_MA,R_MINIMIZE)
        k.calculate()
        k.minimize()

        k.setElementFlag(lastPl, K_ECC, R_MINIMIZE)
        k.setElementFlag(lastPl, K_LOP, R_MINIMIZE)

        k.minimize()

    k.minimize()
    k.calculate()
    #print "Finished Fitting System"

class KList():
    def __init__(self,fn, alloc=True):
        if alloc:
            lib.fopen.restype = c_long
            lib.KL_load.restype = c_long
            lib.KL_load.argtypes = [ c_long, c_int]
            fp = lib.fopen(c_char_p(fn), c_char_p("r"))
            self.ptr = lib.KL_load(fp, c_int(0))
        else:
            self.ptr = None

    def __del__(self):
        lib.KL_free.argtypes = [c_long]
        lib.KL_free(self.ptr)

    def getSize(self):
        lib.KL_getSize.argtypes = [c_long]
        return lib.KL_getSize(self.ptr)

    def getNplanets(self):
        lib.KL_getNplanets.argtypes = [c_long]
        return lib.KL_getNplanets(self.ptr)

    def getElements(self, i, j):
        lib.KL_getElements.restype = c_long
        lib.KL_getElements.argtypes = [c_long, c_int, c_int]
        return vector_to_array(lib.KL_getElements(self.ptr, c_int(i), c_int(j)))

    def getPars(self, i):
        lib.KL_getPars.restype = c_long
        lib.KL_getPars.argtypes = [c_long, c_int]
        return vector_to_array(lib.KL_getPars(self.ptr, c_int(i)))

    # Returns matrix of the statistics for each orbital element
    # the input can be one of (likely with K in front):
    # STAT_MEAN, STAT_MEDIAN, STAT_STDDEV, STAT_MAD
    def getElementsStats(self, i):
        lib.KL_getElementsStats.restype = c_long
        lib.KL_getElementsStats.argtypes = [c_long, c_int]
        return matrix_to_array(lib.KL_getElementsStats(self.ptr, c_int(i)))

    # Same as above but for parameters rather than orbital elements
    def getParsStats(self, i):
        lib.KL_getParsStats.restype = c_long
        lib.KL_getParsStats.argtypes = [c_long, c_int]
        return vector_to_array(lib.KL_getParsStats(self.ptr, c_int(i)))

def KLSave(kl, filename):
    lib.fopen.restype = c_long
    lib.KL_save.argtypes = [c_long, c_long]
    lib.fclose.argtypes = [c_long]
    fp = lib.fopen(c_char_p(filename), c_char_p("w"))
    lib.KL_save(kl.ptr, fp)
    lib.fclose(fp)
    print "Successfully saved Kernel List"


def MCMC_New(k, nChains=2, nTemps=1, skip=1000, discard=200, rStop=1.1):
    lib.MCMC_New.restype = c_long
    args = [c_long, c_int, c_int, c_int, c_int, c_double]
    lib.MCMC_New.argtypes = args
    return lib.MCMC_New(k.ptr, c_int(nChains), c_int(nTemps), c_int(skip), c_int(discard), c_double(rStop))

# MCMC_mult( [k, k] )
# MCMC_single( k )

def MCMC(k, nChains = 2, rStop = 1.1):
    ks = [ k for n in range(0,nChains)]
    return MCMC_mult(ks,nChains=nChains, rStop = rStop)


def MCMC_mult(ks, nChains = 2, nTemps = 1, skip = 1000, discard = 200, rStop = 1.1, params = None):
    kl = KList('', alloc=False)
    lib.K_mcmc_mult.restype = c_long
    lib.K_mcmc_mult.argtypes = [POINTER(c_long), c_int, c_int, c_int, c_int, POINTER(c_double), c_double, c_void_p]
     
    if (len(ks) < nChains):
        print("Fewer kernels than chains, I think this breaks the inner C loop")
        return
 
    ks_ptr = np.array([k.ptr for k in ks], dtype=np.long)
 
    kl_ptr = lib.K_mcmc_mult(ks_ptr.ctypes.data_as(POINTER(c_long)), c_int(nChains), c_int(nTemps), c_int(skip), c_int(discard), None, c_double(rStop), None)

    kl.ptr = kl_ptr
 
    return kl
 
def MCMC_single(k, nSteps=2, skip=1000, discard=200, params=None, cont_list=None, merit=None, tag=0, flag=0):
    kl = KList('', alloc=False)
    lib.K_mcmc_single.restype = c_long
    lib.K_mcmc_single.argtypes = [c_long, c_int, c_int, c_int, POINTER(c_double), c_void_p, c_void_p, c_int, POINTER(c_int)]
     
    kl_ptr = lib.K_mcmc_single(k.ptr, c_int(nSteps), c_int(skip), c_int(discard), None, None, None, c_int(tag), byref(c_int(flag)))
    kl.ptr = kl_ptr
    return kl

def bootstrap(k, trials, warmup, malgo=K_SIMPLEX, miter=5000):
        """Bootstrap error analysis on the kernel K.
        Args:
        trials (int) - Number of bootstrap trials to preform
        warmup (int) - Initial states to discard
        Optional Args:
        malgo=K_SIMPLEX (int) - Minimizer algorithm
        miter=5000      (int) - Number of iterations for minimize function
        """
        kl = KList('', alloc=False)
        lib.K_bootstrap.restype = c_long
        lib.K_bootstrap.argtypes = [c_long, c_int, c_int, c_int, c_int, c_void_p]
        kl_ptr = lib.K_bootstrap(k.ptr, c_int(trials), c_int(warmup), c_int(malgo), c_int(miter), None)
        kl.ptr = kl_ptr
        return kl

def findClosestTransit(state, planet, intMeth, eps, transitType):
    """Find the closest transit to the state of the system passed to this function.
        Args:
        state (pointer) - single snapshot of system elements
        planet (int) - index of planet to find transits for
        intMeth (int) - One of built in integration methods [e.g. K_RK89, K_SWFTRMVS]
        eps (double) - Desired precision to find transits to.
        transitType (int) - Primary transit [K_TDS_PRIMARY] or Secondary tranist [K_TDS_SECONDARY]
    """
    timeout = float()
    error = float()
    lib.ok_find_closest_tranist.restype = c_long
    lib.ok_find_closest_transit.argtypes = [c_long, c_int, c_void_p, c_int, c_double, POINTER(c_double), POINTER(c_int)]
    return lib.ok_find_closest_transit(state, planet, None, intMeth, eps, transitType, timeout, error)


    
# ok_list* K_bootstrap(ok_kernel* k, int trials, int warmup, int malgo, int miter, double mparams[])


# int ok_find_closest_transit(ok_system* sys, const int pidx, ok_integrator_options* options, const int intMethod, const double eps, const int type, double* timeout, int* error)
# double ok_pcalc(double Mcenter, double Mp, double a)
# double ok_acalc(double Mcenter, double Mp, double P)
# ok_list* KL_alloc(const int size, ok_kernel* prototype)
# void KL_free(ok_list* list)
# ok_list* KL_load(FILE* fid, int skip)
# void KL_save(const ok_list* kl, FILE* out)
# void KL_append(ok_list* dest, ok_list* src)
# gsl_vector* KL_getParsStats(const ok_list* kl, const int what)
# gsl_vector* KL_getElements(const ok_list* kl, const int pl, const int el)
# gsl_vector* KL_getPars(const ok_list* kl, const int vo)
# gsl_matrix* KL_getElementsStats(const ok_list* kl, const int what)
# ok_list_item* KL_set(ok_list* kl, const int idx, gsl_matrix* elements, gsl_vector* pars, double merit, int tag)
# int KL_getSize(const ok_list* kl)
# void KL_removeAtIndex(ok_list* kl, const int idx)
# void KL_fprintf(const ok_list* kl, FILE* out, const char* fmt, const char* lfmt)
# void KL_to_ptr(const ok_list* kl, double* out)
# int KL_getNplanets(const ok_list* kl)
# ok_list* K_mcmc_single(ok_kernel* k, unsigned int nsteps, unsigned int skip, unsigned int discard, const double dparams[], ok_list* cont, ok_callback2 merit_function, int tag, int* flag)
# ok_list* K_mcmc_mult(ok_kernel** k, unsigned int nchains, unsigned int ntemps, unsigned int skip, unsigned int discard, const double params[], double Rstop, ok_callback2 merit_function)
# ok_list* K_bootstrap(ok_kernel* k, int trials, int warmup, int malgo, int miter, double mparams[])
# int K_isMstable_coplanar(const gsl_matrix* alle)
# gsl_matrix* gsl_matrix_alloc(size_t n1, size_t n2)
# void gsl_matrix_free(gsl_matrix* m)
# void gsl_vector_free(gsl_vector* v)




    
