import sys
print(sys.version)
print(sys.executable)

# import matplotlib.pyplot as plt

import math as m

# enforce a partition of ONE POSSIBLE region R
# x_values = np.linspace(1,10,10)
# y_values = np.linspace(1,15,10)

# subdivide each of the x and y values into 
# further partitions

# new_array_1 = np.zeros([10,10])

# for X in range(len(x_values)):
    # create partitions of x coordinates
    # of each box
    
# new_array_2 = np.zeros([15,10])

# for Y in range(len(y_values)):
    # create partitions of y coordinates
    # of each box

from math import sqrt
from scipy.stats import norm
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
# import array as arr

# from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
# from pylab import plot, show, grid, axis, xlabel, ylabel, title

def brownian(x0, n, dt, delta, out=None):

    x0 = np.asarray(x0)

    # For each element of x0, generate a sample of n numbers from a
    # normal distribution.
    r = norm.rvs(size=x0.shape + (n,), scale=delta*sqrt(dt))

    # If `out` was not given, create an output array.
    if out is None:
        out = np.empty(r.shape)

    # This computes the Brownian motion by forming the cumulative sum of
    # the random samples. 
    np.cumsum(r, axis=-1, out=out)
   
    # Add the initial condition.
    out += np.expand_dims(x0, axis=-1)

    return out

# The Wiener process parameter.
# delta = 0.25

# Total time.
T = 10.0
# Number of steps.
N = 50
# Time step size
# dt = T/N

Count_array_2 = np.linspace(1,10,100)

Sample_array = { }
    

# Initial values of x and array to 
# store Brownian Bridge samples

def prep_brownian(Count_array,N,T,delta):

    for X in Count_array:

        x = np.empty((2,N+1))
        # declare initial value of the
        # Brownian sampling process.

        # x0 = np.asarray(Count_array[X])
        x[:, 0] = Count_array[X-1]
        
        Z = brownian(x[:,0], N, T/N,delta, out=x[:,1:])
        # print(Z)
        Sample_array[X] = Z
    return Sample_array

Array = prep_brownian(Count_array_2, N,10,0.01)

# take absolute value of sample displacement

def plot_sample_distribution(testArray):

    test_array = { }
    # concat all samples together, also looking to do this with append
    Sample_cat = np.concatenate([testArray[1][0], testArray[2][0], testArray[3][0], testArray[4][0], testArray[5][0], testArray[6][0], testArray[7][0], testArray[8][0], testArray[9][0], testArray[10][0]])
    for I in np.linspace(1,len(Sample_cat), len(Sample_cat)):
        I = I+1
        if I<(len(Sample_cat)):
            test_array[I] = Sample_cat[I] - Sample_cat[I-1]
        else:
            test_array[I] = 0

    # after for loop, get values from dict 
    x = test_array.values()
    # convert list to numpy array
    x_empty = np.zeros(len(x))
    for J in range(len(x)):
        x_empty[J] = x[J]
    
    x_abs = abs(x_empty)

    # next, proceed to plot the distributions
    # of Brownian samples in the signed
    # and unsigned cases

    plt.plot(x_abs, np.linspace(1,x_abs,x_abs))
    plt.plot(x_empty, np.linspace(1,x_empty,x_empty))

    return test_array


Test_Arr = plot_sample_distribution(Array)

def compute_couplings(Test_arry):
    # define the revised couplings
    # as from the post PDF
    Test_arry = Test_arry.values()
    # maximum number of nonzero entries in coupling short range 
    # is 49, corresponding to 50-1 entries using 1 as a reference
    couplings_short = np.empty([(len(Count_array_2)/10) * N, N])
    # on the other hand, there are 
    couplings_long = np.empty([(len(Count_array_2)/10) * 450,450])

    # sep_array = np.empty([10,len(Test_arry)/10])

    temp = np.zeros([ 1 , 50 ])

    # times at which each sample is collected
    time_array_short= np.zeros([N,1])
    time_array_long = np.zeros([((len(Count_array_2)/10)-1)* N, 1])

    for KI in range(len(time_array_short)):
        time_array_short[KI]= 0 + ((T/N) * KI)

    for KJ in range(len(time_array_long)):
        time_array_long[KJ] = 0 + ((T/N) * KJ)

    short_time_sum = sum(time_array_short)
    long_time_sum = sum(time_array_long)

    for I in range(int(len(Count_array_2)/10)):

        # obtain the bridge samples for each increment
        # of 50 steps taken per box (will automate this
        # to allow for arbitrary steps in the future)
        
        temp = Test_arry[(50*I)+1:50*(I+1)]
        # obtain remaining bridge samples rather
        # than those in temp vector above

        tmp_array = Test_arry

        for AB in range(len(np.linspace((50*I)+1 , 50*(I+1), (50*(I+1)- ((50*I)+1)-1)))):
            # modify Test_arry to place a 0 in each
            # entry that is captured in temp array given above
            tmp_array[AB] = tmp_array[AB] * 0

        # obtain indices of nonzero elements in Test_arry_tmp
        
        tmp_array = np.asarray(tmp_array)

        tmp_array = tmp_array(np.nonzero(tmp_array))

        # proceed to compute the couplings from temp and Test_arry_tmp 
        # samples, in addition to the corresponding time_array_short
        # and time_array_long objects defined outside of the loop, resp

        for YI in range(len(couplings_short)):
            for init in range(len(couplings_short[YI])):
                for YK in range(len(couplings_short[YI])):
                    couplings_short[YI][YK] = ( abs(time_array_short[init] - time_array_short[YK]) )/(short_time_sum) * m.exp(-abs(couplings_short[YI][init] - couplings_short[YI][YK]))

        # perform the analagous array assignments for the long range couplings, 
        # which involve the normalisation by the time_array_long sum variable
        
        for YJJ in range(len(couplings_long)):
            for init in range(len(couplings_long[YJJ])):
                for YKJ in range(len(couplings_short[YJJ])):
                    couplings_long[YJJ][YKJ] = ( abs(time_array_long[init] - time_array_long[YKJ]) )/(long_time_sum) * m.exp(-abs(couplings_long[YJJ][init] - couplings_long[YJJ][YKJ]))


    return couplings_short, couplings_long, time_array_short, time_array_long

Couplings_temp = compute_couplings(Test_Arr)

