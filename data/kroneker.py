import time
import random
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import os
os.system("lscpu") # print system information
from jordanchevalley import JCDec
print(JCDec)

random.seed(0)

precision = "mp::1000"
N = 1
init = np.array([[0.9, 0.1],[0.1, 0.9]]) # initial kroneker probabilities

# generate sequence of matrix sizes
n_start = 10
n_step = 1
n_stop = 11
ns = np.arange(n_start, n_stop, n_step)


# Function Definitions
################################################################################

# function to store data
def store_data(meta, times, bitsizes, extra):

    # gather data 
    data = pd.concat( [meta, times, bitsizes, extra], axis = 1)
    
    # name of columns
    names = ["A", "chiA", "muD", "inv", "chev", "D_mat", "N"]
    meta_names = ['n', 'no_iter', 'trivial']
    time_names = [ 'time ' + name for name in names ]
    bitsize_names = [ 'bitsize ' + name for name in names ] 
    extra_names = ['nill', 'comm']

    # store
    print("store data:")
    print(data.tail(5))

    headers = meta_names + time_names + bitsize_names + extra_names
    path = "kroneker_"+str(n_start) + "_" + str(n_step) + "_" + str(n_stop) + "_" + str(N) + ".csv"
    data.to_csv(path, header=headers)

def kroneker_graph(n):
    P = np.eye(1)
    for i in range(0,n):
        P = np.kron(init,P)
    return np.random.binomial(n=1, p=P)

    
def run():
    # colllect data
    meta = pd.DataFrame()
    times = pd.DataFrame()
    bitsizes = pd.DataFrame()
    extra = pd.DataFrame()

    for i in range(0, len(ns)):
        n = ns[i]

        for j in range(0,N):

            start_time = time.time()

            # initialize the defaults
            A = kroneker_graph(n)
            dec = JCDec(A.astype('float64'), verbosity=1)

            # compute the decomposition with defaults
            dec.compute_chiA(precision=precision, round=True)
            print( dec.get_timing() )
            dec.compute_muD(resultant=True)
            print( dec.get_timing() )
            dec.compute_D(mat=True)
            print( dec.get_timing() )
            dec.compute()

            # store data with defaults
            times = pd.concat( [times, pd.DataFrame(dec.get_timings())], ignore_index=True )
            bitsizes = pd.concat( [bitsizes, pd.DataFrame(dec.get_maxBitsizes())], ignore_index=True )

            # store extra
            nill = dec.check_nillpotency()
            comm = dec.check_commutativity()
            extra = pd.concat( [extra, pd.DataFrame( [nill, comm ] ).T ], ignore_index=True )

            # store meta data
            meta = pd.concat( [meta, pd.DataFrame( [dec.__size__, dec.__no_iter__, dec.is_trivial()] ).T ], ignore_index=True)

            print( meta.tail(1).to_string(header=False) 
                + "     time: " + str(time.time() - start_time) 
                + "s     bitsize: " + bitsizes.tail(1).max(axis=1).to_string(header=False, index=False)  )

        store_data(meta, times, bitsizes, extra)
    
    return



# START SCRIPT
################################################################################


run()
