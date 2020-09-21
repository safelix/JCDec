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

precision = "mp::oct"
N = 1
ps = [1/16]

n_start = 5
n_step = 5
n_stop = 150

n_jobs = len(ps)    # gather data for each p on a different core

ns = np.arange(n_start, n_stop, n_step)

# Function Definitions
################################################################################


# function to store data
def store_data(p, meta, times, bitsizes, extra):

    # gather data 
    data = pd.concat( [meta, times, bitsizes, extra], axis = 1)
    
    # name of columns
    names = ["A", "chiA", "muD", "inv", "chev", "D_poly", "N"]
    meta_names = ['n', 'no_iter', 'trivial']
    time_names = [ 'time ' + name for name in names ]
    bitsize_names = [ 'bitsize ' + name for name in names ] 
    extra_names = ['time muD (res)', 'time D_mat', 'diff', 'nill', 'comm']

    # store
    print("store data:")
    print(data.tail(5))

    headers = meta_names + time_names + bitsize_names + extra_names
    path = "full_bernoulli_"+ str(p) + "_" +str(n_start) + "_" + str(n_step) + "_" + str(n_stop) + "_" + str(N) + ".csv"
    data.to_csv(path, header=headers)

def run(p):
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
            A = np.random.choice(a=[0,1], p=[1-p,p], size=(n,n))
            dec = JCDec(A.astype('float64'))

            # compute the decomposition with defaults
            dec.compute_chiA(precision=precision, round=True)
            D = dec.compute().__D__

            # store data with defaults
            times = pd.concat( [times, pd.DataFrame(dec.get_timings())], ignore_index=True )
            bitsizes = pd.concat( [bitsizes, pd.DataFrame(dec.get_maxBitsizes())], ignore_index=True )

            # compute with extra
            dec.compute_muD(resultant=True)
            time1 = dec.get_timing()
            dec.compute_D(mat=True)
            time2 = dec.get_timing()

            # store extra
            diff = np.linalg.norm(D - dec.compute_D(mat=True).__D__)
            nill = dec.check_nillpotency()
            comm = dec.check_commutativity()
            extra = pd.concat( [extra, pd.DataFrame( [time1, time2, diff, nill, comm ] ).T ], ignore_index=True )

            # store meta data
            meta = pd.concat( [meta, pd.DataFrame( [dec.__size__, dec.__no_iter__, dec.is_trivial()] ).T ], ignore_index=True)

            print( meta.tail(1).to_string(header=False) 
                + "     time: " + str(time.time() - start_time) 
                + "s     bitsize: " + bitsizes.tail(1).max(axis=1).to_string(header=False, index=False)  )

        store_data(p, meta, times, bitsizes, extra)



# START SCRIPT in PARALLEL
################################################################################


Parallel(n_jobs=n_jobs)( delayed(run)(proba) for proba in ps )
