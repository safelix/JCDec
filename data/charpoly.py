import time
import random
import numpy as np
import pandas as pd

import os
os.system("lscpu") # print system information
from jordanchevalley import JCDec
print(JCDec)

random.seed(0)

N = 1
p = 0.5  # Probability for Bernoulli matrices
precisions = ["double", "mp::quad", "mp::oct"]

# generate sequence of matrix sizes
n_start = 5
n_step = 5
n_stop = 150
ns = np.arange(n_start, n_stop, n_step)

# function to store data
def store_data(meta, times, bitsizes, extra):

    # gather data 
    data = pd.concat( [meta, times, bitsizes, extra], axis = 1)
    
    # name of columns
    meta_names = ['n']
    time_names = [ 'time ' + name for name in precisions ]
    bitsize_names = [ 'bitsize ' + name for name in precisions ] 
    extra_names = ['err ' + name for name in precisions] + ['err time ' + name for name in precisions]

    # store
    print("store data:")
    print(data.tail(5))

    path = "short_charpoly_"+ str(p) + "_" +str(n_start) + "_" + str(n_step) + "_" + str(n_stop) + "_" + str(N) + ".csv"
    headers = (meta_names + time_names + bitsize_names + extra_names)
    data.to_csv(path, header=headers)


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

        timing = []
        bitsize = []
        err = []
        err_timing = []
        for i in range(0, len(precisions)):

            # compute the decomposition with precision[i]
            dec.compute_chiA(precision=precisions[i], round=True)
            timing.append( dec.get_timing() )
            bitsize.append( dec.get_maxBitsize() )

            err.append( dec.check_cayleyhamilton() )
            err_timing.append( dec.get_timing() )


        # store into dataframe
        meta = pd.concat( [meta, pd.DataFrame( [ n ] ).T ], ignore_index=True)
        times = pd.concat( [times, pd.DataFrame( timing ).T ], ignore_index=True )
        bitsizes = pd.concat( [bitsizes, pd.DataFrame( bitsize ).T ], ignore_index=True )
        extra = pd.concat( [extra, pd.DataFrame( err + err_timing ).T ], ignore_index=True )

        print( meta.tail(1).to_string(header=False) 
            + "     time: " + str(time.time() - start_time))

    store_data(meta, times, bitsizes, extra)


