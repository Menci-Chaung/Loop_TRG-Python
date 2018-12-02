#
#  partit.py
#  Loop_TRG
#  Calculation of the partition function for potts model
#
#  Copyright (C) 2018 Yue Zhengyuan, Liu Honghui and Zhang Wenqian. All rights reserved.
#  Article reference: Phys. Rev. Lett. 118, 110504 (2017)
#

import numpy as np
import filtering as flt
import optimizing as opt
#  this part just means to keep in accordance with the program "main.py" for convenience

beta = 1
q=3 # dimension of the bond in Potts model,which can be changed to other integers

def two(i,j):
    i1 = 2*np.pi*i/float(q)
    j1 = 2*np.pi*j/float(q)
    return beta*np.cos(i1-j1)
#the function defined just to reduce the error from too many multiplication

ts_TA = np.ones((q,q,q,q),dtype=complex)
#  initialize the tensor
for i in range(q):
    for j in range(q):
        for k in range(q):
            for l in range(q):
                ts_TA[i,j,k,l] = np.exp(two(i,j)+two(j,k)+two(k,l)+two(l,i))
ts_TB = ts_TA.copy()

# entanglement filtering

# loop optimize
ts_TA, ts_TB = opt.loop_optimize((ts_TA,ts_TB), 16, 10E-12)
