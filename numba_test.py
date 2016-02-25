from numba import jit
from numpy import arange

import time

# jit decorator tells Numba to compile this function.
# The argument types will be inferred by Numba when function is called.
@jit
def sum2d(arr):
    M, N = arr.shape
    result = 0.0
    for i in range(M):
        for j in range(N):
            result += arr[i,j]
    return result
def sum2d_slow(arr):
    M, N = arr.shape
    result = 0.0
    for i in range(M):
        for j in range(N):
            result += arr[i,j]
    return result

a = arange(64).reshape(8,8)

global_time_start = time.time()
print(sum2d_slow(a))
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)

global_time_start = time.time()
print(sum2d(a))
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)