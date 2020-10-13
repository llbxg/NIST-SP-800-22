import math

import numpy as np
import scipy.special as sc

from tests.src.utils import split_list, __print
from tests.src.rakn import rank

# .5 Binary Matrix Rank Test
def binary_matrix_rank_test(key, n, M=32, Q=32, b_print=True):
    if n < 38912:
        __print(b_print, '{:40} : Error. Need at least 38,912 bits. Got {}.' .format('binary matrix rank test', n))
        return [0], False

    N=n//(M*Q)
    split_key=list(split_list(key,M*Q))

    if len(split_key[-1]) != len(split_key[0]):
        split_key=split_key[0:-1]

    f_split_key= list(map(lambda x : np.reshape(np.array(x),[M,Q]), split_key))

    ranks = list(map(lambda x: rank(list(x)),f_split_key))

    full_rank = M

    FM =ranks.count(full_rank)
    FM1=ranks.count(full_rank-1)
    NFMM1=N - FM -FM1

    chi_squared_obs = (FM-0.2888*N)**2/(0.2888*N) + (FM1-0.5776*N)**2/(0.5776*N)+(NFMM1-0.1336*N)**2/(0.1336*N)
    p=math.e**(-chi_squared_obs/2)

    b = (p >= 0.01)

    __print(b_print, '{:40} : {:.3f} -> {}'.format('binary matrix rank test',p,b))

    return [p],b