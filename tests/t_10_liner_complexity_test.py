import scipy.special as sc

from tests.src.utils import split_list, __print
from tests.src.bma import bma

# .10 Linear Complexity Test
def liner_complexity_test(key, n, M=500, b_print=True):
    if n <1000000:
        __print (b_print, '{:40} : Error. Need at  1,000,000 bits. Got {}.'.format('liner complexity test',n))
        return [0], False

    N=n//M

    K=6

    key = ''.join(list(map(str, key)))
    split_key=list(split_list(key,M))

    if len(split_key[-1]) != len(split_key[0]):
        split_key=split_key[0:-1]

    L=list(map(bma,split_key))

    myu = M/2 +(9+(-1)**(M+1))/36-(M/3+2/9)/2**M

    T = list(map(lambda x : (-1)**M*(x-myu)+2/9,L))

    v=[0]*7

    for i in T:
        if   i <= -2.5 : v[0] += 1
        elif i <= -1.5 : v[1] += 1
        elif i <= -0.5 : v[2] += 1
        elif i <=  0.5 : v[3] += 1
        elif i <=  1.5 : v[4] += 1
        elif i <=  2.5 : v[5] += 1
        else           : v[6] += 1

    pi = [0.010417,0.03125,0.125,0.5,0.25,0.0625,0.020833]

    chi_squared_obs = sum(list(map(lambda x, y : (x-N*y)**2/(N*y),v, pi )))

    p=sc.gammaincc(K/2,chi_squared_obs/2)

    b = (p >= 0.01)

    __print(b_print, '{:40} : {:.3f} -> {} '.format('liner complexity test',p,b))

    return [p], b