import math

import scipy.special as sc

from tests.src.utils import split_list, __print

# .9 Maurer’s “Universal Statistical” Test
def maurers_universal_statical_test(key, n, b_print=True):
    if n <387840:
        __print (b_print, '{:40} : Error. Need at  387,840 bits. Got {}.'.format('maurers universal statical test',n))
        return [0], False

    L=6
    nlist = [904960, 2068480, 4654080, 10342400, 22753280, 49643520, 107560960, 231669760, 496435200, 1059061760]
    num= len(list(filter(lambda x:n > x, nlist)))
    L+=num

    Q=10*2**L

    key = ''.join(list(map(str, key)))
    split_key=list(split_list(key,L))

    if len(split_key[-1]) != len(split_key[0]):
        split_key=split_key[0:-1]

    K=len(split_key)-Q

    split_key=list(map(lambda x : int(x,2),split_key))

    T = [0] * (2**L)

    for i in range(Q):
        contents = split_key[i]
        T[contents]=i

    sums=0.0

    for i in range(Q,Q+K):
        contents=split_key[i]
        bnum=T[contents]

        sums+=math.log2(i-bnum)

        T[contents]=i

    fn=sums/float(K)

    expectedValue_list = [5.2177052,6.1962507,7.1836656,8.1764248,9.1723243,10.170032,11.168765,12.168070,13.167693,14.167488,15.167379]
    variance_list      = [2.954,3.125,3.238,3.311,3.356,3.384,3.401,3.410,3.416,3.419,3.421]

    myu = expectedValue_list[num]
    v   = variance_list[num]

    c   = 0.7-0.8/L+(4+32/L)*(K**(-3/L)/15)
    sig = c*math.sqrt(v/K)

    p=math.erfc(abs((fn-myu)/(math.sqrt(2)*sig)))

    b = (p >= 0.01)

    __print(b_print, '{:40} : {:.3f} -> {} '.format('maurers universal statical test',p,b))

    return [p], b