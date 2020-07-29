import math

import scipy.special as sc

from tests.src.utils import split_list

# .12 Approximate Entropy Test
def approximate_entropy_test(key, n, m=5):
    if 2**m > n:
        print ('{:40} : Error. m is too big .'.format('approximate entropy test'))
        return [0], False

    def compute(s,m):

        if m == 0:
            return 0

        if m == 1: head = ''
        else     : head = s[0:(m-1)]

        s = s + head

        v = [0]*2**m

        for i in range(m):
            ss=s[i:]

            split_key_m=list(split_list(ss,m))
            if len(split_key_m[-1]) != len(split_key_m[0]):
                split_key_m=split_key_m[0:-1]

            split_key_m=list(map(lambda x : int(x,2),split_key_m))
            for i in range(2**m):
                v[i] = v[i]+split_key_m.count(i)

        return v

    c1 = [0]*2**m

    key = ''.join(list(map(str, key)))
    v1 = compute(key,m)

    for i in range(2**m):
        c1[i] = v1[i]/n

    fai1 = sum(list(map(lambda x : x*math.log(x),c1)))

    c2 = [0]*2**(m+1)

    v2 = compute(key,(m+1))

    for i in range(2**(m+1)):
        c2[i] = v2[i]/n

    fai2 = sum(list(map(lambda x : x*math.log(x),c2)))

    ApEn=fai1-fai2

    chi_squared_obs = 2*n*(math.log(2)-ApEn)
    p=sc.gammaincc(2**(m-1),chi_squared_obs/2)

    b = (p >= 0.01)

    print('{:40} : {:.3f} -> {} '.format('approximate entropy test',p,b),end = ' ')
    if m >= math.log2(n)-7:
        print('(m : possess lower reliability)')
    else:
        print('')

    return [p], b