import math

import scipy.special as sc

from tests.src.utils import split_list

# .15 Random Excursions Variant Test
def random_excursion_variant_test(key, n):
    S=[0]*n

    for i in range(n):
        c = key[i]
        if   i == 0 : S[i]=c
        elif c == 1 : S[i] = S[i-1]+1
        else        : S[i] = S[i-1]-1

    S_dash =[0] + S + [0]
    J=S_dash[1:].count(0)

    S_dash= list(map(lambda x : (x+9),S_dash))
    state = [-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9]

    zeta =[0]*18

    for c, i in enumerate(list(map(lambda x : (x+9),state))):
        zeta[c]=S_dash.count(i)

    p=[0]*18
    ps, bs = [], []
    for i in range(len(p)):
        p = math.erfc(abs(zeta[i]-J)/math.sqrt(2*J*(4*abs(state[i])-2)))
        b = (p >= 0.01)

        ps.append(p)
        bs.append(b)

        if   i == 0 : print('{:40} : {:.3f} -> {} '.format('random excursions test (x = {:2})'.format(state[i]),p,b))
        else          : print('{:40} : {:.3f} -> {} '.format('                       (x = {:2})'.format(state[i]),p,b))

    return ps, all(bs)