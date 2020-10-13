import math

import scipy.special as sc

from tests.src.utils import split_list, __print

# .13 Cumulative Sums (Cusum) Test
def cumulative_sums_test(key, n, b_print=True):
    if n < 100:
        __print(b_print, '{:40} : Error. Need at least 100 bits. Got {}.' .format('cumulative sums test', n))

    def compute(s):
        S=[0]*n

        for i in range(n):
            c = s[i]
            if   i == 0 : S[i]=c
            elif c == 1 : S[i] = S[i-1]+1
            else        : S[i] = S[i-1]-1

        abs_S = list(map(abs,S))

        z = max(abs_S)

        def fai(x):
            return 0.5 * math.erfc(-x * math.sqrt(0.5))

        k_start_1 = int(((-n)/z+1)/4)
        k_fin_1   = int((n/z-1)/4)

        num=[]
        for i in range(k_fin_1-k_start_1+1):
            num.append(k_start_1+i)

        pf = 0
        pb = 0

        for k in num:
            pf = pf + fai((4*k+1)*z/math.sqrt(n)) - fai((4*k-1)*z/math.sqrt(n))

        k_start_2 = int(((-n)/z-3)/4)
        k_fin_2   = int((n/z-1)/4)

        num=[]
        for i in range(k_fin_2-k_start_2+1):
            num.append(k_start_2+i)

        for k in num:
            pb = pb + fai((4*k+3)*z/math.sqrt(n)) - fai((4*k+1)*z/math.sqrt(n))

        p = 1 - pf + pb

        return p

    p_forward  = compute(key)
    p_backward = compute(list(reversed(key)))

    b1 = (p_forward >= 0.01)
    b2 = (p_backward >= 0.01)

    __print(b_print, '{:40} : {:.3f} -> {} '.format('cumulative sums test (forward )',p_forward,b1))
    __print(b_print, '{:40} : {:.3f} -> {} '.format('                     (backward)',p_backward,b2))

    return [p_forward, p_backward], all([b1, b2])