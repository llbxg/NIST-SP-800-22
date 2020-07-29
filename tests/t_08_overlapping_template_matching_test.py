import math

import scipy.special as sc

from tests.src.utils import split_list

# Check for a match with a template that consists of all ones.

# .8 Overlapping Template Matching Test
def overlapping_template_matching_test(key, n, m=9, M=1032):
    if n <1000000:
        print ('{:40} : Error. Need at  1,000,000 bits. Got {}.'.format('overlapping template matching test',n))
        return [0], False

    N=n//M
    K=5

    if not(n>=M*N):
        print ('{:40} : Error. "n>=M*N" This conditional expression is not satisfied.'.format('overlapping template matching test'))
        return [0], False

    # lam=(int(M-m+1)/(2**m)==2)

    B = "".join(['1' for x in range(m)])

    split_key=list(split_list(key,M))
    if len(split_key[-1]) != len(split_key[0]):
        split_key=split_key[0:-1]

    def index_multi(txt, x, M):
        ind=0
        ind_total=0
        l=[]
        while ind_total <M:
            if x in txt:
                ind = txt.index(x)
                txt=txt[(ind+1):]
                ind_total+=ind
                l.append(ind_total)
            else:
                return l

    w=[]
    for i in split_key:
        i = ''.join(map(str,i))
        txt=index_multi(i, B, M)
        num=len(txt)
        w.append(num)

    v=[0 for x in range(K+1)]

    for i in range(K+1):
        if i < K: 
            v[i]=w.count(i)
        else :
            v[i]=len([x for x in w if x  >= K ])

    pi = [0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865]  # sample???

    lambd = (M-m+1.0)/(2.0**m)
    eta = lambd/2.0

    pi[0] = math.e**(-eta)
    pi[1] = eta/2*math.e**(-eta)
    pi[2] = eta/8*math.e**(-eta)*(eta+2)
    pi[3] = eta/8*math.e**(-eta)*(eta**2/6+eta+1)
    pi[4] = eta/16*math.e**(-eta)*(eta**3/24+eta**2/2+3*eta/2+1)
    pi[5] = 1- sum(pi[0:5])

    if not (N*min(pi)>5):
        print ('{:40} : Error. "N*min(pi)>5" This conditional expression is not satisfied.'.format('overlapping template matching test'))
        return [0], False

    v=list(map(float,v))

    chi_squared_obs = sum(list(map(lambda x, y : ((x-N*y)**2)/(N*y),v,pi)))

    p=sc.gammaincc(5/2,chi_squared_obs/2)

    b = (p >= 0.01)

    print('{:40} : {:.3f} -> {} '.format('overlapping template matching test',p,b))

    return [p] , b