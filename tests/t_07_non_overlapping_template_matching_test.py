import glob
import random

import scipy.special as sc

from tests.src.utils import split_list

# N has been fixed at 8 in this code.
# I use the templates distributed by Nist's program.

# .7 Non-overlapping Template Matching Test
def non_overlapping_template_matching_test(key, n, m=9, B=None, N=8):
    if n <1000000:
        print ('{:40} : Error. Need at  1,000,000 bits. Got {}.'.format('non overlapping template matching test',n))
        return [0], False

    path_list = glob.glob('tests/templates/*')
    path_list.sort()
    path = path_list[m-2]
    with open(path) as f:
        template = ["".join(s.strip().split()) for s in f.readlines()]

    B = B if B is not None else random.choices(template)[-1]

    M=int(n/N)

    if M <= 0.01*n:
        print ('{:40} : Error. M is NOT selected correctly.'.format('non overlapping template matching test'))
        return [0], False

    myu=(M-m+1)/(2**m)
    sig2=M*(1/(2**m)-(2*m-1)/(2**(2*m)))
    
    split_key=list(split_list(key, M))

    if len(split_key[-1]) != len(split_key[0]):
        split_key=split_key[0:-1]

    def index_multi(txt, x, M):
        ind=0
        ind_total=0
        l=[]
        while ind_total <M:
            if x in txt:
                ind = txt.index(x)
                txt=txt[(ind+m):]
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

    w=list(map(float,w))

    chi_squared_obs = sum(list(map(lambda x : ((x-myu)**2)/sig2,w)))
    p=sc.gammaincc(N/2,chi_squared_obs/2)

    b = (p >= 0.01)

    print('{:40} : {:.3f} -> {} '.format('non overlapping template matching test',p,b))

    return [p], b