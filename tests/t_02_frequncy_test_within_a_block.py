import scipy.special as sc

from tests.src.utils import split_list, __print

# .2 Frequency Test within a Block
def frequency_test_with_block(key, n, M=20000, b_print=True):
    N=n//M
    split_key=list(split_list(key,M))

    if len(split_key[-1]) != len(split_key[0]):
        split_key=split_key[0:-1]

    pi=list(map(lambda x : (x.count(1))/M,split_key))

    chi_squared_obs=4*M*sum(list(map(lambda x : (x-1/2)**2, pi)))
    p=sc.gammaincc(N/2,chi_squared_obs/2)

    b = (p >= 0.01)

    __print(b_print, '{:25}( M = {:<5} )   : {:.3f} -> {}'.format('frequency test with block',M,p,b), end=' ')

    if M >= 20 and M > 0.01*n and N < 100:
        __print(b_print, '( M is selected correctly. )')

    else:
        __print(b_print, '( M is NOT selected correctly. )')

    return [p], b