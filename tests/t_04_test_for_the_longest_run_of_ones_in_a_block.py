import scipy.special as sc

from tests.src.utils import split_list, __print

# How to select N is the biggest problem.
# In the Nist program, |n / M| is adopted.
# On the other hand, in R1a, N is determined by K and M.
# This program adopted |n / M|.

# .4 Test for the Longest Run of Ones in a Block
def test_for_the_longest_run_of_ones_in_a_block(key, n, b_print=True):
    v=[0, 0, 0, 0, 0, 0, 0]

    def search_longest(Mbit):
        N=n//M
        split_key=list(split_list(key,M))

        if len(split_key[-1]) != len(split_key[0]):
            split_key=split_key[0:-1]

        ll=[]

        for i in split_key[:N]:
            longest = 0
            now_longest = 0

            for j in range(M):
                if i[j]==1:
                    now_longest += 1

                    if now_longest > longest:
                        longest = now_longest
                else :
                    now_longest = 0
            ll.append(longest)

        return ll

    def table_pi(M):
        M8      = [0.2148, 0.3672, 0.2305, 0.1875]
        M128    = [0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124]
        M10000  = [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]

        if   (M == 8)    : return M8
        elif (M == 128)  : return M128
        elif (M == 10000): return M10000

    if n >= 128 and n<6272:
        M = 8
        K = 3
        N=n//M
        ll=search_longest(M)
        pi=table_pi(M)

        for i in ll:
            if   i <= 1 : v[0] += 1
            elif i == 2 : v[1] += 1
            elif i == 3 : v[2] += 1
            else        : v[3] += 1

        chi_squared_obs = sum(list(map(lambda x, y : (x-N*y)**2/(N*y),v, pi )))
        p=sc.gammaincc(K/2,chi_squared_obs/2)

    elif n >= 6272 and n < 750000:
        M = 128
        K = 5
        N=n//M
        ll=search_longest(M)
        pi=table_pi(M)

        for i in ll:
            if   i <= 4 : v[0] += 1
            elif i == 5 : v[1] += 1
            elif i == 6 : v[2] += 1
            elif i == 7 : v[3] += 1
            elif i == 8 : v[4] += 1
            else        : v[5] += 1

        chi_squared_obs = sum(list(map(lambda x, y : (x-N*y)**2/(N*y),v, pi )))
        p=sc.gammaincc(K/2,chi_squared_obs/2)

    elif n >= 750000:
        M = 10000
        K = 6
        N=n//M
        ll=search_longest(M)
        pi=table_pi(M)

        for i in ll:
            if   i <= 10 : v[0] += 1
            elif i == 11 : v[1] += 1
            elif i == 12 : v[2] += 1
            elif i == 13 : v[3] += 1
            elif i == 14 : v[4] += 1
            elif i == 15 : v[5] += 1
            else         : v[6] += 1

        chi_squared_obs = sum(list(map(lambda x, y : (x-N*y)**2/(N*y),v, pi )))
        p=sc.gammaincc(K/2,chi_squared_obs/2)

    else :
        __print (b_print, '{:40} : Error. Need at least 128 bits. Got {}.'.format('the longest run of ones in a block', n))
        return [0], False

    b = (p >= 0.01)

    __print(b_print, '{:40} : {:.3f} -> {}'.format('the longest run of ones in a block',p,b))

    return [p], b