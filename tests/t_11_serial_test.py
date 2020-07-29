import scipy.special as sc

from tests.src.utils import split_list

# .11 Serial Test
def serial_test(key, n, m=3):
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

        psi2_m = 2**m/n*(sum(list(map(lambda x : x**2,v)))) - n

        return psi2_m

    key = ''.join(list(map(str, key)))
    psi2_m0=compute(key,m)
    psi2_m1=compute(key,(m-1))
    psi2_m2=compute(key,(m-2))

    d_psi2  = psi2_m0 -   psi2_m1
    d2_psi2 = psi2_m0 - 2*psi2_m1 + psi2_m2

    p1=sc.gammaincc(2**(m-2),d_psi2  / 2)
    p2=sc.gammaincc(2**(m-3),d2_psi2 / 2)

    b1 = (p1 >= 0.01)
    b2 = (p2 >= 0.01)

    print('{:40} : {:.3f} -> {} '.format('serial test',p1,b1))
    print('{:40} : {:.3f} -> {} '.format('',p2,b2))

    return [p1, p2], all([b1,  b2])