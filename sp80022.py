from decimal import *

import numpy as np

import random

from math import *

import matplotlib.pyplot as plt
import time


from scipy.special import gammaincc

from scipy.sparse.linalg import svds

def split_list(l, n):
    for idx in range(0, len(l), n):
        yield l[idx:idx + n]

class nist_800_22():

    def __init__(self,key16,sixteen=True):

        print('hello ! this is nist 800-22 program')
        print('----------------------------------------------------------------------------')

        if sixteen:
            self.key16=key16
            self.key=['0']*len(self.key16)
            for i in range(len(self.key16)):
                k=self.key16[i]
                ind = ['0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'].index(k)
                num = ['0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111']
                self.key[i]=num[ind]
    
            self.key=''.join(self.key)
        else :
            self.key=key16


        self.n=len(self.key)
        self.zeros=self.key.count('0')
        self.ones=self.key.count('1')
        self.b_monobit_test=None
        self.b_monobit_test_with_block=None
        self.b_run_test=None

        self.tau=2/sqrt(self.n)

        if 100 >= self.n:
            print('key size is short')

    def running(self):
        start = time.time()
        self.monobit_test()
        print (round(time.time() - start,1))

        start = time.time()
        self.monobit_test_with_block()
        print (round(time.time() - start,1))

        start = time.time()
        self.run_test()
        print (round(time.time() - start,1))

        start = time.time()
        self.test_for_the_longest_run_of_ones_in_a_block()
        print (round(time.time() - start,1))

        start = time.time()
        self.bmrt()
        print (round(time.time() - start,1))

        #start = time.time()
        #self.discrete_fourier_transform_test_by_NIST()
        #print (round(time.time() - start,1))

        start = time.time()
        self.non_overlapping_template_matching_test()
        print (round(time.time() - start,1))

        start = time.time()
        self.overlapping_template_matching_test()
        print (round(time.time() - start,1))

        start = time.time()
        self.maurers_universal_statical_test()
        print (round(time.time() - start,1))

        start = time.time()
        self.liner_complexity_test()
        print (round(time.time() - start,1))

        start = time.time()
        self.serial_test()
        print (round(time.time() - start,1))

        start = time.time()
        self.approximate_entropy_test(m=10)
        print (round(time.time() - start,1))

        start = time.time()
        self.cumulative_sums_test()
        print (round(time.time() - start,1))

        start = time.time()
        self.random_excursions_test()
        print (round(time.time() - start,1))

        start = time.time()
        self.random_excursion_variant_test()
        print (round(time.time() - start,1))

        start = time.time()

    def monobit_test(self):

        sn=abs(self.ones-self.zeros)
        sobs=sn/sqrt(self.n)
        p=erfc(sobs/sqrt(2))

        if p<0.005 or p>0.995:
            b=False
        else:
            b=True 

        print('{:40} : {:.3f} -> {}'.format('monobit test',p,b))

        self.b_monobit_test=b

        return p, b

    def monobit_test_with_block(self,M=20000):
        N=self.n//M
        split_key=list(split_list(self.key,M))

        if len(split_key[-1]) != len(split_key[0]):
            split_key=split_key[0:-1]

        pi=list(map(lambda x : (x.count('1'))/M,split_key))

        chi_squared_obs=4*M*sum(list(map(lambda x : (x-1/2)**2, pi)))
        p=gammaincc(N/2,chi_squared_obs/2)

        if p<0.01:
            b=False
        else:
            b=True


        print('{:25}( M = {:<5} )   : {:.3f} -> {}'.format('monobit test with block',M,p,b), end=' ')

        if M >= 20 and M > 0.01*self.n and N < 100:
            print('( M is selected correctly. )')

        else:
            print('( M is NOT selected correctly. )')

        self.b_monobit_test_with_block=b

        return p, b

    def run_test(self):
        pi = self.ones/self.n

        if abs(pi-1/2) >= self.tau or not(self.b_monobit_test):
            print('run test Failed')

            self.b_run_test=False
            p=0
            
            return p, self.b_run_test

        vobs=1
        for i in range(self.n-1):
            if self.key[i]==self.key[i+1]:
                vobs+=0
            else:
                vobs+=1
        

        p=erfc(abs(vobs-2*float(self.n)*pi*(1-pi))/(2*sqrt(2*float(self.n))*pi*(1-pi)))
        

        if p<0.01:
            b=False
        else:
            b=True

        print('{:40} : {:.3f} -> {}'.format('run test',p,b))

        self.b_run_test=b


        return p, b

    def test_for_the_longest_run_of_ones_in_a_block(self):
        #Nをどう採用するかの問題あり。NistのProgramを選択。

        v=[0, 0, 0, 0, 0, 0, 0]

        def search_longest(Mbit):
            N=self.n//M
            split_key=list(split_list(self.key,M))

            if len(split_key[-1]) != len(split_key[0]):
                split_key=split_key[0:-1]

            ll=[]

            for i in split_key[:N]:
                longest = 0
                now_longest = 0

                for j in range(M):
                    if i[j]=='1':
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

        if self.n >= 128 and self.n<6272:
            M = 8
            K = 3
            N=self.n//M
            ll=search_longest(M)
            pi=table_pi(M)

            for i in ll:
                if   i <= 1 : v[0] += 1
                elif i == 2 : v[1] += 1
                elif i == 3 : v[2] += 1
                else        : v[3] += 1

            chi_squared_obs = sum(list(map(lambda x, y : (x-N*y)**2/(N*y),v, pi )))
            p=gammaincc(K/2,chi_squared_obs/2)

        elif self.n >= 6272 and self.n < 750000:
            M = 128
            K = 5
            N=self.n//M
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
            p=gammaincc(K/2,chi_squared_obs/2)

        elif self.n >= 750000:
            M = 10000
            K = 6
            N=self.n//M
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
            p=gammaincc(K/2,chi_squared_obs/2)

        else :
            print ('{:40} : Error. Need at least 128 bits. Got {}.'.format('the longest run of ones in a block',self.n))
            return 0, False

        if p < 0.01:
            b = False

        else:
            b = True

        print('{:40} : {:.3f} -> {}'.format('the longest run of ones in a block',p,b))


        return p, b

    def binary_matrix_rank_test(self,M=32,Q=32):
        if self.n < 38912: #N=38
            print('{:40} : Error. Need at least 38,912 bits. Got {}.' .format('binary matrix rank test', self.n))

            return 0, False

        N=self.n//(M*Q)
        split_key=list(split_list(self.key,M*Q))

        if len(split_key[-1]) != len(split_key[0]):
            split_key=split_key[0:-1]


        f_split_key= list(map(lambda x : np.reshape(np.array(list(map(lambda y : int(y), x))),[M,Q]), split_key))
        

        exm=f_split_key[2]
        print(np.triu(exm))
        print(np.any(np.triu(exm) == 1, axis=0))
        print(np.any(np.triu(exm) == 1, axis=1))

        def rank(m,r):
            triu_m=np.triu(m)
            in_col_one=np.any(triu_m == 1, axis=0)
            in_row_one=np.any(triu_m == 1, axis=1)
            cc=list(in_col_one).count(False)
            cr=list(in_row_one).count(False)
            print(cr,cc)

            m=list(triu_m)

            def check(x):
                if 1 in x : return x.index(1)
                else      : return 0
            n = list(map(lambda x : check(list(x)),m))
            #print(n)
            c = 0
            for i in range(32):
                if i in n : c+=1
            print(c)
            print(np.shape(triu_m[:r-cr, cc:]))
            return (triu_m[:r-cr, cc:])
        
        print(rank(exm,32))

        ranks = list(map(lambda x: np.linalg.matrix_rank(x),f_split_key))


        full_rank = M

        length_ranks=len(ranks)

        FM =ranks.count(full_rank)
        FM1=ranks.count(full_rank-1)
        NFMM1=N - FM -FM1

        print( FM, FM1, NFMM1)

        chi_squared_obs = (FM-0.2888*N)**2/(0.2888*N) + (FM1-0.5776*N)**2/(0.5776*N)+(NFMM1-0.1336*N)**2/(0.1336*N)
        p=e**(-chi_squared_obs/2)


        if p < 0.01:
            b = False

        else:
            b = True

        print('{:40} : {:.3f} -> {} (Note - unfinished making this test!)'.format('binary matrix rank test',p,b))

        return p,b

    def non_overlapping_template_matching_test(self,N=8,B='000000001'):
        if self.n <1000000:
            print ('{:40} : Error. Need at  1,000,000 bits. Got {}.'.format('non overlapping template matching test',self.n))
            return 0, False

        ep = self.key[:1000000]
        m=len(B)
        M=int(len(ep)/N)

        myu=(M-m+1)/(2**m)
        sig2=M*(1/(2**m)-(2*m-1)/(2**(2*m)))
    
        split_key=list(split_list(ep,M))

        if len(split_key[-1]) != len(split_key[0]):
            split_key=split_key[0:-1]

        exist = list(map(lambda x : B in x, split_key))

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
            txt=index_multi(i, B, M)
            num=len(txt)
            w.append(num)

        w=list(map(float,w))

        chi_squared_obs = sum(list(map(lambda x : ((x-myu)**2)/sig2,w)))
        
        p=gammaincc(N/2,chi_squared_obs/2)

        if p < 0.01:
            b = False

        else:
            b = True

        print('{:40} : {:.3f} -> {} '.format('non overlapping template matching test',p,b))

        return 0 , False

    def overlapping_template_matching_test(self,N=8,B='111111111',M=1032):
        if self.n <1000000:
            print ('{:40} : Error. Need at  1,000,000 bits. Got {}.'.format('overlapping template matching test',self.n))
            return 0, False

        N=self.n//1032

        ep = self.key[:N*M]
        m=len(B)

        myu=(M-m+1)/(2**m)
        sig2=M*(1/(2**m)-(2*m-1)/(2**(2*m)))

        split_key=list(split_list(ep,M))

        if len(split_key[-1]) != len(split_key[0]):
            split_key=split_key[0:-1]

        exist = list(map(lambda x : B in x, split_key))
        

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
            txt=index_multi(i, B, M)
            num=len(txt)
            w.append(num)

        v=[0,0,0,0,0,0]
        for i in range(6):
            if i < 5: 
                v[i]=w.count(i)
            else :
                v[i]=len([x for x in w if x  >= 5 ])

        pi = [0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865] #sample???

        lambd = (M-m+1.0)/(2.0**m)
        eta = lambd/2.0

        pi[0] = e**(-eta)
        pi[1] = eta/2*e**(-eta)
        pi[2] = eta/8*e**(-eta)*(eta+2)
        pi[3] = eta/8*e**(-eta)*(eta**2/6+eta+1)
        pi[4] = eta/16*e**(-eta)*(eta**3/24+eta**2/2+3*eta/2+1)
        pi[5] = 1- sum(pi[0:5])

        v=list(map(float,v))
        

        chi_squared_obs = sum(list(map(lambda x, y : ((x-N*y)**2)/(N*y),v,pi)))

        p=gammaincc(5/2,chi_squared_obs/2)

        if p < 0.01:
            b = False

        else:
            b = True

        print('{:40} : {:.3f} -> {} '.format('overlapping template matching test',p,b))

        return p , b

    def maurers_universal_statical_test(self):

        if self.n <387840:
            print ('{:40} : Error. Need at  387,840 bits. Got {}.'.format('maurers universal statical test',self.n))
            return 0, False


        L=7
        nlist = [904960, 2068480, 4654080, 10342400, 22753280, 49643520, 107560960, 231669760, 496435200, 1059061760]
        num= len(list(filter(lambda x:self.n > x, nlist)))
        L+=num
        
        Q=10*2**L

        split_key=list(split_list(self.key,L))

        if len(split_key[-1]) != len(split_key[0]):
            split_key=split_key[0:-1]

        K=len(split_key)-Q

        split_key=list(map(lambda x : int(x,2),split_key))

        T = [0] * (2**L)

        for i in range(Q):
            contents = split_key[i]
            T[contents]=i

        sum=0.0

        for i in range(Q,Q+K):
            contents=split_key[i]
            bnum=T[contents]

            sum+=log2(i-bnum)

            T[contents]=i

        fn=sum/float(K)

        expectedValue_list = [5.2177052,6.1962507,7.1836656,8.1764248,9.1723243,10.170032,11.168765,12.168070,13.167693,14.167488,15.167379]
        variance_list      = [2.954,3.125,3.238,3.311,3.356,3.384,3.401,3.410,3.416,3.419,3.421]

        myu = expectedValue_list[num]
        v   = variance_list[num]


        c   = 0.7-0.8/L+(4+32/L)*(K**(-3/L)/15)
        sig = c*sqrt(v/K)

        p=erfc(abs((fn-myu)/(sqrt(2)*sig)))
        print(p)

        if p<0.01 :
            b=False
        else:
            b=True 

        print('{:40} : {:.3f} -> {} '.format('maurers universal statical test',p,b))

        return p, b

    def liner_complexity_test(self,M=500):
        if self.n <1000000:
            print ('{:40} : Error. Need at  1,000,000 bits. Got {}.'.format('liner complexity test',self.n))
            return 0, False

        N=self.n//M

        K=6

        split_key=list(split_list(self.key,M))

        if len(split_key[-1]) != len(split_key[0]):
            split_key=split_key[0:-1]

        #copy ;; i have to rewrite !
        def bma(s):
            n = len(s)
            c = np.zeros(n)
            b = np.zeros(n)
            c[0], b[0] = 1, 1
            l, m, i = 0, -1, 0
            int_s = [int(el) for el in s]
            while i < n:
                v = int_s[(i - l):i]
                v = v[::-1]
                cc = c[1:l + 1]
                d = (int_s[i] + np.dot(v, cc)) % 2
                if d == 1:
                    temp = np.copy(c)
                    p = np.zeros(n)
                    for j in range(0, l):
                        if b[j] == 1:
                            p[j + i - m] = 1
                    c = (c + p) % 2
                    if l <= 0.5 * i:
                        l = i + 1 - l
                        m = i
                        b = temp
                i += 1
            return l

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

        p=gammaincc(K/2,chi_squared_obs/2)

        if p<0.01 :
            b=False
        else:
            b=True 

        print('{:40} : {:.3f} -> {} '.format('liner complexity test',p,b))

        return p, b

    def serial_test(self, m=5):

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

            psi2_m = 2**m/self.n*(sum(list(map(lambda x : x**2,v)))) - self.n

            return psi2_m

        psi2_m0=compute(self.key,m)
        psi2_m1=compute(self.key,(m-1))
        psi2_m2=compute(self.key,(m-2))

        d_psi2  = psi2_m0 -   psi2_m1
        d2_psi2 = psi2_m0 - 2*psi2_m1 + psi2_m2

        p1=gammaincc(2**(m-2),d_psi2  / 2)
        p2=gammaincc(2**(m-3),d2_psi2 / 2)

        if p1<0.01 :
            b1=False
        else:
            b1=True 

        if p2<0.01 :
            b2=False
        else:
            b2=True 

        print('{:40} : {:.3f} -> {} '.format('serial test',p1,b1))
        print('{:40} : {:.3f} -> {} '.format('',p2,b2))

        return p1, b1

    def approximate_entropy_test(self, m=5):

        if 2**m > self.n:
            print ('{:40} : Error. m is too big .'.format('approximate entropy test'))
            return 0, False

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

        v1 = compute(self.key,m)

        for i in range(2**m):
            c1[i] = v1[i]/self.n

        fai1 = sum(list(map(lambda x : x*log(x),c1)))

        c2 = [0]*2**(m+1)

        v2 = compute(self.key,(m+1))

        for i in range(2**(m+1)):
            c2[i] = v2[i]/self.n

        fai2 = sum(list(map(lambda x : x*log(x),c2)))


        ApEn=fai1-fai2

        chi_squared_obs = 2*self.n*(log(2)-ApEn)

        p=gammaincc(2**(m-1),chi_squared_obs/2)

        if p<0.01 :
            b=False
        else:
            b=True 


        print('{:40} : {:.3f} -> {} '.format('approximate entropy test',p,b),end = ' ')
        if m >= log2(self.n)-7: # NIST  log2(self.n)-5
            print('(m : possess lower reliability)')
        else:
            print('')

        return p, b

    def cumulative_sums_test(self):
        if self.n < 100: 
            print('{:40} : Error. Need at least 100 bits. Got {}.' .format('cumulative sums test', self.n))

        def compute(s):
            S=[0]*self.n

            key = list(map(int,s))

            for i in range(self.n):
                c = key[i]
                if   i == 0 : S[i]=c
                elif c == 1 : S[i] = S[i-1]+1
                else        : S[i] = S[i-1]-1

            abs_S = list(map(abs,S))

            z = max(abs_S)

            def fai(x):
                return 0.5 * erfc(-x * sqrt(0.5))

            k_start_1 = int(((-self.n)/z+1)/4)
            k_fin_1   = int((self.n/z-1)/4)

            num=[]
            for i in range(k_fin_1-k_start_1+1):
                num.append(k_start_1+i)

            pf = 0
            pb = 0

            for k in num:
                pf = pf + fai((4*k+1)*z/sqrt(self.n)) - fai((4*k-1)*z/sqrt(self.n)) 

            k_start_2 = int(((-self.n)/z-3)/4)
            k_fin_2   = int((self.n/z-1)/4)


            num=[]
            for i in range(k_fin_2-k_start_2+1):
                num.append(k_start_2+i)

            for k in num:
                pb = pb + fai((4*k+3)*z/sqrt(self.n)) - fai((4*k+1)*z/sqrt(self.n)) 

            p = 1 - pf + pb

            return p

        p_forward  = compute(self.key)
        p_backward = compute(reversed(self.key))

        if p_forward<0.01 :
            b1=False
        else:
            b1=True 

        if p_backward<0.01 :
            b2=False
        else:
            b2=True 

        print('{:40} : {:.3f} -> {} '.format('cumulative sums test (forward )',p_forward,b1))
        print('{:40} : {:.3f} -> {} '.format('                     (backward)',p_backward,b2))

        return p_forward, b1

    def random_excursions_test(self):
        S=[0]*self.n

        key = list(map(int,self.key))

        for i in range(self.n):
            c = key[i]
            if   i == 0 : S[i]=c
            elif c == 1 : S[i] = S[i-1]+1
            else        : S[i] = S[i-1]-1

        #print(S)

        S_dash =[0] + S + [0]

        #print(S_dash[:100])
        #print(S_dash[999900:])
        #print(S_dash[1:])

        zeros = S_dash[1:].count(0)

        #print(zeros)

        c = 0

        while c<0:
            ind = S_dash.index(0)

        def index_multi(txt, x, n):
            ind=0
            ind_total=0
            l=[]
            while (ind_total-1) <n:
                if x in txt:
                    
                    ind = txt.index(x)
                    ind_total=ind_total+len(txt[:(ind+1)])
                    l.append(ind_total-1)
                    txt=txt[(ind+1):] 

                    
                else:
                    return l
            return l

        ind_list=index_multi(S_dash, 0, self.n)
        J = len(ind_list) - 1

        S_dash= list(map(lambda x : str(x+5),S_dash))
        #print(S_dash[:100])

        cyc_list=[]

        for i in range(len(ind_list)):
            if not ((i+1) >= len(ind_list)):
                start = ind_list[i]
                fin   = ind_list[i+1]

                cyc = S_dash[start:fin+1]
                cyc_list.append(cyc)

        #print(cyc_list[:10])

        cycle=[]
        for i in range(len(cyc_list)):
            v=[]
            for j in [1,2,3,4,6,7,8,9]:
                v.append(cyc_list[i].count(str(j)))
            cycle.append(v)

        #print(cycle[:10])

        states =[]
        for i in range(8):
            state=[]
            for j in range(len(cycle)):
                state.append(cycle[j][i])
            states.append(state)
        #print(len(states))
        #print(states)
        
        vs=[]
        #print(vs)

        for i in states:
            vv=[0,0,0,0,0,0]
            vv[0]=len(list(filter(lambda x : x==0, i)))
            vv[1]=len(list(filter(lambda x : x==1, i)))
            vv[2]=len(list(filter(lambda x : x==2, i)))
            vv[3]=len(list(filter(lambda x : x==3, i)))
            vv[4]=len(list(filter(lambda x : x==4, i)))
            vv[5]=len(list(filter(lambda x : x>=5, i)))

            vs.append(vv)

            #print(c)

        #print(vs)

        pis=[
        [0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000],
        [0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625],
        [0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143],
        [0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051]]

        #print(list(map(lambda x : sum(x),vs)))
        #J = sum(vs[0])
        #print(J)

        cnt=0
        for i in [3,2,1,0,0,1,2,3]:
            c = sum(list(map(lambda x,y : (x-J*y)**2/(J*y),vs[cnt],pis[i])))
            #print(c)
            p = gammaincc(5/2,c/2)
            #print(p)
            

            if p<0.01 :
                b=False
            else:
                b=True 
            st = [-4,-3,-2,-1,1,2,3,4]
            if   cnt == 0 : print('{:40} : {:.3f} -> {} '.format('random excursions test (x = {:2})'.format(st[cnt]),p,b))
            else          : print('{:40} : {:.3f} -> {} '.format('                       (x = {:2})'.format(st[cnt]),p,b))
            cnt+=1

    def random_excursion_variant_test(self):


        S=[0]*self.n

        key = list(map(int,self.key))

        for i in range(self.n):
            c = key[i]
            if   i == 0 : S[i]=c
            elif c == 1 : S[i] = S[i-1]+1
            else        : S[i] = S[i-1]-1

        #print(S)

        S_dash =[0] + S + [0]
        J=S_dash[1:].count(0)

        #print(S_dash[:100])
        #print(S_dash[999900:])
        #print(S_dash[1:])

        zeros = S_dash[1:].count(0)

       # print(zeros)

        def index_multi(txt, x, n):
            ind=0
            ind_total=0
            l=[]
            while (ind_total-1) <n:
                if x in txt:
                    
                    ind = txt.index(x)
                    ind_total=ind_total+len(txt[:(ind+1)])
                    l.append(ind_total-1)
                    txt=txt[(ind+1):] 

                    
                else:
                    return l
            return l

        S_dash= list(map(lambda x : (x+9),S_dash))
        state = [-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9]

        zeta =[0]*18
        c=0

        for i in list(map(lambda x : (x+9),state)):

            zeta[c]=S_dash.count(i)

            c+=1

        #print(zeta)

        p=[0]*18

        for i in range(len(p)):
            p=erfc(abs(zeta[i]-J)/sqrt(2*J*(4*abs(state[i])-2)))

            if p<0.01 :
                b=False
            else:
                b=True 
            if   i == 0 : print('{:40} : {:.3f} -> {} '.format('random excursions test (x = {:2})'.format(state[i]),p,b))
            else          : print('{:40} : {:.3f} -> {} '.format('                       (x = {:2})'.format(state[i]),p,b))

    def discrete_fourier_transform_test_by_NIST(self):
        key=[]
        for i in self.key:
            k=1.0
            if i == '0':
                k=-1.0
            key.append(k)

        def fourier_spectrum(x,j):
            jj=float(j)
            n=float(self.n)
            coss=0
            sins=0
            for k in  range(self.n//2):
                kk=float(k)
                coss += x[k]*cos(2*pi*kk*jj/n)
                sins += x[k]*sin(2*pi*kk*jj/n)
            return abs(sqrt(coss**2+sins**2))

        T=sqrt(2.995732274*self.n)

        key=np.array(key)
        S=list(np.abs(np.fft.fft(key)))
        N1=len(list(filter(lambda x:x <T, S[:int(self.n/2)])))

        N0=0.95*self.n/2

        d= (N1-N0)/sqrt(0.95*0.05*self.n/4)

        p = erfc(abs(d)/sqrt(2))

        if p<0.01 :
            b=False
        else:
            b=True 

        print('{:40} : {:.3f} -> {} '.format('discrete fourier transform test by NIST',p,b))

        return p, b

    def neo_discrete_fourier_transform_test(self):
        key=[]
        for i in self.key:
            k=1.0
            if i == '0':
                k=-1.0
            key.append(k)

        key=np.array(key)
        S=list(np.abs(np.fft.fft(key)))

        sumS4=sum(list(map(lambda x : x**4,S[:int(self.n/2)])))

        V=1/sqrt(2*self.n**5)*sumS4-sqrt(self.n/2)

        p = erfc(abs(V)/sqrt(2))

        if p<0.01 :
            b=False
        else:
            b=True 

        print('{:40} : {:.3f} -> {} '.format('neo discrete fourier transform test',p,b))

        return p, b

    def bmrt(self,M=32,Q=32):

        if self.n < 38912: #N=38
            print('{:40} : Error. Need at least 38,912 bits. Got {}.' .format('binary matrix rank test', self.n))

            return 0, False

        #exm=np.array([[0,1,1,1],[1,1,0,1],[1,1,1,1],[1,0,0,0]])

        def exchange_rows(m, i, j):
            if i==j:
                return m
            temp=np.copy(m[i])
            m[i]=m[j]
            m[j]=temp
            return m

        def one_or_zero(m, row, col,row_num,col_num):
            for i in range(row_num):
                if m[row+i][col]==1:
                    #print(row+i, row)
                    return exchange_rows(m, row, row+i), False
            return m, True

        def plus(m, row, v):
            copy_m=m[:]
            m[row]=copy_m[row]+v
            return m%2
         
        def see_you(m,row,col,row_num, col_num):
            v=m[row]
            for i in range(row_num-1):
                if m[row+i+1][col]==1:
                    #print(row+i, row)
                    m=plus(m, row+i+1,v )
            return m

        def compute(m,M,Q):
            now_row=0
            now_col=0
            new_m=np.copy(m)
            for i in range(M):
                #print(i)
                new_m,b =one_or_zero(new_m,i,now_col,M-i,Q-now_col)
                #print(i)
                #print(new_m)
                #print(i)
                new_m=see_you(new_m, i,now_col,M-i,Q-now_col )
                #print(new_m)
                now_col+=1
            return new_m

        def rank(m,M,Q):
            new_m=compute(m,M,Q)

            b=np.all(new_m==0,axis=1)
            zero_num=list(b).count(True)
            #print(new_m)

            part_m =np.copy(new_m[:M-zero_num])
            zero_m = np.copy(new_m[M-zero_num:])
            #print(part_m)

            for i in reversed(range(M-zero_num)):
                if 1 in list(part_m[i]):
                    ind=((list(part_m[i])).index(1))
                    #print(ind)
                    v=part_m[i]
                    #print(v)
                    for j in range(i):
                        if part_m[j][ind]==1:
                            part_m=plus(part_m,j,v)

            #print(part_m)

            #print(list(new_m))
            new_m = np.insert(part_m,M-zero_num,zero_m,axis=0 )
            #print(list(new_m))
            c=0
            for j in range(M-1):
                for i in range(M-1):
                    if (not j == i+1) and (new_m[j]==new_m[i+1]).all():
                        #print(new_m[j],new_m[i+1])
                        c+=1
                        #print('hello')

            b=np.all(new_m==0,axis=1)
            num=list(b).count(False)

            return num-c#np.sum(np.diag(new_m))#new_m#M-c #np.sum(np.diag(new_m))#

        N=self.n//(M*Q)
        split_key=list(split_list(self.key,M*Q))

        if len(split_key[-1]) != len(split_key[0]):
            split_key=split_key[0:-1]

        f_split_key= list(map(lambda x : np.reshape(np.array(list(map(lambda y : int(y), x))),[M,Q]), split_key))
        #print(rank(f_split_key[726],M,Q))

        ranks = list(map(lambda x: rank(x,M,Q),f_split_key))
        #print(ranks)

        full_rank = M

        length_ranks=len(ranks)

        FM =ranks.count(full_rank)
        FM1=ranks.count(full_rank-1)
        NFMM1=N - FM -FM1

        #print( FM, FM1, NFMM1)
        p = 1
        r = 32
        for i in range(r):
            p=p*(1-2**(i-Q))*(1-2**(i-M))/(1-2**(i-r))
        pM=2**(r*(Q+M-r)-M*Q)*p
        
        p = 1
        r = 31
        for i in range(r):
            p=p*(1-2**(i-Q))*(1-2**(i-M))/(1-2**(i-r))
        pM1=2**(r*(Q+M-r)-M*Q)*p

        pM2 = 1 - pM - pM1

        #print(pM,pM1,pM2)

        chi_squared_obs = (FM-pM*N)**2/(pM*N) + (FM1-pM1*N)**2/(pM1*N)+(NFMM1-pM2*N)**2/(pM2*N)
        p=e**(-chi_squared_obs/2)

        if p < 0.01:
            b = False
        
        else:
            b = True

        print('{:40} : {:.3f} -> {}'.format('binary matrix rank test',p,b))

        return p,b