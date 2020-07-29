import scipy.special as sc

# .14 Random Excursions Test
def random_excursions_test(key, n):
    S=[0]*n

    for i in range(n):
        c = key[i]
        if   i == 0 : S[i]=c
        elif c == 1 : S[i] = S[i-1]+1
        else        : S[i] = S[i-1]-1

    S_dash =[0] + S + [0]

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

    ind_list=index_multi(S_dash, 0, n)
    J = len(ind_list) - 1

    S_dash= list(map(lambda x : str(x+5),S_dash))

    cyc_list=[]

    for i in range(len(ind_list)):
        if not ((i+1) >= len(ind_list)):
            start = ind_list[i]
            fin   = ind_list[i+1]

            cyc = S_dash[start:fin+1]
            cyc_list.append(cyc)

    cycle=[]
    for i in range(len(cyc_list)):
        v=[]
        for j in [1,2,3,4,6,7,8,9]:
            v.append(cyc_list[i].count(str(j)))
        cycle.append(v)

    states =[]
    for i in range(8):
        state=[]
        for j in range(len(cycle)):
            state.append(cycle[j][i])
        states.append(state)

    vs=[]
    for i in states:
        vv=[0,0,0,0,0,0]
        vv[0]=len(list(filter(lambda x : x==0, i)))
        vv[1]=len(list(filter(lambda x : x==1, i)))
        vv[2]=len(list(filter(lambda x : x==2, i)))
        vv[3]=len(list(filter(lambda x : x==3, i)))
        vv[4]=len(list(filter(lambda x : x==4, i)))
        vv[5]=len(list(filter(lambda x : x>=5, i)))

        vs.append(vv)

    pis=[
        [0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000],
        [0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625],
        [0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143],
        [0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051]
        ]

    ps, bs = [], []
    for cnt, i in enumerate([3,2,1,0,0,1,2,3]):
        c = sum(list(map(lambda x,y : (x-J*y)**2/(J*y),vs[cnt],pis[i])))
        p = sc.gammaincc(5/2,c/2)
        b = (p >= 0.01)

        ps.append(p)
        bs.append(b)

        st = [-4,-3,-2,-1,1,2,3,4]

        if   cnt == 0 : print('{:40} : {:.3f} -> {} '.format('random excursions test (x = {:2})'.format(st[cnt]),p,b))
        else          : print('{:40} : {:.3f} -> {} '.format('                       (x = {:2})'.format(st[cnt]),p,b))
        cnt+=1

    return ps, all(bs)