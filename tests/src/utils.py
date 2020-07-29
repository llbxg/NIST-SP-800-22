import numpy as np

def open_rand(path, sixteen=False):
    with open(path) as f:
        l_strip = [s.strip() for s in f.readlines()]
    l_strip = ''.join(l_strip)

    if sixteen:
        l_strip=sixteen_to_two(l_strip)

    l_strip = list(l_strip)
    #l = np.array(l_strip)
    #f = np.frompyfunc(int, 1, 1)
    #return f(l)
    return l_strip

def sixteen_to_two(key16, _list=False):
    if not _list:
        key16 = list(key16)
    key=['0']*len(key16)
    for i in range(len(key16)):
        k=key16[i]
        ind = ['0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'].index(k)
        num = ['0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111']
        key[i]=num[ind]
    return ''.join(key)

def split_list(l, n):
    for idx in range(0, len(l), n):
        yield l[idx:idx + n]

def __print(b, txt):
    if b:
        print(txt)