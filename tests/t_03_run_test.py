import math

# .3 Runs Test
def run_test(key, ones, zeros, n, tau):
    pi = ones/n

    vobs=1
    for i in range(n-1):
        if key[i]==key[i+1]:
            vobs+=0
        else:
            vobs+=1

    p=math.erfc(abs(vobs-2*float(n)*pi*(1-pi))/(2*math.sqrt(2*float(n))*pi*(1-pi)))

    b = (p >= 0.01)

    print('{:40} : {:.3f} -> {}'.format('run test',p,b))
    return [p], b