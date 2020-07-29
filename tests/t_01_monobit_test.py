import math

# .1 Frequency (Monobit) Test
def monobit_test(ones, zeros, n):
    sn=abs(ones-zeros)
    sobs=sn/math.sqrt(n)
    p=math.erfc(sobs/math.sqrt(2))

    b = (p >= 0.01)

    print('{:40} : {:.3f} -> {}'.format('monobit test',p,b))
    return [p], b