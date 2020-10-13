import math
from tests.src.utils import __print

# .1 Frequency (Monobit) Test
def monobit_test(ones, zeros, n, b_print=True):
    sn=abs(ones-zeros)
    sobs=sn/math.sqrt(n)
    p=math.erfc(sobs/math.sqrt(2))

    b = (p >= 0.01)

    __print(b_print, '{:40} : {:.3f} -> {}'.format('monobit test',p,b))
    return [p], b