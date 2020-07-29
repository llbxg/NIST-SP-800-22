import numpy as np

from tests.src.utils import open_rand

from sp80022 import Sp80022

path = 'data/rand.txt'
_data=open_rand(path, sixteen=True)
_data=_data[:1000000]
_data = list(map(int, _data))
tst1 = Sp80022(_data)
s = tst1.run()
print(s)

# I use the example (e - Napier's constant) distributed by Nist's program.
path = 'data/e.txt'
_data=open_rand(path, sixteen=False)
_data=_data[:1000000]
_data = list(map(int, _data))
tst2 = Sp80022(_data)
s = tst2.run()
print(s)