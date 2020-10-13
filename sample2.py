import random
from concurrent.futures import ThreadPoolExecutor

from sp80022 import Sp80022
from proportion import proportion

def cal_rand():
    _data = [random.randint(0,1) for _ in range(10**6)]
    tst = Sp80022(_data, console=False)
    tst.run()
    return tst.ps


with ThreadPoolExecutor(max_workers=4, thread_name_prefix="thread") as executor:
    futures = []
    for i in range(4):
        futures.append(executor.submit(cal_rand))

    results = [f.result() for f in futures]

proportion(results)