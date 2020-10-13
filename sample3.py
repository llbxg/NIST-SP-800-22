import random
from multiprocessing import Process, Manager

from sp80022 import Sp80022
from proportion import proportion

def cal_rand(results):
    _data = [random.randint(0,1) for _ in range(10**6)]
    tst = Sp80022(_data, console=False)
    tst.run()
    results.append(tst.ps)

if __name__ == "__main__":
    with Manager() as manager:
        results = manager.list([])

        p_list = []
        for i in range(2):
            p = Process(target=cal_rand, args=(results,))
            p.start()
            p_list.append(p)

        for p in p_list:
            p.join()

        proportion(results)