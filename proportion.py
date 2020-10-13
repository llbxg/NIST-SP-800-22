# Proportion of Sequences Passing a Test

import numpy as np

def proportion(dts, alpha = 0.01):

    l = len(dts)

    p_hat = 1 - alpha
    range_p_min = round(p_hat - (b := np.sqrt((p_hat*(1-p_hat))/l)*3), 3)
    range_p_max = round(p_hat + b, 3)

    sum_p = np.zeros(40)
    for dt in dts:
        sum_p += np.array([1 if(p>=alpha) else 0 for p in dt])

    result = [ True if s>=range_p_min and s<= range_p_max else False for s in sum_p/l]

    msg = f'Congratulations on passing the test.' if all(result) else 'Failed to pass the test.'
    print(f"{msg} : {result.count(True)}/{len(result)} with {l} datas")

    return result