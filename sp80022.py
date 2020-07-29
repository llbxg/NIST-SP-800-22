import numpy as np

from tests import monobit_test, frequency_test_with_block, run_test, test_for_the_longest_run_of_ones_in_a_block, binary_matrix_rank_test, non_overlapping_template_matching_test, overlapping_template_matching_test, maurers_universal_statical_test, liner_complexity_test
from tests import serial_test, approximate_entropy_test, cumulative_sums_test, random_excursions_test, random_excursion_variant_test

import math
#from functools import partial

class Sp80022():

    def __init__(self,key16,sixteen=True, console=False):
        self.key = key16
        self.zeros=self.key.count(0)
        self.ones=self.key.count(1)
        self.n = len(self.key)
        self.tau=2/math.sqrt(self.n)
        self.console = console

    def __repr__(self):
        return "nist sp800-22"

    def run(self):
        bs = []
        print('hello ! this is nist sp800-22 program')
        print('---------------------------------------------------------------------------------------------------------------')
        print('total-bit : {} / 0 : {} / 1 : {}'.format(self.n, self.ones, self.zeros))
        print('---------------------------------------------------------------------------------------------------------------')

        p,b = monobit_test(self.ones, self.zeros, self.n)
        bs.append(b)

        p,b = frequency_test_with_block(self.key, self.n, M=10101)
        bs.append(b)

        p,b = run_test(self.key, self.ones, self.zeros, self.n, self.tau)
        bs.append(b)

        p,b = test_for_the_longest_run_of_ones_in_a_block(self.key, self.n)
        bs.append(b)

        p,b = binary_matrix_rank_test(self.key, self.n)
        bs.append(b)

        p,b = non_overlapping_template_matching_test(self.key, self.n)
        bs.append(b)

        p,b = overlapping_template_matching_test(self.key, self.n)
        bs.append(b)

        p,b = maurers_universal_statical_test(self.key, self.n)
        bs.append(b)

        p,b = liner_complexity_test(self.key, self.n)
        bs.append(b)

        p,b = serial_test(self.key, self.n)
        bs.append(b)

        p,b = approximate_entropy_test(self.key, self.n)
        bs.append(b)

        p,b = cumulative_sums_test(self.key,self.n)
        bs.append(b)

        p,b = random_excursions_test(self.key, self.n)
        bs.append(b)

        p,b = random_excursion_variant_test(self.key, self.n)
        bs.append(b)

        print('---------------------------------------------------------------------------------------------------------------')
        return all(bs)