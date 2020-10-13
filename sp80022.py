import numpy as np

from tests import monobit_test, frequency_test_with_block, run_test, test_for_the_longest_run_of_ones_in_a_block, binary_matrix_rank_test, non_overlapping_template_matching_test, overlapping_template_matching_test, maurers_universal_statical_test, liner_complexity_test
from tests import serial_test, approximate_entropy_test, cumulative_sums_test, random_excursions_test, random_excursion_variant_test

import math
#from functools import partial

class Sp80022():

    def __init__(self,key16,sixteen=True, console=True):
        self.key = key16
        self.zeros=self.key.count(0)
        self.ones=self.key.count(1)
        self.n = len(self.key)
        self.tau=2/math.sqrt(self.n)
        self.console = console

        self.ps=None

    def __repr__(self):
        return "nist sp800-22"

    def __print(self, msg):
        if self.console:
            print(msg)


    def run(self):
        ps, bs = [], []
        self.__print('hello ! this is nist sp800-22 program')
        self.__print('---------------------------------------------------------------------------------------------------------------')
        self.__print(f'total-bit : {self.n} / 0 : {self.ones} / 1 : {self.zeros}')
        self.__print('---------------------------------------------------------------------------------------------------------------')

        p,b = monobit_test(self.ones, self.zeros, self.n, b_print=self.console)
        ps.extend(p)
        bs.append(b)

        p,b = frequency_test_with_block(self.key, self.n, M=10101, b_print=self.console)
        ps.extend(p)
        bs.append(b)

        p,b = run_test(self.key, self.ones, self.zeros, self.n, self.tau, b_print=self.console)
        ps.extend(p)
        bs.append(b)

        f_list = [
            test_for_the_longest_run_of_ones_in_a_block,
            binary_matrix_rank_test,
            non_overlapping_template_matching_test,
            overlapping_template_matching_test,
            maurers_universal_statical_test,
            liner_complexity_test,
            serial_test,
            approximate_entropy_test,
            cumulative_sums_test,
            random_excursions_test,
            random_excursion_variant_test
            ]

        for f in f_list:
            p,b = f(self.key, self.n, b_print=self.console)
            ps.extend(p)
            bs.append(b)
        self.ps = ps

        self.__print('---------------------------------------------------------------------------------------------------------------')
        return all(bs)