import unittest

import numpy as np

from CTIDH.csidh import CSIDH

class TestCSIDH(unittest.TestCase):
    """Showcase for a test with setUp and tearDown methods"""

    CSIDH_instances = []

    def setUp(self) -> None:
        CSIDH_1024_tvelu_ladder_original = CSIDH(
            'p1024_CTIDH', 
            'tvelu',
        )
        self.CSIDH_instances.append(CSIDH_1024_tvelu_ladder_original)

    # NOTE: decrease num_sk in the future
    def test_skgen(self, num_sk=100):
        def test_one_CSIDH_instance(csidh_instance, num_sk = num_sk):
            batch_bound = csidh_instance.prime_info['batch_bound']
            batch_start = csidh_instance.prime_info['batch_start']
            batch_stop = csidh_instance.prime_info['batch_stop']
            batch_num = len(batch_start)
            # num_total_positive = 0
            # num_total_negative = 0
            for _ in range(num_sk):
                sk = csidh_instance.skgen()
                self.assertEqual(len(sk), batch_stop[-1])

                for i in range(batch_num):
                    eis = np.array(sk[batch_start[i]: batch_stop[i]])
                    self.assertLessEqual(
                        np.sum(np.abs(eis)), 
                        batch_bound[i], 
                        f'batch {i}, sum of sk abs greater than batchbound = {batch_bound[i]}, eis = {eis}'
                    ) 
            
        for instance in self.CSIDH_instances:
            test_one_CSIDH_instance(instance)