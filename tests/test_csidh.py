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
        # TODO: Add svelu and p2048_CITDH in the future
        return super().setUp()

    # TODO: reduce num_sk in the future
    # TODO: test the randomness of sk
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

    # TODO: test the correctness of group action
    """
    NOTE: Currently, there's no reliable data... the output of sibc seems doubtful.
        I wrote a simple CSIDH version, which doesn't push points, just calculate coeffs. I think it is reliable.
        This implementation's output is the same as that simple version.
        But in the following example, original CTIDH's code has output different from both sibc's and mine..
    """
    # def test_group_action(self, num_sk=10):
    #     # def test_one_CSIDH_instance(csidh_instance, num_sk = num_sk):
    #     #     for _ in range(num_sk):
    #     #         sk = csidh_instance.skgen()
    #     #         print(f'sk = {sk}')
    #     # for instance in self.CSIDH_instances:
    #     #     test_one_CSIDH_instance(instance)
    #     ska = [-2, 0, 0, 0, -4, -1, -1, 1, 0, 1, -2, 1, 0, 2, 1, 2, 0, 0, 0, 2, 1, 0, 1, 0, 0, 4, 0, 3, 0, 0, 0, -1, 0, 2, -1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #     instance = self.CSIDH_instances[0] # tvelu p1024
    #     a = instance.group_action(0, ska)
    #     print(f'a = {hex(a)}')

    #     skb = [0, 0, 3, 0, -1, 0, 0, 2, 0, 0, 0, 0, 3, -1, 0, -2, -2, 1, 0, -1, 0, 1, 0, 2, 3, 0, 0, -1, -1, 0, 0, -1, 0, -1, -1, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #     b = instance.group_action(0, skb)
    #     print(f'b = {hex(b)}')

    #     # Is sibc's output correct?????
    #     # This result is very doubtful...
    #     # self.assertEqual(
    #     #     a, 
    #     #     0x947e58d0c56c30adf6edd66169bac9ec437e05477a07a67140fd0091796a792d8b282ef489f6f11d17e97a7833813918b79c6773590362eea54b79955568eda09d6f16b4004bb5ac3e2089c1a6e949b1c529b04edf1b4ea257b9a7a099d576074126dab11262c3469f76444c15a6a05264f954d319cad287975b988331bc964
    #     # )

    # TODO: reduce num_protocol in the future
    def test_protocol(self, num_protocols=50):
        def test_one_CSIDH_instance(csidh_instance, num_protocols = num_protocols):
            for _ in range(num_protocols):
                ska, pka = csidh_instance.keygen()
                # print(f'ska = {ska}')
                # print(f'pka = {pka}')                
                skb, pkb = csidh_instance.keygen()
                # print(f'skb = {skb}')
                # print(f'pkb = {pkb}')
                shared_a = csidh_instance.derive(ska, pkb)
                shared_b = csidh_instance.derive(skb, pka)
                # print(f'shared_a = {shared_a}')
                # print(f'shared_b = {shared_b}')
                self.assertEqual(shared_a, shared_b)

        for instance in self.CSIDH_instances:
            test_one_CSIDH_instance(instance)


    def tearDown(self) -> None:
        for ins in self.CSIDH_instances:
            del ins
        # del self.CSIDH_instances

        return super().tearDown()