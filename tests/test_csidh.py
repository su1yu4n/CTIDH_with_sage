import unittest

import numpy as np
from sage.all import EllipticCurve, proof, is_prime, GF

from CTIDH.csidh import CSIDH
    

proof.arithmetic(False)


def isogeny_sage(E, p, l, inverse_action):
    def sample_l_order_point(E, l):
        P = E(0)
        while P == E(0):
            R = E.random_point()
            order_E = p + 1
            P = (order_E // l) * R

        assert l * P == E(0) and P != E(0)
        return P
    if inverse_action == True:
        a = E.a2()
        E = EllipticCurve(E.base_field(), [0, -a, 0, 1, 0])
    P = sample_l_order_point(E, l)
    if l > 200:
        phi = E.isogeny(kernel=P, model="montgomery", algorithm="velusqrt")
    else:
        phi = E.isogeny(kernel=P, model="montgomery")
    E_new = phi.codomain()

    if inverse_action == True:
        a_new = E_new.a2()
        E_new = EllipticCurve(E.base_field(), [0, -a_new, 0, 1, 0])

    return E_new


def group_action_sage(sk: list, L: list, E, p):
    def do_l_e_action(l, e, E):
        if e == 0:
            return E

        inverse_action = True if e < 0 else False
        e = abs(e)
        for _ in range(e):
            E = isogeny_sage(E, p, l, inverse_action)
        return E

    assert len(sk) == len(L)
    E_a = E
    for i in range(len(sk)):
        e_i = sk[i]
        l_i = L[i]
        E_a = do_l_e_action(l_i, e_i, E_a)

    return E_a

def pkgen_sage(sk: list, prime_info: dict):
    p = prime_info["p"]
    L = prime_info["L"]
    assert is_prime(p)
    E0 = EllipticCurve(GF(p), [0, 0, 0, 1, 0])
    pk = group_action_sage(sk, L, E0, p)
    return pk.a2()


class TestCSIDH(unittest.TestCase):
    CSIDH_instances = []

    def setUp(self) -> None:
        # CSIDH_1024_tvelu_ladder_original = CSIDH(
        #     'p1024_CTIDH', 
        #     'tvelu',
        # )        
        # self.CSIDH_instances.append(CSIDH_1024_tvelu_ladder_original)

        # CSIDH_1024_tvelu_SDAC_original_slow_legendre = CSIDH(
        #     'p1024_CTIDH', 
        #     'tvelu',
        #     SDAC=True
        # )
        # self.CSIDH_instances.append(CSIDH_1024_tvelu_SDAC_original_slow_legendre)
        
        CSIDH_1024_tvelu_SDAC_original_fast_kronecker = CSIDH(
            'p1024_CTIDH', 
            'tvelu',
            SDAC=True,
            fast_kronecker=True
        )
        self.CSIDH_instances.append(CSIDH_1024_tvelu_SDAC_original_fast_kronecker)

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


    def test_group_action(self, num_sk=10):
        def test_one_CSIDH_instance(csidh_instance, num_sk = num_sk):
            prime_info = csidh_instance.prime_info
            p = prime_info["p"]
            L = prime_info["L"]

            def test_one_action(a=0):
                print(f'a={a}')
                sk = csidh_instance.skgen()
                print(f'sk = {sk}')
                anew = csidh_instance.group_action(a, sk, debug=False)

                Ea = EllipticCurve(GF(p), [0, a, 0, 1, 0])
                Enew_sage = group_action_sage(sk, L, Ea, p)
                anew_sage = Enew_sage.a2()
                print(f'anew={anew}')
                print(f'anew_sage={anew_sage}\n')
                self.assertEqual(anew, anew_sage)
                return anew

            for _ in range(num_sk // 2):
                anew = test_one_action(a=0)
                test_one_action(a=anew)
                
                
        for instance in self.CSIDH_instances:
            test_one_CSIDH_instance(instance)


    def test_protocol(self, num_protocols=10):
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