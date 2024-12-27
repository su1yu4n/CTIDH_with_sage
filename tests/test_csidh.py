import unittest

import numpy as np
from sage.all import EllipticCurve, proof, is_prime, GF, binomial

from CTIDH.csidh import CSIDH
    

proof.arithmetic(False)


def isogeny_sage(E, p, l, inverse_action):
    def sample_l_order_point(E, l):
        P = E(0)
        while P == E(0):
            R = E.random_point()
            P = (E.order() // l) * R
        # assert l * P == E(0) and P != E(0)
        return P
    
    if inverse_action:
        a = E.a2()
        E = EllipticCurve(E.base_field(), [0, -a, 0, 1, 0])
    E.set_order(p+1)
    P = sample_l_order_point(E, l)
    P.set_order(l)
    if l > 200:
        phi = E.isogeny(kernel=P, model="montgomery", algorithm="velusqrt", check=False)
    else:
        phi = E.isogeny(kernel=P, model="montgomery", degree=l, check=False)
    E_new = phi.codomain()

    if inverse_action:
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
    E0.set_order(p+1)
    pk = group_action_sage(sk, L, E0, p)
    return pk.a2()


class TestCSIDH(unittest.TestCase):
    CSIDH_instances = []

    def setUp(self) -> None:
        CSIDH_1024_tvelu_SDAC_original_fast_kronecker = CSIDH(
            'p1024_CTIDH', 
            'tvelu',
            SDAC=True,
            fast_kronecker=True
        )
        self.CSIDH_instances.append(CSIDH_1024_tvelu_SDAC_original_fast_kronecker)
        
        CSIDH_1024_tvelu_SDAC_original_slow_legendre= CSIDH(
            'p1024_CTIDH', 
            'tvelu',
            SDAC=True,
            fast_kronecker=False
        )
        self.CSIDH_instances.append(CSIDH_1024_tvelu_SDAC_original_slow_legendre)
        
        CSIDH_2048_tvelu_SDAC_original_fast_kronecker = CSIDH(
            'p2048_CTIDH', 
            'tvelu',
            SDAC=True,
            fast_kronecker=True
        )
        self.CSIDH_instances.append(CSIDH_2048_tvelu_SDAC_original_fast_kronecker)
        
        CSIDH_2048_tvelu_SDAC_original_slow_legendre = CSIDH(
            'p2048_CTIDH', 
            'tvelu',
            SDAC=True,
            fast_kronecker=False
        )
        self.CSIDH_instances.append(CSIDH_2048_tvelu_SDAC_original_slow_legendre)
        


        # CSIDH_1024_hvelu_SDAC_original_slow_legendre = CSIDH( # in fact this slow legendre is faster in this python implementation
        #     'p1024_CTIDH', 
        #     'hvelu',
        #     SDAC=True,
        #     fast_kronecker=False
        # )
        # self.CSIDH_instances.append(CSIDH_1024_hvelu_SDAC_original_slow_legendre)

        # CSIDH_2048_hvelu_SDAC_original_slow_legendre = CSIDH(
        #     'p2048_CTIDH', 
        #     'hvelu',
        #     SDAC=True,
        #     fast_kronecker=False
        # )
        # self.CSIDH_instances.append(CSIDH_2048_hvelu_SDAC_original_slow_legendre)
        return super().setUp()

    # NOTE: Randomness of skgen tested in test_randomboundl1
    def test_skgen(self, num_sk=20):
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

           
    # NOTE: this takes about 300s on my laptop
    def test_random_boundedl1(self):
        def phi_x_y(x, y) -> int:
            up = min(x, y)
            result = 0
            for k in range(0, up+1):
                result += binomial(x, k) * 2**k * binomial(y, k)
            return result

        for csidh_instance in self.CSIDH_instances:
            batch_bound = csidh_instance.prime_info['batch_bound']
            batch_start = csidh_instance.prime_info['batch_start']
            batch_stop = csidh_instance.prime_info['batch_stop']

            def test_one_batch(b):
                
                Ni = batch_stop[b] - batch_start[b]
                mi = batch_bound[b]
                keyspace = phi_x_y(Ni, mi)
                # print(f'\nbatch number: {b}')
                # print(f'Ni = {Ni}')
                # print(f'mi = {mi}')
                # print(f'keyspace = {keyspace}')
                if keyspace > 5000: # tooo slow
                    return

                avg_num_per_key = 200 if keyspace < 1000 else 100
                bias = 0.25*avg_num_per_key if keyspace < 1000 else 0.35*avg_num_per_key
                num_ei = avg_num_per_key * keyspace
                
                # print(f'num_ei = {num_ei}')
                ei_stat = {}

                num_out_of_bias = 0
                for _ in range(num_ei):
                    key_str = str(csidh_instance.random_boundedl1(Ni, mi))
                    if key_str not in ei_stat:
                        ei_stat[key_str] = 1
                    else:
                        ei_stat[key_str] += 1
                for _, num_key in ei_stat.items():
                    if abs(num_key - avg_num_per_key) > bias:
                        num_out_of_bias += 1
                
                # print(f'num_out_of_bias = {num_out_of_bias}')
                out_bias_ratio = 0.02 if keyspace < 500 else 0.05
                max_num_out_of_bias = out_bias_ratio * keyspace
                # This should only fail with a negligible probablity
                self.assertLessEqual(
                    num_out_of_bias, 
                    max_num_out_of_bias,
                    msg=f"Too many keys' number out of bias, batch {b}, cardinality of keyspace = {keyspace}"
                )
            
            for b in range(0, 5):
                test_one_batch(b)

    # NOTE: I have run it with num_sk = 50 and that took more than 2h
    def test_group_action(self, num_sk=5, verbose_level=2):
        def test_one_CSIDH_instance(csidh_instance, num_sk = num_sk, name='1024'):
            prime_info = csidh_instance.prime_info
            p = prime_info["p"]
            L = prime_info["L"]

            def test_one_action(a=0, name='1024'):
                # print(f'a={a}')
                sk = csidh_instance.skgen()
                # print(f'sk = {sk}')
                anew = csidh_instance.group_action(a, sk, verbose_level)

                csidh_instance.field.show_runtime("CTIDH-" + name + " GA")
                csidh_instance.field.show_power_invert_time("CTIDH-" + name + " GA")
                csidh_instance.field.reset_runtime()
                csidh_instance.field.reset_power_invert_time()

                Ea = EllipticCurve(GF(p), [0, a, 0, 1, 0])
                Enew_sage = group_action_sage(sk, L, Ea, p)
                anew_sage = Enew_sage.a2()
                # print(f'anew={anew}')
                # print(f'anew_sage={anew_sage}\n')
                self.assertEqual(anew, anew_sage)
                return anew

            for i in range(num_sk):
                print(f'Running action {i}')
                anew = test_one_action(a=0, name=name)
                test_one_action(a=anew, name=name)
                
                
        # for i in range(len(self.CSIDH_instances)):
        #     print(f"Testing instance {i}")
        #     test_one_CSIDH_instance(self.CSIDH_instances[i])
        
        # NOTE: test CTIDH-1024 only
        test_one_CSIDH_instance(self.CSIDH_instances[0]) # fast kro
        # test_one_CSIDH_instance(self.CSIDH_instances[1]) # slow legendre
        # test_one_CSIDH_instance(self.CSIDH_instances[-1], name='2048') # 2048 slow legendre
        # test_one_CSIDH_instance(self.CSIDH_instances[-2], name='2048') # 2048 fast kronecker
        
        


    def test_protocol(self, num_protocols=10):
        def test_one_CSIDH_instance(csidh_instance, num_protocols = num_protocols):
            for i in range(num_protocols):
                print(f'Running protocol {i}')
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

        for i in range(len(self.CSIDH_instances)):
            print(f"Testing instance {i}")
            test_one_CSIDH_instance(self.CSIDH_instances[i])


    def tearDown(self) -> None:
        for ins in self.CSIDH_instances:
            del ins
        # del self.CSIDH_instances

        return super().tearDown()