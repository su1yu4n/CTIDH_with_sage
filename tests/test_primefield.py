import unittest
import json

from sage.all import kronecker_symbol, proof

from CTIDH import PrimeField
from CTIDH.utils import read_prime_info, get_randint

proof.arithmetic(False)

p1024_info = read_prime_info('p1024_CTIDH')
p2048_info = read_prime_info('p2048_CTIDH')

p1024 = p1024_info["p"]
p2048 = p2048_info["p"]

Fp1024 = PrimeField(p1024)
Fp2048 = PrimeField(p2048)


class TestPrimeField(unittest.TestCase):
    def test_init(self):
        a = Fp1024(3)
        self.assertEqual(a.value, 3)
        b = Fp1024(p1024)
        self.assertEqual(b.value, 0)

    def test_eq(self):
        a = Fp1024(3)
        self.assertEqual(a, 3)
        b = Fp1024(p1024)
        self.assertEqual(b, 0)

    def test_add(self):
        Fp1024.reset_runtime()
        a = Fp1024(3)
        b = Fp1024(4)
        self.assertEqual(Fp1024.add_count, 0)
        c = a + b
        self.assertEqual(c, 7)
        self.assertEqual(Fp1024.add_count, 1)
        c = c + 3
        self.assertEqual(c, 10)
        self.assertEqual(Fp1024.add_count, 2)

        Fp2048.reset_runtime()
        a = Fp2048(39)
        b = Fp2048(41)
        self.assertEqual(Fp2048.add_count, 0)
        c = a + b
        self.assertEqual(Fp2048.add_count, 1)
        self.assertEqual(c, 80)
        c = c + 3
        self.assertEqual(Fp2048.add_count, 2)
        self.assertEqual(c, 83)

    def test_radd(self):
        Fp1024.reset_runtime()
        a = Fp1024(3)
        b = Fp1024(4)
        self.assertEqual(Fp1024.add_count, 0)
        c = 3 + a
        self.assertEqual(c, 6)
        self.assertEqual(Fp1024.add_count, 1)
        c = 3 + b
        self.assertEqual(c, 7)
        self.assertEqual(Fp1024.add_count, 2)

    def test_iadd(self):
        Fp1024.reset_runtime()
        a = Fp1024(7)
        b = Fp1024(8)
        self.assertEqual(Fp1024.add_count, 0)
        a += 1
        self.assertEqual(a, 8)
        self.assertEqual(a, b)
        self.assertEqual(Fp1024.add_count, 1)
        b += a
        self.assertEqual(b, 16)
        self.assertEqual(Fp1024.add_count, 2)

    def test_sub(self):
        Fp1024.reset_runtime()
        a = Fp1024(16)
        b = Fp1024(15)
        self.assertEqual(Fp1024.add_count, 0)
        c = a - b
        self.assertEqual(c, 1)
        self.assertEqual(Fp1024.add_count, 1)
        a = a - b
        self.assertEqual(a, 1)
        self.assertEqual(Fp1024.add_count, 2)
        a = a - 1
        self.assertEqual(a, 0)
        self.assertEqual(Fp1024.add_count, 3)

    def test_rsub(self):
        Fp1024.reset_runtime()
        a = Fp1024(p1024)
        b = Fp1024(p1024 - 2)
        self.assertEqual(Fp1024.add_count, 0)
        c = p1024 - b
        self.assertEqual(Fp1024.add_count, 1)
        self.assertEqual(c, 2)
        d = 3 - a
        self.assertEqual(Fp1024.add_count, 2)
        self.assertEqual(d, 3)
        e = 0 - b
        self.assertEqual(Fp1024.add_count, 3)
        self.assertEqual(e, 2)

    def test_isub(self):
        Fp2048.reset_runtime()
        a = Fp2048(90)
        b = Fp2048(42)
        self.assertEqual(Fp2048.add_count, 0)
        a -= 3
        self.assertEqual(a, 87)
        self.assertEqual(Fp2048.add_count, 1)
        a -= b
        self.assertEqual(a, 45)
        self.assertEqual(Fp2048.add_count, 2)

    def test_mul(self):
        Fp2048.reset_runtime()
        a = Fp2048(3)
        b = Fp2048(7)
        self.assertEqual(Fp2048.mul_count, 0)
        c = a * b
        self.assertEqual(c, 21)
        self.assertEqual(Fp2048.mul_count, 1)
        c = p2048 * c
        self.assertEqual(c, 0)
        self.assertEqual(Fp2048.mul_count, 2)
        b = b * (p2048 - 1)
        self.assertEqual(b, -7)
        self.assertEqual(Fp2048.mul_count, 3)

    def test_rmul(self):
        Fp1024.reset_runtime()
        a = Fp1024(p1024-1)
        b = Fp1024(5)
        self.assertEqual(Fp1024.mul_count, 0)
        c = 8*a
        self.assertEqual(c, -8)
        self.assertEqual(Fp1024.mul_count, 1)
        c = a*b
        self.assertEqual(c, -5)
        self.assertEqual(Fp1024.mul_count, 2)
        c = a*(-b)
        self.assertEqual(c, 5)
        self.assertEqual(Fp1024.mul_count, 3)

    def test_imul(self):
        Fp1024.reset_runtime()
        a = Fp1024(p1024-1)
        b = Fp1024(5)
        self.assertEqual(Fp1024.mul_count, 0)
        a *= 3
        self.assertEqual(a, -3)
        self.assertEqual(Fp1024.mul_count, 1)
        b *= a
        self.assertEqual(b, -15)
        self.assertEqual(Fp1024.mul_count, 2)
    
    def test_pow(self):
        Fp2048.reset_runtime()
        Fp2048.reset_power_invert_time()
        a = Fp2048(p2048-2)
        b = 8
        c = 7
        d = 233
        self.assertEqual(Fp2048.mul_count, 0)
        self.assertEqual(Fp2048.sqr_count, 0)
        self.assertEqual(Fp2048.pow_count, 0)

        
        e = a**8  # a -> a^2 -> a^4 -> a^8
        self.assertEqual(e, 256)
        self.assertEqual(Fp2048.mul_count, 0)
        self.assertEqual(Fp2048.sqr_count, 3)
        self.assertEqual(Fp2048.pow_count, 1)
        Fp2048.reset_runtime() 

        e = a**c
        self.assertEqual(Fp2048.mul_count, 2)
        self.assertEqual(Fp2048.sqr_count, 2)
        # reset_runtime won't change pow_count
        self.assertEqual(Fp2048.pow_count, 2)
        Fp2048.reset_runtime() 

        e = a**2
        self.assertEqual(Fp2048.mul_count, 0)
        self.assertEqual(Fp2048.sqr_count, 1)
        self.assertEqual(Fp2048.pow_count, 2) # a^2 is not considered as a power, because it only need a square

        
    def test_ipow(self):
        #TODO
        pass


    def test_safe_pow(self, num_test=100):
        for Fp, p in [(Fp1024, p1024), (Fp2048, p2048)]:
            for _ in range(num_test):
                Fp.reset_runtime()
                Fp.reset_power_invert_time()
                a = Fp(get_randint(1, p-1))
                e = get_randint(3, 2**12)
                e_max = get_randint(3, 2**12)
                if e_max < e:
                    e_max, e = e, e_max
                
                e_maxbitlen = e_max.bit_length()
                result = a.safe_pow(e, e_maxbitlen)
                self.assertEqual(Fp.pow_count, 1)
                self.assertEqual(Fp.add_count, 0)
                self.assertEqual(Fp.mul_count, e_maxbitlen-1)
                self.assertEqual(Fp.sqr_count, e_maxbitlen-1)

                self.assertEqual(result, a**e)


    def test_invert(self):
        #TODO
        pass

    def test_random(self, num_test = 50):
        for _ in range(num_test):
            a = Fp1024.get_random()
            b = Fp1024.get_random()
            self.assertNotEqual(a, b)

    # num_of_test used to be 10**4, and have passed the test several times.
    # Change to 100 to make the test faster
    def test_is_square(self, num_test = 50):
        for _ in range(num_test):
            a = Fp2048.get_random()
            self.assertEqual(a.is_square(), kronecker_symbol(a.value, p2048) == 1)

        for _ in range(num_test):
            a = Fp1024.get_random()
            self.assertEqual(a.is_square(), kronecker_symbol(a.value, p1024) == 1)
