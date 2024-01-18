import unittest
import json

from CTIDH import PrimeField


with open("parameters/p1024_CTIDH", "r") as f:
    p1024_info = json.load(f)
with open("parameters/p2048_CTIDH", "r") as f:
    p2048_info = json.load(f)
p1024 = p1024_info["p"]
p2048 = p2048_info["p"]

Fp1024 = PrimeField(p1024)
Fp2048 = PrimeField(p2048)


class TestPrimeField(unittest.TestCase):
    def test_init(self):
        a = Fp1024(3)
        self.assertEqual(a.x, 3)
        b = Fp1024(p1024)
        self.assertEqual(b.x, 0)

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
        a = Fp2048(p2048-2)
        b = 8
        c = 7
        d = 233
        self.assertEqual(Fp2048.mul_count, 0)
        
        e = a**8  # a -> a^2 -> a^4 -> a^8
        self.assertEqual(e, 256)
        self.assertEqual(Fp2048.mul_count, 0)
        self.assertEqual(Fp2048.sqr_count, 3)
        self.assertEqual(Fp2048.pow_count, 1)
        Fp2048.reset_runtime()

        e = a**c
        self.assertEqual(Fp2048.mul_count, 2)
        self.assertEqual(Fp2048.sqr_count, 2)
        self.assertEqual(Fp2048.pow_count, 2)
        
    def test_ipow(self):
        pass

    def test_invert(self):
        pass

    # TODO: Add tests of pow and invert after implement them
    