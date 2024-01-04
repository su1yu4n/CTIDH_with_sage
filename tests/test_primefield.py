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
    
    # TODO: Add more tests, such as sub rsub isub mul pow invert...



# p = 31258613193250572484686795204590428045325278335756488509417380147989191591013695960523454229439943526174516836804526694803754791806514797454806413422135845959682240605147195662178033945501227284906138519332665421238173115419904786072611649808082328359399330601814276219358255285323504274244891804869091065552318327790140350673759394401692394790149771353631859923204704775535338816817699239207701424948322420803517506192052699284149495781764810879921006421214731560197045604172518737388974607923865313337968911439410704459583975801501147132536789055845508806054526511545228687644365106576314242096626345373821692886509
# Fp = PrimeField(p)
# Fp.reset_runtime()
# a = Fp(3)
# b = Fp(4)
# Fp.show_runtime()
# c = a + b
# Fp.show_runtime()
# c = c + 3
# Fp.show_runtime()
# c = 3 + c
# Fp.show_runtime()
# c += a
# Fp.show_runtime()


# Fp.reset_runtime()
# a -= 5
# Fp.show_runtime()
# a = a-5
# Fp.show_runtime()
# a = 5-a
# Fp.show_runtime()
# a = -a


# Fp.reset_runtime()
# c = a*b
# Fp.show_runtime()
# c = a*3
# Fp.show_runtime()
# c = 3*a
# Fp.show_runtime()
