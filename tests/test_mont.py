import unittest
import json

from sage.all import EllipticCurve, proof, GF

from CTIDH import MontgomeryCurve, PrimeField

proof.arithmetic(False)


with open("parameters/p1024_CTIDH", "r") as f:
    p1024_info = json.load(f)
with open("parameters/p2048_CTIDH", "r") as f:
    p2048_info = json.load(f)
p1024 = p1024_info["p"]
p2048 = p2048_info["p"]

Fp1024 = PrimeField(p1024)
Fp2048 = PrimeField(p2048)

# sage_GFp1024 = GF(p1024)
# sage_GFp2048 = GF(p2048)

# SageFp = type(sage_GFp1024)
# SageEC = type(EllipticCurve(sage_GFp1024, [0, 0, 0, 1, 0]))

MontCurve_p1024 = MontgomeryCurve("p1024_CTIDH")
MontCurve_p2048 = MontgomeryCurve("p2048_CTIDH")


# def get_sage_montgomery_curve(Fp: int, a: int) -> SageEC:
#     return EllipticCurve(Fp, [0, a, 0, 1, 0])

# TODO: Add tests of montgomery curve arithmetics.
class TestMontgomeryCurve(unittest.TestCase):
    def test_elligator(self, num_of_test=1000):
        for field, MontCurve in [(Fp1024, MontCurve_p1024), (Fp2048, MontCurve_p2048)]:
            field.reset_runtime(); field.reset_power_invert_time()

            a = field(0)
            A = (a, field(1)) # 
            T0, T1 = MontCurve.elligator(A)


            T0_x = T0[0] * T0[1]**(-1)
            self.assertEqual(field.is_square(T0_x ** 3 + a * T0_x**2 + T0_x), True)
            T1_x = T1[0] * T1[1]**(-1)
            self.assertEqual(field.is_square(T1_x ** 3 + a * T1_x**2 + T1_x), False)

            for _ in range(num_of_test - 1):
                a = field.get_random()
                # Check if Ea is supersingular.
                # It suffices to check a+2, a-2 are non-zero
                if a == 2 or a == -2:
                    continue
                A = (a, field(1)) # 
                T0, T1 = MontCurve.elligator(A)
                # print(f'a = {a}')
                # print(f'T0 = {T0}')
                # print(f'T1 = {T1}')
                
                T0_x = T0[0] * T0[1]**(-1)
                self.assertEqual(field.is_square(T0_x ** 3 + a * T0_x**2 + T0_x), True)
                T1_x = T1[0] * T1[1]**(-1)
                self.assertEqual(field.is_square(T1_x ** 3 + a * T1_x**2 + T1_x), False)
            
            field.show_runtime(label='{} {} elligators'.format(field, num_of_test))
        
