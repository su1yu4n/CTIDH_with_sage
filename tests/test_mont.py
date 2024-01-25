import unittest
import json

from sage.all import EllipticCurve, proof, GF
from sage.all import *

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

sage_GFp1024 = GF(p1024)
sage_GFp2048 = GF(p2048)



MontCurve_p1024 = MontgomeryCurve("p1024_CTIDH")
MontCurve_p2048 = MontgomeryCurve("p2048_CTIDH")


def get_sage_montgomery_curve(sage_Fp, a: int):
    return EllipticCurve(sage_Fp, [0, a, 0, 1, 0])

def get_affine_from_projective(A: list) -> int:
    """ Given A = (Ax: Az), Az!= 0, compute a = Ax * Az^(-1)

    Args:
        A (list): [Ax, Az] that represents (Ax:Az), where Ax, Az are ZModPrime
    """
    Ax, Az = A[0], A[1]
    assert len(A) == 2
    assert Az != 0

    a = Ax*Az**(-1)
    return int(a.get_int_value())


# TODO: Add tests of montgomery curve arithmetics.
class TestMontgomeryCurve(unittest.TestCase):
    def test_elligator(self, num_of_test=500):
        for field, MontCurve, prime_name in [
            (Fp1024, MontCurve_p1024, 'p1024_CTIDH'), 
            (Fp2048, MontCurve_p2048, 'p2048_CTIDH')
        ]:

            field.reset_runtime(); field.reset_power_invert_time()

            a = field(0)
            A = (a, field(1)) 
            T0, T1 = MontCurve.elligator(A)

            T0_x = get_affine_from_projective(T0)
            self.assertEqual(field.is_square(T0_x ** 3 + a * T0_x**2 + T0_x), True)
            T1_x = get_affine_from_projective(T1)
            self.assertEqual(field.is_square(T1_x ** 3 + a * T1_x**2 + T1_x), False)

            for _ in range(num_of_test - 1):
                a = field.get_random()
                # Ensure Ea is NOT singular.
                # It suffices to check a+2, a-2 are non-zero
                if a == 2 or a == -2:
                    continue

                # TODO: Check the case when A[1] != 1
                A = (a, field(1)) 
                T0, T1 = MontCurve.elligator(A)
                # print(f'a = {a}')
                # print(f'T0 = {T0}')
                # print(f'T1 = {T1}')
                
                T0_x = T0[0] * T0[1]**(-1)
                self.assertEqual(field.is_square(T0_x ** 3 + a * T0_x**2 + T0_x), True)
                T1_x = T1[0] * T1[1]**(-1)
                self.assertEqual(field.is_square(T1_x ** 3 + a * T1_x**2 + T1_x), False)
            

            field.show_runtime(label='{} {} elligators'.format(prime_name, num_of_test))
    

    def test_xdbl(self, num_of_test=500):
        for field, sage_Fp, MontCurve, prime_name in [
            (Fp1024, sage_GFp1024, MontCurve_p1024, 'p1024_CTIDH'), 
            (Fp2048, sage_GFp2048 ,MontCurve_p2048, 'p2048_CTIDH')
        ]:
            
            field.reset_runtime(); field.reset_power_invert_time()

            def test_one_curve(a=field(0), num_point=10):
                A = (a, field(1)) 
                sage_EC = get_sage_montgomery_curve(sage_Fp, a.get_int_value())
                
                T0, T1 = MontCurve.elligator(A)
                T0_x = get_affine_from_projective(T0)

                A24 = (A[0] + 2*A[1], 4*A[1])
                two_T0 = MontCurve.xdbl(T0, A24)
                two_T0_x = get_affine_from_projective(two_T0)

                T0_sage = sage_EC.lift_x(sage_Fp(T0_x))
                two_T0_sage = 2*T0_sage
                self.assertEqual(two_T0_x, two_T0_sage[0])
            
            test_one_curve(a=field(0), num_point=100)
            
            for _ in range(num_of_test - 1):
                a = field.get_random()
                test_one_curve(a)

    def test_xadd(self, num_of_test=500):
        for field, sage_Fp, MontCurve, prime_name in [
            (Fp1024, sage_GFp1024, MontCurve_p1024, 'p1024_CTIDH'), 
            (Fp2048, sage_GFp2048 ,MontCurve_p2048, 'p2048_CTIDH')
        ]:
            
            field.reset_runtime(); field.reset_power_invert_time()

            def test_one_curve(a=field(0), num_point=10):
                A = (a, field(1)) 
                sage_EC = get_sage_montgomery_curve(sage_Fp, a.get_int_value())
                
                P, _ = MontCurve.elligator(A)
                Q, _ = MontCurve.elligator(A)

                Px = get_affine_from_projective(P)
                Qx = get_affine_from_projective(Q)
               
                P_sage = sage_EC.lift_x(sage_Fp(Px))
                Q_sage = sage_EC.lift_x(sage_Fp(Qx))
                P_minus_Q_sage = P_sage - Q_sage
                P_minus_Q = (field(P_minus_Q_sage[0]), field(P_minus_Q_sage[2]))
                
                P_plus_Q = MontCurve.xadd(P, Q, P_minus_Q)
                P_plus_Qx = get_affine_from_projective(P_plus_Q)

                P_plus_Q_sage = P_sage + Q_sage
                self.assertEqual(P_plus_Qx, P_plus_Q_sage[0])
            
            test_one_curve(a=field(0), num_point=100)
            
            for _ in range(num_of_test - 1):
                a = field.get_random()
                test_one_curve(a)

            
