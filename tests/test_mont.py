import unittest
import json

import tqdm
from sage.all import EllipticCurve, proof, GF, kronecker_symbol
# from random import randint
# from sage.all import *

from CTIDH import MontgomeryCurve, PrimeField
from CTIDH.utils import get_randint

proof.arithmetic(False)


with open("data/prime_info/p1024_CTIDH", "r") as f:
    p1024_info = json.load(f)
with open("data/prime_info/p2048_CTIDH", "r") as f:
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


def get_affine_from_projective(A: tuple) -> int:
    """Given A = (Ax: Az), Az!= 0, compute a = Ax * Az^(-1)

    Args:
        A (list): [Ax, Az] that represents (Ax:Az), where Ax, Az are ZModPrime
    """
    Ax, Az = A[0], A[1]
    assert len(A) == 2
    assert Az != 0

    a = Ax * Az ** (-1)
    return int(a.get_int_value())

def load_supersingular_coefficients(prime_name='p1024-CTIDH') -> list:
    coeffs = []
    with open(f'tests/data/supersingular_coeff/{prime_name}', 'r') as file:
        for a in file:
            coeffs.append(int(a))
    return coeffs


# TODO: Create benchmarks (currently these benchmarks with tests are not accurate)
# Should make test and bench as two methods , not in the same time.
# TODO: Refactor these tests since they are in the same pattern.
class TestMontgomeryCurve(unittest.TestCase):
    def test_elligator(self, num_curve=20, num_point=5):
        for p, field, MontCurve, prime_name in [
            (p1024, Fp1024, MontCurve_p1024, "p1024_CTIDH"),
            (p2048, Fp2048, MontCurve_p2048, "p2048_CTIDH"),
        ]:
            field.reset_runtime()
            field.reset_power_invert_time()

            def test_one_curve(a=field(0), num_point=num_point):
                A = (a, field(1))

                for _ in range(num_point):
                    T0, T1 = MontCurve.elligator(A)
                    T0_x = get_affine_from_projective(T0)
                    self.assertEqual(
                        kronecker_symbol(T0_x**3+a.get_int_value()*T0_x**2+T0_x, p), 1
                    )
                    T1_x = get_affine_from_projective(T1)
                    self.assertEqual(
                        kronecker_symbol(T1_x**3+a.get_int_value()*T1_x**2+T1_x, p), -1
                    )

            test_one_curve(a=field(0))

            for _ in range(num_curve - 1):
                a = field.get_random()
                # Ensure Ea is NOT singular.
                # It suffices to check a+2, a-2 are non-zero
                if a == 2 or a == -2:
                    continue
                test_one_curve(a)

            field.show_runtime(
                label="{} {} elligators".format(prime_name, num_curve * num_point)
            )

    # TODO: Check the case when Az is not 1
    def test_xdbl(self, num_curve=20, num_point=5):
        for field, sage_Fp, MontCurve, prime_name in [
            (Fp1024, sage_GFp1024, MontCurve_p1024, "p1024_CTIDH"),
            (Fp2048, sage_GFp2048, MontCurve_p2048, "p2048_CTIDH"),
        ]:
            field.reset_runtime()
            field.reset_power_invert_time()

            def test_one_curve(a=field(0), num_point=num_point):
                A = (a, field(1))
                sage_EC = get_sage_montgomery_curve(sage_Fp, a.get_int_value())

                for _ in range(num_point):
                    T0, _ = MontCurve.elligator(A)
                    T0_x = get_affine_from_projective(T0)

                    A24 = (A[0] + 2 * A[1], 4 * A[1])
                    two_T0 = MontCurve.xdbl(T0, A24)
                    two_T0_x = get_affine_from_projective(two_T0)

                    T0_sage = sage_EC.lift_x(sage_Fp(T0_x))
                    two_T0_sage = 2 * T0_sage

                    self.assertEqual(two_T0_x, two_T0_sage[0])

            test_one_curve(a=field(0))

            for _ in range(num_curve - 1):
                a = field.get_random()
                test_one_curve(a)

            field.show_runtime(
                label="{} {} xdbl + elligator".format(prime_name, num_curve * num_point)
            )

    def test_xadd(self, num_curve=20, num_point=5):
        for field, sage_Fp, MontCurve, prime_name in [
            (Fp1024, sage_GFp1024, MontCurve_p1024, "p1024_CTIDH"),
            (Fp2048, sage_GFp2048, MontCurve_p2048, "p2048_CTIDH"),
        ]:
            field.reset_runtime()
            field.reset_power_invert_time()

            def test_one_curve(a=field(0), num_point=num_point):
                A = (a, field(1))
                sage_EC = get_sage_montgomery_curve(sage_Fp, a.get_int_value())

                for _ in range(num_point):
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

            test_one_curve(a=field(0))

            for _ in range(num_curve - 1):
                a = field.get_random()
                test_one_curve(a)

            field.show_runtime(
                label="{} {} xadd + elligator".format(prime_name, num_curve * num_point)
            )

    # TODO: Check the case when Az is not 1
    def test_xmul_Ladder(self, num_curve=20, num_point=5):
        for field, sage_Fp, MontCurve, prime_name, L in [
            (Fp1024, sage_GFp1024, MontCurve_p1024, "p1024_CTIDH", p1024_info["L"]),
            (Fp2048, sage_GFp2048, MontCurve_p2048, "p2048_CTIDH", p2048_info["L"]),
        ]:
            field.reset_runtime()
            field.reset_power_invert_time()

            def test_one_curve(a=field(0), num_point=num_point):
                A = (a, field(1))
                sage_EC = get_sage_montgomery_curve(sage_Fp, a.get_int_value())
                A24 = (A[0] + 2 * A[1], 4 * A[1])

                for _ in range(num_point):
                    P, _ = MontCurve.elligator(A)
                    Px = get_affine_from_projective(P)
                    sage_P = sage_EC.lift_x(sage_Fp(Px))
                    sage_Q = sage_P
                    i = get_randint(0, len(L) - 1)
                    sage_Q = L[i] * sage_P

                    Q = MontCurve.xmul_Ladder(P, A24, i)
                    Qx = get_affine_from_projective(Q)

                    self.assertEqual(Qx, sage_Q[0])

            test_one_curve(a=field(0))

            for _ in range(num_curve - 1):
                a = field.get_random()
                test_one_curve(a)

            field.show_runtime(
                label="{} {} xmul_Ladder + elligator".format(prime_name, num_curve * num_point)
            )


    def test_issupersingular_original(self, num_randcurve=20):
        print('Testing is_supersingular_original:')
        for field, sage_Fp, MontCurve, prime_name, L in [
            (Fp1024, sage_GFp1024, MontCurve_p1024, "p1024_CTIDH", p1024_info["L"]),
            # (Fp2048, sage_GFp2048, MontCurve_p2048, "p2048_CTIDH", p2048_info["L"]),
        ]:
            def test_one_curve(a=field(0), is_Ea_supersingular=False):
                # A = (a, field(1))
                u = field.get_random()
                A = (a*u, u)
                sage_EC = get_sage_montgomery_curve(sage_Fp, a.get_int_value())
                sage_is_supersingular = sage_EC.is_supersingular(proof=False)
                if is_Ea_supersingular:
                    self.assertEqual(sage_is_supersingular, True)

                self.assertEqual(MontCurve.issupersingular_original(A), sage_is_supersingular)
            
            test_one_curve(a=field(0))
            for _ in range(num_randcurve - 1):
                while True:
                    a = field.get_random()
                    if a != 2 and a != -2: # ensure the curve is non-singular
                        break
                test_one_curve(a)
            
            supersingular_coeffs = load_supersingular_coefficients(prime_name)
            for a in supersingular_coeffs:
                print(f'a = {a}')
                test_one_curve(field(a), True)