import unittest
from random import randint

from sage.all import EllipticCurve, proof, GF

from CTIDH import PrimeField, MontgomeryCurve, MontgomeryIsogeny

from CTIDH.utils import read_prime_info


p1024_info = read_prime_info('p1024_CTIDH')
p2048_info = read_prime_info('p2048_CTIDH')
p1024 = p1024_info["p"]
p2048 = p2048_info["p"]

Fp1024 = PrimeField(p1024)
Fp2048 = PrimeField(p2048)

sage_GFp1024 = GF(p1024)
sage_GFp2048 = GF(p2048)

MontCurve_p1024 = MontgomeryCurve("p1024_CTIDH")
MontCurve_p2048 = MontgomeryCurve("p2048_CTIDH")

isogeny_tvelu_p1024 = MontgomeryIsogeny('tvelu')(MontCurve_p1024)
isogeny_tvelu_p2048 = MontgomeryIsogeny('tvelu')(MontCurve_p2048)



def get_sage_montgomery_curve(sage_Fp, a: int):
    return EllipticCurve(sage_Fp, [0, a, 0, 1, 0])


def get_affine_from_projective(A: list) -> int:
    """Given A = (Ax: Az), Az!= 0, compute a = Ax * Az^(-1).

    Args:
        A (list): [Ax, Az] that represents (Ax:Az), where Ax, Az are ZModPrime
    """
    Ax, Az = A[0], A[1]
    assert len(A) == 2
    assert Az != 0

    a = Ax * Az ** (-1)
    return int(a.get_int_value())

class TestMontgomeryIsogeny(unittest.TestCase):

    def test_kps_t(self, num_curve=10, num_d=5):
        for sage_Fp, field, MontCurve, MontIsogeny in [
            (GF(p1024), Fp1024, MontCurve_p1024, isogeny_tvelu_p1024),
            (GF(p2048), Fp2048, MontCurve_p2048, isogeny_tvelu_p2048)
        ]:
            def test_one_curve(a=field(0), num_d=num_d):
                A = (a, field(1))
                A24 = (a+2, field(4))

                sage_EC = get_sage_montgomery_curve(sage_Fp, a.get_int_value())
                for _ in range(num_d):
                    P, _ = MontCurve.elligator(A)
                    d_fake = (MontCurve.L[randint(0, MontCurve.n - 1)] - 1) // 2
                    Xi_Zis = MontIsogeny.kps_t(d_fake, P, A24) 
                    self.assertEqual(len(Xi_Zis), d_fake)
                    self.assertEqual(Xi_Zis[0], P)

                    k = randint(1, d_fake)
                    
                    Px = get_affine_from_projective(P)
                    P_sage = sage_EC.lift_x(sage_Fp(Px))
                    kP_sage = k * P_sage
                    self.assertEqual(get_affine_from_projective(Xi_Zis[k-1]), kP_sage.xy()[0])
                    
            test_one_curve(a=field(0), num_d=2*num_d)
            for _ in range(num_curve-1):
                test_one_curve(a=field.get_random())

        
    def test_xisog_t(self):
        pass
    def test_xeval_t(self):
        pass
    def test_kps_s(self):
        pass
    def test_xisog_s(self):
        pass
    def test_xeval_s(self):
        pass