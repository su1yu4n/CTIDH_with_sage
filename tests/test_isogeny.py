import unittest
# from random import randint

from sage.all import EllipticCurve, proof, GF

from CTIDH import PrimeField, MontgomeryCurve, MontgomeryIsogeny

from CTIDH.utils import read_prime_info, batchmaxprime_of_Li, get_randint


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

#NOTE: To test velusqrt, Add hvelu in the future. Maybe also add svelu.


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

    # def test_kps_t(self, num_curve=5, num_isogeny=5):
    #     for sage_Fp, field, MontCurve, MontIsogeny in [
    #         (GF(p1024), Fp1024, MontCurve_p1024, isogeny_tvelu_p1024),
    #         (GF(p2048), Fp2048, MontCurve_p2048, isogeny_tvelu_p2048)
    #     ]:
    #         def test_one_curve(a=field(0), num_isogeny=num_isogeny):
    #             A = (a, field(1))
    #             A24 = (a+2, field(4))

    #             sage_EC = get_sage_montgomery_curve(sage_Fp, a.get_int_value())
    #             for _ in range(num_isogeny):
    #                 P, _ = MontCurve.elligator(A)
    #                 d_fake = (MontCurve.L[get_randint(0, MontCurve.n - 1)] - 1) // 2
    #                 Xi_Zis = MontIsogeny.kps_t(d_fake, P, A24) 
    #                 self.assertEqual(len(Xi_Zis), d_fake)
    #                 self.assertEqual(Xi_Zis[0], P)

    #                 k = get_randint(1, d_fake)
                    
    #                 Px = get_affine_from_projective(P)
    #                 P_sage = sage_EC.lift_x(sage_Fp(Px))
    #                 kP_sage = k * P_sage
    #                 self.assertEqual(get_affine_from_projective(Xi_Zis[k-1]), kP_sage.xy()[0])
                    
    #         test_one_curve(a=field(0), num_isogeny=2*num_isogeny)
    #         for _ in range(num_curve-1):
    #             test_one_curve(a=field.get_random())

        
    # def test_xisog_t(self, num_curve=5, num_isogeny=5):
    #     for sage_Fp, field, MontCurve, MontIsogeny in [
    #         (GF(p1024), Fp1024, MontCurve_p1024, isogeny_tvelu_p1024),
    #         (GF(p2048), Fp2048, MontCurve_p2048, isogeny_tvelu_p2048)
    #     ]:           
    #         # test and return a new curve's coefficient, because we need to ensure the curve is supersingular
    #         def test_one_curve(a=field(0), num_isogeny=num_isogeny) -> int:
    #             A = (a, field(1))
    #             A24 = (a+2, field(4))

    #             sage_EC = get_sage_montgomery_curve(sage_Fp, a.get_int_value())
    #             A_new = 1
    #             for _ in range(num_isogeny):
    #                 ind = get_randint(0, MontCurve.n - 1)
    #                 l = MontCurve.L[ind]
    #                 # print(f'ind = {ind}')
    #                 l_fake = batchmaxprime_of_Li(ind, MontCurve.batch_start, MontCurve.batch_stop, MontCurve.L)
    #                 d = (l-1) // 2
    #                 d_fake = (l_fake - 1)//2

    #                 assert d_fake >= d

    #                 while True:
    #                     P, _ = MontCurve.elligator(A)
    #                     P = MontCurve.xdbl(P, A24); P = MontCurve.xdbl(P, A24) # clear cofactor
    #                     for j in range(ind):
    #                         P = MontCurve.xmul_public(P, A24, j)
    #                     for j in range(ind+1, MontCurve.n):
    #                         P = MontCurve.xmul_public(P, A24, j)
    #                     if not MontCurve.isinfinity(P):
    #                         break

    #                 Xi_Zis = MontIsogeny.kps_t(d_fake, P, A24)
    #                 Xi_Zi_hats = [(Xi+Zi, Xi-Zi) for (Xi, Zi) in Xi_Zis]

    #                 A_new = MontIsogeny.xisog_t(d, d_fake, Xi_Zi_hats, A)

    #                 Px = sage_Fp(get_affine_from_projective(P))
    #                 sage_P = sage_EC.lift_x(Px)
    #                 phi = sage_EC.isogeny(kernel=sage_P, model='montgomery')
    #                 # print(f'codomain is {phi.codomain()}')
    #                 sage_A_new = phi.codomain().a2()
    #                 self.assertEqual(sage_Fp(get_affine_from_projective(A_new)), sage_A_new)

    #             return int(get_affine_from_projective(A_new))
            
    #         a_new = test_one_curve()
    #         for _ in range(num_curve-1):
    #             a_new = test_one_curve(field(a_new))       


    # def test_xeval_t(self, num_curve=5, num_isogeny=3):
    #     for sage_Fp, field, MontCurve, MontIsogeny in [
    #         (GF(p1024), Fp1024, MontCurve_p1024, isogeny_tvelu_p1024),
    #         (GF(p2048), Fp2048, MontCurve_p2048, isogeny_tvelu_p2048)
    #     ]:           
    #         # test and return a new curve's coefficient, because we need to ensure the curve is supersingular
    #         def test_one_curve(a=field(0), num_isogeny=num_isogeny) -> int:
    #             A = (a, field(1))
    #             A24 = (a+2, field(4))

    #             sage_EC = get_sage_montgomery_curve(sage_Fp, a.get_int_value())
    #             A_new = 1
    #             for _ in range(num_isogeny):
    #                 ind = get_randint(0, MontCurve.n - 1)
    #                 l = MontCurve.L[ind]
    #                 # print(f'ind = {ind}')
    #                 l_fake = batchmaxprime_of_Li(ind, MontCurve.batch_start, MontCurve.batch_stop, MontCurve.L)
    #                 d = (l-1) // 2
    #                 d_fake = (l_fake - 1)//2

    #                 assert d_fake >= d

    #                 while True:
    #                     T, _ = MontCurve.elligator(A)
    #                     P = T
    #                     P = MontCurve.xdbl(P, A24); P = MontCurve.xdbl(P, A24) # clear cofactor
    #                     for j in range(ind):
    #                         P = MontCurve.xmul_public(P, A24, j)
    #                     for j in range(ind+1, MontCurve.n):
    #                         P = MontCurve.xmul_public(P, A24, j)
    #                     if not MontCurve.isinfinity(P):
    #                         break

    #                 Xi_Zis = MontIsogeny.kps_t(d_fake, P, A24)
    #                 Xi_Zi_hats = [(Xi+Zi, Xi-Zi) for (Xi, Zi) in Xi_Zis]

    #                 A_new = MontIsogeny.xisog_t(d, d_fake, Xi_Zi_hats, A)
    #                 phi_T = MontIsogeny.xeval_t(d, d_fake, Xi_Zi_hats, T)

    #                 Px = sage_Fp(get_affine_from_projective(P))
    #                 sage_P = sage_EC.lift_x(Px)
    #                 self.assertEqual(sage_P.order(), l)
    #                 sage_phi = sage_EC.isogeny(kernel=sage_P, model='montgomery')
    #                 Tx = get_affine_from_projective(T)
    #                 sage_T = sage_EC.lift_x(sage_Fp(Tx))
    #                 sage_phi_Tx = sage_phi(sage_T).xy()[0]

    #                 self.assertEqual(sage_phi_Tx, get_affine_from_projective(phi_T))

    #             return int(get_affine_from_projective(A_new))
           
    #         a_new = test_one_curve()
    #         for _ in range(num_curve-1):
    #             a_new = test_one_curve(field(a_new))


    def test_tvelu(self, num_curve=5, num_isogeny=3):
        for sage_Fp, field, MontCurve, MontIsogeny in [
            (GF(p1024), Fp1024, MontCurve_p1024, isogeny_tvelu_p1024),
            # (GF(p2048), Fp2048, MontCurve_p2048, isogeny_tvelu_p2048) # too slow for a regular test
        ]:           
            # test and return a new curve's coefficient, because we need to ensure the curves we choose are all supersingular
            def test_one_curve(a=field(0), num_isogeny=num_isogeny) -> int:
                A = (a, field(1))
                A24 = (a+2, field(4))

                sage_EC = get_sage_montgomery_curve(sage_Fp, a.get_int_value())
                A_new = 1
                for _ in range(num_isogeny):
                    # ind = get_randint(0, MontCurve.n - 1)
                    ind = get_randint(0, 15)
                    l = MontCurve.L[ind]
                    # print(f'ind = {ind}')
                    l_fake = batchmaxprime_of_Li(ind, MontCurve.batch_start, MontCurve.batch_stop, MontCurve.L)
                    d = (l-1) // 2
                    d_fake = (l_fake - 1)//2

                    assert d_fake >= d

                    while True:
                        T, _ = MontCurve.elligator(A)
                        P = T
                        P = MontCurve.xdbl(P, A24); P = MontCurve.xdbl(P, A24) # clear cofactor
                        for j in range(ind):
                            P = MontCurve.xmul_public(P, A24, j)
                        for j in range(ind+1, MontCurve.n):
                            P = MontCurve.xmul_public(P, A24, j)
                        if not MontCurve.isinfinity(P):
                            break
                    
                    # Do kps_t and test 
                    Xi_Zis = MontIsogeny.kps_t(d_fake, P, A24)
                    Xi_Zi_hats = [(Xi+Zi, Xi-Zi) for (Xi, Zi) in Xi_Zis]
                    
                    self.assertEqual(len(Xi_Zis), d_fake)
                    self.assertEqual(Xi_Zis[0], P)
                    bound = d_fake - 1 if d_fake <= l else l-1
                    k = get_randint(1, bound)
                    Px = get_affine_from_projective(P)
                    P_sage = sage_EC.lift_x(sage_Fp(Px))
                    kP_sage = k * P_sage
                    self.assertEqual(get_affine_from_projective(Xi_Zis[k-1]), kP_sage.xy()[0])

                    # test the correctness of xisog_t, xeval_t 
                    A_new = MontIsogeny.xisog_t(d, d_fake, Xi_Zi_hats, A)
                    phi_T = MontIsogeny.xeval_t(d, d_fake, Xi_Zi_hats, T)

                    Px = sage_Fp(get_affine_from_projective(P))
                    sage_P = sage_EC.lift_x(Px)
                    self.assertEqual(sage_P.order(), l)
                    sage_phi = sage_EC.isogeny(kernel=sage_P, model='montgomery')
                    Tx = get_affine_from_projective(T)
                    sage_T = sage_EC.lift_x(sage_Fp(Tx))
                    sage_phi_Tx = sage_phi(sage_T).xy()[0]

                    sage_A_new = sage_phi.codomain().a2()
                    self.assertEqual(sage_Fp(get_affine_from_projective(A_new)), sage_A_new)

                    self.assertEqual(sage_phi_Tx, get_affine_from_projective(phi_T))

                return int(get_affine_from_projective(A_new))
           
            a_new = test_one_curve()
            for _ in range(num_curve-1):
                a_new = test_one_curve(field(a_new))
    

    def test_matryoshka_isogeny_tvelu(self, num_curve=5, num_isogeny=3):
        pass

    def test_kps_s(self, num_curve=5, num_isogeny=3):
        pass
    
    
    def test_xisog_s(self, num_curve=5, num_isogeny=3):
        pass
    
    
    def test_xeval_s(self, num_curve=5, num_isogeny=3):
        pass

    
    def test_svelu(self, num_curve=5, num_isogeny=3):
        pass


    def test_matryoshka_isogeny_svelu(self, num_curve=5, num_isogeny=3):
        pass
