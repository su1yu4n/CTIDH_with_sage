from typing import List
from math import floor, sqrt

from CTIDH.mont import MontgomeryCurve
from CTIDH.utils import read_velusqrt_steps_info, hamming_weight, bitlength, isequal, batchmaxprime_of_Li, batchminprime_of_Li, batchnumber_of_Li, CMOV, CSWAP
from CTIDH.polymul import PolyMul
from CTIDH.polyredc import PolyRedc
from CTIDH.primefield import PrimeField
import numpy


def doc(s):
    class __doc(object):
        def __init__(self,f):
            self.func = f
            self.desc = s
        def __call__(self,*args,**kwargs):
            return self.func(*args,**kwargs)
        def __repr__(self):
            return self.desc
    return __doc

# Velu and Velusqrt formula that compute small odd prime degree isogeny.
def MontgomeryIsogeny(formula_name='tvelu', uninitialized = False):
    cutoff = 103
    cutoff_string = f' with cutoff ell <= {cutoff}' if formula_name == 'hvelu' else ''
    NAME = 'Isogeny class using the %s Velu\'s formulae%s' % ({'tvelu':'traditional', 'svelu':'square-root', 'hvelu':'hybrid'}[formula_name], cutoff_string)

    @doc(NAME)
    class Formulae:
        def __init__(self, curve, tuned=True, scaled=False):
            """_summary_

            Args:
                curve (MontgomeryCurve): The object returned by MontgomeryCurve() in mont.py
                tuned (bool, optional): True if fine-tuned velusqrt information is presented in data folder. Defaults to True.
                scaled (bool, optional): Use scaled remainder tree or not. 
                If True, it will read and use velusqrt tuned info for scaled version. Defaults to False.
            """
            self.formula_name = formula_name

            if formula_name != 'tvelu':
                if tuned:
                    self.sI_list, self.sJ_list = read_velusqrt_steps_info(curve.prime_name, scaled)
                else:
                    self.sI_list = None
                    self.sJ_list = None

            self.HYBRID_BOUND = {'tvelu':max(curve.L), 'svelu':1, 'hvelu':cutoff}[formula_name]
            
            # Global variables to be used in kps, xisog, and xeval

            # Here, J is a set of cardinality sJ
            self.I = None
            self.J = None
            self.sJ = None

            # Here, ptree_I corresponds with the product tree determined by I, and I is a set of cardinality sJ
            self.ptree_hI = None
            self.sI = (None,)

            # Here, K is a set of cardinality sK
            self.K = None
            self.sK = None

            # An extra global variable which is used in xisog and xeval
            self.XZJ4 = None
            self.XZj_add = None
            self.XZj_sub = None
            self.ADD_SQUARED = None
            self.SUB_SQUARED = None
            self.XZk_add = None
            self.XZk_sub = None
    
            

            self.SCALED_MULTIEVALUATION = scaled
            self.tuned = tuned

            self.prime_name = curve.prime_name
            self.curve = curve
            self.field = self.curve.field
            self.L = self.curve.L
            self.batch_start = self.curve.batch_start
            self.batch_stop = self.curve.batch_stop

            self.poly_mul = PolyMul(self.field)
            self.poly_redc = PolyRedc(self.poly_mul)

            self.c_xeval = list(
                map(self.ceval, self.L)
            )  # list of the costs of each degree-l isogeny evaluation
            self.c_xisog = list(
                map(self.cisog, self.L)
            )  # list of the costs of each degree-l isogeny construction
            
             
            '''
            TODO: Compute the cost.
            NOTE: Currently this is not performed, because it requires finding generators, which need PRAC, cofactor_multiples and more...
                Also it is not clear that whether sibc's routine is adaptable to CTIDH .
            '''
            # # Now, we proceed to store all the correct costs
            # if formula_name != 'tvelu' and uninitialized:
            #     print("// Precomputation cost regarding kps, xisog, and xeval of the velusqrt formulae")
            #     self.velusqrt_cost()
                

        def ceval(self, l: int):
            return numpy.array([2.0 * (l - 1.0), 2.0, (l + 1.0)])


        def cisog(self, l: int):
            return numpy.array(
                [
                    (
                        3.0 * l
                        + 2.0 * hamming_weight(l)
                        - 9.0
                        + isequal[l == 3] * 4.0
                    ),
                    (l + 2.0 * bitlength(l) + 1.0 + isequal[l == 3] * 2.0),
                    (3.0 * l - 7.0 + isequal[l == 3] * 6.0),
                ]
            )


        # TODO: Possibly add a parameter A24 to save a few additions. (not very interesting..)
        def matryoshka_isogeny(self, A: tuple, Ts:List[tuple], Tnewlen: int, P: tuple, i: int):
            """ Computing the L[i]-isogeny phi: EA -> EA'. 
                The computation leverages matryoshka-doll structure to resist timing attack.

            Args:
                A (tuple): fractional representation of affine coefficient a = A[0]/A[1]
                Ts (List[tuple]): Points (possibly) to be pushed.
                Tnewlen (int): 0 -> won't push any point through the isogeny. 1 -> push Ts[0], 2 means push two points. 
                P (tuple): projective x-coordinate of the kernel point P. (i.e. EA' = EA/<P>)
                i (int): the index of the prime

            Returns:
                tuple: Anew, Ts (list of T0, T1). 
                Note that Ts always consists of two points, although some of them may be unchanged (same as the input).
            """
            # print(self.HYBRID_BOUND)
            assert 0 <= Tnewlen <= 2
            assert i < len(self.L)
            # print(self.HYBRID_BOUND)
            l = self.L[i]
            l_fake = batchmaxprime_of_Li(i, self.batch_start, self.batch_stop, self.L)
            # print(f"L:{l} and L_fake:{l_fake}")
            A24 = self.curve.xA24(A)
            # NOTE: Here we use l_fake to decide which formula to use(traditional velu or velusqrt) to avoid timing attack.
            # It may be unsafe if different primes in the same batch use different formulae.
            if l_fake <= self.HYBRID_BOUND:     # Use Velu formula
                d = (l-1)//2
                d_fake = (l_fake-1)//2

                xiP_list = self.kps_t(d_fake, P, A24)
                Xi_Zi_hats = [(Xi+Zi, Xi-Zi) for (Xi, Zi) in xiP_list]
                # for xiP in xiP_list:
                #     Xi, Zi = xiP
                #     Xi_Zi_hat = (Xi + Zi, Xi - Zi)
                #     Xi_Zi_hats.append(Xi_Zi_hat)

                A_new = self.xisog_t(d, d_fake, Xi_Zi_hats, A)

                if Tnewlen > 0:
                    Ts[0] = self.xeval_t(d, d_fake, Xi_Zi_hats, Ts[0])
                if Tnewlen > 1:
                    Ts[1] = self.xeval_t(d, d_fake, Xi_Zi_hats, Ts[1])

                return A_new, Ts
            
            else: 
                # Use Velusqrt formulas
             
                self.sJ = self.sJ_list[i]
                self.sI = self.sI_list[i]
                d = ((l - 2 - 4 * self.sJ * self.sI - 1) // 2) + 1
                d_fake = ((l_fake - 2 - 4 * self.sJ * self.sI - 1) // 2) + 1
                self.sK = d
                self.sK_fake = d_fake
                # print(f"mul_counts:{mul_counts}")
                # print(f"pow_counts:{pow_counts}")
                # else:
                #     if l == 3:
                #         b = 0
                #         c = 0
                #     else:
                #         b = int(floor(sqrt(l - 1) / 2.0))
                #         c = int(floor((l- 1.0) / (4.0 * b)))
                #     d = ((l - 2 - 4 * b * c - 1) // 2) + 1
                #     d_fake = ((l_fake - 2 - 4 * b * c - 1) // 2) + 1
                #     print("%%%%%%%")
                    # self.sJ = b
                    # self.sI = c
                    # self.sK = d
                    # self.sK_fake = d_fake
                    # print("########")    
                # Now sI, sJ, sK are set.
                self.kps_s(P, A24 , i)
                self.SUB_SQUARED = [0 for j in range(0, self.sJ, 1)]  
                self.ADD_SQUARED = [0 for j in range(0, self.sJ, 1)]
                self.XZJ4 = [0 for j in range(0, self.sJ, 1)]
                self.XZj_add = [0 for j in range(0, self.sJ, 1)]
                self.XZj_sub = [0 for j in range(0, self.sJ, 1)]
                self.XZk_add = [0 for k in range(0, self.sK_fake, 1)]
                self.XZk_sub = [0 for k in range(0, self.sK_fake, 1)]
                for j in range(0, self.sJ, 1):
                    self.XZj_add[j] = self.J[j][0] + self.J[j][1]
                    self.XZj_sub[j] = self.J[j][0] - self.J[j][1]
                    self.SUB_SQUARED[j] = self.XZj_sub[j]**2
                    self.ADD_SQUARED[j] = self.XZj_add[j]**2
                    self.XZJ4[j] = self.SUB_SQUARED[j] - self.ADD_SQUARED[j]

                for k in range(0, self.sK_fake, 1):
                    self.XZk_add[k] = self.K[k][0] + self.K[k][1]
                    self.XZk_sub[k] = self.K[k][1] - self.K[k][0]
                   
                P2 = self.curve.xdbl(P, A24)
                if self.sI == 0:
                    self.ptree_hI = None
                else:
                    if self.sI == 1:
                        I = [list(P2)]
                        hI = [
                        list([(0 - P2[0]), P2[1]])
                        ]
                    else:
                        hI = [
                        [(0 - iP[0]), iP[1]] for iP in self.I
                        ]  # we only need to negate x-coordinate of each point
                self.ptree_hI = self.poly_mul.product_tree(
                    hI, self.sI
                )  # product tree of hI
            
                # if not self.SCALED_MULTIEVALUATION:
                #         # Using scaled remainder trees
                #         self.ptree_hI = self.poly_redc.reciprocal_tree(
                #             {'rpoly': [1], 'rdeg': 0, 'fpoly': [1], 'fdeg': 0, 'a': 1},
                #             2 * self.sJ + 1,
                #             self.ptree_hI,
                #             self.sI,
                #         )  # reciprocal of the root is used to compute sons' reciprocals

                # else:
                # Using scaled remainder trees
                self.ptree_hI.extend([None, None, None, None])  # Using scaled remainder trees
                if self.sI < (2 * self.sJ - self.sI + 1) or self.sI == 1:
                    (
                        self.ptree_hI[4],
                        self.ptree_hI[5],
                    ) = self.poly_redc.reciprocal(
                        self.ptree_hI[2][::-1],
                        self.sI + 1,
                        2 * self.sJ - self.sI + 1,
                    )
                    self.ptree_hI[7], self.ptree_hI[6] = (
                        list(self.ptree_hI[4][: self.sI]),
                        self.ptree_hI[5],
                    )

                else:
                    (
                        self.ptree_hI[7],
                        self.ptree_hI[6],
                    ) = self.poly_redc.reciprocal(
                        self.ptree_hI[2][::-1], self.sI + 1, self.sI
                    )
                    self.ptree_hI[4], self.ptree_hI[5] = (
                        list(
                            self.ptree_hI[7][: (2 * self.sJ - self.sI + 1)]
                        ),
                        self.ptree_hI[6],
                    ) 
                
                A_new = self.xisog_s( A24, i) 
                
                if Tnewlen > 0:
                    Ts[0] = self.xeval_s(Ts[0], A24)
                if Tnewlen > 1:
                    Ts[1] = self.xeval_s(Ts[1], A24)  
                    
                return A_new, Ts
                
        
        def kps_t(self, d_fake: int, P: tuple, A24: tuple) -> List[tuple]:
            """Timing attack safe kps for traditional velu formula, 
            used in computing the l-isogeny phi: E -> E/<P>.

            Return the list of x([i]P) for i from 1 to d_fake = (l_fake - 1) // 2

            Args:
                d_fake (int): See above. l_fake is the largest small odd prime in the batch of l
                P (tuple): projective x-coordinate of the kernel point P (i.e. x(P))
                A24 (tuple): A24 = (Ax+2Az : 4Az)
            """
            Xi_Zis = []
            Xi_Zis.append(P)
            if d_fake >= 2:
                Xi_Zis.append(self.curve.xdbl(P, A24))
            for i in range(2, d_fake):
                Xi_Zis.append(self.curve.xadd(Xi_Zis[i-1], P, Xi_Zis[i-2]))        
            
            return Xi_Zis
        

        def xisog_t(self, d: int, d_fake: int, Xi_Zi_hats: List[tuple], A: tuple) -> tuple:
            """Timing attack safe xisog for traditional velu formula,
            Return the fraction representation of quadratic term's coefficient of the l-isogeny's codomain. 
            
            Args:
                d (int): d = (l-1)/2.
                d_fake (int): (l_fake - 1)/2
                Xi_Zi_hats (List[tuple]): list of (Xi+Zi : Xi-Zi) where (Xi : Zi) is x([i]P), P a generator of ker phi. 
                i ranges from 1 to d_fake.
                A (tuple): A = (Ax: Az). Ax and Az must have type ZModPrime(Primefield).
                Ax/Az is the quadratic term's coefficient of domain curve's affine equation. 

            Returns:
                tuple: A' = (Ax': Az'), where Ax'/Az' is the quadratic term's coefficient of the codomain.
                Ax', Az' have type ZModPrime
            """
            assert d_fake >= d
            Ax, Az = A
            l = 2*d + 1; l_maxbitlen = (2*d_fake + 1).bit_length()

            t = Az + Az
            aE = Ax + t; dE = Ax - t
            al = aE.safe_pow(l, l_maxbitlen); dl = dE.safe_pow(l, l_maxbitlen)
            pi_Y = self.field(1); pi_Z = self.field(1)

            # print(f'd = {d}, d_fake = {d_fake}')
            for i in range(d_fake):
                tmp1 = pi_Y * Xi_Zi_hats[i][1]
                tmp2 = pi_Z * Xi_Zi_hats[i][0]
                pi_Y = CMOV(pi_Y, tmp1, i <= d-1)
                pi_Z = CMOV(pi_Z, tmp2, i <= d-1)

            # for i in range(d):
            #     pi_Y *= Xi_Zi_hats[i][1]
            #     pi_Z *= Xi_Zi_hats[i][0]
            # print(f'Correct pi_Y = {pi_Y}, pi_Z = {pi_Z}')
            
            aE_new = al * pi_Z ** 8; dE_new = dl * pi_Y ** 8
            aE_dE = aE_new + dE_new; Ax_new = aE_dE + aE_dE; Az_new = aE_new - dE_new

            return (Ax_new, Az_new)
        

        def xeval_t(self, d: int, d_fake: int, Xi_Zi_hats: List[tuple], T: tuple) -> tuple:
            """Push the point T to the codomain through the l-isogeny phi.

            Args:
                d (int): degree of isogeny
                d_fake (int): the largest prime in the same batch as l.
                Xi_Zi_hats (List[tuple]): (Xi+Zi : Xi-Zi) where (Xi:Zi) x[i]P, i = 1, ..., d_fake
                T (tuple): (projective x-coordinate of) the point to push 

            Returns:
                tuple: (projective x-coordinate of) the image point phi(T)
            """
            X, Z = T
            X_hat, Z_hat = X+Z, X-Z
            X1_hat, Z1_hat = Xi_Zi_hats[0]
            X_prime, Z_prime = self.curve.crisscross(X1_hat, Z1_hat, X_hat, Z_hat)
            for i in range(1, d_fake):
                Xi_hat, Zi_hat = Xi_Zi_hats[i]
                t0, t1 = self.curve.crisscross(Xi_hat, Zi_hat, X_hat, Z_hat)
                tmp1, tmp2 = t0*X_prime, t1*Z_prime
                X_prime = CMOV(X_prime, tmp1, i<=d-1)
                Z_prime = CMOV(Z_prime, tmp2, i<=d-1)
            X_prime, Z_prime = X*X_prime**2, Z*Z_prime**2
            return X_prime, Z_prime


        # NOTE: This functions is used for setting the cardinalities sI, sJ, and sK
        # In sibc it is called by velusqrt_cost()
        # def set_parameters_velu(self, b, c, i):
        #     assert b <= c
        #     # At this step, everythin is correct
        #     self.sJ = b
        #     self.sI = c
        #     # print(f"b:{b}")
        #     # print(f"c:{c}")
        #     d = ((self.L[i] - 2 - 4 * b * c - 1) // 2) + 1
        #     assert d >= 0
        #     self.sK = d
        #     return None

        # TODO: Implement these velusqrt algorithms
        def kps_s(self, P: tuple,A24: tuple, i: int) -> None:
            # print(f"kps_s{self.L[i]}")
            if self.sI == 0:
                self.I = []
                self.J = []
                self.K_fake = [list(P)]
                # print("sI=0")
                # print(f"self.sK:{self.sK}")
                # print(f"self.sK:{self.sI}")
                # print(f"self.sK:{self.sJ}")
                # print(f"self.L[i]:{self.L[i]}")
                # J, b, ptree_hI, c, K, d
                assert self.sJ == 0 and self.sI == 0 and self.sK_fake == 1
                
            elif self.sI == 1:
                assert self.sJ == 1
                # print("sI=1")
                P2 = self.curve.xdbl(P, A24)
                self.J = [list(P)]
                self.I = [list(P2)]
                
                assert self.sK_fake <= 1
                if self.sK_fake == 1:
                    self.K = [list(P2)]
                else:
                    self.K = []
                            
            elif self.sI > 1 :# At this step, sI > 1
                assert self.sI > 1
                if self.sJ == 1:
                    # print("sI>1 and sJ=1")
                    # This branch corresponds with L[i] = 11 and L[i] = 13
                    Q = self.curve.xdbl(P, A24)  # x([2]P)
                    Q2 = self.curve.xdbl(Q, A24)  # x([2]Q)

                    self.J = [list(P)]

                    self.I = [[0, 0]] * self.sI
                    # print(self.sI)
                    self.I[0] = list(Q)  # x(   Q)
                    self.I[1] = self.curve.xadd(Q2, self.I[0], self.I[0])  # x([3]Q)
                    for ii in range(2, self.sI, 1):
                        self.I[ii] = self.curve.xadd(self.I[ii - 1], Q2, self.I[ii - 2])  # x([2**i + 1]Q)
                    # print(len(self.I))    
                    assert self.sK_fake <= 1
                    if self.sK_fake == 1:
                        self.K = [list(Q)]
                    else:
                        self.K = []
                
                else:
                    # print("sJ>1 and sI>1")
                    # Computing [j]P for each j in {1, 3, ..., 2*sJ - 1}
                    self.J = [[0, 0]] * self.sJ
                    self.J[0] = list(P)  
                    P2 = self.curve.xdbl(P, A24) 
                    self.J[1] = self.curve.xadd(P2, self.J[0], self.J[0])  
                    for jj in range(2, self.sJ, 1):
                        self.J[jj] = self.curve.xadd(
                            self.J[jj - 1], P2, self.J[jj - 2]
                        ) 
                    # -------------------------------------------------------
                    # Computing [i]P for i in { (2*sJ) * (2i + 1) : 0 <= i < sI}
                    bhalf_floor = self.sJ // 2
                    bhalf_ceil = self.sJ - bhalf_floor
                    P4 = self.curve.xdbl(P2, A24)  
                    Q = self.curve.xadd(
                    self.J[bhalf_ceil], self.J[bhalf_floor - 1], P4 if self.sJ % 2 else P2
                    )  
                    
                    
                    self.I = [[0, 0]] * self.sI
                    self.I[0] = list(Q)  
                    Q2 = self.curve.xdbl(Q, A24)  
                    self.I[1] = self.curve.xadd(Q2, self.I[0], self.I[0]) 
                    for ii in range(2, self.sI, 1):
                        self.I[ii] = self.curve.xadd(self.I[ii - 1], Q2, self.I[ii - 2])  
                        
                    # --------------------------------------------------------------
                    # Computing [k]P for k in { 4*sJ*sI + 1, ..., l - 6, l - 4, l - 2}
                    self.K = [[0, 0]] * self.sK_fake

                    if self.sK_fake >= 1:
                        self.K[0] = list(P2) 
                    if self.sK_fake >= 2:
                        self.K[1] = list(P4) 

                    for k in range(2, self.sK_fake, 1):
                        self.K[k] = self.curve.xadd(self.K[k - 1], P2, self.K[k - 2])
                        
            # print(f"len:{len(self.I)}")
            # print(self.tuned)
            # # print(f"self.L[i]:{self.L[i]}")
            # print(f"self.sJ{self.sJ}")
            # print(f"self.sI{self.sI}")
            # print(f"self.sK_fake{self.sK_fake}")
            # print(f"self.sK{self.sK}")
            assert len(self.I) == self.sI
            assert len(self.J) == self.sJ
            assert len(self.K) == self.sK_fake
            assert self.sI >= self.sJ
            assert self.sK_fake >= 0
            assert self.sK >= 0
            return None
            
        
        def xisog_s(self, A24: tuple, i: int) -> tuple:
            B = 2*(2*A24[0] - A24[1])
            C = A24[1]
            
            EJ_0 = [[0, 0, 0] for j in range(0, self.sJ, 1)]
            EJ_1 = [[0, 0, 0] for j in range(0, self.sJ, 1)]
            
            
            for j in range(0, self.sJ, 1):
                tadd = (self.ADD_SQUARED[j] * C)
                tsub = (self.SUB_SQUARED[j] * C)
                
                tadd2 = (tadd + tadd)
                tsub2 = (tsub + tsub)
                
                t1 = (self.XZJ4[j] * B)
                linear = (
                    t1 - tadd2
                )
                EJ_0[j] = [tsub, linear, tsub]
                
                linear = (
                    tsub2 - t1
                )
                EJ_1[j] = [tadd, linear, tadd]
                
            poly_EJ_0 = self.poly_mul.product_selfreciprocal_tree(EJ_0, self.sJ)[
                2
            ]
            poly_EJ_1 = self.poly_mul.product_selfreciprocal_tree(EJ_1, self.sJ)[
                2
            ] 
            
             # Approach using scaled remainder trees
            if self.ptree_hI != None:
                poly_EJ_0 = self.poly_redc.poly_redc(
                    poly_EJ_0, 2 * self.sJ + 1, self.ptree_hI
                )
                fg_0 = self.poly_mul.poly_mul_middle(
                    self.ptree_hI[7], self.sI, poly_EJ_0[::-1], self.sI
                )
                remainders_EJ_0 = self.poly_redc.multieval_scaled(
                    fg_0[::-1],
                    self.sI,
                    [1] + [0] * (self.sI - 1),
                    self.sI,
                    self.ptree_hI,
                    self.sI,
                )

                poly_EJ_1 = self.poly_redc.poly_redc(
                    poly_EJ_1, 2 * self.sJ + 1, self.ptree_hI
                )
                fg_1 = self.poly_mul.poly_mul_middle(
                    self.ptree_hI[7], self.sI, poly_EJ_1[::-1], self.sI
                )
                remainders_EJ_1 = self.poly_redc.multieval_scaled(
                    fg_1[::-1],
                    self.sI,
                    [1] + [0] * (self.sI - 1),
                    self.sI,
                    self.ptree_hI,
                    self.sI,
                )
            else:
                remainders_EJ_0 = []
                remainders_EJ_1 = []

            R_0 = self.poly_mul.product(remainders_EJ_0, self.sI)
            R_1 = self.poly_mul.product(remainders_EJ_1, self.sI)
             
            
            # hK_0 = [
            #     [self.XZk_sub[k]]
            #     for k in range(0, self.sK, 1)
            # ]
            # M_0 = self.poly_mul.product(
            #     hK_0, self.sK
            # )  # product of (Zk - Xk) for each k in K
            # Case alpha = -1
            # hK_1 = [
            #     [self.XZk_add[k]]
            #     for k in range(0, self.sK_fake, 1)
            # ]
            # M_1 = self.poly_mul.product(
            #     hK_1, self.sK
            # )  # product of (Zk + Xk) for each k in K
            M_0 = self.field(1); M_1 = self.field(1)
            for k in range(0, self.sK_fake, 1):
                tmp0 = M_0 * self.XZk_sub[k]
                tmp1 = M_1 * self.XZk_add[k]
                M_0 = CMOV(M_0, tmp0, k <= self.sK-1)
                M_1 = CMOV(M_1, tmp1, k <= self.sK-1)
                
            l = self.L[i]
            l_fake = batchmaxprime_of_Li(i,self.batch_start,self.batch_stop,self.L)
            # B_0 = self.field(1); B_1 = self.field(1)
            # for i in range(l_fake):
            #     tmp1 = B_0 * A24[0]
            #     tmp2 = B_1 * (A24[0]- A24[1])
            #     B_0 = CMOV(B_0, tmp1, i <= l-1)
            #     B_1 = CMOV(B_1, tmp2, i <= l-1)
            
            B_0 = self.field(A24[0])
            B_1 = self.field(A24[0]- A24[1])
            B_0 = B_0.safe_pow(l, l_fake)
            B_1 = B_1.safe_pow(l, l_fake)
            
            
            C_1 = (M_1*R_1)**8
            C_0 = (M_0*R_0)**8    
            
            d_0 = B_0*C_1
            d_1 = B_1*C_0
            
            Az_new = d_0 - d_1
            Ax_new = 2*(d_0 + d_1)
            
            return (Ax_new , Az_new)
        
        def xeval_s(self, P: tuple, A24: tuple) -> tuple:
                
            B = 2*(2*A24[0] - A24[1])
            C = A24[1]
            
            XZ_add = P[0] + P[1]
            XZ_sub= P[0] - P[1]
            XZ2 = 2 * P[0] * P[1]
            X2Z2 = P[0]**2 + P[1]**2
            
            EJ_0 = [[0, 0, 0] for j in range(0, self.sJ , 1)]
            for j in range(0, self.sJ, 1):
                t1 = (XZ_sub* self.XZj_add[j])
                t2 = (XZ_add * self.XZj_sub[j])
                # print(f"t1:{t1} t2:{t2}")  
                quadratic = C*((t1 - t2) ** 2)
                constant = C*((t1 + t2) ** 2)  
                linear = 2*((C*X2Z2 * self.XZJ4[j]) - (C * XZ2 *( 2*(self.XZj_add[j])**2 + self.XZJ4[j])) + (B * XZ2 * self.XZJ4[j]))
                # print(f"qua:{quadratic}")
                # print(f"con :{constant}")
                # print(f"linear:{linear}")
                EJ_0[j] = [constant, linear, quadratic]
                
            poly_EJ_0 = self.poly_mul.product_tree(EJ_0, self.sJ)[
                2
            ]  
            poly_EJ_1 = list(
                poly_EJ_0[::-1]
            ) 
            
            # Approach using scaled remainder trees
            if self.ptree_hI != None:
                poly_EJ_0 = self.poly_redc.poly_redc(
                    poly_EJ_0, 2 * self.sJ + 1, self.ptree_hI
                )
                fg_0 = self.poly_mul.poly_mul_middle(
                    self.ptree_hI[7], self.sI, poly_EJ_0[::-1], self.sI
                )
                remainders_EJ_0 = self.poly_redc.multieval_scaled(
                    fg_0[::-1],
                    self.sI,
                    [1] + [0] * (self.sI - 1),
                    self.sI,
                    self.ptree_hI,
                    self.sI,
                )

                poly_EJ_1 = self.poly_redc.poly_redc(
                    poly_EJ_1, 2 * self.sJ + 1, self.ptree_hI
                )
                fg_1 = self.poly_mul.poly_mul_middle(
                    self.ptree_hI[7], self.sI, poly_EJ_1[::-1], self.sI
                )
                remainders_EJ_1 = self.poly_redc.multieval_scaled(
                    fg_1[::-1],
                    self.sI,
                    [1] + [0] * (self.sI - 1),
                    self.sI,
                    self.ptree_hI,
                    self.sI,
                )
            else:
                remainders_EJ_0 = []
                remainders_EJ_1 = []
                # print("*********")
                
            R_0 = self.poly_mul.product(remainders_EJ_0, self.sI)
            R_1 = self.poly_mul.product(remainders_EJ_1, self.sI)

            # print("PPP",R_0)
            # print(R_1)
            # hK_0 = [[0]] * self.sK
            # hK_1 = [[0]] * self.sK
            M_0 = self.field(1); M_1 = self.field(1)
            for k in range(0, self.sK_fake, 1):
                t1 = (XZ_sub* self.XZk_add[k])  
                t2 = (XZ_add * self.XZk_sub[k])  
                tmp0 = M_0 * (t1 + t2)
                tmp1 = M_1 * (t1 - t2)
                M_0 = CMOV(M_0, tmp0, k <= self.sK-1)
                M_1 = CMOV(M_1, tmp1, k <= self.sK-1)
                # hK_0[k] = [(t1 + t2)]
                # hK_1[k] = [(t1 - t2)] 
                
            # M_0 = self.poly_mul.product(
            #     hK_0, self.sK
            # )  
            # M_1 = self.poly_mul.product(
            #     hK_1, self.sK
            # ) 
            # for k in range(0, self.sK_fake, 1):
            #     tmp0 = M_0 * (t1 + t2)
            #     tmp1 = M_1 * (t1 - t2)
            #     M_0 = CMOV(M_0, tmp0, k <= self.sK-1)
            #     M_1 = CMOV(M_1, tmp1, k <= self.sK-1)
            # print(M_0)
            # print(M_1)
            XX = P[0]*((M_1*R_1)**2)         
            ZZ = P[1]*((M_0*R_0)**2)
            
            return [XX, ZZ]
        

    return Formulae
