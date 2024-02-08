from typing import List

from CTIDH.mont import MontgomeryCurve
from CTIDH.utils import read_velusqrt_steps_info, hamming_weight, bitlength, isequal, batchmaxprime_of_Li, batchminprime_of_Li, batchnumber_of_Li
from CTIDH.polymul import PolyMul
from CTIDH.polyredc import PolyRedc

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
    cutoff = 83
    cutoff_string = f' with cutoff ell <= {cutoff}' if formula_name == 'hvelu' else ''
    NAME = 'Isogeny class using the %s Velu\'s formulae%s' % ({'tvelu':'traditional', 'svelu':'square-root', 'hvelu':'hybrid'}[formula_name], cutoff_string)

    @doc(NAME)
    class Formulae:
        def __init__(self, curve, tuned=True, scaled=True):
            """_summary_

            Args:
                curve (MontgomeryCurve): The object returned by MontgomeryCurve() in mont.py
                tuned (bool, optional): True if fine-tuned velusqrt information is presented in data folder. Defaults to True.
                scaled (bool, optional): Use scaled remainder tree or not. 
                If True, it will read and use velusqrt tuned info for scaled version. Defaults to True.
            """
            self.formula_name = formula_name

            if tuned:
                self.sI_list, self.sJ_list = read_velusqrt_steps_info(curve.prime_name, scaled)


            self.HYBRID_BOUND = {'tvelu':max(curve.L), 'svelu':1, 'hvelu':cutoff}[formula_name]

            # Global variables to be used in kps, xisog, and xeval

            # Here, J is a set of cardinality sJ
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
            self.ADD_SQUARED = None
            self.SUB_SQUARED = None

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
            

            # Now, we proceed to store all the correct costs
            if formula_name != 'tvelu' and uninitialized:
                print("// Precomputation cost regarding kps, xisog, and xeval of the velusqrt formulae")
                self.velusqrt_cost()
                

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
            assert 0 <= Tnewlen <= 2
            assert i < len(self.L)
            
            l = self.L[i]
            l_fake = batchmaxprime_of_Li(i, self.batch_start, self.batch_stop, self.L)
            A24 = self.curve.xA24(A)
            # NOTE: Here we use l_fake to decide which formula to use(traditional velu or velusqrt) to avoid timing attack.
            # It may be unsafe if different primes in the same batch use different formulae.
            if l_fake <= self.HYBRID_BOUND:     # Use Velu formula
                xiP_list = self.kps_t(l_fake, P, A24)
                Xi_Zi_hats = []
                for xiP in xiP_list:
                    Xi, Zi = xiP
                    Xi_Zi_hat = (Xi + Zi, Xi - Zi)
                    Xi_Zi_hats.append(Xi_Zi_hat)

                A_new = self.xisog_t(l, l_fake, Xi_Zi_hats, A)

                if Tnewlen > 0:
                    Ts[0] = self.xeval_t(l, l_fake, Xi_Zi_hats, Ts[0])
                if Tnewlen > 1:
                    Ts[1] = self.xeval_t(l, l_fake, Xi_Zi_hats, Ts[1])

                return A_new, Ts
            
            else:    # Use Velusqrt formulas
                self.set_parameters_velu(self.sJ_list[i], self.sI_list[i], i)
                raise NotImplementedError("matryoshka isogeny of velusqrt not implemented yet!")

        
        # TODO: Complete traditional velu formula
        def kps_t(self, l_fake: int, P: tuple, A24: tuple) -> List[tuple]:
            """Timing attack safe kps for traditional velu formula, 
            used in computing the l-isogeny phi: E -> E/<P>.

            Return the list of x([i]P) for i from 1 to (l_fake - 1) // 2

            Args:
                l_fake (int): the largest small odd prime in the batch of l
                P (tuple): projective x-coordinate of the kernel point P (i.e. x(P))
                A24 (tuple): A24 = (Ax+2Az : 4Az)
            """
            raise NotImplementedError
        

        def xisog_t(self, l: int, l_fake: int, Xi_Zi_hats: List[tuple], A: tuple) -> tuple:
            """Timing attack safe xisog for traditional velu formula,
            Return the fraction representation of quadratic term's coefficient of the l-isogeny's codomain. 
            
            Args:
                l (int): See above
                l_fake (int): the largest small odd prime in the batch of l
                Xi_Zi_hats (List[tuple]): list of (Xi+Zi : Xi-Zi) where (Xi : Zi) is x([i]P), P a generator of ker phi. 
                A (tuple): A = (Ax: Az). Ax and Az must have type ZModPrime(Primefield).
                Ax/Az is the quadratic term's coefficient of domain curve's affine equation. 

            Returns:
                tuple: A' = (Ax': Az'), where Ax'/Az' is the quadratic term's coefficient of the codomain.
                Ax', Az' have type ZModPrime
            """
            raise NotImplementedError
        

        def xeval_t(self, l: int, l_fake: int, Xi_Zi_hats: List[tuple], T: tuple) -> tuple:
            """Push the point T to the codomain through the l-isogeny phi.

            Args:
                l (int): degree of isogeny
                l_fake (int): the largest prime in the same batch as l.
                Xi_Zi_hats (List[tuple]): (Xi+Zi : Xi-Zi) where (Xi:Zi) x[i]P, i = 1, ..., d_fake
                T (tuple): (projective x-coordinate of) the point to push 

            Returns:
                tuple: (projective x-coordinate of) the image point phi(T)
            """
            raise NotImplementedError


        # NOTE: This functions is used for setting the cardinalities sI, sJ, and sK
        # In sibc it is called by velusqrt_cost()
        def set_parameters_velu(self, b, c, i):
            assert b <= c
            # At this step, everythin is correct
            self.sJ = b
            self.sI = c
            d = ((self.L[i] - 2 - 4 * b * c - 1) // 2) + 1
            assert d >= 0
            self.sK = d
            return None

        # TODO: Implement these velusqrt algorithms
        def kps_s(self, P: tuple, A: tuple, i: int):
            raise NotImplementedError
        
        def xisog_s(self, A, i):
            raise NotImplementedError
        
        def xeval_s(self, P, A):
            raise NotImplementedError