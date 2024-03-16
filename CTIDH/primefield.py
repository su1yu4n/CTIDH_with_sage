from sage.all import GF, proof, is_prime
from sage.rings.finite_rings.integer_mod import IntegerMod_gmp, IntegerMod_int

from .utils import bitlength, hamming_weight, memoize, CMOV

import copy

proof.arithmetic(False)


@memoize
def PrimeField(p: int):
    if not is_prime(p):
        raise ArithmeticError("Cannot construct Fp: p is not a prime!")

    GFp = GF(p)

    # other can have type ZModPrime or int
    def get_value(other):
        if isinstance(other, ZModPrime):
            return other.value
        elif isinstance(other, int):
            return other
        else:
            raise TypeError(
                "Cannot get the value of (type:{}) {}!".format(type(other), other)
            )

    # TODO: Maybe rewrite get_value and operators with this decorator.
    # def check_other_type_and_return(f):
    #   def wrapper(x, other):
    #       if isinstance(other, ZModPrime):
    #           return f(x.value, other.value)
    #       elif isinstance(other, int):
    #           return f(x.value, other)
    #   return wrapper

    class ZModPrime:
        add_count = 0
        sqr_count = 0
        mul_count = 0

        # Note that computing a^0, a^1, a^2 won't increase pow_count
        # And reset_runtime() won't change pow_count and inv_count
        pow_count = 0
        inv_count = 0
        _p = p
        # p = p

        # self.value always has the type IntegerMod_gmp when p is large or IntegerMod_int
        def __init__(self, elem):
            if isinstance(elem, IntegerMod_gmp) or isinstance(elem, IntegerMod_int):
                self.value = elem
            elif isinstance(elem, ZModPrime):
                self.value = elem.value
            elif isinstance(elem, int):
                self.value = GFp(elem)
            else:
                raise TypeError(
                    "Cannot convert {} type {} to a ZModPrime!".format(type(elem), elem)
                )


        def __add__(self, other):
            ZModPrime.add_count += 1
            other = get_value(other)
            return ZModPrime(self.value + other)

        def __radd__(self, other):
            return self + other

        def __iadd__(self, other):
            ZModPrime.add_count += 1
            other = get_value(other)
            self.value = self.value + other
            return self

        def __sub__(self, other):
            ZModPrime.add_count += 1
            other = get_value(other)
            return ZModPrime(self.value - other)

        def __rsub__(self, other):
            return -self + other

        def __isub__(self, other):
            ZModPrime.add_count += 1
            other = get_value(other)
            self.value -= other
            return self

        def __mul__(self, other):
            ZModPrime.mul_count += 1
            other = get_value(other)
            return ZModPrime(self.value * other)

        def __rmul__(self, other):
            return self * other

        def __imul__(self, other):
            ZModPrime.mul_count += 1
            other = get_value(other)
            self.value = self.value * other
            return self

        # TODO: Maybe implment div rdiv and idiv?
        def __div__(self, other):
            raise NotImplementedError

        def __rdiv__(self, other):
            raise NotImplementedError

        def __idiv__(self, other):
            raise NotImplementedError

        def __pow__(self, e: int):
            # TODO: write a faster constant-time power if it is slow
            # Can we use a nearly optimal addition chain here?
            # See https://en.wikipedia.org/wiki/Exponentiation_by_squaring#Alternatives_and_generalizations
            """
            Exponentiation

            ...
            Parameters
            ----------
                    - self, which is an element of a Prime Field
                    - an integer e
            Returns
            -------
                    - self raised to e
            -----
            Usage:
                    - self.pow(e)
                    - self ** e
            Notes
            -----
                    - The Fp-operations are counted correspond to the cost of left-to-right method.
                    - It allows negative exponents, but any exponent is expected to belong to |[ 0 .. p - 1 ]|
            """
            if e == 0:
                return ZModPrime(1)
            
            elif e > 0: # e > 0
                ZModPrime.sqr_count += bitlength(e) - 1
                ZModPrime.mul_count += hamming_weight(e) - 1
                if e > 2: # e > 2
                    ZModPrime.pow_count += 1
                return ZModPrime(self.value ** e)

            # Seems that this is the only case when e<0 in CSIDH.
            elif e == -1:
                # ~ indicate invert
                return ~self
            
            else:
                raise NotImplementedError('Unexpected behavior: performing a ** e, e < -1 in CTIDH.')

        def __ipow__(self, e: int):
            if e > 0:
                ZModPrime.sqr_count += bitlength(e) - 1
                # hamming weight of 1 and 2 is 1, subtracting 1
                # means we're adding 0 to the counter.
                # by special-casing that here we can avoid
                # ~300k function calls per CSIDH:
                if e > 2:
                    ZModPrime.pow_count += 1
                    ZModPrime.mul_count += hamming_weight(e) - 1
                self.value **= e
                return self
            
            elif e == 0:  # e == 0
                self.value = 1
                return self
            
            else: # e < 0
                raise NotImplementedError('Unexpected behavior: performing a **= e, e < 0 in CTIDH.')

        def safe_pow(self, e:int, e_maxbitlen:int):
            """Timing-safe powmod. Return a ZModPrime with value self.value ** e

            Args:
                e (int): the exponent
                e_maxbitlen (int): the possible max value of e's bitlength. 

            NOTE: To achieve timing-safe property, this function pretends to compute the case exponent equals 2**e_maxbitlen

            """
            ZModPrime.pow_count += 1
            
            a = self
            tmp1 = copy.deepcopy(self) # a ** padded_e
            tmp2 = copy.deepcopy(self) 
            e_bitlen = bitlength(e)

            padded_e = e << (e_maxbitlen - e_bitlen)
            ans = copy.deepcopy(self)
            # tmp and ans have initial value a
            # start from the second MSB
            for i in range(2, e_maxbitlen+1):
                tmp1 = tmp1 ** 2
                tmp2 = tmp1 * a
                should_mov = (padded_e >> (e_maxbitlen - i)) & 1 # if 1 then mov, if 0 then don't move
                tmp1 = CMOV(tmp1, tmp2, should_mov)

                ans = CMOV(ans, tmp1, i <= e_bitlen)
            
            return ans
                
        def __invert__(self):
            # TODO: write a faster constant-time invert.
            # Currently we use self.x**(p-2).
            self.inv_count += 1
            return self ** (ZModPrime._p - 2)

        def __neg__(self):
            return ZModPrime(-self.value)

        def __eq__(self, other):
            other = get_value(other)
            return self.value == other
        
        def __str__(self):
            return 'ZModPrime {} mod {}'.format(self.value, prime_name)
        
        def __repr__(self):
            return str(self)

        def get_int_value(self):
            return int(self.value)

        def copy(self):
            ret = object.__new__(ZModPrime)
            ret.value = self.value
            return ret

        def is_square(self) -> bool:
            legendre_symbol = self ** ( (ZModPrime._p - 1) // 2)
            return True if legendre_symbol == 1 else False

        def is_square_fast(self) -> bool:
            # TODO: Implement this CCS' 23 algorithm.
            raise NotImplementedError('Fast algorithm not implemented yet!')

        @classmethod
        def get_random(cls):
            return ZModPrime(GFp.random_element())


        @classmethod
        def reset_runtime(cls):
            cls.add_count = 0
            cls.sqr_count = 0
            cls.mul_count = 0

        @classmethod
        def reset_power_invert_time(cls):
            cls.pow_count = 0
            cls.inv_count = 0

        @classmethod
        def show_runtime(cls, label: str):
            print(
                "| %s: %7dM + %7dS + %7da"
                % (label, cls.mul_count, cls.sqr_count, cls.add_count),
                end="\n",
            )

        @classmethod
        def show_sqr_pow(cls, label: str):
            print(
                "| %s: %2dP + %2dI" % (label, cls.pow_count, cls.inv_count),
                end="\n",
            )

    prime_name = 'p1024_CTIDH' if p.bit_length() <= 1025 else 'p2048_CTIDH'
    ZModPrime.__name__ = f'ZModPrime with p = {prime_name}'

    return ZModPrime
