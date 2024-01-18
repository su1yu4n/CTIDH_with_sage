from sage.all import GF, proof, is_prime
from sage.rings.finite_rings.integer_mod import IntegerMod_gmp, IntegerMod_int
from utils import bitlength, hamming_weight

proof.arithmetic(False)


def PrimeField(p: int):
    if not is_prime(p):
        raise ArithmeticError("Cannot construct Fp: p is not a prime!")

    GFp = GF(p)

    # other can have type ZModPrime or int
    def get_value(other):
        if isinstance(other, ZModPrime):
            return other.x
        elif isinstance(other, int):
            return other
        else:
            raise TypeError(
                "Cannot get the value of (type:{}) {}!".format(type(other), other)
            )

    # TODO: Maybe rewrite get_value and operators with this decorator. Just for touching fish ×
    # def check_other_type_and_return(f):
    #   def wrapper(x, other):
    #       if isinstance(other, ZModPrime):
    #           return f(x.x, other.x)
    #       elif isinstance(other, int):
    #           return f(x.x, other)
    #   return wrapper

    class ZModPrime:
        add_count = 0
        sqr_count = 0
        mul_count = 0
        pow_count = 0
        inv_count = 0
        p = p

        # self.x always has the type IntegerMod_gmp when p is large or IntegerMod_int
        def __init__(self, x):
            if isinstance(x, IntegerMod_gmp) or isinstance(x, IntegerMod_int):
                self.x = x
            elif isinstance(x, int):
                self.x = GFp(x)
            else:
                raise TypeError(
                    "Cannot convert {} type {} to a ZModPrime!".format(type(x), x)
                )

        def __add__(self, other):
            ZModPrime.add_count += 1
            other = get_value(other)
            return ZModPrime(self.x + other)

        def __radd__(self, other):
            return self + other

        def __iadd__(self, other):
            ZModPrime.add_count += 1
            other = get_value(other)
            self.x = self.x + other
            return self

        def __sub__(self, other):
            ZModPrime.add_count += 1
            other = get_value(other)
            return ZModPrime(self.x - other)

        def __rsub__(self, other):
            return -self + other

        def __isub__(self, other):
            ZModPrime.add_count += 1
            other = get_value(other)
            self.x -= other
            return self

        def __mul__(self, other):
            ZModPrime.mul_count += 1
            other = get_value(other)
            return ZModPrime(self.x * other)

        def __rmul__(self, other):
            return self * other

        def __imul__(self, other):
            ZModPrime.mul_count += 1
            other = get_value(other)
            self.x = self.x * other
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
                    - This is a constant-time implementation by using the left-to-right method
                    - It allows negative exponents, but any exponent is expected to belong to |[ 0 .. p - 1 ]|
            """
            if e == 0:
                return ZModPrime(1)
            
            elif e > 0: # e > 0
                ZModPrime.sqr_count += bitlength(e) - 1
                ZModPrime.mul_count += hamming_weight(e) - 1
                return ZModPrime(self.x ** e)

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
                    ZModPrime.mul_count += hamming_weight(e) - 1
                self.x **= e
                return self
            
            elif e == 0:  # e == 0
                self.x = 1
                return self
            
            else: # e < 0
                raise NotImplementedError('Unexpected behavior: performing a **= e, e < 0 in CTIDH.')

        def __invert__(self):
            # TODO: write a faster constant-time invert.
            # Currently we use self.x**(p-2).
            return self ** (ZModPrime.p - 2)

        def __neg__(self):
            return ZModPrime(-self.x)

        def __eq__(self, other):
            other = get_value(other)
            return self.x == other

        def copy(self):
            ret = object.__new__(ZModPrime)
            ret.x = self.x
            return ret

        @classmethod
        def reset_runtime(cls):
            cls.add_count = 0
            cls.sqr_count = 0
            cls.mul_count = 0

        @classmethod
        def show_runtime(cls, label: str):
            print(
                "| %s: %7dM + %7dS + %7da"
                % (label, cls.mul_count, cls.sqr_count, cls.add_count),
                end="\t",
            )

        @classmethod
        def show_sqr_pow(cls, label: str):
            print(
                "| %s: %2dP + %2dI" % (label, cls.pow_count, cls.inv_count),
                end="\t",
            )

    return ZModPrime
