from sage.all import GF, proof, is_prime
from sage.rings.finite_rings.integer_mod import IntegerMod_gmp, IntegerMod_int

proof.arithmetic(False)

def reset_runtime(field):
	field.add_count = 0
	field.sqr_count = 0
	field.mul_count = 0
    

def PrimeField(p: int):
    if not is_prime(p):
        raise ArithmeticError('Cannot construct Fp: p is not a prime!')
    
    GFp = GF(p)

    # other can have type ZModPrime or int
    def get_value(other):
        if isinstance(other, ZModPrime):
            return other.x
        elif isinstance(other, int):
            return other
        else:
            raise TypeError('Cannot get the value of (type:{}) {}!'.format(type(other), other))


    class ZModPrime:            
        # self.x always has the type IntegerMod_gmp when p is large or IntegerMod_int
        def __init__(self, x):
            if isinstance(x, IntegerMod_gmp) or isinstance(x, IntegerMod_int):
                self.x = x
            elif isinstance(x, int):
                self.x = GFp(x)
            else:
                raise TypeError('Cannot convert {} type {} to a ZModPrime!'.format(type(x), x))

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
        
        # TODO: Maybe add div rdiv and idiv?
        def __div__(self, other):
            raise NotImplementedError
        
        def __rdiv__(self, other):
            raise NotImplementedError
        
        def __idiv__(self, other):
            raise NotImplementedError


        def __pow__(self, other):
            # TODO: write a fast constant-time power if it is slow
            raise NotImplementedError

        def __invert__(self, other):
            # TODO: write a fast constant-time invert if it is slow
            raise NotImplementedError
        
        def __neg__(self): 
            return ZModPrime(-self.x)

        def __eq__(self, other):
            other = get_value(other)
            return self.x == other

        def copy(self):
            ret = object.__new__(ZModPrime)
            ret.x = self.x
            return ret
        

    ZModPrime.add_count = 0
    ZModPrime.sqr_count = 0
    ZModPrime.mul_count = 0

    ZModPrime.pow_count = 0
    ZModPrime.inv_count = 0

    ZModPrime.show_runtime = lambda label='Since last init_runtime': print(
        "| %s: %7dM + %7dS + %7da"
        % (label, ZModPrime.mul_count, ZModPrime.sqr_count, ZModPrime.add_count),
        end="\t",
    )

    ZModPrime.reset_runtime = lambda: reset_runtime(ZModPrime)

    ZModPrime.show_sqr_pow = lambda label='The whole process': print(
        "| %s: %2dP + %2dI"
        % (label, ZModPrime.pow_count, ZModPrime.inv_count),
        end="\t",
    )

    return ZModPrime


