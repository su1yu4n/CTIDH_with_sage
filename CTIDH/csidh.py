from CTIDH.isogeny import MontgomeryIsogeny
from CTIDH.mont import MontgomeryCurve
from CTIDH.utils import read_prime_info

from CTIDH.utils import get_randint

class CSIDH:
    def __init__(
        self,
        prime_name,
        formula_name='hvelu', # Choose tvelu/hvelu/svelu. This is an object returned by MontgomeryIsogeny
        scaled=False, #  use scaled remainder tree if True. always True in sibc
        SDAC=False, # use SDAC scalar mult?
        tuned=True, # velusqrt is tuned? Not needed (always fine-tuned)
        verbose=True,
        uninitialized=False,
        validation='original',
    ):
        # Check parameters
        if formula_name not in ['tvelu', 'svelu', 'hvelu']:
            raise ValueError(f'No such formula: {formula_name}')
        if validation not in ['original', 'doliskani', 'pairing1', 'pairing2']:
            raise ValueError(f'No such validation algorithm: {validation}')

        self.prime_name = prime_name
        self.formula_name = formula_name
        self.SDAC = SDAC
        self.tuned = tuned
        self.scaled = scaled
        self.uninitialized = uninitialized
        self.verbose = verbose
        self.validation=validation

        self.curve = MontgomeryCurve(prime_name, SDAC, validation)
        self.isogeny = MontgomeryIsogeny(formula_name, uninitialized)(
            self.curve,
            self.tuned,
            self.scaled
        )
        self.field = self.curve.field

        self.prime_info = read_prime_info(prime_name)

    def keygen(self) -> tuple:
        """Generate public key and private key.

        Returns:
            tuple: sk(list), pk(int)
        """
        sk = self.skgen()
        pk = self.group_action(0, sk)
        return sk, pk
    
    def skgen(self):
        def random_boundedl1(Ni, mi, b=32):
            if mi == 0:
                return [0] * Ni
            
            rnum = Ni+mi
            while True:
                # get random Ni+mi b-bit ints
                r = [get_randint(-2**(b-1), 2**(b-1) - 1) for _ in range(rnum)]
                # step2 and step3
                for j in range(0, rnum):
                    r[j] &= ~1
                for j in range(0, Ni):
                    r[j] |= 1
                r.sort() # NOTE: use constant-time sort in practice
                # step4: if any adjacent ints are the same outside the bottom bit, start over
                collision = 0
                # NOTE: strangely, in original CTIDH's code the following for loop condition is j<Ni, not rnum.
                for j in range(0, rnum-1):
                    collision |=  ((r[j] ^ r[j+1] & ~1) == 0)
                    # NOTE: ^ if x==0 , then ~( (x>>(b-1)) | (-x>>(b-1))) == -1. otherwise it is 0. x==0 should be implemented like this.
                if collision: # if collision == 1
                    continue
                for j in range(0, rnum):
                    r[j] &= 1
                for j in range(1, rnum):
                    r[j] += r[j-1]
                
                e = [0] * Ni
                for i in range(0, Ni):
                    numi = 0
                    for j in range(0, rnum):
                        numi += (r[j] == i) # this == should also be done in constant time using bit operation
                    e[i] = numi
                
                for i in range(1, Ni):
                    e[i] -= 1
                
                reject = 0
                s = get_randint(-2**(Ni-1), 2**(Ni-1)-1)
                for i in range(0, Ni):
                    if s&1 == 1:
                        e[i] = -e[i]
                        reject |= (e[i] == 0)
                    s >>= 1
                if reject:
                    continue
                return e
        
        batch_bound = self.prime_info['batch_bound']
        batch_start = self.prime_info['batch_start']
        batch_stop = self.prime_info['batch_stop']
        batch_num = len(batch_start)
        sk = []
        for i in range(batch_num):
            batchlen = batch_stop[i] - batch_start[i]
            batch_ei = random_boundedl1(batchlen, batch_bound[i])
            sk += batch_ei
        return sk
    

    def derive(self, sk: list, pk: int) -> int:
        if not self.is_supersingular(pk):
            raise ValueError(f'the public key: {pk} is not supersingular!')
        return self.group_action(pk, sk)
    

    def group_action(self, a:int, exponents:list) -> int:
        raise NotImplementedError
    

    def is_supersingular(self, pk:int) -> bool:
        A = (self.field(pk), self.field(1))
        return self.curve.issupersingular(A)
        
        


