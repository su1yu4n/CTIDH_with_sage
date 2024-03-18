from CTIDH.isogeny import MontgomeryIsogeny
from CTIDH.mont import MontgomeryCurve
from CTIDH.utils import read_prime_info

class CSIDH:
    def __init__(
        self,
        prime_name,
        formula_name='hvelu', # Choose tvelu/hvelu/svelu. This is an object returned by MontgomeryIsogeny
        scaled=False, #  use scaled remainder tree if True. always True in sibc
        use_SDAC=False, # use SDAC scalar mult?
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
        self.use_SDAC = use_SDAC
        self.tuned = tuned
        self.scaled = scaled
        self.uninitialized = uninitialized
        self.verbose = verbose
        self.validation=validation

        self.isogeny = MontgomeryIsogeny(formula_name, uninitialized)(
            self.curve,
            self.tuned,
            self.scaled
        )
        self.curve = MontgomeryCurve(prime_name, use_SDAC, validation)
        self.field = self.curve.field

        self.prime_info = read_prime_info(prime_name)

    def keygen(self) -> tuple:
        # gen private key
        batch_bound = self.prime_info['batch_bound']
        batch_start = self.prime_info['batch_start']
        batch_stop = self.prime_info['batch_stop']
        # TODO
        # gen public key

        raise NotImplementedError
    
    def derive(self, sk: list, pk: int) -> int:
        if not self.is_supersingular(pk):
            raise ValueError(f'the public key: {pk} is not supersingular!')
        return self.group_action(pk, sk)
    

    def group_action(self, a, exponents) -> int:
        raise NotImplementedError
    
    def is_supersingular(self, pk:int) -> bool:
        A = (self.field(pk), self.field(1))
        return self.curve.issupersingular(A)
        
        


