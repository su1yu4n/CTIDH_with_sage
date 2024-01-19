import json


def read_prime_info(prime_name='p2048_CTIDH'):
    """_summary_

    Args:
        prime_name (str, optional): Defaults to 'p2048_CTIDH'.

    Returns:
        prime info : a dict like
    {
        'L': list of small odd primes(SOPs)
        'n': total number of small primes
        'k': see below
        'f': 2**k, which is the cofactor
        'p': value of p. p = 2**k * L[0] * ... * L[n-1] - 1
    }
    """
    with open(f'parameters/{prime_name}') as f:
        return json.load(f)


# Dictionary which provides attribute access to its keys.
class attrdict(dict):
    __getattr__ = dict.__getitem__


# number of bits, use builtin int.bit_length if present:
bitlength = getattr(int, 'bit_length', lambda x: len(bin(x)[2:]))

# python3.10 has builtin popcount aka hamming weight aka int.bit_count:
hamming_weight = getattr(int, 'bit_count', lambda x: bin(x).count(r'1'))
# hamming weight: number of bits equal 1


# conditional MOV, return a
# Note that this PoC implementation does not implement CMOV and CSWAP carefully,
# Maybe it should be implement in another way to make sure it safe.

def CMOV(a, b, control:bool):
    a = b if control else a
    return a

def CSWAP(a, b, control:bool):
    if control:
        a, b = b, a
    return a, b