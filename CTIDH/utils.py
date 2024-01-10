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
    with open(f'../parameters/{prime_name}') as f:
        return json.load(f)


# Dictionary which provides attribute access to its keys.
class attrdict(dict):
    __getattr__ = dict.__getitem__