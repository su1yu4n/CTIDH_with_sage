# from functools import cache

import json

from Crypto.Util.number import sieve_base

from sage.all import product

# @cache
def gen_prime_info(n, k, excluded=[], included=[]):
    """Return a dict like
    {
        'L': list of small odd primes(SOPs)
        'n': total number of small primes
        'k': see below
        'f': 2**k, which is the cofactor
        'p': value of p. p = 2**k * L[0] * ... * L[n-1] - 1
    }

    Args:
        n (int): see n above
        k (int): see p above
        excluded (list): excluded primes
        included (list): extra included primes
    """

    
    li_list = list(sieve_base[1 : n + 2 - len(included)])
    for excluded in excluded:
        li_list.remove(excluded)

    for extra in included:
        li_list.append(extra)

    f = 2**k
    p = f * product(li_list) - 1

    return {"L": li_list, "n": n, "k": k, "f": f, "p": p}


# TODO: Check p1024 and p2048 in original CTIDH
primes = dict(
    p1024_CTIDH= gen_prime_info(n=130, k=2, excluded=[739], included=[983]),
    p2048_CTIDH= gen_prime_info(n=231, k=2, excluded=[5], included=[3413]),
)


# TODO: get batch information and keyspace
# TODO: Compute fine-tuned Velusqrt parmeters and shortest addition chains

def write_primes_to_file(primes:dict):
    for prime_name in primes.keys():
        with open(f'parameters/{prime_name}', 'w') as f:
            json.dump(primes[prime_name], f)

write_primes_to_file(primes=primes)