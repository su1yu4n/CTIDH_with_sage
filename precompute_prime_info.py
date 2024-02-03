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


primes = dict(
    p1024_CTIDH = gen_prime_info(n=130, k=2, excluded=[739], included=[983]),
    p2048_CTIDH = gen_prime_info(n=231, k=2, excluded=[5], included=[3413]),
)


# Add batch, keyspace, sdac max length information 
p1024_CTIDH_batch_info = {
    'batch_start': [0, 2, 5, 10, 14, 20, 26, 32, 38, 44, 51, 58, 65, 71, 78, 85, 90, 96, 101, 111, 114, 124, 129],
    'batch_stop': [2, 5, 10, 14, 20, 26, 32, 38, 44, 51, 58, 65, 71, 78, 85, 90, 96, 101, 111, 114, 124, 129, 130],
    'batch_maxdaclen': [1, 3, 5, 6, 7, 8, 9, 9, 10, 10, 10, 10, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 13],
    'batch_bound': [2, 4, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 3, 6, 2, 6, 2, 0],
}

primes["p1024_CTIDH"].update(p1024_CTIDH_batch_info)

p2048_CTIDH_batch_info = {
    'batch_start': [0, 9, 19, 27, 35, 42, 52, 64, 75, 85, 100, 110, 119, 127, 133, 143, 156, 166, 175, 187, 200, 210, 220, 230],
    'batch_stop': [9, 19, 27, 35, 42, 52, 64, 75, 85, 100, 110, 119, 127, 133, 143, 156, 166, 175, 187, 200, 210, 220, 230, 231],
    'batch_maxdaclen': [5, 7, 8, 9, 9, 10, 10, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 16],
    'batch_bound': [1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 0, 2, 0],
}

primes['p2048_CTIDH'].update(p2048_CTIDH_batch_info)


def write_primes_to_file(primes:dict):
    for prime_name in primes.keys():
        with open(f'data/prime_info/{prime_name}', 'w') as f:
            json.dump(primes[prime_name], f)

write_primes_to_file(primes=primes)