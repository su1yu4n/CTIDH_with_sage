import json
import random
from datetime import datetime
from functools import reduce

# Dictionary which provides attribute access to its keys.
class attrdict(dict):
    __getattr__ = dict.__getitem__


def sign(x):
    return (1, -1)[x < 0]  # Sign of an integer

isequal = {True: 1, False: 0}  # Simulating constant-time integer comparison

# memoize calls to the class constructors for fields
# this helps typechecking by never creating two separate
# instances of a number class.
def memoize(f):
    cache = {}

    def memoizedFunction(*args, **kwargs):
        argTuple = args + tuple(kwargs)
        if argTuple not in cache:
            cache[argTuple] = f(*args, **kwargs)
        return cache[argTuple]

    memoizedFunction.cache = cache
    return memoizedFunction


# number of bits, use builtin int.bit_length if present:
bitlength = getattr(int, "bit_length", lambda x: len(bin(x)[2:]))

# python3.10 has builtin popcount aka hamming weight aka int.bit_count:
hamming_weight = getattr(int, "bit_count", lambda x: bin(x).count(r"1"))
# hamming weight: number of bits equal 1


# Return k0, k1, ..., kn
def binrep(k: int) -> list:
    bin_rep = []
    for _ in range(k.bit_length()):
        bin_rep.append(k & 1)
        k = k >> 1
    return bin_rep


# conditional MOV, return a
# Note that this PoC implementation does not implement CMOV and CSWAP carefully,
# Maybe it should be implement in another way to make sure it safe.
def CMOV(a, b, control: bool):
    return b if control else a


def CSWAP(a, b, control: bool):
    if control:
        a, b = b, a
    return a, b


def read_prime_info(prime_name="p2048_CTIDH"):
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
        'batch_start': the start index of each batch. eg. batch_start[0] = 0
        'batch_stop': the end index + 1 of each batch. eg. batch_stop[-1] = len(L) 
        'batch_maxdaclen': the max sdac length of primes in each batch
        'batch_bound': keyspace info: the max value of absolute values' sum for each batch
    }
    """
    with open(f"data/prime_info/{prime_name}") as f:
        return json.load(f)


def read_velusqrt_steps_info(prime_name="p2048_CTIDH", scaled=True):
    """
    Args:
        prime_name (str, optional): "p1024_CTIDH" or "p2048_CTIDH". Defaults to "p2048_CTIDH".
        scaled (bool, optional): Use scaled remainder tree or not. Defaults to True.

    Returns:
        sI_list, sJ_list: baby steps info, giant steps info resp.
    """
    suffix = '_scaled' if scaled else '_unscaled'

    sJ_list = []
    sI_list = []

    steps_info_path = f"data/prime_info/{prime_name}" + suffix
    with open(steps_info_path, 'r') as f:
        steps = f.readlines()

    if prime_name == 'p1024_CTIDH':
        assert len(steps) == 130
    elif prime_name == 'p2048_CTIDH':
        assert len(steps) == 231
    # else:
    #     raise NotImplementedError()
    
    for step in steps:
        bs, gs = step.split()
        sJ_list.append(int(bs))
        sI_list.append(int(gs))

    return sI_list, sJ_list


def read_SDAC_info(prime_name='p2048_CTIDH'):
    SDAC_info = []
    with open(f"data/sdacs/{prime_name}", 'r') as f:
        for line in f:
            SDAC_info.append(list(map(int, line.split())))
    return SDAC_info


# return the batch number of the batch that L[i] belongs to.
def batchnumber_of_Li(i:int, batch_start: list, batch_stop: list):
    for j in range(len(batch_start)):
        if batch_start[j] <= i < batch_stop[j]:
            return j
        
# return the max prime in the same batch as L[i]
def batchmaxprime_of_Li(i: int, batch_start: list, batch_stop: list, L:list):
    batch_num = batchnumber_of_Li(i, batch_start, batch_stop)
    batchmaxprime_index = batch_stop[batch_num] - 1 # 
    return L[batchmaxprime_index] 


# return the min prime in the same batch as L[i]
def batchminprime_of_Li(i: int, batch_start: list, batch_stop: list, L:list):
    batch_num = batchnumber_of_Li(i, batch_start, batch_stop)
    batchminprime_index = batch_start[batch_num]
    return L[batchminprime_index] 


# for random tests
def get_randint(a: int, b: int):
    random.seed(datetime.now().strftime('%H:%M:%S'))
    return random.randint(a, b)
