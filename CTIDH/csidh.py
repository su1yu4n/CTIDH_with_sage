from copy import deepcopy

from sage.misc.verbose import verbose, set_verbose
from CTIDH.isogeny import MontgomeryIsogeny
from CTIDH.mont import MontgomeryCurve
from CTIDH.utils import read_prime_info, read_prime_info_for_tvelu_test

from CTIDH.utils import get_randint, CMOV, CSWAP, sign


class CSIDH:
    def __init__(
        self,
        prime_name,
        formula_name="hvelu",  # Choose tvelu/hvelu/svelu. This is an object returned by MontgomeryIsogeny
        scaled=False,  #  use scaled remainder tree if True. always True in sibc
        SDAC=False,  # use SDAC scalar mult?
        tuned=True,  # velusqrt is tuned? Not needed (always fine-tuned)
        verbose=True,
        uninitialized=False,
        validation="original",
        fast_kronecker=False,
    ):
        # Check parameters
        if formula_name not in ["tvelu", "svelu", "hvelu"]:
            raise ValueError(f"No such formula: {formula_name}")
        if validation not in ["original", "doliskani", "pairing1", "pairing2"]:
            raise ValueError(f"No such validation algorithm: {validation}")

        self.prime_name = prime_name
        self.formula_name = formula_name
        self.SDAC = SDAC
        self.tuned = tuned
        self.scaled = scaled
        self.uninitialized = uninitialized
        self.verbose = verbose
        self.validation = validation
        self.fast_kronecker = fast_kronecker

        self.curve = MontgomeryCurve(prime_name, SDAC, validation, fast_kronecker)
        self.isogeny = MontgomeryIsogeny(formula_name, uninitialized)(
            self.curve, self.tuned, self.scaled
        )
        self.field = self.curve.field

        # if formula_name == 'tvelu':
        #     # self.prime_info = read_prime_info_for_tvelu_test(prime_name)
        #     self.prime_info = read_prime_info(prime_name)
        # else:
        #     self.prime_info = read_prime_info(prime_name)
        self.prime_info = read_prime_info(prime_name)

    def keygen(self) -> tuple:
        """Generate public key and private key.

        Returns:
            tuple: sk(list), pk(int)
        """
        sk = self.skgen()
        pk = self.group_action(0, sk)
        return sk, pk

    def random_boundedl1(self, Ni, mi, b=32):
        if mi == 0:
            return [0] * Ni

        rnum = Ni + mi
        while True:
            # get random Ni+mi b-bit ints
            r = [get_randint(-(2 ** (b - 1)), 2 ** (b - 1) - 1) for _ in range(rnum)]
            # step2 and step3
            for j in range(0, rnum):
                r[j] &= ~1
            for j in range(0, Ni):
                r[j] |= 1
            r.sort()  # NOTE: use constant-time sort in practice
            # step4: if any adjacent ints are the same outside the bottom bit, start over
            collision = 0
            # NOTE: strangely, in original CTIDH's code the following for loop condition is j<Ni, not rnum.
            for j in range(0, rnum - 1):
                collision |= (r[j] ^ r[j + 1] & ~1) == 0
                # NOTE: ^ if x==0 , then ~( (x>>(b-1)) | (-x>>(b-1))) == -1. otherwise it is 0. x==0 should be implemented like this.
            if collision:  # if collision == 1
                continue
            for j in range(0, rnum):
                r[j] &= 1
            for j in range(1, rnum):
                r[j] += r[j - 1]

            e = [0] * Ni
            for i in range(0, Ni):
                numi = 0
                for j in range(0, rnum):
                    numi += (
                        r[j] == i
                    )  # this == should also be done in constant time using bit operation
                e[i] = numi

            for i in range(1, Ni):
                e[i] -= 1

            reject = 0
            # NOTE: This is the setting of original CTIDH, but I think the bias is a little big
            counter = Ni - mi
            s = get_randint(-(2 ** (Ni - 1)), 2 ** (Ni - 1) - 1)
            for i in range(0, Ni):
                ei_zero_mask = e[i] == 0
                counter -= ei_zero_mask
                sbit = s & 1
                if sbit:
                    e[i] = -e[i]
                reject |= (counter < 0) & ei_zero_mask & sbit
                s >>= 1
            if reject:
                continue
            return e

    def skgen(self):
        batch_bound = self.prime_info["batch_bound"]
        batch_start = self.prime_info["batch_start"]
        batch_stop = self.prime_info["batch_stop"]
        batch_num = len(batch_start)
        sk = []
        for i in range(batch_num):
            batchlen = batch_stop[i] - batch_start[i]
            batch_ei = self.random_boundedl1(batchlen, batch_bound[i])
            sk += batch_ei
        return sk

    def derive(self, sk: list, pk: int) -> int:
        if pk < 0:
            raise ValueError(f"Invalid public key (smaller than 0): {pk} < 0.")
        if self.field(pk) == 2 or self.field(pk) == -2:
            raise ValueError(f"Invalid public key (singular curve): {pk}.")
        if not self.is_supersingular(pk):
            raise ValueError(f"the public key: {pk} is not supersingular!")
        return self.group_action(pk, sk)

    def group_action(self, a: int, e: list, verbose_level=0) -> int:
        set_verbose(verbose_level)

        def PointAccept(P, i: int, j: int) -> bool:
            if self.curve.isinfinity(P):
                return False
            # toss a coin with success probability gamma = (l*(lmin-1)) / (lmin*(l-1))
            r = get_randint(-(2**255), 2**255 - 1)  # get random 256-bit integer
            lmin = L[batch_start[i]]
            l = L[j]
            r %= lmin * (l - 1)  # now r uniformly distributed in [0, lmin*(l-1) - 1]
            assert r >= 0
            if r > l * (lmin - 1):
                return False
            return True

        def clear_public_prime(T, A24, I: list, do_verbose=False):
            if do_verbose:
                verbose(f"In clear_public_prime. Input: {T=}")
            T = self.curve.xdbl(T, A24)
            T = self.curve.xdbl(T, A24)
            for i in range(0, batch_num):
                if i in I:
                    continue
                for j in range(batch_start[i], batch_stop[i]):
                    T = self.curve.xmul_public(
                        T, A24, j
                    )  # public scalar mult： T = L[j] * T
            if do_verbose:
                verbose(f"Output: {T=}")
            return T

        def clear_private_prime(T, A24, I: list, J: list, do_verbose=False):
            if do_verbose:
                verbose(f"In clear_private_prime. Input: {T=}")
            for i in range(len(I)):
                for j in range(batch_start[I[i]], batch_stop[I[i]]):
                    if j == J[i]:
                        continue
                    T = self.curve.xmul_private(T, A24, j)  # private scalar mult
            if do_verbose:
                verbose(f"Output: {T=}")
            return T

        field = self.field
        batch_bound = self.prime_info["batch_bound"]
        batch_start = self.prime_info["batch_start"]
        batch_stop = self.prime_info["batch_stop"]
        n = self.prime_info["n"]  # length of prime table L[]
        L = self.prime_info["L"]
        batch_num = len(batch_start)

        A = (field(a), field(1))
        A24 = (field(a + 2), field(4))
        batchtodo = deepcopy(batch_bound)
        batchtodosum = sum(batchtodo)
        todo = [
            abs(e[i]) for i in range(0, n)
        ]  # NOTE: range(0, n) <=> for(i=0;i<n;i++)

        # NOTE: new keyspace
        batch_bound_max = max(batch_bound)
        assert batch_bound[-2] + batch_bound[-3] == batch_bound_max
        assert batch_bound[-4] + batch_bound[-5] == batch_bound_max
        assert batch_bound[-1] == 0

        while batchtodosum > 0:
            I = []  # 存储当前没做完的batch
            for i in range(0, batch_num - 5):
                if batchtodo[i] != 0:
                    I.append(i)  # 将i插入I的尾部
            if batchtodo[-5] != 0:
                I.append(batch_num - 5)
            elif batchtodo[-4] != 0:
                I.append(batch_num - 4)
            if batchtodo[-3] != 0:
                I.append(batch_num - 3)
            elif batchtodo[-2] != 0:
                I.append(batch_num - 2)
            # NOTE: Currently, batchbound[-1] is always 0
            # if batchtodo[-1] != 0:
            #     I.append(batch_num - 1)
            k = len(I)
            # Now I is in ascending order, eg. 0 1 2 3 4 5
            # Optimize the order of I
            assert k > 0
            tmp_I = I[-2::-1] + I[-1:]
            # Now tmp_I looks like 4 3 2 1 0 5
            if k >= 4:
                I = tmp_I[0:1] + tmp_I[2:-1]
                # Now I is 4 2 1 0
                I += tmp_I[1:2]
                # Now I is 4 2 1 0 3
                I += tmp_I[-1:]
            else:
                I = tmp_I.copy()
            # Now I looks like 4 2 1 0 3 5
            if k >= 4:  # not carefully checked, but 2 0 1 3 -> 1 0 2 3
                I[0], I[-2] = I[-2], I[0]
            # 初始化 CTIDH inner loop的J和epsilon
            J = [
                batch_start[I[i]] for i in range(0, k)
            ]  # Ji将存储第Ii个batch中，本次准备做的那个
            epsilon = []
            for i in range(0, k):
                j = batch_start[I[i]]
                eps_i = 0
                # sign(a)==0, -1, 1, if a==0, <0, >0
                eps_i = CMOV(eps_i, sign(e[j]), todo[j] != 0)
                epsilon.append(eps_i)

            for i in range(0, k):
                for j in range(batch_start[I[i]], batch_stop[I[i]]):
                    if todo[j] != 0:
                        J[i] = j
                        epsilon[i] = sign(e[j])
                        break
                    # NOTE: these should be implemented in constant-time.
                    # eg. use masking technique to avoid break
            verbose("\nCurrent atomic block (CTIDH inner loop):")
            verbose(f"I = {I}")
            verbose(f"J = {J}")
            verbose(f"A = {A}\n")
            # "CTIDH inner loop"
            T0, T1 = self.curve.elligator(A)
            verbose(f"elligator(A) Output:\n\t{T0=}\n\t{T1=}")
            for i in range(0, k):
                verbose(f"Doing batch {I[i]=}, index {J[i]=}, prime {L[J[i]]=}")

                T0, T1 = CSWAP(T0, T1, epsilon[i] < 0)
                if i == 0:
                    # T0 = [r]T0, T1 = [r]T1
                    T0 = clear_public_prime(T0, A24, I)
                    T0 = clear_private_prime(T0, A24, I, J)
                    if k > 1:
                        T1 = clear_public_prime(T1, A24, I)
                        T1 = clear_private_prime(T1, A24, I, J)

                P = deepcopy(T0)
                # 清理掉当前以外的其他待做素数
                for j in range(i + 1, k):
                    P = self.curve.xmul_private(P, A24, J[j])
                verbose(f"After for loop of xmul_private:\n\t{P=}", level=2)
                fi = PointAccept(P, I[i], J[i])
                maskisogeny = (fi == 1) and (epsilon[i] != 0)

                verbose(f"Accept the point? {fi=}", level=2)
                verbose(f"Will do this isogeny? {maskisogeny=}", level=2)

                if i == k - 2 and k > 2:  # 正在做倒数第二个batch
                    T0 = CMOV(T0, T1, (epsilon[i + 1] < 0) ^ (epsilon[i] < 0))

                # 做小同源
                if fi == 1:
                    Tnew = [T0, T1]
                    Tnewlen = 0
                    if i == k - 1:  # 最后一个同源不用推点
                        Tnewlen = 0
                    else:
                        Tnewlen = 2
                    if i == k - 2 and k > 2:  # 倒数第一个同源只需要推一个点
                        Tnewlen = 1
                    # NOTE: checked: matryoshka_isogeny does not modify the input T0 T1
                    Anew, Tnew = self.isogeny.matryoshka_isogeny(
                        A, Tnew, Tnewlen, P, J[i]
                    )
                    verbose(
                        f"matryoshka_isogeny output:\n\t{Anew=}\n\t{Tnew=}", level=2
                    )
                    A = CMOV(A, Anew, maskisogeny)
                    A24 = self.curve.xA24(A)
                    T0 = CMOV(T0, Tnew[0], maskisogeny)
                    T1 = CMOV(T1, Tnew[1], maskisogeny)
                    verbose(f"After 3 CMOVs:\n\t{A=}\n\t{T0=}\n\t{T1=}")
                    if verbose_level >= 3:
                        verbose(f"affine a={A[0] * A[1] ** (-1)}", level=3)

                # 处理T0和T1：做完小同源之后清掉对应小素数，避免未来重复做
                if i == k - 2 and k > 2:
                    T0 = self.curve.xmul_private(T0, A24, J[i])
                    T1 = deepcopy(T0)
                elif i < k - 1:
                    T0, T1 = CSWAP(T0, T1, epsilon[i] < 0)  # T0, T1换回去
                    T0 = self.curve.xmul_private(T0, A24, J[i])
                    T1 = self.curve.xmul_private(T1, A24, J[i])

                verbose(
                    f"Before changing todo, batchtodo, batchtodosum:\n\t{T0=}\n\t{T1=}",
                    level=2,
                )
                # NOTE: 实际中为了安全会用j扫一遍I[i]这个batch，右值改为 maskisogeny & (J[i] == j)
                todo[J[i]] -= maskisogeny
                batchtodo[I[i]] -= fi
                batchtodosum -= fi
                verbose(
                    f"After batchtodosum -= fi:\t{todo[J[i]]=}, {batchtodo[I[i]]=}, {batchtodosum=}",
                    level=2,
                )

                assert batchtodo[I[i]] >= 0
                assert batchtodosum >= 0

        anew = A[0] * A[1] ** (-1)  # anew = Ax * Az^(-1)
        anew = anew.get_int_value()
        verbose(f"Final result: {anew=}")
        set_verbose(0)

        return anew

    def is_supersingular(self, pk: int) -> bool:
        A = (self.field(pk), self.field(1))
        return self.curve.issupersingular(A)
