from copy import deepcopy

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

        if formula_name == 'tvelu':
            self.prime_info = read_prime_info_for_tvelu_test(prime_name)
        else:
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
            r = [
                get_randint(-(2 ** (b - 1)), 2 ** (b - 1) - 1) for _ in range(rnum)
            ]
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
            s = get_randint(-(2 ** (Ni - 1)), 2 ** (Ni - 1) - 1)
            for i in range(0, Ni):
                if s & 1 == 1:
                    e[i] = -e[i]
                    reject |= e[i] == 0
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
        if not self.is_supersingular(pk):
            raise ValueError(f"the public key: {pk} is not supersingular!")
        return self.group_action(pk, sk)


    def group_action(self, a: int, e: list, debug=False) -> int:
        def PointAccept(P, i: int, j: int) -> bool:
            if self.curve.isinfinity(P):
                return False
            # toss a coin with success probability gamma = (l*(lmin-1)) / (lmin*(l-1))
            r = get_randint(-2**255, 2**255-1)  # get random 256-bit integer
            lmin = L[batch_start[i]]
            l = L[j]
            r %= lmin * (l-1) # now r uniformly distributed in [0, lmin*(l-1) - 1]
            assert r >= 0
            if r > l * (lmin-1):
                return False
            return True

        def clear_public_prime(T, A24, I: list):
            T = self.curve.xdbl(T, A24); T = self.curve.xdbl(T, A24)
            for i in range(0, batch_num):
                if i in I:
                    continue
                for j in range(batch_start[i], batch_stop[i]):
                    T = self.curve.xmul_public(
                        T, A24, j
                    )  # public scalar mult： T = L[j] * T
            return T

        def clear_private_prime(T, A24, I: list, J: list):
            for i in range(len(I)):
                for j in range(batch_start[I[i]], batch_stop[I[i]]):
                    if j == J[i]:
                        continue
                    T = self.curve.xmul_private(T, A24, j)  # private scalar mult
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

        while batchtodosum > 0:
            I = []  # 存储当前没做完的batch
            for i in range(0, batch_num):
                if batchtodo[i] != 0:
                    I.append(i)  # 将i插入I的尾部
            k = len(I)

            # 略：按一种特定的方式打乱I，代码中说是一种优化，但应该对效率没影响
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
                epsilon.append( eps_i ) 

            for i in range(0, k):
                for j in range(batch_start[I[i]], batch_stop[I[i]]):
                    if todo[j] != 0:
                        J[i] = j
                        epsilon[i] = sign(e[j])
                        break
                    # NOTE: these should be implemented in constant-time.
                    # eg. use masking technique to avoid break
            
            if debug:
                print(f'I = {I}')
                print(f'J = {J}\n')
            # "CTIDH inner loop"
            T0, T1 = self.curve.elligator(A)
            for i in range(0, k):
                if debug:
                    print(f'Doing batch {I[i]}, index = {J[i]}, prime = {L[J[i]]}')

                T0, T1 = CSWAP(T0, T1, epsilon[i] < 0)
                if i == 0:
                    # T0 = [r]T0
                    T0 = clear_public_prime(T0, A24, I)
                    T0 = clear_private_prime(T0, A24, I, J)

                P = deepcopy(T0)
                # 清理掉当前以外的其他待做素数
                for j in range(i+1, k):
                    P = self.curve.xmul_private(P, A24, J[j])
                fi = PointAccept(P, I[i], J[i])
                maskisogeny = (fi == 1) and (epsilon[i] != 0)

                if debug:
                    # print(f'Accept the point?: {fi}')
                    print(f'Will do this isogeny? {maskisogeny}')

                if i == k - 2 and k > 2:  # 正在做倒数第二个batch
                    T0 = CMOV(T0, T1, (epsilon[i + 1] < 0) ^ (epsilon[i] < 0))

                # 做小同源
                if fi == 1:
                    Tnew = [T0, T1]
                    Tnewlen = 0
                    if i == k - 1:  # 最后一个同源不用推点
                        Tnewlen = 0
                    elif i == 0:  # 第一个同源只需要推一个点
                        Tnewlen = 1
                    else:
                        Tnewlen = 2
                    if i == k - 2 and k > 2:  # 倒数第一个同源只需要推一个点
                        Tnewlen = 1                    
                    # NOTE: checked: matryoshka_isogeny does not modify the input T0 T1
                    Anew, Tnew = self.isogeny.matryoshka_isogeny(
                        A, Tnew, Tnewlen, P, J[i]
                    )
                    A = CMOV(A, Anew, maskisogeny)
                    A24 = self.curve.xA24(A)
                    T0 = CMOV(T0, Tnew[0], maskisogeny)
                    T1 = CMOV(T1, Tnew[1], maskisogeny)

                if debug:
                    a = A[0] * A[1] ** (-1)
                    print(f'After that, a={a}')

                # 处理T0和T1：做完小同源之后清掉对应小素数，避免未来重复做
                if i == 0:
                    # 第一个小同源做完要出一个新的点，保证T0和T1相互独立
                    T_plus, T1 = self.curve.elligator(A)
                    T_plus, T1 = CSWAP(T_plus, T1, epsilon[i] < 0)
                    T1 = clear_public_prime(T1, A24, I)
                    T1 = clear_private_prime(T1, A24, I, J)
                if i == k - 2 and k > 2:
                    T0 = self.curve.xmul_private(T0, A24, J[i])
                    T1 = deepcopy(T0)
                elif i < k - 1:
                    T0, T1 = CSWAP(T0, T1, epsilon[i] < 0) # T0, T1换回去
                    T0 = self.curve.xmul_private(T0, A24, J[i])
                    T1 = self.curve.xmul_private(T1, A24, J[i])

                # NOTE: 实际中为了安全会用j扫一遍I[i]这个batch，右值改为 maskisogeny & (J[i] == j)
                todo[J[i]] -= maskisogeny
                batchtodo[I[i]] -= fi
                batchtodosum -= fi

                assert batchtodo[I[i]] >= 0
                assert batchtodosum >= 0

        anew = A[0] * A[1] ** (-1)  # anew = Ax * Az^(-1)
        anew = anew.get_int_value()
        return anew

    def is_supersingular(self, pk: int) -> bool:
        A = (self.field(pk), self.field(1))
        return self.curve.issupersingular(A)
