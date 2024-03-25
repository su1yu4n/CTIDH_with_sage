import numpy
from sage.all import sqrt
from copy import deepcopy
from .primefield import PrimeField
from .utils import read_prime_info, attrdict, CMOV, CSWAP, memoize, binrep, read_SDAC_info, batchnumber_of_Li

# MontgomeryCurve class determines the family of supersingular elliptic curves over GF(p)


@memoize
def MontgomeryCurve(prime_name="p1024_CTIDH", SDAC=False, validation="original", fast_kronecker=False):
    if validation not in ["original", "doliskani", "pairing1", "pairing2"]:
        raise ValueError

    prime_info = read_prime_info(prime_name)

    L = prime_info["L"]
    n = prime_info["n"]
    k = prime_info["k"]
    f = prime_info["f"]
    p = prime_info["p"]

    batch_start = prime_info["batch_start"]
    batch_stop = prime_info["batch_stop"]
    batch_maxdaclen = prime_info["batch_maxdaclen"]
    # batch_bound = prime_info["batch_bound"]

    field = PrimeField(p, fast_kronecker=fast_kronecker)

    type_field = type(field(2))

    # Shortest Differential Addition Chains (SDACs) for each l_i, used in fast scalar multiplication.
    SDACS = read_SDAC_info(prime_name)
    assert len(SDACS) > 0, f'No precomputed sdac information for {prime_name}'
    SDACS_LENGTH = list(map(len, SDACS))
    SDACS_REVERSED = list(map(lambda x:x[::-1], SDACS))


    def cmul(l: int):
        return numpy.array([4.0 * (SDACS_LENGTH[L.index(l)] + 2), 2.0 * (SDACS_LENGTH[L.index(l)] + 2), 6.0 * (SDACS_LENGTH[L.index(l)] + 2) - 2.0])
    c_xmul = list(map(cmul, L))  # list of the costs of each [l]P

    def measure(x, SQR=1.00, ADD=0.00):
        """
        Field additions, multiplications, and squarings
        SQR = 1.00
        # In F_p, we have SQR_{F_p} = SQR x MUL_{F_p}
        ADD = 0.00
        # In F_p, we have ADD_{F_p} = ADD x MUL_{F_p}
        """
        return x[0] + SQR * x[1] + ADD * x[2]


    def elligator(A: tuple):
        """Elligator from CTIDH original implementation (NOT Elligator 2)

        Args:
            A (tuple): tuple of ZModPrime class objects (Ax, Az), represent
              A = Ax / Az , or (Ax: Az) in P^1

        Returns:
            T+, T- (tuple):  projective x-coordinates of random points on EA(Fp) and EA(Fp^2) respectively.
            T- is the proj x-coord of point (x, iy) correspond to EA^t(Fp)'s point (-x, y)) .
                Each of them is a tuple like (Tx, Tz), and Tx, Tz are ZModPrime class object
        """
        Ax, Az = A

        # Check if Ax and Az are ZModPrime objects
        # Can't write isinstance(Ax, ZModPrime) here because my ZModPrime is inside the func PrimeField().
        # This can cause some invisible bugs with sage's GF element
        # if not hasattr(Ax, "value") or not hasattr(Az, "value"):

        # Now change to
        if not isinstance(Ax, type_field) or not isinstance(Az, type_field):
            raise TypeError("Input must be ZModPrime type tuple!")

        while True:
            one = field(1)

            # TODO: Change to a nice random generator
            u = field.get_random()  # line 1 of my pseudocode
            if u == 0:
                continue
            u2 = u**2
            D = u2 - 1
            if D == 0:
                continue  # line 7 of my pseudocode

            M = u2 * Ax
            T = M * Ax
            ctrl = Ax == 0
            P = Ax
            P = CMOV(P, one, ctrl)  # line 12 of my pseudocode
            M = CMOV(M, one, ctrl)
            T = CMOV(T, one, ctrl)
            D = D * Az
            D2 = D**2
            T = T + D2
            T = T * D
            T = T * P  # line 19 of my pseudocode

            Tplus_x = P
            Tminus_x = -M
            ctrl = not T.is_square()
            Tplus_x, Tminus_x = CSWAP(Tplus_x, Tminus_x, ctrl)

            Tplus_z = D
            Tminus_z = D

            return (Tplus_x, Tplus_z), (Tminus_x, Tminus_z)

    def affine_to_projective(affine) -> tuple:
        """
        affine_to_projective()
        input : the affine Montgomery coefficient A=A'/C with C=1
        output: projective Montgomery constants A24 := A' + 2C and C24 := 4C
                where E : y^2 = x^3 + (A'/C)*x^2 + x
        """
        return (affine + field(2), field(4))

    # def coeff(A24: tuple):
    #     """
    #     ----------------------------------------------------------------------
    #     coeff()
    #     input : projective Montgomery constants A24 := A + 2C and C24 := 4C
    #             where E : y^2 = x^3 + (A/C)*x^2 + x
    #     output: the affine Montgomery coefficient A/C
    #     ----------------------------------------------------------------------
    #     """
    #     A24, C24 = A24
    #     output = A24 + A24  # (2 * A24)
    #     output -= C24  # (2 * A24) - C24
    #     C24_inv = C24 ** (-1)  # 1 / (C24)
    #     output += output  # 4*A = 2[(2 * A24) - C24]
    #     output *= C24_inv  # A/C = 2[(2 * A24) - C24] / C24

    #     return output

    def xA24(A: tuple) -> tuple:
        Ax, Az = A
        if Az == 1: # save two addtion in this case. 
            return (Ax + 2, 4)
        two_Az = Az + Az
        return (Ax + two_Az, two_Az + two_Az)

    def isinfinity(P):
        """isinfinity(P) determines if x(P) := (XP : ZP) = (1 : 0)"""
        return P[1] == 0

    def isequal(P, Q):
        """isequal(P, Q) determines if x(P) = x(Q)"""
        return (P[0] * Q[1]) == (P[1] * Q[0])

    def xdbl(P: tuple, A24: tuple) -> tuple:
        """
        ----------------------------------------------------------------------
        xdbl()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, and
                the  projective Montgomery constants A24:= A + 2C and C24:=4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Montgomery x-coordinate point x([2]P)
        ----------------------------------------------------------------------
        """
        XP, ZP = P
        # NOTE: in fact XP and ZP can be zero during the protocol
        # assert XP != 0 and ZP != 0 

        V1 = XP + ZP  # line 1 of my pseudo code
        V1 **= 2
        V2 = XP - ZP
        V2 **= 2
        Z2P = A24[1] * V2
        X2P = Z2P * V1  # line 6 of my pseudo code

        V1 -= V2
        Z2P += A24[0] * V1
        Z2P *= V1

        return (X2P, Z2P)

    def xadd(P: tuple, Q: tuple, PQ: tuple) -> tuple:
        """
        ----------------------------------------------------------------------
        xadd()
        input : the projective Montgomery x-coordinate points x(P) := XP/ZP,
                x(Q) := XQ/ZQ, and x(P-Q) := XPQ/ZPQ
        output: the projective Montgomery x-coordinate point x(P+Q)
        ----------------------------------------------------------------------
        """
        XP, ZP = P
        XQ, ZQ = Q
        XPQ, ZPQ = PQ
        # assert XPQ != 0

        V0 = XP + ZP
        V1 = XQ - ZQ
        V1 = V1 * V0
        V0 = XP - ZP
        V2 = XQ + ZQ
        V2 = V2 * V0
        V3 = V1 + V2

        V3 **= 2
        V4 = V1 - V2
        V4 **= 2
        X_plus = ZPQ * V3
        Z_plus = XPQ * V4
        # if ZPQ == 0:
        #     assert X_plus == 0 and Z_plus == 0 
            # X_plus = 1; Z_plus = 0
        return (X_plus, Z_plus)
    

    def crisscross(alpha, beta, gamma, delta):
        """ crisscross() computes a*c + b*d, and a*c - b*d """
        t_1 = (alpha * delta)
        t_2 = (beta * gamma)
        #return (t_1 + t_2), (t_1 - t_2)
        # shave off a FF allocation: ##
        t_3 = t_1.copy() # object.__new__(t_1.__class__); t_3.x = t_1.x #      ## copy(t_1)
        t_1 += t_2                    ##
        t_3 -= t_2                   ##
        return t_1, t_3 # (t_1 + t_2), (t_1 - t_2)


    def xmul_Ladder(P: tuple, A24: tuple, j: int) -> tuple:
        """
        ----------------------------------------------------------------------
        xmul_Ladder():  Constant-time Montgomery Ladder
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
                projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and an positive integer j
        output: the projective Montgomery x-coordinate point x([L[j]]P)
        ----------------------------------------------------------------------
        """
        XP, ZP = P
        # NOTE: in fact XP and ZP can be zero during the protocol
        # assert XP != 0 and ZP != 0
        if ZP == 0:
            return (XP, field(0))
        kbits = binrep(L[j])
        kbitlen = len(kbits)

        x0, x1 = xdbl(P, A24), P
        for i in reversed(range(kbitlen-1)):
            x0, x1 = CSWAP(x0, x1, kbits[i+1] ^ kbits[i])
            x0, x1 = xdbl(x0, A24), xadd(x0, x1, P)
        x0, x1 = CSWAP(x0, x1, kbits[0])
        
        return x0


    def xmul_SDAC(P: tuple, A24: tuple, j: int) -> tuple:
        """
        ----------------------------------------------------------------------
        Scalar mult for PUBLIC primes that use Shortest Differential Addition Chain (SDAC)
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
                projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and an positive integer j
        output: the projective Montgomery x-coordinate point x([L[j]]P)
        ----------------------------------------------------------------------
        """
        P2 = xdbl(P, A24)
        R = [P, P2, xadd(P2, P, P)]
        for sdac in SDACS_REVERSED[j]:
            if isinfinity(R[sdac]): # if isinfinity(R[sdac]):
                R[:] = R[sdac ^ 1], R[2], xdbl(R[2], A24)
            else:
                R[:] = R[sdac ^ 1], R[2], xadd(R[2], R[sdac ^ 1], R[sdac])
        return R[2]


    def xmul_SDAC_safe(P: tuple, A24: tuple, j: int) -> tuple:
        """
        ----------------------------------------------------------------------
        Timing attack safe scalar mult for PRIVATE primes that use Shortest Differential Addition Chain (SDAC).
        This algorithm consider each batch's max dac length to resist timing attack.
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
                projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and an positive integer j
        output: the projective Montgomery x-coordinate point x([L[j]]P)
        ----------------------------------------------------------------------
        """
        P2 = xdbl(P, A24)
        R = [P, P2, xadd(P2, P, P)]

        b = batchnumber_of_Li(j, batch_start, batch_stop)
        maxdac_len = batch_maxdaclen[b]
        daclen = SDACS_LENGTH[j]
        SDACS_padded = SDACS_REVERSED[j] + [0] * (maxdac_len - daclen)

        for sdac in SDACS_padded:
            want = (daclen > 0)

            if isinfinity(R[sdac]):  # if isinfinity(R[sdac]):
                TMP = [R[sdac ^ 1], R[2], xadd(R[0],R[2],R[1])]
            else:
                TMP = [R[sdac ^ 1], R[2], xadd(R[sdac ^ 1],R[2], R[sdac])]
            
            R = CMOV(R, TMP, want)
            # del TMP
            daclen -= 1

        return R[2]


    xmul_public = xmul_SDAC if SDAC else xmul_Ladder
    xmul_private = xmul_SDAC_safe if SDAC else xmul_Ladder

    # TODO: Add more useful things such as PRAC, eucild2d, cofactor_multiples, crisscross...
    # Read papers and see sibc...


    def issupersingular_original(A: tuple):
        def order_rec(A24: tuple, Q: tuple, lower:int, upper:int, order:int):
            # print(f'Calling order_rec with lower = {lower}, upper = {upper}')
            if upper - lower == 1:
                if Q[1] == 0:
                    return order
                
                if not criticaltestdone:
                    Q = xmul_public(Q, A24, lower) # scalarmult by l_upper
                    if Q[1] == 0:
                        return L[lower] * order
                    else: # (p+1) * P != O, not supersingular!
                        return 0
            
            mid = lower + (upper - lower + 1) // 2

            Left = deepcopy(Q)
            for i in range(lower, mid):
                Left = xmul_public(Left, A24, i)

            order = order_rec(A24, Left, mid, upper, order)
            if order > int(4 * sqrt(p)):
                return order
            
            Right = deepcopy(Q)
            for i in range(mid, upper):
                Right = xmul_public(Right, A24, i)

            order = order_rec(A24, Right, lower, mid, order)    
            return order
        
        A24 = xA24(A)
        u = field.get_random()
        P = (u, 1); P = xdbl(P, A24); P = xdbl(P, A24)
        criticaltestdone = False # global wrt order_rec

        order = order_rec(A24, P, 0, n, 1)

        if order == 0:
            return False
        elif order > int(4 * sqrt(p)):
            return True
        else:
            print('Original validation algorithm failed. This should almost never happen. Retry now.')
            issupersingular_original(A)


    def issupersingular_doliskani(A: tuple):
        raise NotImplementedError

    def issupersingular_pairing1(A: tuple):
        raise NotImplementedError

    def issupersingular_pairing2(A: tuple):
        raise NotImplementedError

    validation_options = {
        "original": issupersingular_original,
        "doliskani": issupersingular_doliskani,
        "pairing1": issupersingular_pairing1,
        "pairing2": issupersingular_pairing2,
    }

    issupersingular = validation_options[validation]

    return attrdict(**locals())
