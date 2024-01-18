from primefield import PrimeField
from utils import read_prime_info, attrdict


# MontgomeryCurve class determines the family of supersingular elliptic curves over GF(p)
def MontgomeryCurve(prime_name="p2048_CTIDH", SDAC=False, validation_algo="origin"):
    if validation_algo not in ["origin", "doliskani", "pairing1", "pairing2"]:
        raise ValueError

    prime_info = read_prime_info(prime_name)

    L = prime_info["L"]
    n = prime_info["n"]
    k = prime_info["k"]
    f = prime_info["f"]
    p = prime_info["p"]

    field = PrimeField(p)

    # TODO: read SDAC info

    def measure(x, SQR=1.00, ADD=0.00):
        """
        Field additions, multiplications, and squarings
        SQR = 1.00
        # In F_p, we have SQR_{F_p} = SQR x MUL_{F_p}
        ADD = 0.00
        # In F_p, we have ADD_{F_p} = ADD x MUL_{F_p}
        """
        return x[0] + SQR * x[1] + ADD * x[2]

    def elligator(A):
        # TODO elligator
        raise NotImplementedError

    def affine_to_projective(affine):
        """
        affine_to_projective()
        input : the affine Montgomery coefficient A=A'/C with C=1
        output: projective Montgomery constants A24 := A' + 2C and C24 := 4C
                where E : y^2 = x^3 + (A'/C)*x^2 + x
        """
        return [affine + field(2), field(4)]

    def coeff(A):
        """
        ----------------------------------------------------------------------
        coeff()
        input : projective Montgomery constants A24 := A + 2C and C24 := 4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the affine Montgomery coefficient A/C
        ----------------------------------------------------------------------
        """
        output = A[0] + A[0]  # (2 * A24)
        output -= A[1]  # (2 * A24) - C24
        C24_inv = A[1] ** -1  # 1 / (C24)
        output += output  # 4*A = 2[(2 * A24) - C24]
        output *= C24_inv  # A/C = 2[(2 * A24) - C24] / C24

        return output

    # TODO: Check if this is the same as CTIDH
    def get_A(P, Q, PQ):
        """
        ----------------------------------------------------------------------
        coeff()
        input : the affine Montgomery x-coordinate points x(P) := (XP : YP),
                x(Q) := (XQ : ZQ), and x(P - Q) := (XPQ : ZPQ) on the curve
                E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Montgomery coefficient (A + 2C : 4C)
        ----------------------------------------------------------------------
        """
        t0 = P[0] + P[1]  # XP + ZP
        t1 = Q[0] + Q[1]  # XQ + ZQ

        t = t0 * t1  # (XP + ZP) * (XQ + ZQ)
        XPXQ = P[0] * Q[0]  # XP * XQ
        ZPZQ = P[1] * Q[1]  # ZP * ZQ

        t -= XPXQ
        t -= ZPZQ  # XPZQ + ZPXQ
        s = XPXQ - ZPZQ  # XPXQ - ZPZQ

        t0 = t * PQ[0]  # (XPZQ + ZPXQ) * XPQ
        t1 = s * PQ[1]  # (XPXQ - ZPZQ) * ZPQ
        t0 += t1  # (XPZQ + ZPXQ) * XPQ + (XPXQ - ZPZQ) * ZPQ
        t0 **= 2  # [(XPZQ + ZPXQ) * XPQ + (XPXQ - ZPZQ) * ZPQ] ^ 2

        t1 = t * PQ[1]  # (XPZQ + ZPXQ) * ZPQ
        s = ZPZQ * PQ[0]  # ZPZQ * XPQ
        t1 += s  # (XPZQ + ZPXQ) * ZPQ + ZPZQ * XPQ
        s = XPXQ * PQ[0]  # (XPXQ) * XPQ
        s += s  # 2 * [(XPXQ) * XPQ]
        s += s  # 4 * [(XPXQ) * XPQ]
        t1 *= s  # [(XPZQ + ZPXQ) * ZPQ + ZPZQ * XPQ] * (4 * [(XPXQ) * XPQ])

        t = ZPZQ * PQ[1]  # ZPZQ * ZPQ

        XPXQ = (
            t0 - t1
        )  # [(XPZQ + ZPXQ) * XPQ + (XPXQ - ZPZQ) * ZPQ] ^ 2 - [(XPZQ + ZPXQ) * ZPQ + ZPZQ * XPQ] * (4 * [(XPXQ) * XPQ])
        ZPZQ = s * t  # (4 * [(XPXQ) * XPQ]) * (ZPZQ * ZPQ)

        # ---
        B1 = ZPZQ + ZPZQ  # 2C
        B0 = XPXQ + B1  # A + 2C
        B1 += B1  # 4C
        return [B0, B1]

    def isinfinity(P):
        """isinfinity(P) determines if x(P) := (XP : ZP) = (1 : 0)"""
        return P[1] == 0

    def isequal(P, Q):
        """isequal(P, Q) determines if x(P) = x(Q)"""
        return (P[0] * Q[1]) == (P[1] * Q[0])

    def xdbl(P, A):
        """
        ----------------------------------------------------------------------
        xdbl()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, and
                the  projective Montgomery constants A24:= A + 2C and C24:=4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Montgomery x-coordinate point x([2]P)
        ----------------------------------------------------------------------
        """
        XP, ZP = P[0], P[1]
        V1 = XP + ZP
        V1 **= 2
        V2 = XP - ZP
        V2 **= V2

        raise NotImplementedError

    def xadd(P, Q, P_minus_Q):
        """
        ----------------------------------------------------------------------
        xadd()
        input : the projective Montgomery x-coordinate points x(P) := XP/ZP,
                x(Q) := XQ/ZQ, and x(P-Q) := XPQ/ZPQ
        output: the projective Montgomery x-coordinate point x(P+Q)
        ----------------------------------------------------------------------
        """
        XP, ZP = P[0], P[1]
        XQ, ZQ = Q[0], Q[1]
        X_minus, Z_minus = P_minus_Q[0], P_minus_Q[1]

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
        X_plus = Z_minus * V3 
        Z_plus = X_minus * V4

        return [X_plus, Z_plus]

    def xmul_Ladder(P, A, j):
        """
        ----------------------------------------------------------------------
        xmul_Ladder():  Constant-time Montgomery Ladder
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
                projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and an positive integer j
        output: the projective Montgomery x-coordinate point x([L[j]]P)
        ----------------------------------------------------------------------
        """
        raise NotImplementedError

    def xmul_SDAC(P, n, A):
        """
        ----------------------------------------------------------------------
        xmul(): scalar mult that use Shortest Differential Addition Chain (SDAC)
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
                projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and an positive integer j
        output: the projective Montgomery x-coordinate point x([L[j]]P)
        ----------------------------------------------------------------------
        """
        raise NotImplementedError

    xmul = xmul_SDAC if SDAC else xmul_Ladder

    # TODO: Add more useful things such as PRAC, eucild2d, cofactor_multiples, crisscross...
    # Read papers and see sibc...

    def issupersingular_origin(A):
        raise NotImplementedError

    def issupersingular_doliskani(A):
        raise NotImplementedError

    def issupersingular_pairing1(A):
        raise NotImplementedError

    def issupersingular_pairing2(A):
        raise NotImplementedError

    validation_algo_options = {
        "origin": issupersingular_origin,
        "doliskani": issupersingular_doliskani,
        "pairing1": issupersingular_pairing1,
        "pairing2": issupersingular_pairing2,
    }

    issupersingular = validation_algo_options[validation_algo]

    return attrdict(**locals())
