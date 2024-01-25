from random import randint

from .primefield import PrimeField
from .utils import read_prime_info, attrdict, CMOV, CSWAP, memoize


# MontgomeryCurve class determines the family of supersingular elliptic curves over GF(p)


@memoize
def MontgomeryCurve(prime_name="p1024_CTIDH", SDAC=False, validation="origin"):
    if validation not in ["origin", "doliskani", "pairing1", "pairing2"]:
        raise ValueError

    prime_info = read_prime_info(prime_name)

    L = prime_info["L"]
    n = prime_info["n"]
    k = prime_info["k"]
    f = prime_info["f"]
    p = prime_info["p"]

    field = PrimeField(p)

    type_field = type(field(2))


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

    def coeff(A: tuple):
        """
        ----------------------------------------------------------------------
        coeff()
        input : projective Montgomery constants A24 := A + 2C and C24 := 4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the affine Montgomery coefficient A/C
        ----------------------------------------------------------------------
        """
        A24, C24 = A
        output = A24 + A24  # (2 * A24)
        output -= C24  # (2 * A24) - C24
        C24_inv = C24 ** (-1)  # 1 / (C24)
        output += output  # 4*A = 2[(2 * A24) - C24]
        output *= C24_inv  # A/C = 2[(2 * A24) - C24] / C24

        return output

    # TODO: Check if this is the same as CTIDH
    def get_A24(P: tuple, Q: tuple, PQ: tuple) -> tuple:
        """
        ----------------------------------------------------------------------
        coeff()
        input : the affine Montgomery x-coordinate points x(P) := (XP : YP),
                x(Q) := (XQ : ZQ), and x(P - Q) := (XPQ : ZPQ) on the curve
                E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Montgomery coefficient (A + 2C : 4C)
        ----------------------------------------------------------------------
        """

        XP, ZP = P; XQ, ZQ = Q; XPQ, ZPQ = PQ

        t0 =XP +ZP  # XP + ZP
        t1 = XQ + ZQ  # XQ + ZQ

        t = t0 * t1  # (XP + ZP) * (XQ + ZQ)
        XPXQ =XP * XQ  # XP * XQ
        ZPZQ = ZP * ZQ  # ZP * ZQ

        t -= XPXQ
        t -= ZPZQ  # XPZQ + ZPXQ
        s = XPXQ - ZPZQ  # XPXQ - ZPZQ

        t0 = t * XPQ  # (XPZQ + ZPXQ) * XPQ
        t1 = s * ZPQ  # (XPXQ - ZPZQ) * ZPQ
        t0 += t1  # (XPZQ + ZPXQ) * XPQ + (XPXQ - ZPZQ) * ZPQ
        t0 **= 2  # [(XPZQ + ZPXQ) * XPQ + (XPXQ - ZPZQ) * ZPQ] ^ 2

        t1 = t * ZPQ  # (XPZQ + ZPXQ) * ZPQ
        s = ZPZQ * XPQ  # ZPZQ * XPQ
        t1 += s  # (XPZQ + ZPXQ) * ZPQ + ZPZQ * XPQ
        s = XPXQ * XPQ  # (XPXQ) * XPQ
        s += s  # 2 * [(XPXQ) * XPQ]
        s += s  # 4 * [(XPXQ) * XPQ]
        t1 *= s  # [(XPZQ + ZPXQ) * ZPQ + ZPZQ * XPQ] * (4 * [(XPXQ) * XPQ])

        t = ZPZQ * ZPQ  # ZPZQ * ZPQ

        XPXQ = (
            t0 - t1
        )  # [(XPZQ + ZPXQ) * XPQ + (XPXQ - ZPZQ) * ZPQ] ^ 2 - [(XPZQ + ZPXQ) * ZPQ + ZPZQ * XPQ] * (4 * [(XPXQ) * XPQ])
        ZPZQ = s * t  # (4 * [(XPXQ) * XPQ]) * (ZPZQ * ZPQ)

        # ---
        B1 = ZPZQ + ZPZQ  # 2C
        B0 = XPXQ + B1  # A + 2C
        B1 += B1  # 4C
        return (B0, B1)

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

    def xadd(P: tuple, Q: tuple, PQ: tuple):
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

        return [X_plus, Z_plus]

    def xmul_Ladder(P: tuple, A24: tuple, j) -> tuple:
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

    def xmul_SDAC(P: tuple, n, A24: tuple) -> tuple:
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

    # TODO: Decide whether these verification algorithms should use A24 or affine A
    def issupersingular_origin(A: tuple):
        raise NotImplementedError

    def issupersingular_doliskani(A: tuple):
        raise NotImplementedError

    def issupersingular_pairing1(A: tuple):
        raise NotImplementedError

    def issupersingular_pairing2(A: tuple):
        raise NotImplementedError

    validation_options = {
        "origin": issupersingular_origin,
        "doliskani": issupersingular_doliskani,
        "pairing1": issupersingular_pairing1,
        "pairing2": issupersingular_pairing2,
    }

    issupersingular = validation_options[validation]

    return attrdict(**locals())
