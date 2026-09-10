"""Exact real arithmetic in Q(sqrt(2),sqrt(3),sqrt(5),sqrt(11)).

Basis index m denotes the product of positive square roots whose generator
bits are set in m. Canonical integer numerators/common denominators provide
exact equality, hashing, arithmetic and sign; floats never decide equality.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from math import gcd, isqrt, lcm, prod, sqrt

RADICANDS = (2, 3, 5, 11)
DEGREE = 1 << len(RADICANDS)
BASIS_RADICANDS = tuple(prod(p for i, p in enumerate(RADICANDS) if m & (1 << i))
                       for m in range(DEGREE))
_FACTORS = tuple(tuple(BASIS_RADICANDS[a & b] for b in range(DEGREE))
                 for a in range(DEGREE))
_FLOAT_BASIS = tuple(sqrt(p) for p in BASIS_RADICANDS)


@dataclass(frozen=True, slots=True, init=False)
class MQ:
    coefficients: tuple[int, ...]
    denominator: int

    def __init__(self, coefficients=0, denominator=1):
        if isinstance(coefficients, MQ):
            numerator = coefficients.coefficients
            denominator *= coefficients.denominator
        elif isinstance(coefficients, (int, Fraction)):
            q = Fraction(coefficients, denominator)
            numerator = (q.numerator,) + (0,) * (DEGREE-1)
            denominator = q.denominator
        else:
            numerator = tuple(coefficients)
            if len(numerator) != DEGREE or any(type(x) is not int for x in numerator):
                raise ValueError("Expected 16 integer coefficients")
        if type(denominator) is not int or not denominator:
            raise ValueError("Expected nonzero integer denominator")
        factor = reduce(gcd, numerator, abs(denominator))
        if denominator < 0:
            factor = -factor
        if factor != 1:
            numerator = tuple(x // factor for x in numerator)
            denominator //= factor
        object.__setattr__(self, "coefficients", numerator)
        object.__setattr__(self, "denominator", denominator)

    @classmethod
    def from_json(cls, values):
        if len(values) != DEGREE:
            raise ValueError("Expected exactly 16 rational coefficients")
        fractions = tuple(Fraction(x) for x in values)
        denominator = lcm(*(x.denominator for x in fractions))
        return cls(tuple(x.numerator * (denominator // x.denominator)
                         for x in fractions), denominator)

    def to_json(self):
        return [str(Fraction(x, self.denominator)) for x in self.coefficients]

    @staticmethod
    def _coerce(x):
        return x if isinstance(x, MQ) else MQ(x)

    def __add__(self, other):
        other = self._coerce(other)
        g = gcd(self.denominator, other.denominator)
        a, b = other.denominator // g, self.denominator // g
        return MQ(tuple(x*a + y*b for x, y in zip(self.coefficients, other.coefficients)),
                  self.denominator*a)

    __radd__ = __add__

    def __neg__(self):
        return MQ(tuple(-x for x in self.coefficients), self.denominator)

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) - self

    def __mul__(self, other):
        if isinstance(other, (int, Fraction)):
            q = Fraction(other)
            return MQ(tuple(x*q.numerator for x in self.coefficients),
                      self.denominator*q.denominator)
        other = self._coerce(other)
        out = [0] * DEGREE
        left = [(i, x) for i, x in enumerate(self.coefficients) if x]
        right = [(i, x) for i, x in enumerate(other.coefficients) if x]
        for a, x in left:
            for b, y in right:
                out[a ^ b] += x*y*_FACTORS[a][b]
        return MQ(tuple(out), self.denominator*other.denominator)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, (int, Fraction)):
            return self * (1/Fraction(other))
        return self * self._coerce(other).inverse()

    def __rtruediv__(self, other):
        return self._coerce(other) * self.inverse()

    def flip(self, bit):
        """Galois automorphism negating one generator square root."""
        return MQ(tuple(-x if m & (1 << bit) else x
                        for m, x in enumerate(self.coefficients)), self.denominator)

    def inverse(self):
        if self == ZERO:
            raise ZeroDivisionError
        denominator = self
        numerator = ONE
        for bit in range(len(RADICANDS)):
            conjugate = denominator.flip(bit)
            numerator *= conjugate
            denominator *= conjugate
        assert not any(denominator.coefficients[1:])
        return numerator * Fraction(denominator.denominator, denominator.coefficients[0])

    def __pow__(self, exponent):
        if type(exponent) is not int:
            raise TypeError("Integer exponent required")
        if exponent < 0:
            return self.inverse() ** -exponent
        result, base = ONE, self
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent >>= 1
        return result

    def sign(self):
        """Determine sign by successively tighter *rational* root intervals."""
        if self == ZERO:
            return 0
        bits = 8
        while True:
            lo = hi = 0
            for m, coefficient in enumerate(self.coefficients):
                if not coefficient:
                    continue
                square = BASIS_RADICANDS[m] << (2*bits)
                lower = isqrt(square)
                upper = lower if lower*lower == square else lower+1
                if coefficient > 0:
                    lo += coefficient*lower
                    hi += coefficient*upper
                else:
                    lo += coefficient*upper
                    hi += coefficient*lower
            if lo > 0:
                return 1
            if hi < 0:
                return -1
            bits *= 2

    def __lt__(self, other):
        return (self-other).sign() < 0

    def __le__(self, other):
        return (self-other).sign() <= 0

    def approximate(self):
        return sum(x*y for x, y in zip(self.coefficients, _FLOAT_BASIS)) / self.denominator


ZERO = MQ(0)
ONE = MQ(1)
SQRT2, SQRT3, SQRT5, SQRT11 = tuple(
    MQ(tuple(int(m == (1 << bit)) for m in range(DEGREE))) for bit in range(4))


@dataclass(frozen=True, slots=True)
class Point:
    """An exact Euclidean point, also used as a complex rotation multiplier."""
    x: MQ = ZERO
    y: MQ = ZERO

    def __post_init__(self):
        object.__setattr__(self, "x", MQ._coerce(self.x))
        object.__setattr__(self, "y", MQ._coerce(self.y))

    def __add__(self, other):
        if not isinstance(other, Point):
            other = Point(other)
        return Point(self.x+other.x, self.y+other.y)

    __radd__ = __add__

    def __neg__(self):
        return Point(-self.x, -self.y)

    def __sub__(self, other):
        return self + (-other if isinstance(other, Point) else Point(-other))

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        if not isinstance(other, Point):
            return Point(self.x*other, self.y*other)
        return Point(self.x*other.x-self.y*other.y, self.x*other.y+self.y*other.x)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return Point(self.x/other, self.y/other)

    def conjugate(self):
        return Point(self.x, -self.y)

    def norm_squared(self):
        return self.x*self.x+self.y*self.y

    def is_unit(self):
        return self.norm_squared() == ONE

    def __pow__(self, exponent):
        if type(exponent) is not int:
            raise TypeError("Integer exponent required")
        if exponent < 0:
            return (self.conjugate()/self.norm_squared()) ** -exponent
        result, base = Point(ONE), self
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent >>= 1
        return result

    def approximate(self):
        return complex(self.x.approximate(), self.y.approximate())

    def to_json(self):
        return {"x": self.x.to_json(), "y": self.y.to_json()}

    @classmethod
    def from_json(cls, obj):
        return cls(MQ.from_json(obj["x"]), MQ.from_json(obj["y"]))


ORIGIN = Point()
EAST = Point(ONE)
