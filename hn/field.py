"""Exact complex arithmetic in Q(exp(pi*i/21))[sqrt(5)].

An element is represented as (a(t) + sqrt(5)*b(t))/d, deg a,b < 12,
with integer coefficients and a positive common denominator.  The selected
complex embedding is t=exp(pi*i/21), sqrt(5)>0.  Floats are display/search
helpers only; equality, hashing, conjugation and norms use integers.
"""

from __future__ import annotations

import cmath
from fractions import Fraction
from functools import reduce
from math import gcd, lcm, sqrt
from dataclasses import dataclass

ORDER = 42
DEGREE = 12
RADICAND = 5
# Low degree first, including the monic leading coefficient.
MODULUS = (1, 1, 0, -1, -1, 0, 1, 0, -1, -1, 0, 1, 1)
_REDUCTION = tuple((j, -v) for j, v in enumerate(MODULUS[:-1]) if v)
_ZEROP = (0,) * DEGREE


def _reduce(p):
    p = list(p)
    for k in range(len(p) - 1, DEGREE - 1, -1):
        v = p[k]
        if v:
            for j, coefficient in _REDUCTION:
                p[k - DEGREE + j] += v * coefficient
    return tuple((p + [0] * DEGREE)[:DEGREE])


def _mul(a, b):
    if not any(a) or not any(b):
        return _ZEROP
    p = [0] * (2 * DEGREE - 1)
    aa = [(j, v) for j, v in enumerate(a) if v]
    bb = [(j, v) for j, v in enumerate(b) if v]
    for j, u in aa:
        for k, v in bb:
            p[j + k] += u * v
    return _reduce(p)


def _add(a, b):
    return tuple(x + y for x, y in zip(a, b))


_POWERS = tuple(_reduce([0] * j + [1]) for j in range(ORDER))
_CONJUGATE_BASIS = tuple(_POWERS[(-j) % ORDER] for j in range(DEGREE))
_FLOAT_BASIS = tuple(cmath.exp(cmath.pi * 1j * j / 21) for j in range(DEGREE))


def _conjugate(a):
    out = [0] * DEGREE
    for coefficient, basis in zip(a, _CONJUGATE_BASIS):
        if coefficient:
            for j, v in enumerate(basis):
                out[j] += coefficient * v
    return tuple(out)


@dataclass(frozen=True, slots=True, init=False)
class Element:
    """Canonical and hashable exact complex algebraic number."""

    coefficients: tuple[int, ...]
    denominator: int

    def __init__(self, coefficients=0, denominator=1):
        if isinstance(coefficients, Element):
            numerator = coefficients.coefficients
            denominator *= coefficients.denominator
        elif isinstance(coefficients, (int, Fraction)):
            q = Fraction(coefficients, denominator)
            numerator = (q.numerator,) + (0,) * (2 * DEGREE - 1)
            denominator = q.denominator
        else:
            numerator = tuple(coefficients)
            if len(numerator) == DEGREE:
                numerator += _ZEROP
            if len(numerator) != 2 * DEGREE:
                raise ValueError("Expected 12 or 24 integer coefficients")
            if any(not isinstance(v, int) for v in numerator):
                raise TypeError("Internal coefficients must be integers")
        if not isinstance(denominator, int) or denominator == 0:
            raise ValueError("Expected a nonzero integer denominator")
        common = reduce(gcd, numerator, abs(denominator))
        if denominator < 0:
            common = -common
        if common != 1:
            numerator = tuple(v // common for v in numerator)
            denominator //= common
        object.__setattr__(self, "coefficients", numerator)
        object.__setattr__(self, "denominator", denominator)

    @classmethod
    def from_json(cls, obj):
        fractions = tuple(Fraction(v) for v in obj["a"] + obj["b"])
        denominator = lcm(*(v.denominator for v in fractions))
        return cls(tuple(v.numerator * (denominator // v.denominator)
                         for v in fractions), denominator)

    def to_json(self):
        values = [str(Fraction(v, self.denominator)) for v in self.coefficients]
        return {"a": values[:DEGREE], "b": values[DEGREE:]}

    @staticmethod
    def _coerce(other):
        return other if isinstance(other, Element) else Element(other)

    def __add__(self, other):
        other = self._coerce(other)
        d = gcd(self.denominator, other.denominator)
        a = other.denominator // d
        b = self.denominator // d
        return Element(tuple(x * a + y * b for x, y in
                             zip(self.coefficients, other.coefficients)),
                       self.denominator * a)

    __radd__ = __add__

    def __neg__(self):
        return Element(tuple(-v for v in self.coefficients), self.denominator)

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) - self

    def __mul__(self, other):
        if isinstance(other, (int, Fraction)):
            q = Fraction(other)
            return Element(tuple(v * q.numerator for v in self.coefficients),
                           self.denominator * q.denominator)
        other = self._coerce(other)
        a, b = self.coefficients[:DEGREE], self.coefficients[DEGREE:]
        c, d = other.coefficients[:DEGREE], other.coefficients[DEGREE:]
        ac, bd = _mul(a, c), _mul(b, d)
        ad, bc = _mul(a, d), _mul(b, c)
        real = tuple(x + RADICAND * y for x, y in zip(ac, bd))
        radical = _add(ad, bc)
        return Element(real + radical, self.denominator * other.denominator)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, (int, Fraction)):
            return self * (1 / Fraction(other))
        return self * self._coerce(other).inverse()

    def __pow__(self, power):
        if not isinstance(power, int):
            raise TypeError("Exponent must be an integer")
        if power < 0:
            return self.inverse() ** (-power)
        value = ONE
        base = self
        while power:
            if power & 1:
                value = value * base
            base = base * base
            power >>= 1
        return value

    def conjugate(self):
        return Element(_conjugate(self.coefficients[:DEGREE]) +
                       _conjugate(self.coefficients[DEGREE:]), self.denominator)

    def norm_squared(self):
        return self * self.conjugate()

    def is_unit(self):
        return self.norm_squared() == ONE

    def inverse(self):
        norm = self.norm_squared()
        if norm == ONE:
            return self.conjugate()
        if all(v == 0 for v in norm.coefficients[1:]):
            return self.conjugate() * Fraction(norm.denominator, norm.coefficients[0])
        # Inverting a generic element is intentionally a slow, exact fallback.
        import sympy as sp
        x = sp.Symbol("t")
        f = sum(v * x**j for j, v in enumerate(MODULUS))
        a = sum(v * x**j for j, v in enumerate(self.coefficients[:DEGREE]))
        b = sum(v * x**j for j, v in enumerate(self.coefficients[DEGREE:]))
        inverse_norm = sp.invert(a*a - RADICAND*b*b, f, x)
        aa = sp.Poly(sp.rem(a * inverse_norm, f, x), x)
        bb = sp.Poly(sp.rem(-b * inverse_norm, f, x), x)
        return Element.from_json({
            "a": [str(aa.nth(j) * self.denominator) for j in range(DEGREE)],
            "b": [str(bb.nth(j) * self.denominator) for j in range(DEGREE)],
        })

    def approximate(self):
        c = self.coefficients
        return sum((c[j] + sqrt(RADICAND) * c[j+DEGREE]) * _FLOAT_BASIS[j]
                   for j in range(DEGREE)) / self.denominator


ZERO = Element(0)
ONE = Element(1)
ZETA = Element((0, 1) + (0,) * (DEGREE - 2))
SQRT5 = Element(_ZEROP + (1,) + (0,) * (DEGREE - 1))
FIELD_METADATA = {
    "cyclotomic_order": ORDER,
    "radicand": RADICAND,
    "embedding": "zeta=exp(pi*i/21), sqrt(5)>0",
    "modulus_low_to_high": list(MODULUS),
}
