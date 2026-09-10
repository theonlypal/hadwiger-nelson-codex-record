"""Exact complex arithmetic in Q(zeta_42)(sqrt(5),sqrt(33)).

The embedding is zeta_42=exp(pi*i/21) with positive real radicals. Parts
have masks 0,1,2,3, denoting 1,sqrt(5),sqrt(33),sqrt(165). Each part has
12 coefficients in the cyclotomic basis; a positive common denominator
and gcd normalization give canonical equality and hashing. Approximation
is used only for display and candidate enumeration.
"""
from __future__ import annotations

import cmath
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from math import gcd, lcm, sqrt

ORDER = 42
DEGREE = 12
RADICANDS = (5, 33)
PART_COUNT = 4
MODULUS = (1, 1, 0, -1, -1, 0, 1, 0, -1, -1, 0, 1, 1)
_REDUCTION = tuple((j, -v) for j, v in enumerate(MODULUS[:-1]) if v)
_ZERO_PART = (0,)*DEGREE
_MASK_FACTORS = (1, 5, 33, 165)


def _reduce_polynomial(p):
    p = list(p)
    for k in range(len(p)-1, DEGREE-1, -1):
        value = p[k]
        if value:
            for j, coefficient in _REDUCTION:
                p[k-DEGREE+j] += value*coefficient
    return tuple((p+[0]*DEGREE)[:DEGREE])


def _multiply_polynomials(a, b):
    aa = [(j, v) for j, v in enumerate(a) if v]
    bb = [(j, v) for j, v in enumerate(b) if v]
    if not aa or not bb:
        return _ZERO_PART
    out = [0]*(2*DEGREE-1)
    for j, u in aa:
        for k, v in bb:
            out[j+k] += u*v
    return _reduce_polynomial(out)


_POWERS = tuple(_reduce_polynomial([0]*j+[1]) for j in range(ORDER))
_CONJUGATE_BASIS = tuple(_POWERS[-j % ORDER] for j in range(DEGREE))
_FLOAT_BASIS = tuple(cmath.exp(cmath.pi*1j*j/21) for j in range(DEGREE))
_FLOAT_FACTORS = tuple(sqrt(n) for n in _MASK_FACTORS)


def _conjugate_polynomial(a):
    out = [0]*DEGREE
    for coefficient, basis in zip(a, _CONJUGATE_BASIS):
        if coefficient:
            for j, value in enumerate(basis):
                out[j] += coefficient*value
    return tuple(out)


@dataclass(frozen=True, slots=True, init=False)
class Element:
    coefficients: tuple[int, ...]
    denominator: int

    def __init__(self, coefficients=0, denominator=1):
        if not isinstance(denominator, int) or isinstance(denominator, bool) or denominator == 0:
            raise ValueError('Expected a nonzero integer denominator')
        if isinstance(coefficients, Element):
            numerator = coefficients.coefficients
            denominator *= coefficients.denominator
        elif isinstance(coefficients, (int, Fraction)):
            value = Fraction(coefficients, denominator)
            numerator = (value.numerator,)+(0,)*(DEGREE*PART_COUNT-1)
            denominator = value.denominator
        else:
            numerator = tuple(coefficients)
            if len(numerator) not in (DEGREE, 2*DEGREE, PART_COUNT*DEGREE):
                raise ValueError('Expected 12, 24 or 48 integer coefficients')
            if any(type(v) is not int for v in numerator):
                raise TypeError('Internal coefficients must be integers')
            numerator += (0,)*(PART_COUNT*DEGREE-len(numerator))
        common = reduce(gcd, numerator, abs(denominator))
        if denominator < 0:
            common = -common
        if common != 1:
            numerator = tuple(v//common for v in numerator)
            denominator //= common
        object.__setattr__(self, 'coefficients', numerator)
        object.__setattr__(self, 'denominator', denominator)

    @property
    def parts(self):
        return tuple(self.coefficients[j*DEGREE:(j+1)*DEGREE] for j in range(PART_COUNT))

    @classmethod
    def from_legacy(cls, element):
        if len(element.coefficients) != 2*DEGREE:
            raise ValueError('Expected legacy 24-coefficient field element')
        return cls(element.coefficients, element.denominator)

    @classmethod
    def from_json(cls, obj):
        parts = obj['parts']
        if len(parts) != PART_COUNT or any(len(p) != DEGREE for p in parts):
            raise ValueError('Expected four parts, each with 12 rational strings')
        fractions = tuple(Fraction(v) for part in parts for v in part)
        denominator = lcm(*(v.denominator for v in fractions))
        return cls(tuple(v.numerator*(denominator//v.denominator) for v in fractions), denominator)

    def to_json(self):
        return {'parts': [[str(Fraction(v, self.denominator)) for v in part] for part in self.parts]}

    @staticmethod
    def _coerce(other):
        return other if isinstance(other, Element) else Element(other)

    def __add__(self, other):
        other = self._coerce(other)
        common = gcd(self.denominator, other.denominator)
        left, right = other.denominator//common, self.denominator//common
        return Element(tuple(x*left+y*right for x, y in zip(self.coefficients, other.coefficients)),
                       self.denominator*left)

    __radd__ = __add__

    def __neg__(self):
        return Element(tuple(-v for v in self.coefficients), self.denominator)

    def __sub__(self, other):
        return self+-self._coerce(other)

    def __rsub__(self, other):
        return self._coerce(other)-self

    def __mul__(self, other):
        if isinstance(other, (int, Fraction)):
            value = Fraction(other)
            return Element(tuple(v*value.numerator for v in self.coefficients),
                           self.denominator*value.denominator)
        other = self._coerce(other)
        out = [[0]*DEGREE for _ in range(PART_COUNT)]
        left = [(mask, part) for mask, part in enumerate(self.parts) if any(part)]
        right = [(mask, part) for mask, part in enumerate(other.parts) if any(part)]
        for a, aa in left:
            for b, bb in right:
                product = _multiply_polynomials(aa, bb)
                factor = _MASK_FACTORS[a & b]
                for j, value in enumerate(product):
                    out[a ^ b][j] += factor*value
        return Element(tuple(v for part in out for v in part), self.denominator*other.denominator)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, (int, Fraction)):
            return self*(1/Fraction(other))
        return self*self._coerce(other).inverse()

    def __rtruediv__(self, other):
        return self._coerce(other)*self.inverse()

    def __pow__(self, power):
        if type(power) is not int:
            raise TypeError('Exponent must be an integer')
        if power < 0:
            return self.inverse()**(-power)
        out, base = ONE, self
        while power:
            if power & 1:
                out *= base
            power >>= 1
            if power:
                base *= base
        return out

    def conjugate(self):
        return Element(tuple(v for part in self.parts for v in _conjugate_polynomial(part)),
                       self.denominator)

    def flip_radicals(self, mask):
        if type(mask) is not int or not 0 <= mask < PART_COUNT:
            raise ValueError('Radical-flip mask must be 0, 1, 2 or 3')
        return Element(tuple((-1 if (part_mask & mask).bit_count() % 2 else 1)*v
                             for part_mask, part in enumerate(self.parts) for v in part),
                       self.denominator)

    def norm_squared(self):
        return self*self.conjugate()

    def is_unit(self):
        return self.norm_squared() == ONE

    def inverse(self):
        if self == ZERO:
            raise ZeroDivisionError('Zero has no inverse')
        conjugate = self.conjugate()
        norm = self*conjugate
        if norm == ONE:
            return conjugate
        if not any(norm.coefficients[1:]):
            return conjugate*Fraction(norm.denominator, norm.coefficients[0])
        # First descend the two quadratic extensions, then invert one
        # polynomial modulo Phi_42. This fallback is exact, but deliberately
        # optimized for correctness rather than bulk graph construction.
        adjugate = self.flip_radicals(1)*self.flip_radicals(2)*self.flip_radicals(3)
        relative_norm = self*adjugate
        if any(relative_norm.coefficients[DEGREE:]):
            raise ArithmeticError('Radical norm did not descend to cyclotomic base')
        import sympy as sp
        t = sp.Symbol('t')
        modulus = sum(v*t**j for j, v in enumerate(MODULUS))
        polynomial = sum(v*t**j for j, v in enumerate(relative_norm.coefficients[:DEGREE]))
        inverse = sp.Poly(sp.invert(polynomial, modulus, t), t)
        values = [Fraction(inverse.nth(j))*relative_norm.denominator for j in range(DEGREE)]
        denominator = lcm(*(v.denominator for v in values))
        base_inverse = Element(tuple(v.numerator*(denominator//v.denominator) for v in values), denominator)
        return adjugate*base_inverse

    def approximate(self):
        return sum(factor*sum(v*b for v, b in zip(part, _FLOAT_BASIS))
                   for part, factor in zip(self.parts, _FLOAT_FACTORS))/self.denominator


ZERO = Element(0)
ONE = Element(1)
ZETA = Element((0, 1)+(0,)*(DEGREE-2))
SQRT5 = Element(_ZERO_PART+(1,)+(0,)*(3*DEGREE-1))
SQRT33 = Element((0,)*(2*DEGREE)+(1,)+(0,)*(2*DEGREE-1))
FIELD_METADATA = {
    'cyclotomic_order': ORDER,
    'radicands': list(RADICANDS),
    'embedding': 'zeta=exp(pi*i/21),sqrt(5)>0,sqrt(33)>0',
    'modulus_low_to_high': list(MODULUS),
    'part_basis': ['1', 'sqrt(5)', 'sqrt(33)', 'sqrt(165)'],
}
