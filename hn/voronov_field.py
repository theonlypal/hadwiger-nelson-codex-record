"""Exact arithmetic in Q(sqrt2,sqrt3)(sqrt(h1),sqrt(h7)).

Physical embedding: sqrt2,sqrt3,A=sqrt(h1),B=sqrt(h7) are positive.
The 16 coefficients are ordered by extension part [1,A,B,AB], then base
mask [1,sqrt2,sqrt3,sqrt6]. Coordinates use integer numerators and one
canonical positive denominator; floats only provide candidate locations.

Base sign flips are allowed only after eliminating A and B. In particular,
flipping sqrt2 in an arbitrary tower element while holding A fixed would
not preserve its defining equation and is deliberately not implemented.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from math import gcd, lcm, sqrt

DEGREE = 16
SCHEMA = 'hn-voronov-field-v1'
BASE_RADICANDS = (2,3)
EXTENSION_SQUARES = ((5,3,-1,3), (5,3,1,-3))
FIELD_METADATA = {'base_radicands':[2,3],
                  'extension_squares':[list(p) for p in EXTENSION_SQUARES]}
_BASE_ROOTS = (1,2,3,6)
_BASE_FACTORS = tuple(tuple(_BASE_ROOTS[i&j] for j in range(4)) for i in range(4))


def _base_product(left, right):
    result = [0]*4
    for i,x in enumerate(left):
        for j,y in enumerate(right):
            result[i^j] += x*y*_BASE_FACTORS[i][j]
    return tuple(result)


_PRODUCT_SQUARE = _base_product(*EXTENSION_SQUARES)
if _PRODUCT_SQUARE != (-14,48,0,0):
    raise ArithmeticError('Unexpected product of extension radicands')
_PART_SQUARES = ((1,0,0,0), *EXTENSION_SQUARES, _PRODUCT_SQUARE)


def _basis_product(i,j):
    pi,bi = divmod(i,4); pj,bj = divmod(j,4)
    base_mask = bi^bj; part = (pi^pj)*4
    common_factor = _BASE_FACTORS[bi][bj]
    return tuple((part+(base_mask^r), common_factor*c*_BASE_FACTORS[base_mask][r])
                 for r,c in enumerate(_PART_SQUARES[pi&pj]) if c)


_PRODUCTS = tuple(tuple(_basis_product(i,j) for j in range(DEGREE)) for i in range(DEGREE))
_BASE_FLOATS = tuple(sqrt(x) for x in _BASE_ROOTS)
_A_FLOAT, _B_FLOAT = (sqrt(sum(c*r for c,r in zip(h,_BASE_FLOATS))) for h in EXTENSION_SQUARES)
_FLOAT_BASIS = tuple(part*base for part in (1.0,_A_FLOAT,_B_FLOAT,_A_FLOAT*_B_FLOAT) for base in _BASE_FLOATS)


@dataclass(frozen=True, slots=True, init=False)
class RealElement:
    coefficients: tuple[int,...]
    denominator: int

    def __init__(self, coefficients=0, denominator=1):
        if type(denominator) is not int or denominator==0:
            raise ValueError('Expected a nonzero integer denominator')
        if isinstance(coefficients, RealElement):
            numerator = coefficients.coefficients
            denominator *= coefficients.denominator
        elif isinstance(coefficients,(int,Fraction)):
            q = Fraction(coefficients,denominator)
            numerator = (q.numerator,)+(0,)*(DEGREE-1)
            denominator = q.denominator
        else:
            numerator = tuple(coefficients)
            if len(numerator)!=DEGREE or any(type(c) is not int for c in numerator):
                raise ValueError('Expected 16 integer coefficients')
        if type(denominator) is not int or denominator==0:
            raise ValueError('Expected a nonzero integer denominator')
        factor = reduce(gcd,numerator,abs(denominator))
        if denominator<0: factor = -factor
        if factor != 1:
            numerator = tuple(c//factor for c in numerator); denominator //= factor
        object.__setattr__(self,'coefficients',numerator)
        object.__setattr__(self,'denominator',denominator)

    @staticmethod
    def _coerce(value):
        if isinstance(value,RealElement): return value
        if isinstance(value,(int,Fraction)): return RealElement(value)
        return NotImplemented

    @classmethod
    def from_json(cls, values):
        if len(values)!=DEGREE: raise ValueError('Expected 16 rational coefficients')
        fractions = tuple(Fraction(c) for c in values)
        denominator = lcm(*(q.denominator for q in fractions))
        return cls(tuple(q.numerator*(denominator//q.denominator) for q in fractions),denominator)

    def to_json(self):
        return [str(Fraction(c,self.denominator)) for c in self.coefficients]

    @classmethod
    def from_base_mq(cls, value):
        """Lift only the Q(sqrt2,sqrt3) subfield of the existing MQ basis."""
        from .multiquadratic import MQ
        if not isinstance(value,MQ): raise TypeError('Expected an MQ element')
        if any(value.coefficients[4:]):
            raise ValueError('A base lift cannot contain sqrt5 or sqrt11')
        return cls(value.coefficients[:4]+(0,)*12,value.denominator)

    def __add__(self,other):
        other = self._coerce(other)
        if other is NotImplemented: return NotImplemented
        g = gcd(self.denominator,other.denominator)
        a,b = other.denominator//g,self.denominator//g
        return RealElement(tuple(a*x+b*y for x,y in zip(self.coefficients,other.coefficients)),self.denominator*a)

    __radd__ = __add__

    def __neg__(self):
        return RealElement(tuple(-c for c in self.coefficients),self.denominator)

    def __sub__(self,other):
        other = self._coerce(other)
        return NotImplemented if other is NotImplemented else self+(-other)

    def __rsub__(self,other):
        other = self._coerce(other)
        return NotImplemented if other is NotImplemented else other-self

    def __mul__(self,other):
        if isinstance(other,(int,Fraction)):
            q = Fraction(other)
            return RealElement(tuple(c*q.numerator for c in self.coefficients),self.denominator*q.denominator)
        other = self._coerce(other)
        if other is NotImplemented: return NotImplemented
        left = [(i,c) for i,c in enumerate(self.coefficients) if c]
        right = [(i,c) for i,c in enumerate(other.coefficients) if c]
        result = [0]*DEGREE
        for i,x in left:
            products = _PRODUCTS[i]
            for j,y in right:
                xy = x*y
                for target,factor in products[j]: result[target] += xy*factor
        return RealElement(tuple(result),self.denominator*other.denominator)

    __rmul__ = __mul__

    def __truediv__(self,other):
        if isinstance(other,(int,Fraction)): return self*(1/Fraction(other))
        other = self._coerce(other)
        return NotImplemented if other is NotImplemented else self*other.inverse()

    def __rtruediv__(self,other):
        other = self._coerce(other)
        return NotImplemented if other is NotImplemented else other*self.inverse()

    def flip_extension(self,bit):
        if bit not in (0,1): raise ValueError('Extension bit is 0 for A or 1 for B')
        mask = 1<<bit
        return RealElement(tuple(-c if (i//4)&mask else c for i,c in enumerate(self.coefficients)),self.denominator)

    def flip_base(self,bit):
        if bit not in (0,1): raise ValueError('Base bit is 0 for sqrt2 or 1 for sqrt3')
        if any(self.coefficients[4:]):
            raise ValueError('Base conjugation is allowed only after eliminating A and B')
        return RealElement(tuple(-c if i&(1<<bit) else c for i,c in enumerate(self.coefficients)),self.denominator)

    def inverse(self):
        if self==ZERO: raise ZeroDivisionError('Zero has no inverse')
        numerator,denominator = ONE,self
        for bit in (1,0):
            conjugate = denominator.flip_extension(bit)
            numerator *= conjugate; denominator *= conjugate
        if any(denominator.coefficients[4:]): raise ArithmeticError('Extension norm did not reach the base')
        for bit in (1,0):
            conjugate = denominator.flip_base(bit)
            numerator *= conjugate; denominator *= conjugate
        if any(denominator.coefficients[1:]) or denominator.coefficients[0]==0:
            raise ArithmeticError('Nonzero tower element has an invalid rational norm')
        return numerator*Fraction(denominator.denominator,denominator.coefficients[0])

    def field_norm(self):
        """The rational field norm down the quadratic tower."""
        value = self
        for bit in (1,0): value *= value.flip_extension(bit)
        for bit in (1,0): value *= value.flip_base(bit)
        if any(value.coefficients[1:]): raise ArithmeticError('Norm is not rational')
        return Fraction(value.coefficients[0],value.denominator)

    norm = field_norm

    def conjugate(self):
        """Ordinary complex conjugation fixes every real element."""
        return self

    conj = conjugate

    def __pow__(self,exponent):
        if type(exponent) is not int: raise TypeError('Integer exponent required')
        if exponent<0: return self.inverse()**(-exponent)
        result,base = ONE,self
        while exponent:
            if exponent&1: result *= base
            exponent >>= 1
            if exponent: base *= base
        return result

    def approximate(self):
        return sum(c*r for c,r in zip(self.coefficients,_FLOAT_BASIS))/self.denominator


ZERO = RealElement(0)
ONE = RealElement(1)
SQRT2 = RealElement((0,1)+(0,)*14)
SQRT3 = RealElement((0,0,1)+(0,)*13)
A = SQRT_H1 = RealElement((0,)*4+(1,)+(0,)*11)
B = SQRT_H7 = RealElement((0,)*8+(1,)+(0,)*7)
H1 = RealElement(EXTENSION_SQUARES[0]+(0,)*12)
H7 = RealElement(EXTENSION_SQUARES[1]+(0,)*12)


@dataclass(frozen=True, slots=True)
class Point:
    """An exact Euclidean point and complex rotation multiplier."""
    x: RealElement = ZERO
    y: RealElement = ZERO

    def __post_init__(self):
        x,y = RealElement._coerce(self.x),RealElement._coerce(self.y)
        if x is NotImplemented or y is NotImplemented:
            raise TypeError('Point coordinates must be RealElement or rational')
        object.__setattr__(self,'x',x);object.__setattr__(self,'y',y)

    @staticmethod
    def _coerce(value):
        if isinstance(value,Point): return value
        coordinate = RealElement._coerce(value)
        return NotImplemented if coordinate is NotImplemented else Point(coordinate)

    def __add__(self,other):
        other = self._coerce(other)
        return NotImplemented if other is NotImplemented else Point(self.x+other.x,self.y+other.y)

    __radd__ = __add__

    def __neg__(self): return Point(-self.x,-self.y)

    def __sub__(self,other):
        other = self._coerce(other)
        return NotImplemented if other is NotImplemented else self+(-other)

    def __rsub__(self,other):
        other = self._coerce(other)
        return NotImplemented if other is NotImplemented else other-self

    def __mul__(self,other):
        if isinstance(other,Point):
            return Point(self.x*other.x-self.y*other.y,self.x*other.y+self.y*other.x)
        other = RealElement._coerce(other)
        return NotImplemented if other is NotImplemented else Point(self.x*other,self.y*other)

    __rmul__ = __mul__

    def __truediv__(self,other):
        if isinstance(other,Point): return self*other.inverse()
        other = RealElement._coerce(other)
        return NotImplemented if other is NotImplemented else Point(self.x/other,self.y/other)

    def __rtruediv__(self,other):
        other = self._coerce(other)
        return NotImplemented if other is NotImplemented else other*self.inverse()

    def conjugate(self): return Point(self.x,-self.y)
    conj = conjugate

    def norm_squared(self): return self.x*self.x+self.y*self.y
    norm = norm_squared

    def is_unit(self): return self.norm_squared()==ONE

    def inverse(self):
        squared = self.norm_squared()
        if squared==ZERO: raise ZeroDivisionError('Origin has no complex inverse')
        return self.conjugate()/squared

    def __pow__(self,exponent):
        if type(exponent) is not int: raise TypeError('Integer exponent required')
        if exponent<0: return self.inverse()**(-exponent)
        result,base = EAST,self
        while exponent:
            if exponent&1: result *= base
            exponent >>= 1
            if exponent: base *= base
        return result

    def approximate(self): return complex(self.x.approximate(),self.y.approximate())

    def to_json(self): return {'x':self.x.to_json(),'y':self.y.to_json()}

    @classmethod
    def from_json(cls,value):
        return cls(RealElement.from_json(value['x']),RealElement.from_json(value['y']))

    @classmethod
    def from_base_mq(cls,value):
        from .multiquadratic import Point as MQPoint
        if not isinstance(value,MQPoint): raise TypeError('Expected an MQ Point')
        return cls(RealElement.from_base_mq(value.x),RealElement.from_base_mq(value.y))


ORIGIN = Point()
EAST = Point(ONE)
I = Point(ZERO,ONE)


def lift_base_point(point):
    return Point.from_base_mq(point)


PHI0 = Point((SQRT2*SQRT3+SQRT2)/4,(SQRT2*SQRT3-SQRT2)/4)
PHI1 = Point(SQRT2*SQRT3/3,SQRT3/3)


def psi_rotation(j):
    """Primary-source unit rotations for j=1 or j=7, positive radicals."""
    if j==1:
        wx,wy = (4+SQRT2+SQRT2*SQRT3)/4,(SQRT2*SQRT3-SQRT2)/4
        s,radical = (4+SQRT2+SQRT2*SQRT3)/2,A
    elif j==7:
        wx,wy = (4+SQRT2-SQRT2*SQRT3)/4,(SQRT2+SQRT2*SQRT3)/4
        s,radical = (4+SQRT2-SQRT2*SQRT3)/2,B
    else:
        raise ValueError('This field supports the specified j=1 and j=7 rotations')
    inverse_s = -s**3+8*s**2-20*s+16
    if s*inverse_s!=ONE: raise ArithmeticError('Invalid exact reciprocal formula')
    rotation = Point(((s+3)*wx+radical*wy)*inverse_s/4,
                     (radical*wx-(s+3)*wy)*inverse_s/4)
    if not rotation.is_unit() or not (Point(2)-rotation*Point(wx,wy)).is_unit():
        raise ArithmeticError('The primary rotation identities failed')
    return rotation


PSI1 = psi_rotation(1)
PSI7 = psi_rotation(7)
