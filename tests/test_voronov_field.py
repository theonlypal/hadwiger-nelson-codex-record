from fractions import Fraction
import random
import unittest

from hn import multiquadratic as mq
from hn.voronov_field import (A,B,EAST,H1,H7,I,ONE,ORIGIN,PHI0,PHI1,PSI1,PSI7,
                              SQRT2,SQRT3,ZERO,Point,RealElement,lift_base_point,
                              psi_rotation)


class VoronovFieldTests(unittest.TestCase):
    def test_defining_radicals_and_product(self):
        self.assertEqual(SQRT2**2,RealElement(2))
        self.assertEqual(SQRT3**2,RealElement(3))
        self.assertEqual(A*A,5+3*SQRT2-SQRT3+3*SQRT2*SQRT3)
        self.assertEqual(B*B,5+3*SQRT2+SQRT3-3*SQRT2*SQRT3)
        self.assertEqual((A*B)**2,-14+48*SQRT2)

    def test_scalar_normalization_hash_and_serialization(self):
        values=tuple(range(-8,8))
        value=RealElement(values,-6)
        same=RealElement(tuple(-2*c for c in values),12)
        self.assertEqual(value,same)
        self.assertEqual(hash(value),hash(same))
        self.assertGreater(value.denominator,0)
        self.assertEqual(RealElement.from_json(value.to_json()),value)
        self.assertEqual(RealElement((0,)*16,-23),ZERO)
        self.assertEqual(RealElement(Fraction(5,7),3),RealElement(Fraction(5,21)))

    def test_random_ring_laws_and_inverse(self):
        rng=random.Random(724)
        for _ in range(10):
            a,b,c=[RealElement(tuple(rng.randrange(-2,3) for _ in range(16)),rng.randrange(1,8)) for _ in range(3)]
            self.assertEqual((a*b)*c,a*(b*c))
            self.assertEqual(a*(b+c),a*b+a*c)
            self.assertEqual(a*b,b*a)
            self.assertEqual(a*a.inverse(),ONE)
            self.assertEqual((a/b)*b,a)
            self.assertEqual(a**-2,(a*a).inverse())

    def test_quadratic_norm_descent(self):
        value=1+SQRT2+SQRT3+A+2*B+A*B
        no_b=value*value.flip_extension(1)
        self.assertFalse(any(no_b.coefficients[8:]))
        base=no_b*no_b.flip_extension(0)
        self.assertFalse(any(base.coefficients[4:]))
        self.assertEqual(RealElement(2).field_norm(),Fraction(2**16))
        self.assertEqual(SQRT2.field_norm(),Fraction(2**8))
        self.assertEqual((A*B).field_norm(),A.field_norm()*B.field_norm())

    def test_base_flips_refuse_uneliminated_extension(self):
        with self.assertRaises(ValueError):
            (1+A).flip_base(0)
        with self.assertRaises(ValueError):
            (2+B).flip_base(1)
        self.assertEqual((1+SQRT2+SQRT3).flip_base(0),1-SQRT2+SQRT3)
        self.assertEqual((1+SQRT2+SQRT3).flip_base(1),1+SQRT2-SQRT3)

    def test_base_lift_preserves_existing_exact_arithmetic(self):
        x=(mq.SQRT2+3*mq.SQRT3+mq.SQRT2*mq.SQRT3)/7
        y=(2-mq.SQRT2+mq.SQRT3)/3
        lift=RealElement.from_base_mq
        self.assertEqual(lift(x*y),lift(x)*lift(y))
        self.assertEqual(lift(x/y),lift(x)/lift(y))
        point=mq.Point(x,y)
        self.assertEqual(lift_base_point(point),Point(lift(x),lift(y)))
        with self.assertRaises(ValueError): lift(mq.SQRT5)
        with self.assertRaises(ValueError): lift_base_point(mq.Point(mq.SQRT11,0))

    def test_rotation_units_and_period(self):
        for rotation in (PHI0,PHI1,PSI1,PSI7):
            self.assertTrue(rotation.is_unit())
            self.assertEqual(rotation.inverse(),rotation.conjugate())
            self.assertEqual(rotation**-1,rotation.conjugate())
        self.assertEqual(PHI0**24,EAST)
        self.assertNotEqual(PHI0**12,EAST)
        self.assertEqual(PHI0**6,I)

    def test_primary_spindle_edges(self):
        for sign,root,rotation in ((1,A,PSI1),(-1,B,PSI7)):
            wx=(4+SQRT2+sign*SQRT2*SQRT3)/4
            wy=(SQRT2*SQRT3-SQRT2)/4 if sign==1 else (SQRT2+SQRT2*SQRT3)/4
            s=(4+SQRT2+sign*SQRT2*SQRT3)/2
            self.assertEqual(Point(wx,wy).norm_squared(),s)
            self.assertEqual(root*root,(s-1)*(9-s))
            self.assertTrue((Point(2)-rotation*Point(wx,wy)).is_unit())

    def test_complex_division_and_json(self):
        p=Point(A+SQRT2,B-SQRT3)
        q=Point(2+B,1+A*B)
        self.assertEqual((p/q)*q,p)
        self.assertEqual(Point.from_json(p.to_json()),p)
        self.assertEqual((p*p.conjugate()).y,ZERO)
        self.assertEqual((p*p.conjugate()).x,p.norm_squared())
        self.assertEqual(A*Point(1,2),Point(A,2*A))

    def test_invalid_zero_and_exponents(self):
        with self.assertRaises(ZeroDivisionError): ZERO.inverse()
        with self.assertRaises(ZeroDivisionError): ORIGIN.inverse()
        with self.assertRaises(ZeroDivisionError): ONE/0
        with self.assertRaises(ValueError): RealElement((0,)*15)
        with self.assertRaises(ValueError): RealElement(1,0)
        with self.assertRaises(TypeError): A**Fraction(1,2)
        with self.assertRaises(ValueError): psi_rotation(3)


if __name__=='__main__': unittest.main()
