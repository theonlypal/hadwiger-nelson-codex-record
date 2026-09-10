"""Strict arithmetic-only parser for Heule's primary Mathematica .vtx data."""
from __future__ import annotations

import ast
from fractions import Fraction
from pathlib import Path

from .multiquadratic import MQ,Point,RADICANDS


def rational_square_root(value):
    if any(value.coefficients[1:]):
        raise ValueError('Only rational radicands are supported')
    q=Fraction(value.coefficients[0],value.denominator)
    if q<0:raise ValueError('Negative real radicand')
    if q==0:return MQ(0)
    integer=q.numerator*q.denominator
    outside=1;mask=0
    for i,p in enumerate(RADICANDS):
        exponent=0
        while integer%p==0:
            integer//=p;exponent+=1
        outside*=p**(exponent//2)
        if exponent%2:mask|=1<<i
    if integer!=1:raise ValueError('Radicand outside the specified multiquadratic field')
    coefficients=[0]*16;coefficients[mask]=outside
    return MQ(tuple(coefficients),q.denominator)


def evaluate(node):
    if isinstance(node,ast.Constant) and type(node.value) is int:return MQ(node.value)
    if isinstance(node,ast.UnaryOp):
        if isinstance(node.op,ast.USub):return -evaluate(node.operand)
        if isinstance(node.op,ast.UAdd):return evaluate(node.operand)
    if isinstance(node,ast.BinOp):
        a,b=evaluate(node.left),evaluate(node.right)
        if isinstance(node.op,ast.Add):return a+b
        if isinstance(node.op,ast.Sub):return a-b
        if isinstance(node.op,ast.Mult):return a*b
        if isinstance(node.op,ast.Div):return a/b
    if (isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='sqrt'
        and len(node.args)==1 and not node.keywords):
        return rational_square_root(evaluate(node.args[0]))
    raise ValueError('Unsupported syntax in coordinate file')


def load_vtx(path):
    points=[]
    for line in Path(path).read_text().splitlines():
        if not line.strip():continue
        expression=line.replace('Sqrt[','sqrt(').replace(']',')').replace('{','(').replace('}',')')
        node=ast.parse(expression,mode='eval').body
        if not isinstance(node,ast.Tuple) or len(node.elts)!=2:
            raise ValueError('Expected a coordinate pair')
        points.append(Point(evaluate(node.elts[0]),evaluate(node.elts[1])))
    if len(set(points))!=len(points):raise ValueError('Duplicate exact source coordinates')
    return tuple(points)
