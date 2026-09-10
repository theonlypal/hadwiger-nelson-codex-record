#!/usr/bin/env python3
"""Exact rotation, square-class, residue, and primary-expression checks.

No Euclidean graph is constructed. Symbolic A+B*sqrt(h) arithmetic is used
only for the two rotation identities, independently of the new field module.
"""
from __future__ import annotations
import ast
import hashlib
from itertools import product
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from hn.multiquadratic import Point,ONE,ZERO,SQRT2,SQRT3

ORIGIN=Point(ZERO,ZERO);EAST=Point(ONE,ZERO)
R2,R3=SQRT2,SQRT3;R6=R2*R3
PHI0=Point((R6+R2)/4,(R6-R2)/4)
PHI1=Point(R6/3,R3/3)


def finite_add(a,b):return ((a[0]+b[0])%5,(a[1]+b[1])%5)
def finite_mul(a,b):return ((a[0]*b[0]+2*a[1]*b[1])%5,(a[0]*b[1]+a[1]*b[0])%5)
def finite_scale(n,a):return (n*a[0]%5,n*a[1]%5)


def source_check(line,h,base,radical):
    """Parse the author's arithmetic expression, then cross-multiply exactly."""
    def add(a,b):return(a[0]+b[0],a[1]+b[1])
    def mul(a,b):return(a[0]*b[0]+(a[1]*b[1])*h,a[0]*b[1]+a[1]*b[0])
    def ev(node):
        if isinstance(node,ast.Constant) and type(node.value)is int:return(Point(ONE*node.value,ZERO),ORIGIN)
        if isinstance(node,ast.Name) and node.id=='I':return(Point(ZERO,ONE),ORIGIN)
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):
            a=ev(node.operand);return(-a[0],-a[1])
        if isinstance(node,ast.BinOp):
            a,b=ev(node.left),ev(node.right)
            if isinstance(node.op,ast.Add):return add(a,b)
            if isinstance(node.op,ast.Sub):return add(a,(-b[0],-b[1]))
            if isinstance(node.op,ast.Mult):return mul(a,b)
            if isinstance(node.op,ast.Div):
                if b[1]!=ORIGIN:raise ValueError('Unexpected inner radical denominator')
                inverse=Point(b[0].x,-b[0].y)/b[0].norm_squared()
                return(a[0]*inverse,a[1]*inverse)
        if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='Sqrt':
            if len(node.args)!=1 or node.keywords:raise ValueError('Unexpected Sqrt syntax')
            value=ev(node.args[0])
            if value[1]!=ORIGIN or value[0].y!=ZERO:raise ValueError('Unexpected nested extension')
            for root in (ONE,R2,R3,R6,ONE/R2):
                if value[0].x==root*root:return(Point(root,ZERO),ORIGIN)
                if value[0].x==h*root*root:return(ORIGIN,Point(root,ZERO))
            raise ValueError('Unrecognized positive radical')
        raise ValueError('Unrecognized primary arithmetic data syntax')
    tree=ast.parse(line.strip()[1:-1].replace('[','(').replace(']',')'),mode='eval').body
    if not isinstance(tree,ast.BinOp) or not isinstance(tree.op,ast.Div):raise ValueError('Expected quotient')
    numerator,denominator=ev(tree.left),ev(tree.right)
    matches=[]
    for k in range(24):
        rotation=PHI0**k
        if mul((base*rotation,radical*rotation),denominator)==numerator:matches.append(k)
    if len(matches)!=1:raise ArithmeticError('Primary rotation identity did not have one match')
    return matches[0]


def main():
    h1=5-R3+3*R6+3*R2;h7=5+R3+3*R2-3*R6
    assert h1*h7==-14+48*R2
    assert PHI0.is_unit() and PHI1.is_unit() and PHI0**24==EAST
    assert not any(PHI0**k==EAST for k in (1,2,3,4,6,8,12))
    notebook_generator=Point(ONE/2+R6/6,R2/2-R3/6)
    assert notebook_generator==PHI0**4*Point(PHI1.x,-PHI1.y)
    report={'scope':'exact construction diagnostics only; no chromatic lower bound',
            'field':{'base_radicands':[2,3],'extension_squares':[[5,3,-1,3],[5,3,1,-3]],'degree':16},
            'program_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'positive_root_convention':True,'h1_h7_product':(h1*h7).to_json()[:4],
            'notebook_generator_equals_phi0_power4_times_conjugate_phi1':True,
            'embeddings':[],'finite_images':[],'rotations':[],'success_certificate':False}
    for signs in product((1,-1),repeat=2):
        images=[]
        for value in (h1,h7,h1*h7):
            for bit,sign in enumerate(signs):
                if sign<0:value=value.flip(bit)
            images.append(value.sign())
        report['embeddings'].append({'sqrt2_sign':signs[0],'sqrt3_sign':signs[1],
                                      'signs_h1_h7_product':images})
    for value in (h1,h7):assert value**4-20*value**3+1432*value-4412==ZERO
    for e2,e3 in product((1,-1),repeat=2):
        r2=(0,e2%5);r3=(0,2*e3%5);r6=finite_mul(r2,r3)
        a=finite_add(finite_add(finite_scale(-1,r3),finite_scale(3,r6)),finite_scale(3,r2))
        b=finite_add(finite_add(r3,finite_scale(3,r2)),finite_scale(-3,r6))
        roots_a=[x for x in product(range(5),repeat=2) if finite_mul(x,x)==a]
        roots_b=[x for x in product(range(5),repeat=2) if finite_mul(x,x)==b]
        ab=finite_mul(a,b);norm=(ab[0]**2-2*ab[1]**2)%5
        assert bool(roots_a)!=bool(roots_b) and norm==3
        report['finite_images'].append({'sqrt2':r2,'sqrt3':r3,'h1':a,'h7':b,
                     'h1_square_roots':roots_a,'h7_square_roots':roots_b,
                     'product':ab,'product_norm_to_F5':norm,'joint_extension_possible':False})
    source=Path('sources/voronov-dist-graphs/plane/series 2/vtx/rotation')
    lines=source.read_text().splitlines();report['primary_rotation_file_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
    for j,h,mask,row in ((1,h1,1,2),(7,h7,2,3)):
        w=EAST+PHI0**j;s=w.norm_squared();c=s+3
        assert h==(s-1)*(9-s) and h.sign()==1 and (s-1).sign()==1
        assert s**4-8*s**3+20*s*s-16*s+1==ZERO
        assert s*(-s**3+8*s*s-20*s+16)==ONE
        base=Point(c*w.x/(4*s),-c*w.y/(4*s));radical=Point(w.y/(4*s),w.x/(4*s))
        assert base.norm_squared()+h*radical.norm_squared()==ONE
        assert 2*(base.x*radical.x+base.y*radical.y)==ZERO
        edge_base=2*EAST-base*w;edge_radical=-(radical*w)
        assert edge_base.norm_squared()+h*edge_radical.norm_squared()==ONE
        assert 2*(edge_base.x*edge_radical.x+edge_base.y*edge_radical.y)==ZERO
        expected={}
        for axis in ('x','y'):
            values=['0']*16;values[:4]=getattr(base,axis).to_json()[:4]
            values[4*mask:4*mask+4]=getattr(radical,axis).to_json()[:4]
            expected[axis]=values
        match=source_check(lines[row-1],h,base,radical)
        assert match==(13 if j==1 else 19)
        report['rotations'].append({'j':j,'w':w.to_json(),'s':s.to_json()[:4],
                   'h':h.to_json()[:4],'psi_expected_coordinates':expected,
                   'unit_norm_identity_checked':True,'unit_edge_2_to_psi_w_checked':True,
                   'primary_rotation_row':row,'primary_equals_phi0_power_times_psi':match})
    Path('research/voronov_rotation_algebra.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rotations','embeddings','finite_images')}),flush=True)


if __name__=='__main__':main()
