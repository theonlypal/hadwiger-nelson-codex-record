# A separate multiquadratic coordinate family

This implementation uses exact real coordinates in

\[
K=\mathbb Q(\sqrt2,\sqrt3,\sqrt5,\sqrt{11}).
\]

Index `m` in a 16-entry coefficient array means
`product(sqrt(r[i]) for each set bit i of m)`, where `r=[2,3,5,11]`.
All square roots are positive. Multiplication of basis entries `a,b` gives
the entry `a XOR b`, multiplied by the product of radicands whose bits occur
in `a AND b`. Distinct products of these four primes represent independent
square classes, so these 16 elements form a rational basis of the field.
Consequently canonical rational arrays establish point equality exactly.

An exported vertex has separate arrays `x,y`. Each retained edge satisfies
`(x1-x2)^2+(y1-y2)^2=1` by rational coefficient equality. Vertex selection
inside a circle uses exact sign evaluation: rational intervals enclosing
each positive square root are refined until the sign is determined. No
floating-point comparison decides an algebraic equality or clipping bound.

## Primary-source construction

The starting construction is from Marijn J. H. Heule,
[Computing Small Unit-Distance Graphs with Chromatic Number 5,
arXiv:1805.12181v1](https://arxiv.org/html/1805.12181v1), Sections 2.1 and 4.1.
Using complex multiplication on real coordinate pairs, set

\[
h=(1/2,\sqrt3/2),\qquad q=(\sqrt{33}/6,\sqrt3/6).
\]

The 30 directions are `h^j*q^k` for `0<=j<6` and `-2<=k<=2`, and adjoining
the origin yields V31. Form V151 by keeping the points of V31+V31 whose
squared radius is at most one, then form V1939=V31+V151. The exact computed
cardinalities are 31, 151, and 1,939. The unit rotation

\[
\theta_4=(7/8,\sqrt{15}/8)
\]

produces the baseline union `V1939 union theta4*V1939`, which has 3,877
vertices and 26,814 retained exact unit edges. Its construction is a
baseline, not a claim of non-6-colorability.

For a smaller geometry check, the diamond
`{0,(1,0),h,(1,0)+h}` united with its copy rotated by
`theta3=(5/6,sqrt(11)/6)` gives the 7-vertex Moser spindle.

## Additional rotations

The explicitly checked norm-one multipliers

\[
r_1=(\sqrt2/2,\sqrt2/2),\qquad r_2=(1/3,2\sqrt2/3)
\]

introduce the independent square root of two. The code supports unions of
rotated copies and sums with a rotated six-direction unit star. These are
search families; entering a larger coordinate field does not establish an
increase in chromatic number.

Numerical squared-distance screening is permitted only to propose edges.
Every retained edge is checked exactly, and graph files explicitly decline
to claim induced-graph completeness. A verified non-6-colorable subgraph
would satisfy the target even if candidate screening omitted other edges.

Implementation: `hn/multiquadratic.py`, `hn/mq_construct.py`.
