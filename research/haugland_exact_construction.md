# Exact algebraic reconstruction

Use the complex coordinate `z=x+iy` and let

\[
t=e^{\pi i/21},\qquad
\Phi_{42}(X)=X^{12}+X^{11}-X^9-X^8+X^6-X^4-X^3+X+1.
\]

The field is `Q(t)[sqrt(5)]`, with the positive square root. A vertex is the
pair of rational coefficient arrays `(a,b)` representing
`a(t)+sqrt(5)*b(t)`, each polynomial of degree below 12. Complex conjugation
replaces `t` by `t^-1` and fixes `sqrt(5)`. Thus an edge is valid exactly when
`(z-w)*conjugate(z-w)=1` in this field. The implementation performs these
operations with integer numerators and a common denominator.

The cyclotomic polynomial is irreducible. The quadratic extension is proper:
the cyclotomic field of conductor 42 is unramified at 5, whereas `Q(sqrt(5))`
is ramified at 5, so it cannot be a subfield. Consequently these normal forms
also establish that distinct exported coordinates designate distinct points
in the selected complex embedding.

Haugland's complex heptagon coordinates can be written directly as

\[
P_j=\frac{t^{6j}}{t^{24}-t^{18}},\quad
Q_j=\frac{t^{6j+14}}{t^6-t^{36}},\quad
R_j=\frac{t^{6j+28}}{t^{12}-t^{30}},\qquad 0\le j<7.
\]

Reducing these fractions modulo `Phi42` gives the coefficient constants in
`hn/construct.py`. The ordered 84 unit directions are

\[
u_{2j}=t^j,\quad u_{2j+1}=v t^j,\qquad 0\le j<42,
\]

where

\[
7v=4t^{11}-4t^{10}-6t^9+t^8+t^7+5t^6+t^5-2t^4+3t^2+t+2.
\]

Here `v=Q0-R0`; all 84 norm identities are checked exactly. The distinguished
path endpoint is `B=2*t^7-1=i*sqrt(3)`. All 231 paths in the source Appendix A
end exactly at this value, with 740 distinct vertices and 3,985 accepted
edges. The two rotated copies have 1,066 vertices. Applying the exact unit
rotation

\[
\rho=\frac{7+(2t^7-1)\sqrt5}{8}=\frac{7+i\sqrt{15}}8
\]

produces 2,131 vertices and 12,530 accepted edges, matching the source.

The independently derived `t_points(n)` enumerates unordered half paths with
at most three steps. Exact complementary sums yield every step multiset
with total endpoint `B`; every submultiset sum is a vertex on an ordering of
that path. This gives exactly 1,042 vertices for T5 and 12,856 for T6. Starting
with T5, adding T6 vertices with at least seven T5 neighbors gives 1,294
vertices; taking its 7-core gives exactly the Appendix A vertex set and all
3,985 edges. No numerical separation lemma is used for vertex construction
or identification.

Replacing the 740-vertex core by all of T6 yields the search candidate in
`runs/t6-spindle`: the first union has 16,121 vertices and the spindle has
32,241 vertices and 229,108 accepted edges. The saved proper 6-coloring
rejects this candidate as a solution to the requested target.

Source: Jan Kristian Haugland, [arXiv:2608.04542v1](https://arxiv.org/html/2608.04542v1),
CC BY 4.0. Geometric edge acceptance is exact, but the numerical edge
candidate enumeration does not certify induced-graph completeness. An
edge-omitting subgraph would still suffice if proved non-6-colorable.
