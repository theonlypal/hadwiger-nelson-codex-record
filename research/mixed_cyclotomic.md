# Exact Haugland–Heule Minkowski sums

The candidate field is

```
F = Q(t)(sqrt(5),sqrt(33)),   t=exp(pi*i/21).
```

The defining cyclotomic polynomial, in ascending powers, is
`Phi_42(t)=1+t-t³-t⁴+t⁶-t⁸-t⁹+t¹¹+t¹²`. Its degree is 12.
Only 3 and 7 ramify in `Q(t)`. At a prime above 5 the valuation of 5 is
therefore 1, so 5 cannot be a square in `Q(t)`. Adjoining `sqrt(5)` doubles
the degree; the resulting compositum is unramified at 11. The valuation
of 33 at a prime above 11 is consequently 1, so 33 is not a square in that
compositum. Thus `[F:Q]=48`, with the independent basis

```
t^j * sqrt(5)^a * sqrt(33)^b,  0<=j<12, a,b in{0,1}.
```

This ramification argument ensures that coefficient equality in the
implementation is equality in a field, rather than equality in a reducible
formal quotient. Both radicals select their positive real roots, and complex
conjugation sends `t` to `t^-1` while fixing them.

`hn/cyclotomic_tower.py` implements canonical integer coefficient vectors
and a common positive denominator, exact multiplication and conjugation,
and exact inverses. The JSON protocol is:

```
{
  "schema": "hn-cyclotomic-tower-v1",
  "field": {"cyclotomic_order": 42, "radicands": [5,33], ...},
  "vertices": [{"parts": [[12 rational strings], [12 rational strings],
                          [12 rational strings], [12 rational strings]]}],
  "edges": [[u,v], ...],
  "metadata": {...}
}
```

Part masks 0,1,2,3 correspond respectively to `1,sqrt5,sqrt33,sqrt165`.
Every retained edge is accepted only if exact multiplication gives
`(z_u-z_v)*conjugate(z_u-z_v)=1`. Floating-point arithmetic proposes edge
candidates, and does not validate any edge or prove enumeration completeness.

## Heule directions and exact coupling

Put `B=2t^7-1=i*sqrt3` and define

```
hex = t^7,
half_theta3 = (sqrt33+B)/6,
theta3 = (5+B*sqrt33/3)/6,
theta4 = (7+B*sqrt5)/8.
```

The exact identities are `B²=-3`, `half_theta3²=theta3`, and every displayed
rotation has complex norm 1. In particular,
`B*sqrt33/3=i*sqrt11`, giving Heule's usual `theta3=(5+i*sqrt11)/6`.
The 30 unit directions `hex^j * half_theta3^k` for `0<=j<6`, `-2<=k<=2`,
together with the origin, give V31. This matches the construction from
[Heule, *Computing Small Unit-Distance Graphs with Chromatic Number 5*](https://arxiv.org/html/1805.12181).

The existing exact Haugland G1 and G3 coordinates lift by filling parts 2
and 3 with zero. Their construction originates in
[Haugland's paper](https://arxiv.org/html/2608.04542v4); the archived v1 HTML
is retained as the particular coordinate source used by `hn.construct`.
The mixed candidates are the distinct sums `G1+V31` and `G3+V31`. For each
V31 vertex the sum contains a translated Haugland seed, and for each
Haugland vertex it contains a translated V31. The summands therefore
produce unit edges in both directions across the point set. All further
proposed edges still undergo the same exact norm test.

The existing F25 ring map cannot extend to this cyclotomic ring: a root of
`Phi_42` in characteristic 5 has order 42, but `F25*` has order 24. This
is only an obstruction to that particular five-coloring method. It supplies
no lower bound on the chromatic number of these candidates.

## Entry points and verification

`hn.mixed_cyclotomic` exports `lift_haugland(seed)`, `v31_points()`,
`mixed_points(seed,max_vertices=66000)`, `graph_from_points(...)` and
`export(...)`. No truncation is performed if the vertex cap is exceeded.

```
.venv/bin/python -m hn.mixed_cyclotomic --seed g1 --output OUTPUT.json
.venv/bin/python -m hn.mixed_cyclotomic --seed g3 --output OUTPUT.json
```

`--lift-only` and `--v31-only` provide controls. The constructors verify
their exact rotation identities and distinct direction counts. Initial
arithmetic checks cover root order, radical squares, half-angle identity,
all rotation norms, generic inversion, and JSON round trips. The separate
verifier `tools/verify_tower.py` does not import this arithmetic module.
No graph in this note is claimed non-six-colorable; a solver result and
independently checked certificate are required for any such conclusion.
