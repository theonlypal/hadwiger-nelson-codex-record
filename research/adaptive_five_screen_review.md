# Five-color screen design review

This review concerns the adaptive workflow before a five-color screen was
integrated. The preserved run program is
`runs/adaptive-mixed-tower-v2/program.py`; the inspected functional search is
`tools/direction_color_screen.py`.

For a rational coefficient vector whose denominators are coprime to five,
an additive functional `L` over `F5` gives a proper coloring for every
translate of every unit displacement `d` satisfying `L(d) != 0`. Equality of
the colors of a candidate center's old neighbors under some different
six-coloring cannot invalidate this explicit five-coloring.

The safe gate order is:

1. Try each already known functional on **every edge of the newly exported
   graph**. A direction list used to propose centers may omit additional
   unit edges found during export.
2. If a functional passes, emit its complete five-coloring and independently
   verify exact geometry and every edge before rejecting that graph as a
   target candidate.
3. If all known functionals fail, retain the offending exact unit
   displacements. These displacements lie in the corresponding functional
   kernels and identify what actually invalidated the old certificates.
4. Search another additive functional against all current directions. A
   verified coloring is another rejection certificate. UNSAT or timeout in
   this restricted functional search establishes no chromatic lower bound;
   only then is a full graph-coloring SAT diagnostic potentially useful.

Two integration details matter. A denominator divisible by five means this
coefficient map is unavailable, and an edge with zero residue makes every
linear residue coloring fail. Neither condition should abort the broader
exact graph search. They should produce explicit screen outcomes with the
relevant coefficient or edge retained.

The functional search fixes its first coefficient to one. This normalization
is sound when the real unit vector is an observed direction, because that
coefficient must then be nonzero and can be scaled. A general caller should
check that precondition or omit/change the normalization. Independent
verification must use explicit checks rather than relying on Python
assertions being enabled.

Finally, an old functional does not color the whole coordinate field merely
because it colors the current finite direction set. Enlarging a point set
can create previously unobserved unit displacements. The candidate graph,
its complete exported edges, and the actual coloring witness remain the
objects that determine whether it has been eliminated.
