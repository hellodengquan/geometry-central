# BFF: turning boundary curvature constraints into interior coordinates

An evaluation of `BFF` in `src/surface/boundary_first_flattening.cpp`
(header: `include/geometrycentral/surface/boundary_first_flattening.h`),
focused on four things the user asked about — topological premise, boundary
closure, scale degrees of freedom, and singular linear systems — and on one
distinction that the API blurs: **local solvability is not the same thing as a
flip-free map, and neither of them implies global injectivity.**

Every number below was measured against the repository's own source, compiled
with Eigen `SimplicialLDLT` (no SuiteSparse). Where an early hypothesis of mine
did not survive re-measurement it has been removed rather than reported; only
results that reproduced under a clean protocol appear here. Reproduction
commands are in [Appendix A](#appendix-a-reproducing-the-measurements).

---

## 1. What the code actually computes

Constructor (`:23-74`):

1. `GC_SAFETY_ASSERT(mesh.eulerCharacteristic() == 2 && mesh.nBoundaryLoops() == 1, ...)` (`:25`)
2. builds the cotan Laplacian (`:34`)
3. splits vertices into interior `A` / boundary `B` and block-decomposes (`:55`)
4. factorises only the interior block `Lii` (`:62`)
5. stores angle defects `Omega_i = 2pi - sum(angles)` (interior) and
   `Omega_b = pi - sum(angles)` (boundary) (`:67-73`).

Both public input paths funnel into `flattenFromBoth` (`:109-129`):

* `flattenFromScaleFactors(u)`: `k = dirichletToNeumann(u)` (`:87`) then `flattenFromBoth(u, k)`.
* `flattenFromExteriorAngles(k)`: `u = neumannToDirichlet(k)` (`:103`) then `flattenFromBoth(u, k)`.

`flattenFromBoth` does two things:

1. **Boundary layout.** `computeBoundaryPositions(u, k)` (`:147-197`) gives each
   boundary edge a direction from accumulated turning and a target length
   `l_e * exp((u_a+u_b)/2)`, then projects the length vector onto the
   exact-closure subspace with a 2x2 normal-equation solve.
2. **Interior.** Two scalar Laplace solves, `Liisolver->solve(-Lib * boundaryX/Y)`
   (`:114-115`) — a harmonic extension of the boundary.

So the interior coordinates are a Dirichlet solve, not an energy minimisation over
maps. Conformality, orientation and global validity all inherit from step 1.

---

## 2. Topological premise

```cpp
GC_SAFETY_ASSERT(mesh.eulerCharacteristic() == 2 && mesh.nBoundaryLoops() == 1,
                 "Input to BFF must be a topological disk");
```

`eulerCharacteristic()` (`src/surface/manifold_surface_mesh.cpp:515`) returns
`V - E + F + nBoundaryLoops`, and `genus()` divides the same bookkeeping number,
so for a genuine disk both read `2` / `0` and for a sphere too. Combined with
`nBoundaryLoops() == 1`, the pair rejects closed surfaces and higher-genus
surfaces with one loop. **The gap is that the test is per-mesh, not
per-component.**

A mesh made of a **disk plus a closed torus** has `chi = 2 + 0 = 2` and one
boundary loop, and is **accepted**:

```
torus+disk: declared chi=2  boundaryLoops=1  components=2  genus=0
    ACCEPTED by the topology assert; finite output = 1
```

The torus has no boundary data, so its vertices are filled by solving a Laplace
problem with zero right-hand side — an arbitrary constant field. The result is
finite and plausible-looking, and nothing warns. `genus()` likewise reports `0`.

Other unvalidated premises:

* **Connectivity** is never checked (same point).
* **Degenerate inputs.** A single-triangle "disk" (`nInterior == 0`) constructs
  and flattens; `Liisolver` is factorised on a `0x0` matrix, which Eigen reports
  as success. It survives only because the `nInterior == 0` path never touches
  `Lii`. A path that did would index empty vectors.
* **Index-order coupling.** `bIdx`/`iIdx` are built by counting over
  `mesh.vertices()`, and the boundary layouts are indexed by those counts. This
  matches `VertexData::toVector()` (`utilities/mesh_data.ipp:151`), so it is
  correct — but it is an undocumented assumption tying the algorithm to that
  container's iteration order.

**Verdict.** The gate catches the common error but is not a classification:
disconnected multi-component input is silently accepted and given arbitrary
values.

---

## 3. Boundary closure

### 3.1 It is an exact projection, applied unconditionally

```
l = target - N^{-1} T^T (T N^{-1} T^T)^{-1} (T target)
```

with `N = diag(|e_i|)` — the code fills `N` with the **length**, not its inverse
(`:168`). As written this is a legitimate projection of `target` onto
`{ l : sum l_i T_i = 0 }` under the `N^{-1}` inner product. Closure is therefore
exact to machine precision (`~1e-16`) for every input.

The point that matters: **the projection is not a small numerical fix, and it is
not conditional.** Whenever the requested metric does not already close — which
the discrete DtN does not guarantee — the projection silently moves the boundary
lengths, so the achieved scale factor is not the requested one.

### 3.2 Measured size of the correction

Equilateral 96-gon (`u` smooth in arc length, `k = DtN(u)`, mode 2):

| request | max per-edge \|achieved/target − 1\| |
|---|---|
| `amp = 0.05` | `~1e-11` |
| `amp = 0.15` | `~1e-11` |
| `amp = 0.30` | `~1e-11` |

i.e. for this symmetric case the projection is a no-op and achieved == target
exactly. But on a square grid, where boundary edges are not equilateral, the
same kind of request is no longer consistent with exact closure:

| mesh / request | achieved/target ratio | projection change / mean length |
|---|---|---|
| square 13, mode 2, amp 0.15 | `[1.000, 1.000]` | `~1e-15` |
| square 13, mode 3, amp 0.05 | `[0.977, 1.022]` | `0.023` |
| square 13, mode 3, amp 0.15 | `[0.925, 1.072]` | `0.077` |
| square 13, mode 3, amp 0.30 | `[0.843, 1.162]` | `0.162` |

And for a request that is a *genuine conformal map* — `u = log|phi'|` for the disk
automorphism `phi(z) = (z-a)/(1-az)`, `a = 0.4` — the distortion is severe:

```
Möbius test, n = 128 boundary verts
  achieved/target edge-length ratio in [0.4696, 2.0811]
  output boundary radius about its centroid: min=0.23734 max=0.41306  ratio=1.74
  (phi maps the unit circle to a circle, so a correct implementation gives ratio ~1)
```

The true image of the circle is a circle; the produced boundary has radius ratio
1.74. So **the requested scale factors are not achieved**, and the deviation
grows with how much the request departs from exact discrete closure. The
correction and the achieved-vs-target deviation track each other, which identifies
the projection as the mechanism.

Note the two behaviours must be kept apart: when the discrete target metric
already closes (equilateral mesh, symmetric modes) the projection is a no-op and
the request is met to `1e-11`; the distortion appears specifically for requests
that are inconsistent with exact discrete closure. The Möbius witness is the
clearest case because a genuine conformal map of the flat disk *is* such a request.

### 3.3 Direction convention

`phi` advances by `kBdy(tip)` (`:170`), and the direction assigned to edge `i` is
`phi` at its tail vertex — a half-step offset from the continuum (edge midpoint).
Changing this convention moves the closure defect by well under 1% in the tests,
so it is a discretisation detail, not the source of the projection's work. It
does mean the projection is compensating for two effects at once, so it cannot be
read as "fixing floating-point closure error".

### 3.4 Degenerate closure matrix

`(T N^{-1} T^T)` is 2x2. If all turning is concentrated at one vertex, all `T_i`
coincide and the matrix is singular. `.inverse()` on a singular `2x2` Eigen matrix
yields inf/nan with no check. Measured (96-gon, fraction of `2pi` given to `v0`):

```
frac=0.900  max|P|=3.112
frac=0.990  max|P|=4.138
frac=0.999  max|P|=4.263
frac=1.000  THREW: checkFinite() failure. Non-finite row vector entry [1] = nan
```

The failure surfaces later, inside the interior solve, as an opaque
`checkFinite` exception, not at the boundary projection where the degeneracy
actually is.

---

## 4. Scale degrees of freedom

### 4.1 A constant scale factor is a pure gauge

`dirichletToNeumann` is invariant to adding a constant to `u`; measured on a
96-gon:

| `u` offset | `max |k - k(0)|` | area | area / `e^{2c}` |
|---|---|---|---|
| 0.0 | — | 0.95 | 0.95 |
| 0.5 | identical k | 2.58 | 0.95 |
| 1.0 | identical k | 7.02 | 0.95 |

`k` is invariant across offsets (`max |k - k0| ≤ 2e-12`); only the global scale of
the output changes, exactly as `e^{2c}` in area. This is the gauge the BFF paper
uses for normalisation. This implementation **never documents it and never
normalises the output**: no centring, no unit area, no rotation. Two unrelated
calls produce maps at unrelated scale and position.

### 4.2 The two entry points disagree on that same gauge

`neumannToDirichlet` removes the gauge by force:

```cpp
double uMean = uBdy.mean(); // Ensure that u has mean 0
for (int i = 0; i < uBdy.size(); i++) uBdy(i) -= uMean;
```

Measured with `mean(u) = 0.05`:

```
mean(u)                 = 0.050000
mean(NtD(DtN(u)))       = 0.000000
scale difference between the two entry points = e^{0.05} = 1.0513
```

So `flattenFromScaleFactors(u)` and `flattenFromExteriorAngles(DtN(u))` — which
should describe the same map — differ by a global factor `e^{0.05}`. One path
keeps the caller's gauge; the other silently normalises.

### 4.3 `neumannToDirichlet` returns the *negative* preimage, up to the gauge

This is a concrete, measured inconsistency separate from the gauge. The leading
minus at `:139`:

```cpp
Vector<double> fullSolution = -Lsolver->solve(rhs);
```

negates the preimage. Measured on the curved disk, for several `u` (the third
column is the test: it is the spread of `NtD(DtN(u)) + u` about its own mean):

```
mean(u)=+0.000  mean(NtD+u)=+0.000  max|(NtD+u) - mean| = 5.8e-11
mean(u)=+0.050  mean(NtD+u)=+0.050  max|(NtD+u) - mean| = 6.7e-11
mean(u)=-0.200  mean(NtD+u)=-0.200  max|(NtD+u) - mean| = 3.7e-11
```

So the exact relation is

```
NtD(DtN(u)) = -u + mean(u)
```

(up to `1e-10`), i.e. the negative preimage with the gauge re-fixed to zero mean.
Decisively, on the round trip:

```
amp=0.05  |DtN( NtD(k)) - k| = 1.794e-01   |DtN(-NtD(k)) - k| = 1.839e-12
amp=0.20  |DtN( NtD(k)) - k| = 7.175e-01   |DtN(-NtD(k)) - k| = 1.689e-12
```

So `NtD` inverts `DtN` only after negating its own output.

Consequences:

* `flattenFromExteriorAngles` realises the scale field `−u` (relative to the
  scale-factor path) and its literal contract still holds — measured
  `max|achieved exterior angle − requested k| = 1.1e-10` — but it is not the
  inverse of `flattenFromScaleFactors`, and the two cannot be used
  interchangeably.
* Because it also re-centres the mean (§4.2), the two entry points additionally
  differ by a global factor `e^{mean(u)}`.

Together these mean `flattenFromScaleFactors(u)` and
`flattenFromExteriorAngles(DtN(u))` are not two routes to the same map: they
differ by a reflection of the scale profile and a global scale.

---

## 5. Singular and near-singular linear systems

There are three solves.

### 5.1 `Lii` — the interior block (reliable in practice)

`Lii` is the Dirichlet Laplacian on interior vertices. It is SPD for the tested
family even on strongly non-Delaunay meshes:

| shear grid | negative cotan weights | min eigenvalue of `Lii` |
|---|---|---|
| 0.5 | 64 / 208 | 0.338 |
| 2.0 | 64 / 208 | 0.872 |
| 5.0 | 64 / 208 | 4.075 |
| 20.0 | 64 / 208 | 61.33 |

The constructor shifts the *local copy* of `L` by `1e-12` (`:53`) before
decomposing, and `Liisolver` is built from that shifted `Lii`. This is safe here
because the cotan Laplacian is **dimensionless** (weights are cotangents, with no
length factor), so the shift does not depend on mesh scale — verified by uniformly
rescaling a disk by `1e-3` and `1e3`: the cotan weight range and the whole `Lii`
spectrum are unchanged (`min/max = 0.304482 / 7.695518` in all three cases), and
the normalised output shape is identical to `1e-14`. The absolute shift is thus
harmless in practice. (It is still conceptually odd to shift a matrix that is
merely semi-definite rather than pin one value, but it introduces no scale
dependence.)

### 5.2 `L` — the full Laplacian used by `Lsolver` (structurally singular)

`ensureHaveLSolver` builds the solver on the **member** `L`, which is the
*unshifted* cotan Laplacian:

```cpp
Lsolver.reset(new PositiveDefiniteSolver<double>(L));
```

`L` has a constant nullspace. It is only positive **semi**definite, and Eigen's
`SimplicialLDLT` reports success with a zero pivot:

```
||L*1||_inf = 0                    (to machine precision)
SimplicialLDLT(full L).info() = 0  (Success, no warning)
solve(rhs with sum != 0): finite=1  max|u| = 2.25e+15
```

The large magnitude appears only for an inconsistent right-hand side. The code
survives because `sum(kBdy) == 2pi` makes the RHS consistent and because the
result is re-centred at `:142-143`. Empirically even an inconsistent RHS of size
`1e-3` does not blow up the final map (measured `max|P| ≈ 1.514` for
`sum(k)-2pi ∈ {0, 1e-8, 1e-6, 1e-4, 1e-3}`), because the re-centring absorbs the
nullspace component. So this is **not** an observable failure today — but it is a
structural fragility: correctness rests on a consistency accident, a loose
absolute tolerance, and a downstream projection. Building `Lsolver` from the
shifted matrix (as `:53` already does for the local copy) would remove the
reliance entirely.

Also note `Liisolver` and `Lsolver` solve **different matrices** (shifted vs not),
so DtN and NtD are not exact inverses even before the sign issue in §4.3.

### 5.3 Input checks around the solves

* `abs(kBdy.sum() - 2*M_PI) < 1e-3` (`:98`) is an **absolute** tolerance with no
  resolution dependence. Summed turning is dimensionless and always equals `2pi`
  for a simple polygon, so an absolute `1e-3` says nothing about *relative*
  error: on a mesh with `n` boundary vertices each `k ≈ 2pi/n`, so `1e-3` is a
  loose check for fine meshes and a strict one for coarse ones. A `k` summing to
  `4pi` is correctly rejected (verified); one summing to just under `2pi + 1e-3`
  is accepted. A relative tolerance (`1e-6 * max(1, |k|_1)`) would be the right
  shape.
* `exp((u_a+u_b)/2)` (`:163`) overflows to `inf` for `u ≳ 709`, after which
  `checkFinite` throws with no context.
* There is **no check** on `edgeCotanWeights`, so a mesh that produces an
  indefinite `Lii` would fail with a bare "matrix is not positive definite" and no
  hint about which mesh property is responsible.

---

## 6. Local solvability vs. no flips vs. global injectivity

This is the distinction the API collapses.

* **Local solvability** — the linear systems have a solution. Guaranteed by the
  topology check plus successful factorisation.
* **No flips** — every output triangle has positive signed area.
* **Global injectivity / no overlap** — the map is one-to-one onto its image; no
  overlapping triangles, boundary simple.

The code provides the first and checks neither of the others. There is no
orientation test, no overlap test, and no status returned — `flatten*` returns a
plain `VertexData<Vector2>`.

### 6.1 Flips appear while everything still "works"

`flattenFromBoth` is a harmonic extension. A harmonic extension of a **convex**
boundary polygon is injective (Rado–Kneser–Choquet) — that is the guarantee the
paper leans on. The implementation never checks that the layout it produces is
convex. Measured on an equilateral 96-gon with `u = amp * sin(2 * mode)`:

| `amp` | flips | boundary convex? | boundary simple? | crossings | min interior angle |
|---|---|---|---|---|---|
| 0.2 | 0 | yes | yes | 0 | 3.040 |
| 0.4 | 0 | **NO** | yes | 0 | 3.003 |
| 0.5 | 0 | **NO** | yes | 0 | 2.985 |
| 0.6 | **34** | NO | yes | 0 | 2.967 |
| 0.8 | **132** | NO | **NO** | 2 | 3.062 |
| 1.0 | 32 | NO | NO | 2 | 3.025 |

The ordering is instructive. At `amp = 0.4–0.5` the boundary is already
non-convex but the interior is still fold-free. At `amp = 0.6` the layout is
non-convex and **34 triangles flip while the boundary is still a simple
polygon**. So "the boundary looks fine" is *not* evidence of an injective
interior, and flips occur with both solves returning finite answers and `Lii`
still SPD. Local solvability held; no flips did not.

### 6.2 Reflex (concave) target angles are handled correctly

A reflex vertex needs `k < 0`. Measured with orientation-corrected interior
angles (96-gon, `k` given at `v0`):

```
req k0        expected interior at v0    achieved    kind
 0.5                 2.6416               2.6416    convex
 0.0                 3.1416               3.1416    (straight)
-0.3                 3.4416               3.4416    reflex
-0.5                 3.6416               3.6416    reflex
-1.5                 4.6416               4.6416    reflex
```

and the total interior-angle sum equals `(n-2)pi` in every case. So the sign
bookkeeping at `:170` (`phi += kBdy`) is correct for `k < 0`; the earlier worry
that reflex targets are inverted does **not** hold. (Naive unsigned-angle probes
can look like an inversion because the boundary is traversed with the interior on
the right; it is not a defect.)

`k` values outside `(-pi, pi)` are still accepted without check, though a single
turning greater than `pi` is impossible for a simple polygon.

### 6.3 Global no-overlap is not guaranteed

From the table above, at `amp = 0.8` the boundary self-intersects (2 crossings)
and 132 triangles flip. Nothing reports this. A caller that does not
independently test triangle orientation will consume a folded, non-injective map
without any warning. Because the boundary crosses *after* flips begin, checking
only the boundary is insufficient too.

**Verdict.** The implementation guarantees exact boundary closure and the
existence of the harmonic extension. It does **not** guarantee the requested
scale factors (see §3), does not guarantee a simple boundary, and does not
guarantee a flip-free or globally injective interior. The gap between "the solver
returned" and "the map is valid" is wide and entirely undocumented.

---

## 7. Summary

| Property | Implied by the API | Actually guaranteed | Evidence |
|---|---|---|---|
| Input is a disk | yes (assert) | no — disk+torus accepted | `chi=2, loops=1, comps=2` |
| Boundary closes | yes | yes — exact projection | closure `~1e-16` |
| Interior solve exists | yes | yes for valid disks | `Lii` SPD, minEig > 0 |
| Requested scale factor achieved | yes | **no** — projection changes it | Möbius: ratio `[0.47, 2.08]` |
| Reflex `k` honoured | yes | **yes** | interior angles match exactly |
| Two entry points describe the same map | implied | **no** | `NtD∘DtN(u) = −u + mean(u)`; plus a `e^{mean(u)}` scale gap |
| Global scale normalised | no | only on the exterior-angle path | `mean(NtD)=0`, `mean(u)=0.05` |
| No triangle flips | implied | **no** | 34 flips at `amp=0.6`, simple boundary |
| Global injectivity | implied | **no** | 132 flips + 2 crossings at `amp=0.8` |
| Degenerate `k` handled | no | throws late, `NaN` | all turning at one vertex |
| Singular full solve | — | survives by accident | `‖L·1‖=0`, re-centred |

---

## 8. Recommendations (ordered by severity)

1. **Make `Lsolver` use the shifted matrix** (or pin a boundary vertex) as `:53`
   already does locally. Removes reliance on RHS consistency and the absolute
   `1e-3` tolerance.
2. **Fix or document the sign in `neumannToDirichlet`** (`:139`). As written, the
   DtN/NtD pair is not mutually inverse, and the two public entry points produce
   maps differing by a global scale and a sign. Either negate correctly or state
   the convention explicitly.
3. **Validate `Lii` positivity** and report the responsible mesh property
   (negative cotan weights) instead of a bare "matrix is not positive definite".
4. **Replace the topology assert** with a per-component test —
   `nConnectedComponents() == 1` in addition to the current conditions.
5. **Make the closure tolerance relative** (e.g. `1e-6 * max(1, |k|₁)`) instead
   of the fixed `1e-3`.
6. **Document the projection** in `docs/docs/surface/algorithms/parameterization.md`
   (`:31-33`, `:70-73`): target lengths are projected, so the achieved `u` can
   differ from the request; and state that the exterior-angle path normalises the
   global scale while the scale-factor path does not.
7. **Expose validity.** Returning a flip/injectivity status (or a minimum-area
   metric) would let callers avoid silently consuming folded maps. Nothing in the
   file currently makes that possible.

---

## Appendix A: reproducing the measurements

The probes compile against the workspace's own `src/surface/boundary_first_flattening.cpp`.
The library was built with Eigen 3.4 (`SimplicialLDLT`), SuiteSparse off:

```
cmake -S <repo> -B /tmp/gc_build \
      -DGC_EIGEN_LOCATION=<eigen-3.4> -DSUITESPARSE=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/gc_build -j8
```

Test disks:

* **square grid**, `m x m`, triangulated, flat or lifted to a dome
  `z = h sin(pi x) sin(pi y)` — has unequal boundary edge lengths at the corners;
* **equilateral polar disk** (centre + two rings, outer ring at radius 0.55) —
  all boundary edges equal, no negative cotan weights. Used to isolate effects
  from mesh anisotropy.

Important measurement detail: `computeBoundaryPositions` walks the boundary in
**halfedge order**, while `uB`/`kB` are indexed in **vertex-index order**
(`bIdx`). The two orders coincide for the square grid only partially and differ
for the polar disk, so any comparison must map through `bIdx`. Several
inconsistencies during this evaluation came from mixing those orders; the numbers
above were taken after fixing that.

## Appendix B: notes on numbers that did *not* reproduce cleanly

For honesty about the process: an initial reconstruction suggested that the raw
pre-projection closure defect was `O(0.1)` and resolution-independent for all
inputs, and that `flattenFromScaleFactors` never achieves its request. Both were
partly artifacts. The corrected picture is §3.2: for inputs whose discrete metric
already closes (equilateral mesh, symmetric modes) the projection is a no-op and
the request is met exactly; the deviation is real but appears specifically when
the requested metric is inconsistent with exact discrete closure (square grids,
strongly varying requests, and the Möbius witness). The Möbius numbers in §3.2 are
the clearest single witness and are the ones to trust.
