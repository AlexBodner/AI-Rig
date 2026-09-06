# Line A — Coercive Lyapunov functionals (regularity via a monotone quantity)

Status legend (one tag per claim): **PROVED-HERE** (complete proof written out), **KNOWN** (cited, not reproved), **COMPUTED** (script plus printed output below), **CONJECTURE**, **FAILED** (attempted, with the exact step where it breaks). Citation tags: [checked via WebSearch] / [unverified] as in `00_problem.md`.

Files: this note and the script `line_a_lyapunov.py` (same directory; all tables below are pasted from its output, total runtime about 52 s on 4 CPUs). Everything concerns statement (B) of `00_problem.md` (unit torus `T^3 = R^3/Z^3`, `ν = 1`, `f = 0`); remarks about `R^3` are marked as such.

## Goal

Decide, as rigorously as possible, whether global regularity could follow from a *coercive Lyapunov functional*: a functional `Q` of the velocity field at a fixed time that (i) is nonincreasing along every smooth Navier–Stokes solution and (ii) controls a critical or subcritical norm (`L^3`, `Ḣ^{1/2}`, enstrophy, ...), so that the ESS criterion (no-go item 2.4) or Leray's `H^1` continuation criterion (2.5) would give regularity.

Outcome. For every functional in a large explicit class (Fréchet differentiable or convex, with an amplitude-homogeneous top part, or energy plus a lower-order part with amplitude asymptotics) monotonicity forces the functional to be *conserved by the 3D Euler equations* (Lemma A1 and its variants A1b–A1d, A1c'). In the classes where Euler invariants can be classified — radial Fourier-multiplier quadratic forms of Sobolev type (A2) and zeroth-order local integrals `∫ φ(u) dx` (A6) — the only invariants are energy-type, hence supercritical and not coercive. Every standard candidate (enstrophy, `L^3`, `Ḣ^{1/2}`, `∫|ω|^{3/2}`, `L^4`, `‖u‖_∞`, `‖ω‖_∞`, helicity of either sign, localized energy, `∫|u|^2 log(1+|u|)`) is refuted by an explicit trigonometric polynomial at which its Navier–Stokes time derivative is positive at `ν = 1` (A3, A5; exact rational arithmetic where the candidate is polynomial). The lemma passes the positive controls (2D, swirl-free axisymmetric, Burgers: their Lyapunov functionals are inviscid invariants, A1e). The line stops (Breakpoint) at the classification of *nonlocal* (in particular cubic) and *first-order local* Euler invariants, and it says nothing about non-monotone a priori bounds or trajectory functionals, which is where ESS-type results live.

## Setting and notation

- Domain `Ω = T^3 = R^3/Z^3` (statement (B)); `R^3` where stated. Fourier modes on `T^3`: `u(x) = Σ_{n ∈ Z^3} û_n e^{2πi n·x}`, so `∂_j` has symbol `2πi n_j` and `∫_{T^3} e^{2πi n·x} dx = δ_{n,0}`. Because of this convention, an integral of a product of `m` derivatives of trigonometric polynomials with Gaussian-integer coefficients is a rational number times `(2π)^m`; the script tracks the rational part exactly (Python `Fraction`) and the power of `2π` separately.
- `V := H^∞_σ(Ω) = {u ∈ C^∞(Ω; R^3) : div u = 0, u ∈ H^s(Ω) for every s ≥ 0}`. On `T^3` this is the set of all smooth periodic divergence-free fields (nonzero mean allowed). On `R^3` it is the natural class: the Schwartz class is *not* preserved by Navier–Stokes because the nonlocal pressure creates algebraic tails `|u(x,t)| ~ |x|^{-4}` (I recall this from work of Brandolese and of Dobrokhotov–Shafarevich [unverified]), whereas membership in every `H^s` is preserved.
- Navier–Stokes with `ν = 1`, `f = 0`, written as an evolution equation:

```text
∂_t u = Δu + B(u,u),      B(u,v) := −P((u·∇)v),      P = Leray projector onto divergence-free fields.
```

Euler is `∂_t u = B(u,u)`; the heat flow is `∂_t u = Δu`. The pressure of a field `u ∈ V` is `p = p(u) := −Δ^{-1} ∂_i ∂_j (u_i u_j)` (on `T^3` normalized to mean zero), so that `B(u,u) = −(u·∇)u − ∇p(u)`; in Fourier variables `p̂_n = −(n_i n_j / |n|^2) (u_i u_j)^_n` for `n ≠ 0`. The pressure depends on `u` only through `∇u` (it is unchanged by adding a constant vector to `u`).

- Vorticity `ω = curl u`, strain `S = (∇u + ∇u^T)/2` with `(∇u)_{ji} = ∂_j u_i`. For a functional `Q` with derivative `DQ(u)`, the **Euler flux** at `u` is `DQ(u)[B(u,u)]` and the **heat flux** is `DQ(u)[Δu]`; the Navier–Stokes derivative of `Q` at `t = 0` is their sum (Lemma A1, step 1). For a homogeneous `Q` of degree `k` these scale as `A^{k+1} F(u)` and `A^k G(u)` under `u ↦ Au`, so with `F(u) > 0` the **threshold amplitude** is `A* := −G(u)/F(u)`: the Navier–Stokes derivative of `Q` at data `Au` is positive iff `A > A*`. The product `A* ‖u‖_2` is invariant under rescaling `u`, and is the `L^2` norm of the data at which `Q` starts to increase.

- KNOWN inputs used without proof:

  - **K1 (local theory for Navier–Stokes in `V`).** For every `u^0 ∈ V` there are `T = T(u^0) > 0` and a unique `u` with `u ∈ C([0,T); H^s) ∩ C^1([0,T); H^{s−2})` for every `s ≥ 2`, solving `∂_t u = Δu + B(u,u)` classically, `u(0) = u^0`, `u(t) ∈ V` for `t ∈ [0,T)`. Sources I recall: Fujita–Kato 1964, Kato 1972 (`H^s(R^3)`), Temam's book for the periodic case; persistence of every `H^s` norm on the common lifespan follows from Leray's `H^1` continuation criterion (2.5) [unverified exact references; statement standard].
  - **K2 (local theory for Euler in `V`).** Same statement for `∂_t u = B(u,u)` with `C^1([0,T); H^{s−1})`, with a lifespan independent of `s` by Beale–Kato–Majda (2.7) (Kato 1972, Ebin–Marsden 1970 [unverified]).
  - **K3 (heat flow).** For `u^0 ∈ V`, `t ↦ e^{tΔ}u^0` is `C^1([0,∞); H^s)` for every `s`, with derivative `Δe^{tΔ}u^0`; as `t → ∞`, `e^{tΔ}u^0 → 0` in every `H^s` on `R^3`, and `e^{tΔ}u^0 → ū := ∫_{T^3} u^0 dx` (the mean, a constant field) in every `H^s` on `T^3`. (Elementary in Fourier variables; dominated convergence.)
  - **K4 (Danskin's theorem; directional derivatives of Lipschitz functions).** If `f(x,s)` is continuous on `K × (−ε,ε)` with `K` compact and `∂_s f` continuous, then `m(s) := max_{x∈K} f(x,s)` has right derivative `m'_+(0) = max_{x ∈ M} ∂_s f(x,0)`, `M := argmax_x f(x,0)`. For a locally Lipschitz function, existence of one-sided directional derivatives implies Hadamard directional differentiability. [unverified references; standard, e.g. Bonnans–Shapiro's book.]
  - **K5 (energy identity).** `⟨B(u,u), u⟩ = 0` for `u ∈ V`: `∫ ((u·∇)u)·u = ∫ (u·∇)(|u|^2/2) = −∫ (div u)|u|^2/2 = 0` and `∫ ∇p·u = −∫ p div u = 0`. (Used only in A1c' and A2.)

## Results

Summary of statuses: A1, A1b, A1c, A1c', A1d, A1f PROVED-HERE; A1e PROVED-HERE (by hand) + COMPUTED; A2 PROVED-HERE modulo an exact-rational certificate (COMPUTED, exact); A3 COMPUTED (exact where stated) with the helicity example PROVED-HERE; A4 discussion (KNOWN + PROVED-HERE parts); A5 PROVED-HERE modulo one COMPUTED nonzero constant, plus COMPUTED direct verification; A6 PROVED-HERE modulo an exact-rational certificate.

### A1. Amplitude-scaling lemma — **PROVED-HERE**

**Hypotheses.**

- **(H1) Differentiability.** There is `s_0 ≥ 0` such that for every `u ∈ V ∖ {0}` there is a linear functional `DQ(u): V → R`, bounded with respect to `‖·‖_{H^{s_0}}`, with `Q(u+h) − Q(u) − DQ(u)[h] = o(‖h‖_{H^{s_0}})` as `h → 0` in `H^{s_0}`, `h ∈ V`.
- **(H2) Initial monotonicity.** For every `u^0 ∈ V` there is `ε > 0` such that the K1 solution satisfies `Q(u(t)) ≤ Q(u^0)` for `0 < t < ε`. (This is weaker than the task's "`d/dt Q(u(t)) ≤ 0` at every time along every smooth solution", which implies it; the lemma uses nothing else about the time evolution.)
- **(H3) Absolute homogeneity.** There is `k > 0` with `Q(λu) = |λ|^k Q(u)` for all `λ ∈ R ∖ {0}` and all `u ∈ V`. (Every power of a norm, the enstrophy, the helicity, `∫|u|^p`, `∫|ω|^q` satisfy this.)

**Statement.** Under (H1)–(H3), for every `u ∈ V ∖ {0}`:

- (a) `DQ(u)[B(u,u)] = 0` (the Euler flux vanishes identically);
- (b) `DQ(u)[Δu] ≤ 0` (the heat flux is nonpositive);
- (c) `Q` is constant along every Euler solution of class K2, and nonincreasing along the heat flow.

**Proof.**

*Step 1 (flux inequality at `t = 0`).* Fix `u^0 ∈ V ∖ {0}` and let `u` be the K1 solution. Taking `s = s_0 + 2` in K1, `t ↦ u(t)` is `C^1([0,T); H^{s_0})` with `u'(0) = Δu^0 + B(u^0,u^0) =: h^0 ∈ V`. Put `h_t := u(t) − u^0 ∈ V`; then `h_t → 0` and `h_t / t → h^0` in `H^{s_0}` as `t ↓ 0`. By (H1),

```text
(Q(u(t)) − Q(u^0))/t = DQ(u^0)[h_t / t] + o(‖h_t‖_{H^{s_0}})/t  →  DQ(u^0)[h^0]      (t ↓ 0),
```

using boundedness of `DQ(u^0)` for the first term and `‖h_t‖_{H^{s_0}}/t ≤ C` for the second. By (H2) the left side is `≤ 0` for small `t > 0`, hence

```text
(F)   DQ(u^0)[Δu^0 + B(u^0,u^0)] ≤ 0     for every u^0 ∈ V ∖ {0}.
```

*Step 2 (homogeneity of the derivative).* Let `λ ∈ R ∖ {0}`, `u ∈ V ∖ {0}`, `h ∈ V`. For real `s → 0`, (H3) and (H1) at `u` give `Q(λu + sh) = |λ|^k Q(u + (s/λ)h) = |λ|^k Q(u) + |λ|^k (s/λ) DQ(u)[h] + o(|s|)`, while (H1) at `λu` gives `Q(λu + sh) = Q(λu) + s DQ(λu)[h] + o(|s|)`. Since `Q(λu) = |λ|^k Q(u)`, subtracting, dividing by `s` and letting `s → 0`:

```text
DQ(λu)[h] = |λ|^k λ^{-1} DQ(u)[h];   in particular  DQ(λu) = λ^{k−1} DQ(u) for λ > 0,  and  DQ(−u) = −DQ(u).
```

*Step 3 (amplitude scaling).* Apply (F) to `u^0 = λu`, `λ > 0`. Since `B` is bilinear, `B(λu,λu) = λ^2 B(u,u)`, and `Δ(λu) = λΔu`; by step 2,

```text
0 ≥ DQ(λu)[λΔu + λ^2 B(u,u)] = λ^k DQ(u)[Δu] + λ^{k+1} DQ(u)[B(u,u)]     for all λ > 0.
```

Dividing by `λ^{k+1}` and letting `λ → ∞` gives `DQ(u)[B(u,u)] ≤ 0`; dividing by `λ^k` and letting `λ ↓ 0` gives `DQ(u)[Δu] ≤ 0`, which is (b).

*Step 4 (sign symmetry).* Apply the first conclusion of step 3 to `−u ∈ V ∖ {0}`: `DQ(−u)[B(−u,−u)] ≤ 0`. But `B(−u,−u) = B(u,u)` (bilinearity) and `DQ(−u) = −DQ(u)` (step 2), so `−DQ(u)[B(u,u)] ≤ 0`. With step 3, `DQ(u)[B(u,u)] = 0`, which is (a).

*Step 5 (conservation).* Let `u` be an Euler solution of class K2 with `u(t) ∈ V ∖ {0}`. As in step 1 (now two-sided, since `t ↦ u(t)` is `C^1` into `H^{s_0}` on the open interval), `d/dt Q(u(t)) = DQ(u(t))[B(u(t),u(t))] = 0` by (a), so `Q(u(t))` is constant. (If `u(t_0) = 0` for some `t_0` then `u ≡ 0` by uniqueness in K2 and there is nothing to prove.) For the heat flow the same computation with K3 gives `d/dt Q(e^{tΔ}u^0) = DQ(e^{tΔ}u^0)[Δ e^{tΔ}u^0] ≤ 0` by (b), as long as `e^{tΔ}u^0 ≠ 0`. ∎

**What the proof used.** Only: `B` is bilinear (so `B(λu,λu) = λ^2 B(u,u)` and `B(−u,−u) = B(u,u)`), local existence of solutions that are `C^1` in time for all data in `V` (K1), and (H1)–(H3). It did **not** use the energy identity `⟨B(u,u),u⟩ = 0`, incompressibility, the Leray projector, or any estimate on `B`; it also holds verbatim with `Δ` replaced by `−(−Δ)^α` (hyperdissipation, no-go 2.10). This is why it transfers to Tao's averaged equation, to dyadic models, to Burgers and to the hyperdissipative equation (see Transfer check), and why it is a filter rather than a route to regularity.

**Reading.** Amplitude `λ` at viscosity `1` is the same as amplitude `1` at viscosity `1/λ`: `(F)` at `λu` is `λ^{k+1} DQ(u)[λ^{-1}Δu + B(u,u)] ≤ 0`. So "monotone at all amplitudes" is "monotone uniformly down to zero viscosity", and the lemma is the inviscid limit of the flux inequality. Since the Clay problem fixes `ν` but quantifies over all data, this is a genuine constraint on any Lyapunov functional for (B).

### A1b. Odd-degree functionals vanish — **PROVED-HERE**

Replace (H3) by sign homogeneity of odd integer degree: `Q(λu) = λ^k Q(u)` for all real `λ ≠ 0`, `k ≥ 1` odd; assume (H1) also at `u = 0` (so `Q` is continuous at `0` in `H^{s_0}`) and (H2). Then `Q(u) = Q(ū)` for every `u ∈ V` on `T^3` (`ū` = mean of `u`), and `Q ≡ 0` on `R^3`. In particular no odd-degree homogeneous functional can be monotone on mean-zero fields unless it vanishes.

*Proof.* Step 2 now gives `DQ(λu) = λ^{k−1} DQ(u)` for all `λ ≠ 0` (differentiate `Q(λ(u + (s/λ)h)) = λ^k Q(u + (s/λ)h)`). Step 3 with `λ < 0`: `0 ≥ λ^k DQ(u)[Δu] + λ^{k+1} DQ(u)[B(u,u)]`; as `λ ↑ 0`, `λ^k < 0`, so dividing by `λ^k` reverses the inequality: `DQ(u)[Δu] + λ DQ(u)[B(u,u)] ≥ 0`, hence `DQ(u)[Δu] ≥ 0`. With step 3 (`λ ↓ 0`), `DQ(u)[Δu] = 0` for all `u ∈ V ∖ {0}`. By step 5, `Q(e^{tΔ}u)` is constant in `t` (if `e^{tΔ}u = 0` for some `t` then `u = 0`, trivially in Fourier variables). By K3, `e^{tΔ}u → ū` (`T^3`) or `→ 0` (`R^3`) in `H^{s_0}`, and by continuity of `Q` at the limit (H1 at `ū` or at `0`), `Q(u) = Q(ū)` resp. `Q(u) = Q(0) = 0` (`Q(0) = λ^k Q(0)` forces `Q(0) = 0`). ∎

Example: on `T^3` the momentum `∫ u dx` (degree `1`, conserved by Navier–Stokes) has exactly the form `Q(u) = Q(ū)`; the corollary says these are the only monotone odd-degree functionals.

### A1c. Non-homogeneous functionals: the top-degree component — **PROVED-HERE**

Replace (H3) by:

- **(H3')** there are `k > 0` and, for each `u ∈ V ∖ {0}`, a linear functional `L(u): V → R` with `lim_{λ → +∞} λ^{1−k} DQ(λu)[h] = L(u)[h]` for every `h ∈ V`.

Then under (H1), (H2), (H3'): `L(u)[B(u,u)] ≤ 0` for all `u ∈ V ∖ {0}`; if moreover `L(−u) = −L(u)` (true whenever `Q` is even), then `L(u)[B(u,u)] = 0`.

*Proof.* (F) at `λu`, multiplied by `λ^{−k−1}`: `0 ≥ λ^{-1}·(λ^{1−k} DQ(λu)[Δu]) + λ^{1−k} DQ(λu)[B(u,u)] → 0·L(u)[Δu] + L(u)[B(u,u)]` as `λ → ∞`. The sign step is step 4 with `L` in place of `DQ`. ∎

*Examples.* (i) `Q = Σ_{j=1}^m Q_j` with each `Q_j` satisfying (H1) and (H3) with degrees `k_1 < ... < k_m`: then `DQ(λu) = Σ_j λ^{k_j − 1} DQ_j(u)`, so (H3') holds with `k = k_m`, `L = DQ_m`, and the top component `Q_m` must be Euler-conserved. E.g. `‖u‖_2^2 + ‖u‖_3^3` (top component `L^3`) and any `‖u‖_2^2 + ε‖u‖_{Ḣ^1}^2` (homogeneous of degree 2, Euler flux `ε` times the enstrophy flux) are excluded by A2/A3 below. (ii) `Q(u) = ∫ |u|^2 log(1 + |u|) dx` has `DQ(λu)[h] = λ ∫ g(λ|u|) u·h` with `g(r) = 2 log(1+r) + r/(1+r)`, so `λ^{1−k} DQ(λu)[h]` tends to `0` for `k = 2` and diverges for `k < 2`: there is no nonzero top-degree component and A1c is silent. A1c' handles it.

*Why the hypothesis is on `DQ(λu)` and not on `Q(λu)`.* `Q(λu) = λ^k Q_k(u) + o(λ^k)` does not imply (H3'): for `Q(u) = ‖u‖_2^2 + ‖u‖_2 sin(‖u‖_2^2)` one has `Q(λu)/λ^2 → ‖u‖_2^2` but `λ^{-1} DQ(λu)[h] = 2⟨u,h⟩ + 2λ‖u‖_2 ⟨u,h⟩ cos(λ^2‖u‖_2^2) + O(λ^{-1})`, which does not converge. (This `Q` is a function of the energy, hence Euler-conserved; the example only shows the hypothesis cannot be weakened to the value level.)

### A1c'. Energy plus a lower-order part with amplitude asymptotics — **PROVED-HERE** (uses K5)

Replace (H3) by:

- **(H3'')** there are `k > 0`, a function `a: (0,∞) → R` with `a(λ) = o(λ^k)` as `λ → ∞`, and for each `u ∈ V ∖ {0}` a linear functional `L(u): V → R` such that `lim_{λ→∞} λ^{1−k} ( DQ(λu)[h] − a(λ) ⟨u,h⟩ ) = L(u)[h]` for every `h ∈ V`.

Then under (H1), (H2), (H3''): `L(u)[B(u,u)] ≤ 0` for all `u ∈ V ∖ {0}`, and `L(u)[B(u,u)] = 0` if `L(−u) = −L(u)`.

*Proof.* Write `R_λ(u)[h] := λ^{1−k}(DQ(λu)[h] − a(λ)⟨u,h⟩) − L(u)[h]`, so `R_λ(u)[h] → 0` for each fixed `h ∈ V`, and `R_λ(u)` is linear in `h`. Then (F) at `λu` reads

```text
0 ≥ DQ(λu)[λΔu + λ^2 B] = a(λ)( λ⟨u,Δu⟩ + λ^2 ⟨u,B⟩ ) + λ^{k−1}( L(u)[λΔu + λ^2 B] + λ R_λ(u)[Δu] + λ^2 R_λ(u)[B] ),   B := B(u,u).
```

By K5, `⟨u,B⟩ = 0`, and `⟨u,Δu⟩ = −‖∇u‖_2^2`. Dividing by `λ^{k+1}`: `0 ≥ −a(λ)λ^{−k}‖∇u‖_2^2 + λ^{−1}( L(u)[Δu] + R_λ(u)[Δu] ) + L(u)[B] + R_λ(u)[B]`. As `λ → ∞` every term but `L(u)[B]` tends to `0` (`a(λ)λ^{-k} → 0`, `R_λ(u)[Δu] → 0`, `R_λ(u)[B] → 0`), so `L(u)[B(u,u)] ≤ 0`; the sign step is step 4. ∎

This is the first place where the energy identity enters (only to kill `a(λ) λ^2 ⟨u,B⟩`); A1c is the case `a ≡ 0`. It is applied in A5 to `∫|u|^2 log(1+|u|)`, where `a(λ) = λ(2 log λ + 1)`, `k = 2`.

### A1d. Convex functionals (norms, including `L^∞` norms) — **PROVED-HERE**

Replace (H1) by:

- **(H1')** `Q: V → R` is convex, even (`Q(−u) = Q(u)`), absolutely homogeneous of some degree `k ≥ 1`, and Hadamard directionally differentiable at every `u ∈ V ∖ {0}` with respect to `H^{s_0}`: for every `h ∈ V` the limit `Q'(u;h) := lim (Q(u + s h') − Q(u))/s` over `s ↓ 0`, `h' ∈ V`, `h' → h` in `H^{s_0}`, exists.

Then under (H1'), (H2): `Q'(u; B(u,u)) = Q'(u; −B(u,u)) = 0` and `Q'(u; Δu) ≤ 0` for every `u ∈ V ∖ {0}`.

*Proof.* Facts about one-sided derivatives of a convex `Q` (all elementary): (i) `Q'(u; a h) = a Q'(u;h)` for `a > 0`, `Q'(u;0) = 0`; (ii) sublinearity `Q'(u; h_1 + h_2) ≤ Q'(u;h_1) + Q'(u;h_2)` (from `Q(u + s(h_1+h_2)) ≤ (1/2)Q(u + 2sh_1) + (1/2)Q(u + 2sh_2)`); (iii) `Q'(λu; h) = λ^{k−1} Q'(u; h)` for `λ > 0` (from `Q(λu + sh') = λ^k Q(u + (s/λ)h')`); (iv) `Q'(−u; h) = Q'(u; −h)` (evenness). Chain rule: with `h_t = (u(t) − u^0)/t → h^0` in `H^{s_0}` (K1), Hadamard differentiability gives `(Q(u(t)) − Q(u^0))/t = (Q(u^0 + t h_t) − Q(u^0))/t → Q'(u^0; h^0)`, so (H2) yields `(F')`: `Q'(u; Δu + B(u,u)) ≤ 0` for all `u ∈ V ∖ {0}`. Apply `(F')` at `λu` and use (iii), (i): `λ^k Q'(u; Δu + λB(u,u)) ≤ 0`. By (ii), `λ Q'(u; B) = Q'(u; λB) ≤ Q'(u; Δu + λB) + Q'(u; −Δu) ≤ Q'(u; −Δu)` for all `λ > 0`, so `Q'(u;B) ≤ 0`. At `−u`, using `B(−u,−u) = B(u,u)` and (iv): `Q'(u; −B) ≤ 0`. Then `0 = Q'(u;0) ≤ Q'(u;B) + Q'(u;−B) ≤ 0`, forcing both to vanish. Finally `Q'(u;Δu) ≤ Q'(u; Δu + λB) + Q'(u; −λB) ≤ 0 + λ Q'(u; −B) = 0`. ∎

*Consequence for sup norms (via K4).* For `Q(u) = ‖u‖_{L^∞(T^3)}` (convex, even, degree `1`, Lipschitz on `C^0 ⊃ H^2`, hence Hadamard differentiable by K4): if `|u|` has a unique maximizer `x*`, then `Q'(u;h) = u(x*)·h(x*)/|u(x*)|`, and since `∇|u|^2(x*) = 0`, `u·(u·∇)u = (u·∇)|u|^2/2 = 0` there, so the Euler flux is `−u(x*)·∇p(x*)/|u(x*)|`. A1d says a monotone `‖u‖_∞` would need `u·∇p(u) = 0` at the speed maximum of every field. Likewise for the Beale–Kato–Majda quantity `Q(u) = ‖ω‖_{L^∞}` (Lipschitz on `C^1 ⊃ H^3`): at a unique maximizer `x*` of `|ω|`, `∇|ω| = 0` kills the transport term and the Euler flux is the stretching rate `ω·Sω/|ω|` at `x*`. Both are refuted numerically in A3.3.

### A1e. Positive controls: the lemma does not exclude the known Lyapunov functionals — **PROVED-HERE** (by hand) and **COMPUTED** (2D case)

A necessary condition is only useful if the systems that *do* have coercive Lyapunov functionals satisfy it. They do, and in each case the inviscid flux vanishes for the structural reason that is destroyed in 3D by vortex stretching or by pressure:

- 2D Navier–Stokes (no-go 2.8), `Q = ∫|ω|^p`, `ω` scalar: inviscid flux `= −p ∫ |ω|^{p−2} ω (u·∇)ω = −∫ (u·∇)|ω|^p = ∫ |ω|^p div u = 0`. For `Q = ‖ω‖_∞` (A1d): at a maximizer `∇|ω| = 0`, flux `0`. The script checks this numerically: five random planar fields fed to the 3D enstrophy code give Euler flux at most `4.4e-11` (round-off), because `ω·Sω = 0` for planar flows (A3.5); the exact planar triad certificate in A2 shows the mechanism at the level of a single triad.
- Axisymmetric swirl-free (2.9), `Q = ‖ω_θ/r‖_∞`: `f := ω_θ/r` is transported, `(∂_t + u·∇)f = 0` inviscidly, so at a maximizer of `|f|` the flux vanishes; likewise `∫|f|^p dx` (transport of a scalar by a divergence-free field).
- Viscous Burgers (T2 list), `Q = ‖u_i‖_∞`: inviscid flux at a maximizer `x*` of `|u_i|` is `−(u·∇)u_i(x*) · sign(u_i(x*)) = 0` since `∇u_i(x*) = 0`. There is no pressure to spoil it; in Navier–Stokes the term `−u·∇p` survives (A1d and row `Linf_u` of A3.3).

So the lemma is validated as a filter: every known coercive Lyapunov functional is inviscid-conserved, and every 3D candidate fails precisely at the stretching term `∫ ω·Sω` or at the pressure term. This is the T2 line of this attack.

### A1f. Two-parameter scaling on `R^3` (Newton polygon) — **PROVED-HERE**

On `R^3` the family `u_{λ,μ}(x) := λ u(μx)` (`λ, μ > 0`) lies in `V` whenever `u` does. Suppose `Q = Σ_{j=1}^m Q_j` where each `Q_j` satisfies (H1) and is bi-homogeneous: `DQ_j(λ u(μ·))[λ h(μ·)] = λ^{k_j} μ^{m_j} DQ_j(u)[h]` for all `λ, μ > 0` (every `∫|D^a u|^p`, `Ḣ^s` seminorm and the helicity are of this form). Write `G_j(u) := DQ_j(u)[Δu]`, `F_j(u) := DQ_j(u)[B(u,u)]`. Then, because `Δ(λu(μ·)) = λμ^2 (Δu)(μ·)` and `B(λu(μ·), λu(μ·)) = λ^2 μ (B(u,u))(μ·)`,

```text
DQ(u_{λ,μ})[Δu_{λ,μ} + B(u_{λ,μ},u_{λ,μ})] = Σ_j ( λ^{k_j} μ^{m_j + 2} G_j(u) + λ^{k_j + 1} μ^{m_j + 1} F_j(u) ).
```

If (H2) holds, this is `≤ 0` for all `λ, μ > 0`. Let `N ⊂ R^2` be the set of exponent pairs `(k_j, m_j + 2)` and `(k_j + 1, m_j + 1)`, and let `(a,b)` be a vertex of the convex hull of `N`, i.e. a point at which some linear functional `(α,β)` is uniquely maximized over `N`. Substituting `λ = t^α`, `μ = t^β` and letting `t → ∞`, the term with exponent `(a,b)` dominates, so its coefficient is `≤ 0`. If that term is an Euler-flux term `F_j`, then since `F_j(−u) = −F_j(u)` (step 2, `Q_j` even) we get `F_j ≡ 0`: **every Euler-flux term sitting at a vertex of the Newton polygon must be an Euler invariant.** A1 is the special case of the vertex with maximal `k`. (On `T^3` only `μ ∈ N` is available, which is why A5 uses A1c' instead.)

### A2. No `Ḣ^s` or `H^s` norm with `s ≠ 0` is an Euler invariant — **PROVED-HERE** modulo an exact-rational certificate (**COMPUTED**, exact)

*Setup.* For a real symbol `m: Z^3 ∖ {0} → R` put `Q_m(u) := Σ_{n≠0} m(n) |û_n|^2` (so `Ḣ^s` is `m(n) = (2π|n|)^{2s}` and `H^s` is `m(n) = (1 + (2π|n|)^2)^s`). Then `DQ_m(u)[h] = 2 Σ_n m(n) Re(conj(û_n)·ĥ_n)`, and the Euler flux is `F_m(u) = 2 Σ_n m(n) T_n(u)`, `T_n(u) := Re(conj(û_n)·(B(u,u))^_n)`, the energy transfer into wavevector `n`; for real fields `T_{−n} = T_n`. The Euler flux of the energy (`m ≡ 1`) is `0` for every `u` (K5).

*Triad.* Let `k = (1,0,0)`, `p = (0,1,1)`, `q = (−1,−1,−1)` (so `k + p + q = 0`, `|k|^2 = 1`, `|p|^2 = 2`, `|q|^2 = 3`) and let `u` be a real trigonometric polynomial supported on `{±k, ±p, ±q}`. The script enumerates all triples from this set and confirms that the only zero-sum triples are `±(k,p,q)` (no degenerate interactions such as `2k = −p`). Hence `(B(u,u))^_n` for `n` in the support is produced only by the interaction of the other two modes, the energy identity K5 restricted to this single triad reads `T_k + T_p + T_q = 0` (detailed conservation; printed exactly as `sum = 0`), and the Euler flux of `Q_m` is (summing over the six modes `±k, ±p, ±q`)

```text
F_m(u) = 4 ( m(k) T_k + m(p) T_p + m(q) T_q ) = 4 ( (m(p) − m(k)) T_p + (m(q) − m(k)) T_q ).
```

*Certificate (script, Part 2).* Two amplitude choices with Gaussian-integer coefficients (`u = 2 Re Σ a_n e^{2πi n·x}`):

```text
choice 1: a_k = (0, 1, i),     a_p = (1, i, −i),   a_q = (1, −1, 0):    (T_k, T_p, T_q) = (−1, 1, 0)·2π
choice 2: a_k = (0, 1+i, −1),  a_p = (1−i, 1, −1), a_q = (1, 1, −2):    (T_k, T_p, T_q) = (1, −3, 2)·2π
det [ (T_p, T_q)(choice 1) ; (T_p, T_q)(choice 2) ] = det [1 0; −3 2] = 2 ≠ 0.
```

(Each `a_n` is orthogonal to `n`, so the fields are divergence-free; the script asserts this.) Therefore the vectors `(T_p, T_q)(u)` span `R^2`, and for every symbol with `(m(p) − m(k), m(q) − m(k)) ≠ (0,0)` some `u` on the triad has `F_m(u) ≠ 0`. For `Ḣ^s` this pair is `(2π)^{2s}(2^s − 1, 3^s − 1)`, nonzero iff `s ≠ 0`; for `H^s` it is `((1+8π^2)^s − (1+4π^2)^s, (1+12π^2)^s − (1+4π^2)^s)`, nonzero iff `s ≠ 0`. By A1, **no power of an `Ḣ^s` or `H^s` norm with `s ≠ 0` — in particular neither the enstrophy (`s = 1`) nor the critical `Ḣ^{1/2}` norm — satisfies (H2)**, and neither does `a‖u‖_2^2 + b Q_m` for a single such norm `Q_m` and `b ≠ 0` (its Euler flux is `b F_m`). ∎ (Combinations of two or more distinct Sobolev orders can satisfy `m(k) = m(p) = m(q)` on this one triad, so the single-triad certificate is silent about them; the general statement "a radial multiplier form is Euler-invariant iff `m` is constant" needs triads connecting all `|n|` and is left to What would be needed next. Hand check of the certificate: `(B(u,u))^_n = −P_n Σ_{n_1+n_2=n} i2π (a_{n_1}·n_2) a_{n_2}` over the six modes, `T_n = Re(conj(a_n)·(B(u,u))^_n)`, `P_n v = v − n(n·v)/|n|^2`; each `T_n` is a sum of two terms.)

Explicit values for choice 1 (script, exact): energy `∫|u|^2 = 14`, enstrophy `∫|ω|^2 = 28 (2π)^2`, enstrophy Euler flux `2⟨−Δu, B(u,u)⟩ = 4 (2π)^3` — computed two independent ways, via the Leray projection and via `2∫ω·(ω·∇)u`, which agree exactly — enstrophy heat flux `2⟨−Δu, Δu⟩ = −128 (2π)^4`; `Ḣ^{1/2}` Euler flux `(2π)^2 (−4 + 4√2) ≈ 65.41`.

*2D positive control.* For the planar triad `k = (1,0,0)`, `p = (0,1,0)`, `q = (−1,−1,0)` with planar fields, the same computation gives `(T_k,T_p,T_q) = (1,−1,0)·2π` and `(7,−7,0)·2π` for two choices: `T_k + T_p + T_q = 0` **and** `|k|^2 T_k + |p|^2 T_p + |q|^2 T_q = 0` (enstrophy conservation in 2D), the `(T_p, T_q)` vectors are confined to a line and the corresponding determinant is `0` (printed). This is exactly where the 2D/3D difference sits: with a third dimension the enstrophy constraint is lost.

### A3. Explicit fields at which each candidate increases at `ν = 1` — **COMPUTED** (exact where stated); helicity **PROVED-HERE**

#### A3.1 Exact certificate on the triad (Part 2 of the script; exact rational arithmetic)

**Proposition (explicit enstrophy increase).** Let `u` be choice 1 of A2, i.e. with `θ_p := 2π(x_2 + x_3)`, `θ_q := 2π(x_1 + x_2 + x_3)`,

```text
u(x) = ( 2 cos θ_p + 2 cos θ_q ,   2 cos 2πx_1 − 2 sin θ_p − 2 cos θ_q ,   −2 sin 2πx_1 + 2 sin θ_p ).
```

Then `u ∈ V` (`div u = −4π sin θ_q − 4π cos θ_p + 4π sin θ_q + 4π cos θ_p = 0`), `∫|u|^2 = 14`, `∫|ω|^2 = 28(2π)^2`, and for the Navier–Stokes solution with data `A u` (`A > 0`, `ν = 1`) on the unit torus

```text
d/dt ∫|ω(t)|^2 dx |_{t=0} = −128 (2π)^4 A^2 + 4 (2π)^3 A^3 = 4 (2π)^3 A^2 (A − 64π).
```

So the enstrophy increases initially for every `A > 64π ≈ 201.06` (`A*‖u‖_2 = 64π√14 ≈ 752`). Status: COMPUTED (exact); the numbers are checkable by hand from the formulas in A2 and from `d/dt ∫|ω|^2 = −2∫|∇ω|^2 + 2∫ω·(ω·∇)u` (identities I1–I5 below).

Same field, other candidates (exact): `∫|u|^4 = 286`, `L^4` Euler flux `4∫|u|^2 u·B(u,u) = −(32/9)(2π)` (agrees exactly with `−4∫|u|^2 u·∇p`), heat flux `4∫|u|^2 u·Δu = −2320(2π)^2`; hence for the field `−u` (which flips the sign of the Euler flux) `d/dt ∫|u|^4 |_{t=0} = 2π A^4 ((32/9) A − 4640π) > 0` for `A > 1305π ≈ 4099.8`. For choice 2: enstrophy flux `4(2π)^3`, heat flux `−292(2π)^4`, threshold `146π`; `L^4` for `−u`: `A > (12987/160)·2π ≈ 510.0`.

#### A3.2 Candidate table on Gaussian-integer trigonometric polynomials (Part 3)

Fields: `u = Σ_{n ∈ M} 2 Re[ (α_n b_n^1 + β_n b_n^2) e^{2πi n·x} ]` with `M` the 13 wavevectors below (one representative per `±` pair, `|n|_∞ ≤ 1`), `b_n^1, b_n^2` the fixed integer basis of `n^⊥` printed by the script, and `α_n, β_n ∈ {−1,0,1} + i{−1,0,1}` drawn with seed `20260906` (60 fields). All polynomial quantities are evaluated exactly on a `16^3` grid (the integrands have degree `≤ 5` in each variable); the non-polynomial candidates (`L^3`, `∫|ω|^{3/2}`, `Flog`) use `N^3` quadrature with `N = 32, 64, 96` (convergence shown in A3.4). Candidate names: `enstrophy = ∫|ω|^2`, `L3 = ∫|u|^3`, `H12 = Σ 2π|n| |û_n|^2` (the `Ḣ^{1/2}` norm squared), `omega32 = ∫|ω|^{3/2}`, `L4 = ∫|u|^4`, `local_energy = ∫|u|^2 (1 + cos 2πx_1)`, `Flog = ∫|u|^2 log(1+|u|)`. `G` = heat flux, `F` = Euler flux at unit amplitude; the sign of the field is chosen so that `F > 0` (the Euler flux is odd in `u`); "best" minimizes `A*‖u‖_2`.

| candidate      | fields with nonzero Euler flux | best field    | `Q(u)`  | heat flux `G` | Euler flux `F` | `‖u‖_2` | `A* = −G/F` | `A* ‖u‖_2` |
| -------------- | ------------------------------ | ------------- | ------- | ------------- | -------------- | ------- | ----------- | ---------- |
| enstrophy      | 58/60                          | #24 (sign +1) | 8211.51 | -1.54608e+06  | 35719.2        | 9.79796 | 43.2842     | 424.096    |
| L3             | 60/60                          | #40 (sign -1) | 1945.69 | -502352       | 7306.91        | 11.4891 | 68.7504     | 789.882    |
| H12            | 60/60                          | #19 (sign -1) | 836.971 | -149283       | 1869.37        | 9.59166 | 79.8578     | 765.969    |
| omega32        | 60/60                          | #24 (sign +1) | 816.495 | -115224       | 2645.8         | 9.79796 | 43.5499     | 426.701    |
| L4             | 60/60                          | #40 (sign -1) | 32308   | -1.10499e+07  | 312725         | 11.4891 | 35.3341     | 405.958    |
| `local_energy` | 60/60                          | #20 (sign -1) | 94      | -15870.3      | 501.817        | 9.59166 | 31.6257     | 303.343    |
| Flog           | 60/60                          | #40 (sign -1) | 355.931 | -72185.5      | 322.331        | 11.4891 | 223.949     | 2572.97    |

(For `Flog` the column `A*` is meaningless because `Q` is not homogeneous; see A5 for the correct amplitude analysis.) For Gaussian-integer fields on these modes the enstrophy Euler flux is an integer multiple of `(2π)^3` (exact arithmetic; values up to `184`), and for two of the sixty fields it happens to be exactly `0` — an integer coincidence with no symmetry behind it (checked: neither field is even or odd under `x ↦ −x`); those two fields have nonzero `L^3` flux.

Explicit coefficients (`(α_n, β_n)` per mode, Gaussian integers written `(re, im)`), copied from the script output:

```text
modes M: [(1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1), (0,1,1), (1,-1,0), (1,0,-1), (0,1,-1), (1,1,1), (1,1,-1), (1,-1,1), (-1,1,1)]
bases (b^1, b^2): [((0,0,1),(0,-1,0)), ((0,0,-1),(1,0,0)), ((0,1,0),(-1,0,0)), ((0,0,-1),(1,-1,0)), ((0,1,0),(-1,0,1)), ((0,1,-1),(-1,0,0)), ((0,0,1),(-1,-1,0)), ((0,-1,0),(1,0,1)), ((0,-1,-1),(1,0,0)), ((0,1,-1),(-1,0,1)), ((0,-1,-1),(1,0,1)), ((0,1,1),(-1,0,1)), ((0,1,-1),(-1,0,-1))]
field #24 (enstrophy, omega32; sign +1): ((1,-1),(-1,0)) ((1,-1),(-1,0)) ((-1,-1),(0,0)) ((1,-1),(0,0)) ((1,0),(-1,-1)) ((1,0),(0,-1)) ((1,1),(1,-1)) ((1,0),(-1,-1)) ((0,-1),(0,-1)) ((-1,-1),(0,0)) ((0,0),(0,-1)) ((1,1),(1,-1)) ((0,-1),(0,0))
field #40 (L3, L4, Flog; sign -1): ((1,0),(-1,1)) ((1,1),(-1,1)) ((-1,1),(1,1)) ((-1,0),(1,1)) ((1,1),(1,0)) ((-1,1),(1,1)) ((1,-1),(1,1)) ((1,1),(-1,-1)) ((0,1),(1,1)) ((0,-1),(1,-1)) ((-1,1),(-1,-1)) ((0,1),(1,1)) ((0,0),(-1,-1))
field #19 (H12; sign -1): ((1,0),(-1,0)) ((-1,0),(-1,0)) ((1,0),(-1,-1)) ((-1,0),(-1,1)) ((1,0),(1,-1)) ((0,1),(1,0)) ((1,-1),(0,-1)) ((0,-1),(-1,-1)) ((-1,0),(-1,0)) ((-1,0),(-1,1)) ((0,0),(0,0)) ((0,-1),(0,-1)) ((-1,0),(0,1))
field #20 (local_energy; sign -1): ((0,0),(-1,-1)) ((1,1),(1,1)) ((1,-1),(-1,0)) ((1,0),(1,0)) ((0,-1),(1,0)) ((1,-1),(0,0)) ((-1,-1),(1,0)) ((0,-1),(-1,0)) ((0,1),(-1,1)) ((-1,1),(-1,1)) ((0,1),(-1,0)) ((-1,1),(0,-1)) ((-1,0),(0,-1))
field #11 (Linf_u; sign -1): ((-1,0),(1,-1)) ((-1,-1),(1,0)) ((1,0),(-1,1)) ((0,0),(0,0)) ((0,1),(0,1)) ((0,0),(0,-1)) ((1,-1),(0,0)) ((-1,-1),(1,-1)) ((1,-1),(0,0)) ((0,0),(0,1)) ((0,1),(1,1)) ((0,-1),(-1,0)) ((-1,1),(0,1))
field #8 (Linf_omega; sign +1): ((-1,-1),(0,1)) ((1,1),(0,1)) ((1,1),(-1,0)) ((1,0),(0,1)) ((-1,1),(0,-1)) ((1,1),(1,0)) ((-1,0),(0,0)) ((0,-1),(1,-1)) ((1,-1),(0,0)) ((1,-1),(1,1)) ((1,-1),(1,-1)) ((1,-1),(1,1)) ((0,1),(0,0))
```

Direct verification at `A = 2A*` (the Navier–Stokes derivative of `Q` at data `2A* u`, computed from scratch, must be positive; it equals `−(2A*)^k G > 0` for a degree-`k` homogeneous `Q`), with the quadrature resolution for the non-polynomial candidates:

| candidate      | field | `N` | `G`        | `F`       | flux at `A = 2A*` (should be `> 0`) |
| -------------- | ----- | --- | ---------- | --------- | ----------------------------------- |
| enstrophy      | #24   | 16  | -1546077.1 | 35719.231 | 1.1586419e+10                       |
| L3             | #40   | 32  | -502352.45 | 7306.9062 | 1.3059408e+12                       |
| L3             | #40   | 64  | -502352.53 | 7306.8331 | 1.3059808e+12                       |
| L3             | #40   | 96  | -502352.52 | 7306.8363 | 1.3059791e+12                       |
| H12            | #19   | 16  | -149283.38 | 1869.3658 | 3.8080771e+09                       |
| omega32        | #24   | 32  | -115224.45 | 2645.8004 | 93663621                            |
| omega32        | #24   | 64  | -115224.33 | 2645.6286 | 93672506                            |
| omega32        | #24   | 96  | -115224.35 | 2645.6013 | 93673993                            |
| L4             | #40   | 16  | -11049851  | 312725.02 | 2.7558247e+14                       |
| `local_energy` | #20   | 16  | -15870.324 | 501.81707 | 63493096                            |

Cross-checks printed by the script: the FFT path reproduces the exact triad values of A3.1 to all printed digits (`G = −199493.818... = −128(2π)^4`, `F = 992.2008... = 4(2π)^3`, `L^4` `G = −91589.93 = −2320(2π)^2`, `F = −22.34 = −(32/9)2π`, `Ḣ^{1/2}` `F = 65.40998`). The `L^3` and `∫|ω|^{3/2}` quadratures are stable to 5–6 digits between `N = 32` and `N = 96` (the integrands are Lipschitz, resp. bounded with `|ω|^{-1/2}` singularities at the isolated zeros of `ω`).

Smallest threshold found: a Nelder–Mead descent on the 52 real coefficients of the enstrophy candidate lowers `A*‖u‖_2` from `424.1` to `63.5` (optimized field: `‖u‖_2 = 18.63`, `G = −1.879e6`, `F = 5.514e5`, `A* = 3.41`), i.e. there are fields with `|n|_∞ ≤ 1` and `L^2` norm `63.5` whose enstrophy increases at `ν = 1`. (Not optimized further; the value is not a bound.)

#### A3.3 Sup norms via Danskin (A1d): `‖u‖_∞` and `‖ω‖_∞`

The maximizer `x*` of `|u|` (resp. `|ω|`) is found on a `48^3` grid and refined by BFGS with the analytic gradient of `|u|^2`; the gradient norm at `x*` and the second-largest grid local maximum are printed (the latter being clearly smaller shows the maximizer is unique up to a margin far larger than the grid error, which is `O(h^2 |D^2|u|^2|)`, about `0.1` in `|u|`). `G = u·Δu/|u|` at `x*` is the heat flux (always `≤ 0` at a maximum of `|u|^2`, because `Δ|u|^2 = 2u·Δu + 2|∇u|^2 ≤ 0` there); `F = −u·∇p/|u|` (for `‖u‖_∞`) resp. `F = ω·Sω/|ω|` (for `‖ω‖_∞`) is the Euler flux; the one-sided Navier–Stokes derivative at data `Au` is `A G + A^2 F`.

| candidate    | field         | `x*` (maximizer)               | max     | second grid local max | gradient norm at `x*` | `G`      | `F`     | `A*`    |
| ------------ | ------------- | ------------------------------ | ------- | --------------------- | --------------------- | -------- | ------- | ------- |
| `Linf_u`     | #11 (sign -1) | (0.873467, 0.271629, 0.269145) | 20.6466 | 18.8483               | 2.2e-07               | -1768.3  | 86.2941 | 20.4916 |
| `Linf_omega` | #8 (sign +1)  | (0.783314, 0.657194, 0.663356) | 198.212 | 180.758               | 5.0e-09               | -18616.5 | 9207.56 | 2.02187 |

So `‖u‖_∞` increases initially for `A > 20.5` (pressure gradient at the speed maximum), and `‖ω‖_∞` increases for `A > 2.02` (vortex stretching at the vorticity maximum). Neither is a Lyapunov functional, consistent with BKM (2.7) being a criterion rather than a bound.

#### A3.4 Helicity — **PROVED-HERE** (by hand) and confirmed numerically

Helicity `H(u) = ∫ u·ω dx` is an Euler invariant (KNOWN; here: `DH(u)[h] = 2∫ω·h`, and `2∫ω·B(u,u) = −2∫ω·((u·∇)u) = −2∫ω·(∇(|u|^2/2) − u×ω) = 0` by identities I7, I8 below), so it passes A1(a); it fails A1(b). Take the ABC field on the unit torus,

```text
u(x) = ( sin 2πx_3 + cos 2πx_2 ,  sin 2πx_1 + cos 2πx_3 ,  sin 2πx_2 + cos 2πx_1 ),     curl u = 2π u,     ∫|u|^2 dx = 3.
```

Then `B(u,u) = 0` (`(u·∇)u = ∇(|u|^2/2) − u × 2πu = ∇(|u|^2/2)` is a gradient, so its Leray projection vanishes), `Δu = −(2π)^2 u`, and along the Navier–Stokes solution (which is `e^{−(2π)^2 t} u`) `d/dt H = 2∫ω·Δu = −2(2π)^3 ∫|u|^2 = −48π^3`. The mirror image `v(x) := u(−x)` satisfies `curl v = −2π v`, `H(v) = −6π`, and `d/dt H(v) = +48π^3`. So neither `H` nor `−H` is nonincreasing, and helicity is indefinite (`H(u) = 6π`, `H(v) = −6π`). Script: helicity `18.849556 = 6π`, heat flux `−1488.301281 = −48π^3`, Euler flux `0.00e+00`; mirrored: `+1488.301281`; on 30 random fields the helicity Euler flux is `≤ 1e-12` and the heat flux is positive for 15 of 30.

#### A3.5 Positive control and local energy (Part 3)

- 2D control: five random planar fields (modes with `n_3 = 0`, `u_3 = 0`) give enstrophy Euler flux `≤ 4.4e-11` (zero up to round-off), while the `L^3` and `Ḣ^{1/2}` fluxes of the same planar fields are of order `10^2` — in 2D the vorticity `L^p` norms, not the velocity norms, are the invariants.
- Localized energy `∫|u|^2 χ`, `χ = 1 + cos 2πx_1` (a homogeneous degree-2 functional, the torus analogue of the quantity whose evolution is the CKN local energy inequality): Euler flux `∫(|u|^2 + 2p) u·∇χ ≠ 0` — field #20 (sign `−1`) gives `F = 501.8`, `A* = 31.6`. The local energy inequality is an inequality for a space-time distribution with an indefinite flux term, not a monotone functional.

#### A3.6 Sympy verification of the pointwise identities (Part 1)

All ten identities are verified for generic smooth `u(x,y,z)` (no divergence-free assumption; the `div u` terms are kept explicitly): I1 `ω·(u·∇)ω = (1/2)u·∇|ω|^2`; I2 `(1/2)u·∇|ω|^2 = div(u|ω|^2/2) − (|ω|^2/2) div u`; I3 `curl((u·∇)u) = (u·∇)ω − (ω·∇)u + ω div u` (the vorticity equation); I4 `ω·(ω·∇)u = ω·Sω`; I5 `ω·Δω = Δ(|ω|^2/2) − |∇ω|^2`; I6 `div((u·∇)u) = ∂_i u_j ∂_j u_i + u·∇(div u)` (pressure equation); I7 `(u·∇)u = ∇(|u|^2/2) − u×ω`; I8 `ω·(u×ω) = 0`; I9 for `φ(v) = F(|v|)` the traceless part of the Hessian is `(F'' − F'/r)(v v^T/r^2 − I/3)`; I10 `∇(F(|u|))·(u·∇)u = (u·∇)F(|u|)`. Together they give the enstrophy budget `d/dt ∫|ω|^2 = −2∫|∇ω|^2 + 2∫ω·Sω` used in A3.1 (script output: ten lines `OK`).

### A4. Why ESS-type `L^3` control cannot come from a monotone quantity — discussion (KNOWN + PROVED-HERE parts)

- KNOWN (2.4): if `sup_{t<T*} ‖u(t)‖_{L^3} < ∞` then the solution is smooth past `T*`; at a first singular time `‖u(t)‖_{L^3} → ∞` (Seregin), with a triple-logarithmic lower rate (Tao). So a functional `Q` on `V` with (H2) and `Q(u) ≥ Φ(‖u‖_{L^3})`, `Φ → ∞`, would prove (B). This is the cleanest imaginable regularity proof, and this note shows it cannot exist within the following classes:
  - `Q` a power of the `L^3` norm, or any function of it: the `L^3` norm is not monotone (A3.2, field #40: the Navier–Stokes derivative of `∫|u|^3` at data `A u` is `A^3 G + A^4 F` with `F = 7306.8 > 0`, positive for `A > 68.75`).
  - `Q` with a top-degree component (A1c) or convex (A1d): the top component must be an Euler invariant; the `L^3`-type candidates (`∫|u|^3`, `∫|u|^4`, `‖u‖_∞`, `∫|ω|^{3/2}`, `Ḣ^{1/2}`) are not (A2, A3), and neither is any `∫φ(u)` with non-quadratic `φ` (A6).
  - `Q = ‖u‖_2^2 + (lower-order part)`: A1c' applies whenever the lower-order part has amplitude asymptotics; A5 shows how the log-corrected energy fails.
- PROVED-HERE (shrinking bumps): the Euler invariants available in the covered classes are energy-type, and no function of the energy controls `‖u‖_{L^3}` (or the enstrophy) on `T^3`. Take `ψ ∈ C_c^∞(B(0,1/4); R^3)`, `u_0 := curl ψ` (divergence-free, compactly supported, not identically zero), and for `λ ≥ 1` let `u_λ(x) := λ u_0(λ(x − x_0))`, periodized; its support lies in a ball of radius `1/(4λ)` inside one cell, so `u_λ ∈ V(T^3)`. Substituting `y = λ(x − x_0)`: `∫_{T^3}|u_λ|^3 dx = ∫_{R^3}|u_0|^3 dy` (fixed), `∫_{T^3}|u_λ|^2 dx = λ^{-1}∫_{R^3}|u_0|^2 dy → 0`, and `∫_{T^3}|∇u_λ|^2 dx = λ ∫_{R^3}|∇u_0|^2 dy → ∞`. (The naive torus dilation `λu(λx)`, `λ ∈ N`, does not work: the cell has fixed volume, so it scales `∫|u|^3` by `λ^3`.) Helicity is scale-invariant (like `L^3`) but indefinite, bounded by `‖u‖_{Ḣ^{1/2}}^2` (Fourier: `|∫u·ω| ≤ Σ 2π|n||û_n|^2`), vanishes on all mirror-symmetric fields, and is not heat-monotone (A3.4). So within the classes of A1–A6 there is no coercive monotone quantity above `L^3`, and any such `Q` outside these classes would have to be a new Euler invariant of degree `3` (in the amplitude) that is not a function of energy and helicity — see Breakpoint.
- Trajectory functionals: `sup_{s ≤ t} ‖u(s)‖_{L^3}` and `∫_0^t ‖∇u(s)‖_2^2 ds` are monotone in `t` trivially; the content of such a quantity is its *boundedness*, which is the problem itself (for the dissipation integral it is the supercritical energy bound). ESS is a theorem about the trajectory obtained by rescaling, compactness and backward uniqueness for the heat operator (Carleman inequalities); it is not the sublevel set of a state functional, and if a coercive Lyapunov functional existed the ESS argument would be superfluous. The quantitative version (Tao 2019, triple-log growth) is likewise a trajectory statement.
- Non-monotone a priori bounds are outside the reach of A1: an inequality `sup_t Q(u(t)) ≤ C(Q(u^0))` for all data at `ν = 1` corresponds, after amplitude rescaling, to a bound at viscosity `1/λ` with a constant `C(λ^k Q(u^0))` that may degenerate as `λ → ∞`, so no inviscid conclusion follows. This is exactly the shape of the hyperdissipative regularity proofs (2.10): for `α ≥ 5/4` the enstrophy is still not monotone (A1 applies verbatim with `−(−Δ)^α`, and the Euler flux of A3.1 is unchanged), regularity there comes from a bounded, non-monotone Gronwall argument. The same holds for the small-data theory (2.15). Hence no known 3D regularity result proceeds via a coercive Lyapunov functional; the only ones that do (2D, swirl-free axisymmetric, Burgers) are exactly those with a transported quantity (A1e).

### A5. `∫ F(|u|)` with slowly growing `F`, viscosity-dependent and local-energy candidates — **PROVED-HERE** modulo one **COMPUTED** constant, plus direct **COMPUTED** verification

**Setting.** `Q_F(u) := ∫_{T^3} F(|u|) dx` with `F(r) = r^2 log(1+r)` — the natural "energy plus a logarithm" candidate: it grows slower than every `∫|u|^{2+ε}`, so A1c is silent (A1c, example (ii)), and it is not in the class of A1f on the torus. `Q_F` is Fréchet differentiable on `H^2 ⊂ L^∞` with `DQ_F(u)[h] = ∫ g(|u|) u·h`, `g(r) := F'(r)/r = 2 log(1+r) + r/(1+r)`, which is increasing and satisfies `λ g'(λ r) → 2/r` as `λ → ∞`.

**Claim A5.1 (structure of the fluxes).** For `u ∈ V ∖ {0}` and `λ > 0`,

```text
DQ_F(λu)[B(λu,λu)] = λ^3 ∫ g(λ|u|) u·B(u,u) dx = λ^4 ∫ p(u) g'(λ|u|) u·∇|u| dx,
DQ_F(λu)[Δ(λu)]    = −λ^2 ∫ g(λ|u|) |∇u|^2 dx − λ^3 ∫ g'(λ|u|) |u| |∇|u||^2 dx   ≤ 0.
```

*Proof.* `∫ g(|v|) v·(v·∇)v = ∫ (v·∇)F(|v|) = 0` (I10 and `div v = 0`) and `−∫ g(|v|) v·∇p = ∫ p div(g(|v|)v) = ∫ p (v·∇)g(|v|)` give the first line at `v = λu`. The second line is `∫ g(|v|)v·Δv = −∫ ∂_j(g(|v|)v_i)∂_j v_i` with `∂_j(g(|v|)v_i) = g ∂_j v_i + g'(|v|)(∂_j|v|) v_i` and `v_i ∂_j v_i = |v|∂_j|v|`. (All integrands are bounded: `|v·∇|v|| ≤ |v||∇v|`, and `|v| ∈ W^{1,∞}`.) ∎

**Claim A5.2 (amplitude asymptotics; A1c' applies).** For `u, h ∈ V`,

```text
λ^{-1} DQ_F(λu)[h] − (2 log λ + 1)⟨u,h⟩  →  L(u)[h] := 2 ∫ log|u| u·h dx      (λ → ∞),
```

so (H3'') holds with `k = 2`, `a(λ) = λ(2 log λ + 1) = o(λ^2)`. Moreover `L(u)[B(u,u)] = 2Π(u)`, where

```text
Π(u) := ∫ p(u) u·∇ log|u| dx = ∫ p(u) u_j u_i ∂_j u_i / |u|^2 dx,     Π(−u) = −Π(u).
```

*Proof.* `g(λ r) − 2 log λ − 2 log r − 1 = 2 log(1 + 1/(λr)) − 1/(1+λr) =: ρ_λ(r) → 0` for `r > 0`, and `|ρ_λ(r)| ≤ 2 log(1 + 1/r) + 1` for `λ ≥ 1`, so `|ρ_λ(|u|) u·h| ≤ (2|u| log(1 + 1/|u|) + |u|)|h| ≤ (2 + |u|)|h| ∈ L^1` (using `r log(1+1/r) ≤ 1`); dominated convergence gives the limit, with `log|u| u ∈ L^∞`. For `h = B(u,u) = −(u·∇)u − ∇p`: `∫ log|u| u·(u·∇)u = ∫ (u·∇)[ (|u|^2/2) log|u| − |u|^2/4 ] = 0` (check: `d/ds[(s^2/2) log s − s^2/4] = s log s`), and `−∫ log|u| u·∇p = ∫ p div(u log|u|) = ∫ p u·∇ log|u|`, where `u·∇ log|u| = u_j u_i ∂_j u_i/|u|^2` is bounded (it vanishes quadratically at zeros of `u`). Oddness of `Π` is clear (`p` is even in `u`). ∎

**Claim A5.3 (conclusion).** By A1c' (with `L(−u) = −L(u)`), if `Q_F` satisfied (H2) then `Π(u) = 0` for every `u ∈ V`. The script computes, for the explicit field `w :=` field #40 with sign `−1` of A3.2,

| `N` | `‖∇w‖_2^2` | `Π(w)`  | Euler flux `/ A^3` at `A = 1e3` | at `A = 1e5` | heat flux `/ (A^2 log A)` at `A = 1e5` |
| --- | ---------- | ------- | ------------------------------- | ------------ | -------------------------------------- |
| 32  | 11448.7    | 167.435 | 335.369                         | 335.382      | -29088.5                               |
| 64  | 11448.7    | 167.741 | 335.31                          | 335.323      | -29088.5                               |
| 96  | 11448.7    | 167.659 | 335.315                         | 335.328      | -29088.5                               |

so `Π(w) ≈ 167.66 ≠ 0` (COMPUTED; stable to three digits under refinement, and the large-amplitude Euler flux divided by `A^3` converges to `2Π = 335.3` as A5.2 predicts). Hence **`Q_F = ∫|u|^2 log(1+|u|)` is not a Lyapunov functional on `T^3`** (PROVED-HERE modulo the COMPUTED value `Π(w) ≠ 0`). The heat flux only grows like `A^2 log A` (second row of A5.1 with `g(A r) ~ 2 log A`), so for `A` large the cubic Euler term wins; direct evaluation (`N = 96`; prediction = the four-term asymptotic formula `2ΠA^3 − A^2[(2 log A + 1)‖∇w‖_2^2 + 2C_1 + 2C_2]` with `C_1 = ∫ log|w| |∇w|^2 = 25911.8`, `C_2 = ∫|∇|w||^2 = 4002.46`, obtained from A5.1–A5.2 by the same dominated convergence, error `o(A^2)`):

| `A`  | heat flux    | Euler flux  | total flux   | prediction   |
| ---- | ------------ | ----------- | ------------ | ------------ |
| 1    | -72185.5     | 322.296     | -71863.2     | -70941.9     |
| 10   | -1.24092e+07 | 334014      | -1.20752e+07 | -1.20647e+07 |
| 100  | -1.76733e+09 | 3.35197e+08 | -1.43213e+09 | -1.43192e+09 |
| 200  | -7.704e+09   | 2.6821e+09  | -5.0219e+09  | -5.02127e+09 |
| 300  | -1.81694e+10 | 9.05268e+09 | -9.11675e+09 | -9.11558e+09 |
| 400  | -3.3355e+10  | 2.14589e+10 | -1.18961e+10 | -1.18943e+10 |
| 500  | -5.33945e+10 | 4.19127e+10 | -1.14818e+10 | -1.14793e+10 |
| 700  | -1.08428e+11 | 1.15011e+11 | 6.58304e+09  | 6.58652e+09  |
| 1000 | -2.29448e+11 | 3.35315e+11 | 1.05867e+11  | 1.0587e+11   |
| 3000 | -2.29143e+12 | 9.05374e+12 | 6.76231e+12  | 6.76215e+12  |

The Navier–Stokes derivative of `Q_F` at data `A w` is negative up to `A = 500` and positive from `A = 700` on (the asymptotic formula puts the crossover near `A ≈ 655`); at `A = 1000` the value `1.0587e11` is stable to five digits for `N = 48, 64, 96, 128` (printed), and the asymptotic prediction agrees to four digits. (COMPUTED.) At small amplitude the heat term dominates (the `A = 1..100` rows), which is why a scan limited to moderate amplitudes would wrongly suggest monotonicity.

**Remarks (what escapes and what does not).**

- Any `F ∈ C^2` that is not quadratic gives a `Q_F` that is not an Euler invariant (A6 with `φ(v) = F(|v|)`: the traceless Hessian is `(F'' − F'/r)(v̂v̂^T − I/3)`, I9). Whether the viscous term can nevertheless dominate at all amplitudes is a case-by-case question: settled here for `F` with a top-degree component of degree `≠ 2` (A1c), for `F = r^2 log(1+r)` (A1c'), and — by the same computation, provided the dominated-convergence bounds hold, i.e. `λ g'(λ r) ≤ C/r` uniformly in `λ ≥ 1` and `g(λ r) − g(λ) = O(1 + |log r|)` — for every `F` with `g = F'/r` increasing, `g(λ) = o(λ)` and `λ g'(λ r) → c/r`, `c ≠ 0` (for `F = r^2 log(1+r)`: `λ g'(λr) ≤ 3λ/(1+λr) ≤ 3/r`). Not settled: `F` whose growth oscillates between regimes (no amplitude asymptotics), where none of A1c, A1c', A1f applies — these are outside every lemma here.
- Viscosity-dependent candidates: with `ν = 1` fixed, "depends on `ν`" is vacuous; what matters is amplitude dependence, and the amplitude is not at the disposal of the functional (the Clay problem quantifies over all data). A functional that adapts its cutoff to the solution (e.g. `Q(u) = ∫ F_{‖u‖_2}(|u|)`) is still a functional on `V` and is covered by A1/A1c' as soon as it is differentiable with amplitude asymptotics of its derivative.
- Time-dependent weights `Q(t,u)`, functionals of `(u,p)` treated as independent, functionals of the whole trajectory, and non-monotone bounded quantities are not covered (A4, last item). Non-differentiable non-convex functionals (e.g. `‖u‖_{BMO^{-1}}` is a norm hence convex and covered; `Ḃ^{-1}_{∞,∞}` likewise) are not covered when convexity also fails.
- The `R^3` dilation route (A1f) is not needed on the torus but gives, on `R^3`, a stronger statement: the two-parameter family `λ u(μ·)` makes every vertex Euler-flux term an invariant.

### A6. Zeroth-order local functionals: `∫ φ(u) dx` is an Euler invariant only if `φ` is quadratic — **PROVED-HERE** modulo an exact certificate (**COMPUTED**, exact)

**Theorem.** Let `φ ∈ C^2(R^3; R)` and `Q(u) := ∫_{T^3} φ(u(x)) dx` on `V` (all smooth periodic divergence-free fields, mean allowed). If `DQ(u)[B(u,u)] = 0` for every `u ∈ V`, then `φ(v) = a|v|^2 + b·v + c` for constants `a ∈ R`, `b ∈ R^3`, `c ∈ R`. Conversely these are Euler (and Navier–Stokes) invariants up to the energy dissipation.

*Proof.* (1) *Flux formula.* `DQ(u)[h] = ∫ ∇φ(u)·h`. With `h = B(u,u) = −(u·∇)u − ∇p`: `−∫ ∂_iφ(u) u_j ∂_j u_i = −∫ (u·∇)(φ(u)) = ∫ φ(u) div u = 0`, and `−∫ ∂_iφ(u) ∂_i p = ∫ p ∂_i(∂_iφ(u)) = ∫ p ∂_{ij}φ(u) ∂_i u_j` (integration by parts on the torus, `φ ∈ C^2`). Hence

```text
DQ(u)[B(u,u)] = ∫_{T^3} p(u) tr( D^2φ(u) ∇u ) dx,     (∇u)_{ji} = ∂_j u_i.
```

(2) *Only the traceless Hessian matters.* Write `D^2φ(v) = (Δφ(v)/3) I + H_0(v)` with `H_0` symmetric and traceless. Then `tr((Δφ/3) I ∇u) = (Δφ(u)/3) div u = 0`, so the flux is `∫ p(u) tr(H_0(u) ∇u)`.

(3) *Perturbation around a constant field.* Suppose `H_0(v_0) ≠ 0` for some `v_0 ∈ R^3`. Let `w ∈ V` be mean-zero and put `u_ε := v_0 + εw` (a constant vector is in `V` on the torus). Since `p` depends only on `∇u` and `∇u_ε = ε∇w`, `p(u_ε) = ε^2 p(w)` (explicitly: `∂_i∂_j(v_{0,i} w_j) = v_{0,i} ∂_i div w = 0`). Hence `DQ(u_ε)[B(u_ε,u_ε)] = ε^3 ∫ p(w) tr(H_0(v_0 + εw) ∇w) dx`. As `ε → 0`, `H_0(v_0 + εw(x)) → H_0(v_0)` uniformly in `x` (`H_0` is continuous, `w` bounded), so

```text
ε^{-3} DQ(u_ε)[B(u_ε,u_ε)]  →  tr( H_0(v_0) M(w) ),      M(w)_{ij} := ∫_{T^3} p(w) ∂_j w_i dx.
```

Since `H_0(v_0)` is symmetric and traceless, only `N(w) := (M(w) + M(w)^T)/2` enters, and `tr N(w) = ∫ p div w = 0`, so `N(w)` lies in the 5-dimensional space `Sym_0` of traceless symmetric matrices. If the flux vanishes identically then `tr(H_0(v_0) N(w)) = 0` for every `w`.

(4) *Certificate.* The script exhibits five mean-zero Gaussian-integer fields `w^{(1)}, ..., w^{(5)}` (coefficients printed, modes `(1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1), (0,1,1), (1,−1,0), (1,0,−1), (0,1,−1)`) and computes `N(w^{(m)})` exactly, in units of `2π`, with coordinates `(N_11, N_22, N_12, N_13, N_23)`:

```text
N(w^(1)) = (-6, 4, 5, 6, 5)
N(w^(2)) = (-10, 4, -1, 4, -1)
N(w^(3)) = (2, 4, -2, -6, -10)
N(w^(4)) = (-2, -2, -1, -2, -7)
N(w^(5)) = (-6, -10, 12, 10, 14)
det = 7936 ≠ 0.
```

So the `N(w^{(m)})` span `Sym_0`, and `tr(H_0(v_0) N) = 0` for all `N ∈ Sym_0` forces `H_0(v_0) = 0` (the trace pairing is nondegenerate on `Sym_0`), a contradiction. Hence `H_0 ≡ 0`.

(5) *Integration.* `D^2φ(v) = c(v) I` for all `v` means `∂_1∂_2 φ = ∂_1∂_3 φ = ∂_2∂_3 φ = 0` and `∂_1^2 φ = ∂_2^2 φ = ∂_3^2 φ`. From the mixed derivatives, `φ(v) = f_1(v_1) + f_2(v_2) + f_3(v_3)` (write `φ = α(v_1,v_3) + β(v_2,v_3)` from `∂_1∂_2 φ = 0`, then split each by the remaining conditions); then `f_1''(v_1) = f_2''(v_2) = f_3''(v_3)` for all arguments forces a common constant `2a`, so `φ(v) = a|v|^2 + b·v + c`. The converse is K5 plus conservation of the mean. ∎

**Corollary (no coercive local Lyapunov functional).** If `Q(u) = ∫φ(u)` with `φ ∈ C^2` satisfies (H1)–(H2) and either (H3) or (H3') or (H3''), then by A1/A1c/A1c' its top part `Q_k` is an Euler invariant; if `Q_k` is itself of the form `∫φ_k(u)` then `φ_k` is quadratic, so `Q_k ≤ C(1 + ‖u‖_2^2)`, and by the shrinking-bump argument of A4 no such `Q` controls `‖u‖_{L^3}` or the enstrophy. In particular `∫|u|^p` (`p ≠ 2`), `∫ F(|u|)` (non-quadratic `F` with a top-degree component), and every `∫φ(u)` with a non-quadratic top-degree part are excluded.

*Remark (`R^3`).* Constants are not in `V(R^3)`; the same proof works there whenever `H_0(0) ≠ 0` (take `v_0 = 0`), which covers every `φ` that is not quadratic near `0`. The general `R^3` case would need a divergence-free field equal to `v_0` on a large ball (e.g. `curl(χ_R (v_0 × x)/2)`) and control of the pressure of the background; not done here.

## Breakpoint

The exact point where this line stops short of the full problem:

1. A1 (with A1b–A1d, A1c') **converts** "coercive Lyapunov functional for Navier–Stokes" into "coercive Euler invariant in the class (H1)–(H3)/(H3')/(H3'')/(H1')". This step is complete (PROVED-HERE) and uses nothing specific to the exact nonlinearity.
2. The **classification of Euler invariants** is then needed, and it is available here only in two classes: radial-multiplier quadratic forms of Sobolev type (A2: only the energy) and zeroth-order local densities `∫φ(u)` (A6: only energy and momentum). For first-order local densities `∫φ(u,∇u)` I recall Serre's 1984 classification (only energy, helicity, momentum, angular momentum — in 2D also `∫f(ω)`) [unverified, not reproved here]. For **nonlocal** functionals — in particular cubic forms `∫∫∫ K(x,y,z) u(x)u(y)u(z)`, which is where any candidate of `L^3` scaling would have to live — there is no classification and no method in this note. The precise unproved statement is:

```text
(E)  Every functional Q on V of class (H1) with (H3) (degree k ≥ 3) that is conserved by 3D Euler and satisfies Q(u) ≥ c ‖u‖_{L^3}^3 vanishes.
```

I state (E) as **CONJECTURE**; it is a rigidity statement about 3D Euler (no "critical Casimir"), and proving it would need a triad/Lie-algebraic analysis of all Euler-invariant polynomial functionals, of which A2 is the quadratic radial case.

3. Even if (E) were proved, the Lyapunov route would be closed but the problem not solved: ESS-type control is a statement about trajectories, and non-monotone a priori bounds are untouched by amplitude scaling (A4). Thus this line establishes that a regularity proof cannot be of the form "exhibit a coercive monotone state functional" within the classes above, and it produces the explicit fields (A3, A5) at which each named candidate increases; it does not produce any regularity statement, nor any blowup statement.

## Transfer check

Everything in A1 uses only bilinearity of `B`, the symmetry `u ↦ −u`, and local existence; A1c' and A2 additionally use the energy identity; A2 (detailed triad conservation) and A6 (pressure formula) use the exact Leray structure. Consequences and the table required by `00_problem.md` §4:

| item                                              | verdict                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| ------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 2.1 supercriticality (Tao 2007)                   | consistent and sharpened: the only monotone quantities in the covered classes are Euler invariants; the energy is the only one that is coercive for anything, and it is supercritical (1.6).                                                                                                                                                                                                                                                                                         |
| 2.2 averaged Navier–Stokes (Tao 2016)             | A1, A1b–A1d, A1f apply verbatim to `∂_t u = Δu + B̃(u,u)` (they never use the form of `B`); A1c' applies because `⟨B̃(u,u),u⟩ = 0`. So the averaged equation has no coercive Lyapunov functional in these classes either — as it must, since it blows up. A2 and A6 do not transfer (they use triad conservation with the Leray projector and the pressure formula), but they are negative results, so no conflict. T1 is passed trivially: this line contains no regularity argument. |
| 2.3 NRŠ/Tsai self-similar                         | not applicable (no blowup ansatz).                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| 2.4 ESS / Seregin / Tao 2019                      | consistent; A4 explains why `L^3` control is a trajectory statement and cannot be the sublevel set of a monotone functional in the covered classes.                                                                                                                                                                                                                                                                                                                                  |
| 2.5 LPS / Leray lower bounds                      | not applicable; no LPS-type bound is produced (T6 not triggered).                                                                                                                                                                                                                                                                                                                                                                                                                    |
| 2.6 CKN partial regularity                        | consistent: the localized energy `∫ χ (u·u) dx` is not monotone (A3.5); the local energy inequality is a distributional inequality with an indefinite flux, not a Lyapunov functional.                                                                                                                                                                                                                                                                                               |
| 2.7 BKM                                           | consistent: `‖ω‖_∞` increases for the explicit field of A3.3 at amplitude `> 2.02`; BKM is a criterion, not a bound.                                                                                                                                                                                                                                                                                                                                                                 |
| 2.8 2D regularity                                 | positive control passes (A1e, A2 planar certificate, A3.5): the 2D invariants `∫f(ω)` have zero inviscid flux; the 3D stretching term is exactly what breaks it.                                                                                                                                                                                                                                                                                                                     |
| 2.9 axisymmetric swirl-free                       | positive control passes by hand (A1e): `ω_θ/r` is transported.                                                                                                                                                                                                                                                                                                                                                                                                                       |
| 2.10 hyperdissipation (Lions, Tao 2009, BMR 2014) | A1 holds verbatim with `−(−Δ)^α`; the enstrophy is not monotone for any `α` (same Euler flux, A3.1). So Lions' regularity for `α ≥ 5/4` is not via a Lyapunov functional but via a bounded non-monotone quantity — consistent with A4.                                                                                                                                                                                                                                               |
| 2.11 non-uniqueness (BV, ABC, Hou et al.)         | not applicable: no uniqueness argument here. A1 uses smooth solutions (K1) only.                                                                                                                                                                                                                                                                                                                                                                                                     |
| 2.12 Cheskidov–Dai–Palasek                        | not applicable.                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| 2.13 Euler blowup (Elgindi, Chen–Hou)             | consistent: a coercive Euler invariant continuous on the relevant class would contradict Elgindi's `C^{1,α}` blowup for norms it controls; our results only concern smooth fields and say the known invariants (energy, helicity) are not coercive. T4 not applicable (no argument uniform in `ν`).                                                                                                                                                                                  |
| 2.14 PINN singularities                           | not applicable.                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| 2.15 standing constraints                         | consistent: no critical-space a priori bound is produced; the small-data theory is a bounded, non-monotone statement (A4).                                                                                                                                                                                                                                                                                                                                                           |
| T1                                                | passes (negative result); the lemma is explicitly `B`-independent.                                                                                                                                                                                                                                                                                                                                                                                                                   |
| T2                                                | passes: the three regular comparison systems satisfy the necessary condition (A1e), and the mechanism by which 3D candidates fail is the stretching term or the pressure term.                                                                                                                                                                                                                                                                                                       |
| T3                                                | not applicable (no blowup scenario).                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| T4                                                | not applicable.                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| T5                                                | not applicable; (H2) is a hypothesis on smooth solutions, and K1 uniqueness is used only to define "the" solution.                                                                                                                                                                                                                                                                                                                                                                   |
| T6                                                | not triggered: no estimate is closed; all conclusions are negative.                                                                                                                                                                                                                                                                                                                                                                                                                  |

What the transfer says about the value of the lemma: because A1 cannot distinguish Navier–Stokes from the averaged equation, it can never be part of a regularity proof; its value is as a filter that (i) eliminates the "monotone quantity" strategy in the classes above, (ii) redirects any Lyapunov-type attempt to the Euler-invariant question (E), and (iii) shows that the known regularity mechanisms in 3D (hyperdissipation, small data, ESS) are of the bounded-non-monotone type.

## What would be needed next

- Prove or refute (E) for **cubic forms**: classify all trilinear symmetric forms `C(u,u,u)` on `V` (Fourier kernels `K(n_1,n_2,n_3)` supported on `n_1+n_2+n_3 = 0`) with `C'(u)[B(u,u)] = 0`. This is a finite linear-algebra problem on each "quartet" of wavevectors and could be attacked with the exact triad machinery of the script; a negative answer would close the Lyapunov route for every polynomial functional of `L^3` scaling.
- Extend A2 from Sobolev symbols to all bounded translation-invariant `L` (matrix symbols, not necessarily radial): expected answer `L = aI + b·curl` (energy and helicity); needs triads connecting all wavevectors.
- Reproduce Serre's first-order classification [unverified] with exact certificates in the style of A6 (perturbation around linear background fields `v_0 + M x` is not periodic; the whole-space version would be needed).
- Give up monotonicity: the useful direction is *bounded* critical quantities along trajectories (ESS-type), i.e. quantitative versions of backward uniqueness and unique continuation, or "almost monotone" quantities of the form `Q(u(t)) + ∫_0^t (dissipation) ≤ Q(u(0)) + (controlled error)` — A1 says the error term cannot vanish for any state functional in the covered classes, so its precise form (which must use the exact nonlinearity, by T1) is the object to study.

## How to reproduce

```bash
cd research/navier-stokes
ruff check line_a_lyapunov.py && ruff format --check line_a_lyapunov.py
python3 line_a_lyapunov.py
```

Fixed seed `20260906`; no files are written; runtime about 52 s on 4 CPUs (sympy 2 s, exact certificates 0.2 s, FFT part 50 s). The printed Markdown tables are the ones pasted above.

## References

- T. Tao, "Why global regularity for Navier–Stokes is hard", blog post, 2007 \[checked via WebSearch in `00_problem.md`; content from memory\].
- T. Tao, "Finite time blowup for an averaged three-dimensional Navier–Stokes equation", J. Amer. Math. Soc. 29 (2016) \[checked via WebSearch in `00_problem.md`\].
- L. Escauriaza, G. Seregin, V. Šverák, Russian Math. Surveys 58 (2003); G. Seregin, Comm. Math. Phys. 312 (2012); T. Tao, arXiv:1908.04958 (2019) \[checked via WebSearch in `00_problem.md`\].
- J. Leray, Acta Math. 63 (1934); H. Fujita, T. Kato, Arch. Rational Mech. Anal. 16 (1964); T. Kato, J. Funct. Anal. 9 (1972) (local `H^s` theory); R. Temam, "Navier–Stokes Equations", North-Holland 1977/1984 (periodic case) [unverified; used for K1].
- D. Ebin, J. Marsden, Ann. of Math. 92 (1970); T. Kato 1972; J. T. Beale, T. Kato, A. Majda, Comm. Math. Phys. 94 (1984) \[BKM checked via WebSearch in `00_problem.md`; others unverified; used for K2\].
- D. Serre, "Les invariants du premier ordre de l'équation d'Euler en dimension trois", Physica D 13 (1984) 105–136 — I recall the statement that the only first-order local conserved densities of 3D Euler are linear combinations of energy, helicity, momentum and angular momentum, with `∫f(ω)` appearing only in 2D [unverified; not used in any proof here]. <!-- codespell:ignore les -->
- V. I. Arnold, B. A. Khesin, "Topological Methods in Hydrodynamics", Springer 1998 — helicity as the Casimir of the 3D Euler Lie–Poisson structure [unverified; background only].
- L. Brandolese, and S. Dobrokhotov, A. Shafarevich — algebraic `|x|^{-4}` spatial decay of Navier–Stokes solutions on `R^3` (non-preservation of the Schwartz class) \[unverified; used only to motivate `V = H^∞_σ`\].
- J. M. Danskin, "The Theory of Max-Min", Springer 1967; J. F. Bonnans, A. Shapiro, "Perturbation Analysis of Optimization Problems", Springer 2000 (Hadamard directional derivatives of Lipschitz functions) [unverified; used for K4].
- Beltrami/ABC flows: V. I. Arnold 1965; T. Dombre et al., J. Fluid Mech. 167 (1986) \[unverified; the identity `curl u = 2πu` is verified in the script and by hand in A3.4\].
