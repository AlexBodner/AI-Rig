# 00 — Problem statement, no-go checklist, transfer tests

Status legend used throughout this project (every claim carries exactly one): **PROVED-HERE** (complete proof written out, every step checkable), **KNOWN** (cited result, not reproved), **COMPUTED** (script plus printed output in the notes), **CONJECTURE**, **FAILED** (attempted, with the exact step where it breaks).

Citation tags: **[checked via WebSearch]** means a WebSearch result snippet confirmed the item named in the tag (usually title, authors, venue, year; sometimes the headline statement). **[unverified]** means the item is from memory and could not be checked in this environment (arXiv, Wikipedia and most sites are blocked). When a tag says "venue checked; statement from memory" the bibliographic data was confirmed but the precise mathematical statement was not.

Notation: `u = (u_1,u_2,u_3)` velocity, `p` pressure, `ν > 0` viscosity, `f` body force, `ω = curl u` vorticity, `∂_j = ∂/∂x_j`, `Δ = Σ_j ∂_j^2`, `div u = Σ_j ∂_j u_j`, `T^3 = R^3/Z^3` the unit torus, `T*` a putative first blowup time, `L^p`, `H^s`, `Ḃ^s_{p,q}` the usual Lebesgue, Sobolev and homogeneous Besov spaces, `L^q_t L^p_x` the mixed space-time norm.

## 1. The Clay Millennium problem as posed by Fefferman (2000)

Source: C. L. Fefferman, "Existence and smoothness of the Navier–Stokes equation", Clay Mathematics Institute official problem description, 2000 (reprinted in "The Millennium Prize Problems", CMI/AMS 2006) [checked via WebSearch: document exists; the wording below is reconstructed from memory and is a faithful paraphrase, not a verbatim quotation — unverified]. Status of this whole section: KNOWN (it is a problem statement, not a theorem).

### 1.1 Equations (Fefferman's (1)–(3), `n = 3`, `ν > 0` fixed)

```text
(1)  ∂_t u_i + Σ_{j=1}^{3} u_j ∂_j u_i = ν Δ u_i − ∂_i p + f_i(x,t)      (x ∈ R^3, t ≥ 0, i = 1,2,3)
(2)  div u = Σ_{j=1}^{3} ∂_j u_j = 0                                   (x ∈ R^3, t ≥ 0)
(3)  u(x,0) = u^0(x)                                                   (x ∈ R^3)
```

Here `u^0` is a given `C^∞` divergence-free vector field on `R^3`, `f` is a given `C^∞` vector field on `R^3 × [0,∞)`, and the unknowns are `u` and `p`.

### 1.2 Hypotheses on the data (Fefferman's (4)–(5)) and admissible solutions (Fefferman's (6)–(7)), whole-space case

```text
(4)  |∂^α u^0(x)| ≤ C_{α,K} (1 + |x|)^{-K}              on R^3, for every multi-index α and every K > 0
(5)  |∂^α_x ∂^m_t f(x,t)| ≤ C_{α,m,K} (1 + |x| + t)^{-K}  on R^3 × [0,∞), for every α, m, K
(6)  p, u ∈ C^∞(R^3 × [0,∞))
(7)  ∫_{R^3} |u(x,t)|^2 dx < C   for all t ≥ 0          ("bounded energy")
```

Fefferman calls a pair `(p,u)` satisfying (1),(2),(3),(6),(7) a "physically reasonable" solution. Two remarks that are part of the specification:

- Condition (7) is what rules out the trivial family `u(x,t) = a(t)`, `p(x,t) = −a'(t)·x` (smooth for any smooth `a`, solves (1)–(2), and shows that without (7) the problem is empty). It also normalizes the pressure only weakly: `p` is determined by `u` up to an additive function of `t` once one imposes decay, but (6)–(7) alone say nothing about the growth of `p` at spatial infinity. See the end of section 1.5.
- Uniqueness is not part of the question in (A)/(B); existence in the class (6)–(7) is.

### 1.3 Periodic case (Fefferman's (8)–(11))

```text
(8)   u^0(x + e_j) = u^0(x)                     for j = 1,2,3 (e_j the unit coordinate vectors)
(9)   f(x + e_j, t) = f(x, t)                   for j = 1,2,3
(10)  |∂^α_x ∂^m_t f(x,t)| ≤ C_{α,m,K} (1 + |t|)^{-K}   for every α, m, K
(11)  u(x + e_j, t) = u(x, t),  p(x + e_j, t) = p(x, t),   and  p, u ∈ C^∞(R^3 × [0,∞))
```

On the torus the energy is automatically bounded on any compact time interval for a smooth solution, so (7) is replaced by the periodicity requirement (11); `u^0` is `C^∞` and divergence-free as before.

### 1.4 The four statements (a proof of any one of them is asked for; (A) and (B) are not mutually exclusive)

- **(A) Existence and smoothness on `R^3`.** Take `ν > 0`, `n = 3`, `f ≡ 0`. For every `u^0` that is `C^∞`, divergence-free and satisfies (4), there exist `p, u` on `R^3 × [0,∞)` satisfying (1),(2),(3),(6),(7).
- **(B) Existence and smoothness on `T^3`.** Take `ν > 0`, `n = 3`, `f ≡ 0`. For every `u^0` that is `C^∞`, divergence-free and satisfies (8), there exist `p, u` on `R^3 × [0,∞)` satisfying (1),(2),(3),(11).
- **(C) Breakdown on `R^3`.** Take `ν > 0`, `n = 3`. There exist a `C^∞` divergence-free `u^0` satisfying (4) and a `C^∞` force `f` satisfying (5) such that no `p, u` on `R^3 × [0,∞)` satisfy (1),(2),(3),(6),(7).
- **(D) Breakdown on `T^3`.** Take `ν > 0`, `n = 3`. There exist a `C^∞` divergence-free `u^0` satisfying (8) and a `C^∞` force `f` satisfying (9),(10) such that no `p, u` on `R^3 × [0,∞)` satisfy (1),(2),(3),(11).

Reading notes (KNOWN, from the problem description):

- Forcing is allowed only in the breakdown statements (C),(D); the regularity statements (A),(B) are unforced. A breakdown proof may choose `f`, but `f` must remain smooth and decaying (in `t`, and in `x` on `R^3`) for all `t ∈ [0,∞)`, in particular through and beyond any blowup time.
- For any smooth divergence-free `v(x,t)` satisfying the decay hypotheses, the pair `(v, f := ∂_t v + v·∇v − νΔv + ∇q)` solves (1)–(3) for any smooth `q`. So with forcing, "some smooth solution exists" is trivially arrangeable; (C)/(D) assert the opposite — that for a specific `(u^0,f)` **no** smooth bounded-energy solution exists on the whole half-line. A blowup at `T* < ∞` of the solution one constructs is not enough by itself: one must also exclude every other smooth solution with the same data (see 1.5 for why this is automatic on `T^3`).
- The four statements are not known to be equivalent to each other. Tao, "Localisation and compactness properties of the Navier–Stokes global regularity problem", Analysis & PDE 6 (2013) [unverified] proves a number of implications among variants (Schwartz-class data, `H^1` data, periodic data, smooth forcing, finite-time versions) and discusses normalizations of the pressure that make the formulations well posed; I recall in particular that with mild extra normalization the periodic and whole-space regularity problems are shown to be equivalent to their `H^1`-data versions, but I do not recall the full implication diagram and do not rely on it below.

### 1.5 PROVED-HERE: uniqueness of smooth periodic solutions, and the reduction of (D) to a finite maximal time

**Proposition 1.** Let `u, v` be two solutions of (1),(2),(3),(11) on `R^3 × [0,T]` (same `ν`, same `u^0`, same periodic smooth `f`), with pressures `p_u, p_v` periodic and smooth. Then `u = v` on `R^3 × [0,T]`.

Proof. Set `w = u − v`, `q = p_u − p_v`. Subtracting the equations,

```text
∂_t w + (u·∇) w + (w·∇) v = ν Δ w − ∇ q,     div w = 0,     w(·,0) = 0.
```

All functions are smooth and periodic, so every integral below is over one period cell `Q = [0,1)^3` and every integration by parts is justified with no boundary terms. Take the `L^2(Q)` inner product with `w`:

- `∫_Q ∂_t w · w dx = (1/2) d/dt ∫_Q |w|^2 dx` (differentiation under the integral, `w` smooth on the compact set `Q × [0,T]`).
- `∫_Q ((u·∇) w)·w dx = (1/2) ∫_Q u·∇(|w|^2) dx = −(1/2) ∫_Q (div u) |w|^2 dx = 0` since `div u = 0`.
- `∫_Q ∇q · w dx = −∫_Q q div w dx = 0`.
- `∫_Q ν Δ w · w dx = −ν ∫_Q |∇ w|^2 dx ≤ 0`.
- `|∫_Q ((w·∇) v)·w dx| ≤ ‖∇v‖_{L^∞(Q×[0,T])} ∫_Q |w|^2 dx`, and `M := ‖∇v‖_{L^∞(Q×[0,T])} < ∞` because `∇v` is continuous on the compact set `Q × [0,T]` (compactness is where periodicity plus smoothness is used; nothing else about `v` enters).

Hence `E(t) := ∫_Q |w(x,t)|^2 dx` satisfies `E'(t) ≤ 2M E(t)`, `E(0) = 0`, `E ≥ 0`. Then `(e^{−2Mt} E(t))' ≤ 0`, so `e^{−2Mt}E(t) ≤ E(0) = 0`, so `E ≡ 0` on `[0,T]`. Since `w` is continuous, `w ≡ 0`. ∎

**Corollary 2 (reduction of (D)).** KNOWN input: for smooth periodic divergence-free `u^0` and smooth periodic `f` there is a unique smooth periodic solution on some interval `[0,T_0)`, `T_0 > 0`, and if the maximal such interval `[0,T*)` has `T* < ∞` then `‖u(t)‖_{H^1(Q)} → ∞` as `t → T*` (Leray 1934 for the whole space; the periodic local theory and the `H^1` continuation criterion are standard, e.g. Fujita–Kato 1964, Kato 1984, Temam's book [unverified]). Then: statement (D) holds if and only if there exist admissible `(u^0, f)` on `T^3` whose unique local smooth solution has finite maximal time `T* < ∞`. Proof. (⇐) Suppose the local solution `u` of some admissible `(u^0,f)` has `T* < ∞`. If a global smooth periodic solution `v` with the same data existed, Proposition 1 on `[0,T]` for every `T < T*` would give `v = u` on `[0,T*)`, so `‖u(t)‖_{H^1}` would stay bounded as `t → T*` (as `v` is smooth on the compact set `Q × [0,T*]`), contradicting the continuation criterion. So no global smooth solution exists for this `(u^0,f)`, which is (D). (⇒) Suppose (D) holds with witness `(u^0,f)`, and let `u` be its local solution with maximal time `T*`. If `T* = ∞` then `u` itself is a smooth periodic solution on `R^3 × [0,∞)`, contradicting (D). So `T* < ∞`. ∎

**What is NOT proved here.** The analogue of Proposition 1 on `R^3` inside Fefferman's class (6)–(7). The class only gives `∫|u|^2 dx < C`; it gives neither `∇u ∈ L^2` nor `u ∈ L^∞` nor decay of `p`, so the integrations by parts above are not justified as written. With additional decay (e.g. `u ∈ L^∞_t H^1_x` or Schwartz-class in `x`) the same proof goes through; Tao 2013 [unverified] discusses exactly these normalizations. For this project we will always work with solutions that are Schwartz in `x` on `R^3` (or periodic), which is strictly stronger than (6)–(7), so that uniqueness is available; any claimed (C) must at the end be checked against the weaker Clay class.

### 1.6 PROVED-HERE: scaling and the supercriticality of the energy

Let `(u,p)` solve (1)–(2) with `f ≡ 0`. For `λ > 0` set `u_λ(x,t) = λ u(λx, λ^2 t)`, `p_λ(x,t) = λ^2 p(λx, λ^2 t)`. Then `∂_t u_λ = λ^3 (∂_t u)(λx,λ^2 t)`, `(u_λ·∇)u_λ = λ^3 ((u·∇)u)(λx,λ^2 t)`, `Δ u_λ = λ^3 (Δu)(λx,λ^2 t)`, `∇ p_λ = λ^3 (∇p)(λx,λ^2 t)`, so `(u_λ, p_λ)` solves (1)–(2) with the same `ν`. Substituting `y = λx` (`dx = λ^{-3} dy`):

```text
∫ |u_λ(x,t)|^2 dx        = λ^{2−3}   ∫ |u(y, λ^2 t)|^2 dy    = λ^{-1} ∫ |u(·, λ^2 t)|^2       (energy: supercritical)
∫ |∇u_λ(x,t)|^2 dx       = λ^{4−3}   ∫ |∇u(y, λ^2 t)|^2 dy   = λ^{+1} ∫ |∇u(·, λ^2 t)|^2      (enstrophy: subcritical)
∫ |u_λ(x,t)|^3 dx        = λ^{3−3}   ∫ |u(y, λ^2 t)|^3 dy    = ∫ |u(·, λ^2 t)|^3              (L^3: critical)
∫_0^∞ ∫ |∇u_λ|^2 dx dt   = λ^{4−3−2} ∫_0^∞ ∫ |∇u|^2 dy ds    = λ^{-1} ∫∫ |∇u|^2              (dissipation: supercritical)
```

Reading: as `λ → ∞` (zooming into fine scales, where blowup would live) the energy and the total dissipation of the rescaled solution go to zero, so a uniform bound on them constrains fine scales less and less; a bound on a critical norm (`L^3`, `H^{1/2}`, `L^∞_t L^3_x`, `BMO^{-1}`) is scale-invariant; a bound on a subcritical norm (enstrophy, `L^∞`) blows up under zooming and therefore forces regularity (KNOWN: an a priori `H^1` bound implies global regularity, Leray 1934 [unverified]).

Hyperdissipation. Replace `νΔ` by `−ν(−Δ)^α`. Then `u_λ(x,t) = λ^{2α−1} u(λx, λ^{2α} t)` is the symmetry (each of the four terms then scales as `λ^{4α−1}`), and

```text
∫ |u_λ|^2 dx = λ^{2(2α−1) − 3} ∫ |u|^2 dy = λ^{4α−5} ∫ |u|^2
```

which is scale-invariant exactly at `α = 5/4`, subcritical for `α > 5/4`, supercritical for `α < 5/4`. This is the arithmetic behind the Lions exponent (no-go item 10).

## 2. No-go checklist

Every proposed argument produced in this project must be run against each item. "Kills" states the concrete test the item imposes. Unless stated otherwise the item is KNOWN.

### 2.1 Supercriticality barrier (Tao 2007)

- Source: T. Tao, "Why global regularity for Navier–Stokes is hard", blog post, March 2007 [checked via WebSearch: post exists; statement paraphrased from memory].
- Statement: the only globally controlled coercive quantities known for 3D Navier–Stokes are the energy `∫|u(t)|^2 dx`, the cumulative dissipation `∫_0^t ∫ |∇u|^2 dx ds` (both from the energy identity `(1/2) d/dt ∫|u|^2 = −ν∫|∇u|^2 + ∫ f·u`), and their localized version (the local energy inequality of CKN). By 1.6 all of these are supercritical. Every quantity that would imply regularity (enstrophy, any `H^s` norm with `s ≥ 1/2`, `L^p` with `p ≥ 3`, etc.) is critical or subcritical and is not known to be controlled; the enstrophy identity `d/dt ∫|ω|^2 dx = −2ν∫|∇ω|^2 dx + 2∫ ω·(ω·∇)u dx` has the vortex-stretching term with no sign. Tao's post classifies the three known strategies (exact/explicit solutions, perturbative/small-data or short-time arguments, abstract/compactness arguments producing weak limits) and explains why each is blocked by supercriticality: strong convergence is needed in a critical or subcritical topology while the available bounds give it only in supercritical ones.
- Hypotheses: none beyond the equations; this is a barrier, not a theorem.
- Kills: any argument whose only a priori inputs are the energy identity/inequality and the local energy inequality, combined with scaling and interpolation. Concretely: if the proposed proof would go through with the enstrophy production term `∫ ω·(ω·∇)u dx` replaced by any other term with the same scaling and no sign, it is wrong. A regularity proof must exhibit a new monotone or bounded quantity that is critical or subcritical, or a new structural mechanism.

### 2.2 Averaged Navier–Stokes blowup (Tao 2016)

- Source: T. Tao, "Finite time blowup for an averaged three-dimensional Navier–Stokes equation", J. Amer. Math. Soc. 29 (2016), 601–674, arXiv:1402.0290 [checked via WebSearch: title, venue, pages; abstract confirms the averaging description].
- Statement: write Navier–Stokes as `∂_t u = Δu + B(u,u)` with `B(u,v) = −P((u·∇)v)`, `P` the Leray projector. Tao constructs a bilinear operator `B̃` which is an average of `B` over a compact family of rotations and order-zero Fourier multipliers, so that `B̃` obeys the same cancellation `⟨B̃(u,u),u⟩ = 0` (hence the same energy identity `(1/2) d/dt ‖u‖_2^2 = −‖∇u‖_2^2`), the same scaling symmetry, and the same harmonic-analysis estimates (all the bilinear bounds used in the local well-posedness theory in Sobolev, Besov, Lorentz spaces, etc.), yet `∂_t u = Δu + B̃(u,u)` has smooth Schwartz-class data whose solution blows up in finite time. The blowup is engineered as a cascade in which a quadratic interaction transfers the energy to the next dyadic scale in a time that shrinks geometrically, via a "self-replicating machine" built from a small number of localized frequency modes at each scale.
- Hypotheses: 3D, unforced, whole space; the Laplacian is not modified; the averaging is over symmetries under which `B` is not invariant but the estimates are.
- Kills: (Transfer test T1.) Any proof of (A)/(B) that uses only (i) the energy identity, (ii) scaling, (iii) any inequality about `B` of the form `|⟨B(u,v),w⟩| ≤ C ‖u‖_X ‖v‖_Y ‖w‖_Z` with `X, Y, Z` translation- and rotation-invariant function spaces, (iv) parabolic smoothing of `e^{tΔ}`. The proof must at some step use a property of the exact Euler nonlinearity that is destroyed by averaging: vorticity transport `∂_t ω + (u·∇)ω = (ω·∇)u + νΔω`, Kelvin circulation / Lagrangian structure, the Biot–Savart geometry (direction of vorticity, Constantin–Fefferman 1993 [unverified]), the pointwise pressure formula `p = Σ R_i R_j (u_i u_j)`, or a maximum principle. The test to apply to any argument: name that property and the line where it is used; if no such line exists the argument is wrong.

### 2.3 No nontrivial backward self-similar solutions (Nečas–Růžička–Šverák 1996; Tsai 1998)

- Sources: J. Nečas, M. Růžička, V. Šverák, "On Leray's self-similar solutions of the Navier–Stokes equations", Acta Math. 176 (1996), 283–294 [checked via WebSearch]; T.-P. Tsai, "On Leray's self-similar solutions of the Navier–Stokes equations satisfying local energy estimates", Arch. Rational Mech. Anal. 143 (1998), 29–51 [checked via WebSearch: venue and headline statement].
- Statement: Leray (1934) asked whether there are blowup solutions of the form `u(x,t) = (T*−t)^{-1/2} U(x/(T*−t)^{1/2})` with a profile `U`. NRŠ: if `U ∈ L^3(R^3)` solves the profile equation `−ΔU − (1/2)U − (1/2)(y·∇)U + (U·∇)U + ∇P = 0`, `div U = 0`, then `U ≡ 0`. Tsai: the same conclusion if `U` satisfies local energy estimates (`U ∈ L^2_loc` type profile bounds appropriate to a suitable weak solution, and in particular if `U ∈ L^p` for some `p > 3`), with the exception of the trivial case; more precisely he shows that any self-similar solution in that class has zero profile or infinite local energy.
- Hypotheses: exact backward self-similarity with the Navier–Stokes scaling; `U` in `L^3` (NRŠ) or with local energy estimates (Tsai). Not covered: discretely self-similar or asymptotically self-similar blowup (not excluded in general), and forward self-similar solutions (which exist, Jia–Šverák 2014 [unverified]).
- Kills: any (C)/(D) attempt in which the singularity is exactly backward self-similar with a finite-energy or `L^3` profile. Any ansatz of the form above must be checked for `U ∉ L^3` and infinite local energy — i.e., outside the Clay class anyway.

### 2.4 Critical `L^∞_t L^3_x` regularity (Escauriaza–Seregin–Šverák 2003; Seregin 2012; Tao 2019)

- Sources: L. Escauriaza, G. Seregin, V. Šverák, "`L_{3,∞}`-solutions of the Navier–Stokes equations and backward uniqueness", Uspekhi Mat. Nauk 58 (2003), 3–44; Russian Math. Surveys 58 (2003), 211–250 [checked via WebSearch: venue, pages, headline statement]. G. Seregin, "A certain necessary condition of potential blow up for Navier–Stokes equations", Comm. Math. Phys. 312 (2012), 833–845, arXiv:1104.3615 [checked via WebSearch: venue, headline statement]. T. Tao, "Quantitative bounds for critically bounded solutions to the Navier–Stokes equations", arXiv:1908.04958 (2019), published in the Simons/"Nine Mathematical Challenges" proceedings [checked via WebSearch: arXiv id and headline statement; publication venue unverified].
- Statements: (ESS) if a Leray–Hopf (suitable) weak solution on `R^3 × (0,T)` satisfies `sup_{t<T} ‖u(t)‖_{L^3(R^3)} < ∞` then it is smooth on `(0,T]`; equivalently, at a first blowup time `T*` one has `limsup_{t→T*} ‖u(t)‖_{L^3} = ∞`. The proof uses a blowup/rescaling argument, CKN partial regularity, and a backward uniqueness theorem for the heat operator with lower-order terms (Carleman inequalities). (Seregin 2012) sharpens `limsup` to `lim_{t→T*} ‖u(t)‖_{L^3} = ∞`. (Tao 2019) makes ESS quantitative: higher norms are bounded by a triple-exponential function of the `L^3` bound, and at a blowup time `‖u(t_n)‖_{L^3} ≥ c (log log log (1/(T*−t_n)))^c` along some sequence `t_n → T*`. Improvements to `log log` in some critical spaces: S. Palasek, "Improved quantitative regularity for the Navier–Stokes equations in a scale of critical spaces", Arch. Rational Mech. Anal. (2021) [checked via WebSearch: title only; author and content from memory, unverified].
- Hypotheses: 3D, whole space (versions on the torus and with forcing exist [unverified]); solution in the Leray–Hopf class; the `L^3` norm is the critical Lebesgue norm. Smooth solutions with Schwartz data are in this class.
- Kills: (Transfer test T3.) Any blowup scenario in which `‖u(t)‖_{L^3}` stays bounded up to `T*`, or tends to infinity only along a subsequence, or grows slower than `(log log log 1/(T*−t))^c`. In particular a scenario in which vorticity blows up while the velocity stays in `L^3` is impossible. Any blowup in which `u` is bounded in a critical space contained in `L^3` (`H^{1/2} ⊂ L^3`, the Lorentz spaces `L^{3,q}` with `q ≤ 3`) is excluded; for the larger critical Besov spaces `Ḃ^{-1+3/p}_{p,q}` with `3 < p, q < ∞` (which contain `L^3`) the same conclusion is due to Gallagher–Koch–Planchon (2013–2016) [unverified], and the endpoint `Ḃ^{-1}_{∞,∞}` is not covered at all (Bourgain–Pavlović 2008 norm inflation shows ill-posedness there [unverified]).

### 2.5 Ladyzhenskaya–Prodi–Serrin conditions and Leray's lower bounds

- Sources: G. Prodi 1959, J. Serrin 1962, O. Ladyzhenskaya 1967; endpoint `p = 3` is 2.4 [unverified as to exact references; statements standard]. J. Leray, "Sur le mouvement d'un liquide visqueux emplissant l'espace", Acta Math. 63 (1934), 193–248 [unverified]. <!-- codespell:ignore mouvement -->
- Statement (LPS): if a Leray–Hopf solution satisfies `u ∈ L^q_t L^p_x` on `(0,T)` with `2/q + 3/p = 1`, `3 < p ≤ ∞`, then it is smooth on `(0,T]`; with `2/q + 3/p < 1` the condition is subcritical (e.g. `L^∞_t L^∞_x` suffices). Weak–strong uniqueness holds under the same condition. Statement (Leray): if `T*` is the first blowup time then `‖u(t)‖_{L^p} ≥ c_p (T*−t)^{−(p−3)/(2p)}` for `3 < p ≤ ∞` and `‖∇u(t)‖_{L^2} ≥ c (T*−t)^{−1/4}` (the exponents are forced by the scaling in 1.6; the constants depend on `ν` in a way I do not recall precisely and do not use), and the set of singular times of a Leray–Hopf solution has Hausdorff dimension at most `1/2` (Leray 1934; Scheffer 1976 [unverified]).
- Hypotheses: Leray–Hopf class; whole space or torus.
- Kills: any blowup scenario with `‖u(t)‖_{L^∞} = o((T*−t)^{−1/2})` (this is the Type I threshold: Type I means `‖u(t)‖_∞ ≤ C(T*−t)^{−1/2}`, and Leray's bound says blowup cannot be slower than Type I), with `‖∇u(t)‖_2 = o((T*−t)^{−1/4})`, or in which singular times accumulate on a set of dimension `> 1/2`. Also any regularity argument that proves an a priori bound in any LPS space is already a full solution of (A)/(B) — so if such a bound comes out easily, it is wrong.

### 2.6 Partial regularity (Caffarelli–Kohn–Nirenberg 1982)

- Source: L. Caffarelli, R. Kohn, L. Nirenberg, "Partial regularity of suitable weak solutions of the Navier–Stokes equations", Comm. Pure Appl. Math. 35 (1982), 771–831 [checked via WebSearch: venue, pages, headline statement]. Simplified proof: F. Lin, CPAM 1998 [unverified]; precursor: V. Scheffer 1976–1980 [unverified].
- Statement: for a suitable weak solution (Leray–Hopf plus the local energy inequality `∂_t (|u|^2/2) + div((|u|^2/2 + p)u) ≤ νΔ(|u|^2/2) − ν|∇u|^2` in the distributional sense, with `p ∈ L^{3/2}_{loc}` or `L^{5/3}_{loc}`), the singular set `S ⊂ R^3 × (0,∞)` (points near which `u` is not essentially bounded) has one-dimensional parabolic Hausdorff measure zero: `P^1(S) = 0`. Key lemma (`ε`-regularity): there is an absolute `ε_0 > 0` such that if `limsup_{r→0} r^{-1} ∫∫_{Q_r(x,t)} |∇u|^2 < ε_0` then `(x,t)` is regular; also a small-`L^3` version, `r^{-2}∫∫_{Q_r} |u|^3 + |p|^{3/2} < ε_0` implies regularity in `Q_{r/2}`.
- Hypotheses: suitability (local energy inequality); no smallness of data; whole space or bounded domains.
- Kills: (Transfer test T3.) Any blowup scenario whose singular set at time `T*` has positive length (a vortex filament that is singular along its whole length, a singular sheet, or a singularity persisting along a space-time curve of positive `P^1` measure). Any scenario in which the local dissipation `r^{-1}∫∫_{Q_r}|∇u|^2` stays below `ε_0` near the putative singular point. At a first singular time the singular set is compact, of `P^1` measure zero; a blowup must concentrate at least `ε_0`-much dissipation per unit scale at each singular point.

### 2.7 Beale–Kato–Majda 1984

- Source: J. T. Beale, T. Kato, A. Majda, "Remarks on the breakdown of smooth solutions for the 3-D Euler equations", Comm. Math. Phys. 94 (1984), 61–66 [checked via WebSearch: venue, pages, headline statement].
- Statement: for 3D Euler (and the same proof gives Navier–Stokes), a smooth solution on `[0,T*)` extends past `T*` if and only if `∫_0^{T*} ‖ω(t)‖_{L^∞} dt < ∞`. Hence at blowup `∫_0^{T*} ‖ω(t)‖_∞ dt = ∞`, in particular `limsup ‖ω(t)‖_∞ = ∞`. Refinements: Kozono–Taniuchi (BMO), Constantin–Fefferman 1993 (regularity if the direction `ω/|ω|` is Lipschitz where `|ω|` is large) [unverified].
- Hypotheses: smooth (`H^s`, `s > 5/2`) solutions, whole space or torus, viscosity irrelevant.
- Kills: any blowup scenario in which vorticity stays bounded, or is only `L^1`-in-time-integrable in sup norm (e.g. `‖ω(t)‖_∞ ~ (T*−t)^{−β}` with `β < 1` is impossible; by scaling the natural rate is `β = 1`). Together with 2.4: velocity must leave `L^3` and vorticity must leave `L^1_t L^∞_x` simultaneously.

### 2.8 Two-dimensional global regularity

- Sources: Leray 1933 (bounded domains), Ladyzhenskaya 1959, Lions–Prodi 1959 [unverified as to exact references; statement standard].
- Statement: in `R^2` and `T^2` the Navier–Stokes equations with smooth decaying data have a unique global smooth solution. Mechanism: the 2D vorticity is a scalar satisfying `∂_t ω + (u·∇)ω = νΔω` with no stretching term, so `‖ω(t)‖_{L^p}` is nonincreasing for every `p` (maximum principle); in particular the enstrophy `∫|ω|^2 = ∫|∇u|^2` is bounded for all time, and in 2D that is a subcritical bound (redo 1.6 with `d = 2`: `∫|u_λ|^2 = λ^{2−2}∫|u|^2` is scale-invariant, so the energy is critical, while `∫|∇u_λ|^2 = λ^{4−2}∫|∇u|^2 = λ^2 ∫|∇u|^2` is subcritical), so the `H^1` bound bootstraps to regularity.
- Kills: (Transfer test T2.) Any blowup mechanism that does not use the vortex-stretching term `(ω·∇)u` (or, equivalently, that survives the restriction to planar flows `u = (u_1(x_1,x_2), u_2(x_1,x_2), 0)`) is wrong. Note that 2D Euler allows double-exponential growth of `‖∇ω‖_∞` (Kiselev–Šverák 2014, on the disc [unverified]) — growth is not blowup.

### 2.9 Axisymmetric flows without swirl

- Sources: O. A. Ladyzhenskaya 1968; M. R. Ukhovskii and V. I. Yudovich 1968 [checked via WebSearch: the attribution and the statement are confirmed by several snippets; exact journal references unverified].
- Statement: for axisymmetric data `u = u_r(r,z) e_r + u_z(r,z) e_z` (zero swirl `u_θ ≡ 0`) with finite energy, the solution is globally smooth. Mechanism: `ω = ω_θ e_θ` and `ω_θ/r` is transported and diffused with a favorable sign: `(∂_t + u·∇)(ω_θ/r) = ν(Δ + (2/r)∂_r)(ω_θ/r)`, so `‖ω_θ/r‖_{L^p}` is nonincreasing, which gives a subcritical bound. With swirl the problem is open; known there: no Type I blowup (Chen–Strain–Tsai–Yau 2008/2009, Koch–Nadirashvili–Seregin–Šverák 2009 [unverified]). <!-- codespell:ignore yau -->
- Kills: (Transfer test T2.) Any blowup mechanism that survives the swirl-free axisymmetric reduction is wrong. Any candidate axisymmetric blowup must have nonzero swirl and cannot be Type I.

### 2.10 Hyperdissipative regularity at and beyond the Lions exponent

- Sources: J.-L. Lions, "Quelques méthodes de résolution des problèmes aux limites non linéaires", Dunod/Gauthier-Villars, 1969 \[checked via WebSearch: book and the `α ≥ 5/4` attribution\]; T. Tao, "Global regularity for a logarithmically supercritical hyperdissipative Navier–Stokes equation", Analysis & PDE 2 (2009), 361–366, arXiv:0906.3070 [checked via WebSearch]; D. Barbato, F. Morandin, M. Romito, "Global regularity for a slightly supercritical hyperdissipative Navier–Stokes system", Analysis & PDE 7 (2014), 2009–2027, arXiv:1407.6734 [checked via WebSearch: venue, pages, and that it proves Tao's conjectured optimal condition].
- Statement: with dissipation `−ν(−Δ)^α` in 3D, smooth data give global smooth solutions for `α ≥ 5/4` (Lions; elementary proof by Mattingly–Sinai [unverified]). Tao 2009: the same holds for the weaker dissipation `−D^2 u`, `D` a Fourier multiplier with symbol `|ξ|^{5/4}/g(|ξ|)` (so the dissipation symbol is `|ξ|^{5/2}/g(|ξ|)^2`), for any nondecreasing `g ≥ 1` with `∫^∞ ds/(s g(s)^4) = ∞`, e.g. `g(s) = log^{1/4}(2+s)`; Barbato–Morandin–Romito 2014 relax this to the optimal condition `∫^∞ ds/(s g(s)^2) = ∞`, e.g. `g(s) = log^{1/2}(2+s)` [both conditions from memory, unverified]. Mechanism: at `α = 5/4` the energy is critical (1.6), and an enstrophy-type estimate closes with a logarithmic loss that is absorbed by the extra dissipation.
- Hypotheses: 3D, whole space or torus, dissipation strictly stronger than the critical one by a divergent-integral amount.
- Kills: (Transfer test T2.) Any blowup mechanism insensitive to the strength of the dissipation at high frequency — in particular any mechanism that does not degrade when `νΔ` is replaced by `−ν(−Δ)^{5/4}` — is wrong. Conversely, any regularity proof that works for `α = 1` must be checked for what it uses beyond what the `α ≥ 5/4` proofs use: if its estimates would also close for `α = 1` in the averaged equation of 2.2, it is wrong.

### 2.11 Non-uniqueness of weak solutions (Buckmaster–Vicol 2019) and of Leray–Hopf solutions (Albritton–Brué–Colombo 2022; Hou et al. 2025)

- Sources: T. Buckmaster, V. Vicol, "Nonuniqueness of weak solutions to the Navier–Stokes equation", Ann. of Math. 189 (2019), 101–144 [checked via WebSearch: venue, pages, headline statement]. D. Albritton, E. Brué, M. Colombo, "Non-uniqueness of Leray solutions of the forced Navier–Stokes equations", Ann. of Math. 196 (2022), 415–455 [checked via WebSearch: venue, pages, headline statement]. T. Y. Hou and two coauthors, "Nonuniqueness of Leray–Hopf solutions to the unforced incompressible 3D Navier–Stokes equation", arXiv:2509.25116 (2025), computer-assisted [checked via WebSearch: title, arXiv id, headline statement; I do not recall the coauthors and cannot verify the proof's status — unverified beyond the snippet].
- Statements: (BV) there are infinitely many weak solutions in `C^0_t L^2_x` (not Leray–Hopf: they need not satisfy the energy inequality, and the energy profile can be prescribed), obtained by convex integration with intermittent building blocks. (ABC) there is a smooth-away-from-`t=0` force `f` such that two distinct Leray–Hopf solutions with `u^0 = 0` exist; the mechanism is a linearly unstable self-similar vortex-ring background in similarity variables, with the instability producing a second solution; the force is self-similar and singular at `t = 0`, so it does not satisfy Fefferman's (5), and neither solution is a candidate for (C). (Hou et al. 2025) non-uniqueness of Leray–Hopf solutions with `f = 0`, again via an unstable self-similar profile, established with rigorous numerics (linearized operator decomposed as coercive plus compact, with computer-assisted bounds). Related: non-uniqueness of smooth solutions from critical (non-smooth, `L^3`-type) data, Inventiones 2025, arXiv:2503.14699, I recall the authors as Coiculescu–Palasek [checked via WebSearch: title and venue only; authors unverified].
- Hypotheses: BV works in a class much weaker than Leray–Hopf; ABC needs a singular force; Hou et al. is unforced but works "in the self-similar setting", which I read as data that are `(−1)`-homogeneous and hence singular at the origin, not smooth (inferred from the abstract; not verified beyond it) — none of these produce non-uniqueness from smooth data, which is consistent with weak–strong uniqueness while a smooth solution exists.
- Kills: (Transfer test T5.) Any regularity or uniqueness argument that applies to every Leray–Hopf weak solution, using only the (local) energy inequality, is wrong: it would prove uniqueness in the Leray–Hopf class, contradicting ABC (forced) or Hou et al. (unforced). Any argument that applies to all `C^0_t L^2_x` weak solutions is wrong by BV. So a regularity proof must use smoothness of the data or of the solution on `[0,T*)` in an essential way (as ESS does via the blowup-and-backward-uniqueness route).

### 2.12 Instantaneous Type I blowup and non-uniqueness of smooth solutions (Cheskidov–Dai–Palasek 2025)

- Source: A. Cheskidov, M. Dai, S. Palasek, "Instantaneous Type I blow-up and non-uniqueness of smooth solutions of the Navier–Stokes equations", arXiv:2511.09556 (Nov 2025, revised Jan 2026) [checked via WebSearch: authors, title, arXiv id, and the abstract statements paraphrased below; nothing beyond the abstract is verified].
- Statement (from the abstract): for any smooth divergence-free data on `T^d`, `d ≥ 2`, there is a solution which is smooth in space and time on `T^d × ([0,T] \ {T*})`, exhibits Type I blowup of `‖u‖_{L^∞}` at `T*`, and bifurcates from the classical solution by an instantaneous injection of energy from infinite wavenumber (a complete inverse energy cascade), producing an infinite family of spatially smooth solutions with the same data.
- Hypotheses/scope as far as the abstract states: periodic, any dimension including `d = 2`; the solutions are not smooth at `T*` and are not in the Leray–Hopf class (energy enters at `T*`, so the energy inequality fails there).
- Kills / does not kill: this does NOT resolve (D): the constructed solutions fail (11) at `T*`, and the existence of a bad solution says nothing about whether the classical solution is global (Corollary 2 needs the classical solution to blow up). It does kill: any reading of (B)/(D) in which "solution" means merely "smooth in `x` for each `t` and smooth for `t ≠ T*`"; and any uniqueness argument that does not use smoothness across every `t`, or the energy inequality, in an essential way. It also shows (since it works for `d = 2`) that any "blowup" phenomenon that survives in 2D with this weak solution concept is an artifact of the solution class, not a Navier–Stokes singularity — a further instance of T2.

### 2.13 Euler blowup: Elgindi 2021 and Chen–Hou 2022–2023

- Sources: T. M. Elgindi, "Finite-time singularity formation for `C^{1,α}` solutions to the incompressible Euler equations on `R^3`", Ann. of Math. 194 (2021), 647–727 [checked via WebSearch: venue, pages]. J. Chen, T. Y. Hou, "Stable nearly self-similar blowup of the 2D Boussinesq and 3D Euler equations with smooth data I: Analysis" (arXiv:2210.07191) and "II: Rigorous numerics" (arXiv:2305.05660) [checked via WebSearch: titles, arXiv ids, headline statement]. Numerical precursor: G. Luo, T. Y. Hou, PNAS 2014 [unverified].
- Statements: (Elgindi) there are axisymmetric, swirl-free, `C^{1,α}` (small `α`), finite-energy velocity fields on `R^3` whose Euler solutions blow up in finite time, via a self-similar profile for a leading-order model of the vorticity equation; the blowup is unstable in the class of smooth data but is stabilized by symmetry in `C^{1,α}`. (Chen–Hou) with smooth data, in a cylinder with boundary, the axisymmetric Euler equations (with swirl) blow up in finite time, via a nearly self-similar profile whose nonlinear stability is established with computer assistance (weighted `L^∞` and `C^{1/2}` norms, rigorous error control of the linearized operator). Smooth Euler blowup without boundary on `R^3` or `T^3` remains open (recent progress on `C^{1,α}`-type asymptotically self-similar blowup with better regularity appears in snippets, unverified).
- Hypotheses: `ν = 0`; low regularity (Elgindi) or boundary (Chen–Hou). For Navier–Stokes with `ν > 0`: Chen–Hou's mechanism relies on the boundary; Elgindi's `C^{1,α}` data are not in the Clay class; whether viscosity destroys either scenario is open (Hou has numerical evidence of "potentially singular behavior" for 3D Navier–Stokes with the Luo–Hou-type scenario: T. Y. Hou, Found. Comput. Math. 2023 [checked via WebSearch: title and venue; content from memory, unverified]).
- Kills: (Transfer test T4.) Any regularity argument for Navier–Stokes that is uniform in `ν → 0` in a `C^{1,α}` class, or that would apply verbatim to axisymmetric Euler with swirl in a cylinder, is wrong. Conversely, no proposed Navier–Stokes blowup mechanism can be justified by Euler blowup alone: the viscous term `νΔω` is of the same scaling order as the stretching term at the natural blowup rate `|ω| ~ (T*−t)^{-1}`, `length ~ (T*−t)^{1/2}` (Navier–Stokes scaling), whereas the Euler self-similar blowups of Elgindi and Chen–Hou have a different spatial scaling (`length ~ (T*−t)^{1/λ}` with `λ ≠ 2`), so under Navier–Stokes scaling the viscous term dominates at small scales in those profiles (CONJECTURE as a general principle; the specific comparison is a scaling computation to be done per scenario).

### 2.14 Unstable self-similar singularities found with physics-informed neural networks (2025)

- Source: Y. Wang, M. Bennani, J. Martens, S. Racanière, S. Blackwell, A. Matthews, S. Nikolov, G. Cao-Labora, D. S. Park, M. Arjovsky, D. Worrall, C. Qin, F. Alet, B. Kozlovskii, N. Tomašev, A. Davies, P. Kohli, T. Buckmaster, B. Georgiev, J. Gómez-Serrano, R. Jiang, C.-Y. Lai, "Discovery of unstable singularities", arXiv:2509.14185 (2025) [checked via WebSearch: author list, title, arXiv id, headline statement]. Follow-up: "Resolving sharp gradients of unstable singularities to machine precision via neural networks", arXiv:2511.22819 [checked via WebSearch: title only]. Precursor: Wang, Lai, Gómez-Serrano, Buckmaster, "Asymptotic self-similar blow-up profile for three-dimensional axisymmetric Euler equations using neural networks", Phys. Rev. Lett. 130 (2023) [checked via WebSearch: title and venue].
- Statement: using PINNs trained to near machine precision (with a second-order optimizer, from memory, unverified), the authors find new families of unstable self-similar blowup profiles for the incompressible porous media (IPM) equation, the 2D Boussinesq equation (equivalently 3D axisymmetric Euler in the Luo–Hou scenario, with boundary), and the Córdoba–Córdoba–Fontelos-type 1D model [the third equation is from memory, unverified], indexed by the number of unstable directions, with the self-similar exponent `λ` depending near-linearly on the instability index. Unstable singularities require infinitely precisely tuned data and are invisible to standard time-stepping DNS. Some of the profiles are intended as inputs to computer-assisted proofs in the style of Chen–Hou.
- Hypotheses: numerical evidence, not theorems, except where a subsequent computer-assisted proof exists; none of the equations is 3D Navier–Stokes; all three are inviscid or have boundary.
- Kills: the heuristic "if 3D Navier–Stokes blew up, DNS would have seen it" — an unstable blowup would not be seen. It also sets the methodological bar for a (C)/(D) attempt: a candidate profile must come with a self-similar (or discretely self-similar) formulation, an unstable-manifold count, and a plan for a rigorous stability estimate of the linearized operator; a profile without those is not a candidate. For Navier–Stokes specifically, 2.3 says an exactly self-similar profile in `L^3` cannot exist, so any PINN-style search must target discretely or asymptotically self-similar profiles, or profiles outside `L^3` (hence outside the Clay class), or forced problems in the sense of (C)/(D).

### 2.15 Additional standing constraints (brief)

- Local well-posedness and the critical spaces: Fujita–Kato 1964 (`H^{1/2}`), Kato 1984 (`L^3`), Koch–Tataru 2001 (`BMO^{-1}`), all KNOWN [unverified as to exact references]. Small data in any of these spaces give global smooth solutions; hence any blowup needs data (or the solution at some time) large in `BMO^{-1}`, and any argument that gives an a priori bound in one of these spaces on `[0,T*)` proves regularity.
- Leray–Hopf weak solutions exist globally for `L^2` data (Leray 1934 on `R^3`, Hopf 1951 on domains [unverified]); the open problem is exactly whether they stay smooth. Weak–strong uniqueness: a Leray–Hopf solution coincides with the smooth solution as long as the latter exists (Prodi–Serrin [unverified]).
- Dyadic/shell models: Katz–Pavlović 2005 blowup for a dyadic model [unverified]; Cheskidov 2008 blowup of Leray–Hopf-type dyadic solutions [unverified]. These share the energy identity and scaling with Navier–Stokes and are the precursors of 2.2; the same T1 test applies.
- Ill-posedness at the endpoint: Bourgain–Pavlović 2008 norm inflation in `Ḃ^{-1}_{∞,∞}` [unverified]. Any regularity proof yielding continuity of the data-to-solution map in `Ḃ^{-1}_{∞,∞}` is wrong.
- Hou's DNS of potentially singular Navier–Stokes behavior (2.13) and the Kerr, Pelz, Hou–Li, Boratav–Pelz scenarios: all earlier candidate scenarios were found to be either non-singular or unstable under refinement [unverified, from memory]. A blowup claim should say which known scenario it is closest to and why the earlier depletion mechanism does not apply.

## 3. Transfer tests

Every argument produced in this project must pass all of the following before it is written up as anything other than FAILED. The tests are mechanical: apply the argument, line by line, to the comparison system, and either locate the line that fails or conclude that the argument is wrong.

- **T1 (averaging test).** If a regularity argument for (A)/(B) applies verbatim to Tao's averaged equation `∂_t u = Δu + B̃(u,u)` of 2.2 — i.e., if it uses only the energy identity, scaling, parabolic smoothing, and bilinear estimates on `B` valid for every order-zero average of `B` — it is wrong. The write-up must point to the line where a property specific to the exact nonlinearity is used (vorticity transport, Kelvin/Lagrangian structure, Biot–Savart geometry, pressure formula, maximum principle) and explain why the averaged operator lacks it. The same test applies to arguments that use only the (local) energy inequality plus abstract compactness (2.1).
- **T2 (globally regular comparison systems).** If a blowup mechanism for (C)/(D) applies verbatim to any system known to be globally regular, it is wrong. The comparison list: 2D Navier–Stokes (2.8); axisymmetric swirl-free 3D Navier–Stokes (2.9); hyperdissipative 3D Navier–Stokes with `α ≥ 5/4`, including the logarithmically supercritical versions (2.10); the viscous Burgers system `∂_t u + (u·∇)u = νΔu` in any dimension (no pressure, no divergence constraint; globally regular by the componentwise maximum principle `‖u_i(t)‖_∞ ≤ ‖u_i(0)‖_∞`, which is a subcritical bound — KNOWN, standard [unverified as to a reference]); 3D Navier–Stokes with small data in `BMO^{-1}` (2.15); and any 1D or 2D model in which the mechanism is claimed to be visible. The Burgers item forces any mechanism to use pressure/incompressibility in an essential way; the 2D and swirl-free items force it to use vortex stretching; the hyperdissipative item forces it to be sensitive to the strength of dissipation at high frequency.
- **T3 (consistency with the blowup profile constraints).** Any (C)/(D) blowup must satisfy simultaneously: `lim_{t→T*} ‖u(t)‖_{L^3} = ∞` with at least `(log log log 1/(T*−t))^c` growth along a sequence (2.4); `∫_0^{T*} ‖ω(t)‖_∞ dt = ∞` (2.7); Leray's lower bounds `‖u(t)‖_∞ ≥ c(T*−t)^{−1/2}`, `‖∇u(t)‖_2 ≥ c(T*−t)^{−1/4}` (2.5); a singular set at `T*` of `P^1` measure zero with at least `ε_0` dissipation per unit scale at each singular point (2.6); no exact backward self-similarity with an `L^3` or finite-local-energy profile (2.3); and, if axisymmetric, nonzero swirl and no Type I rate (2.9). A scenario violating any one of these is FAILED, and the violated constraint is the recorded breakpoint.
- **T4 (inviscid-limit test).** A regularity argument for Navier–Stokes that is uniform in `ν → 0` on a class containing Elgindi's `C^{1,α}` data, or that applies verbatim to axisymmetric Euler with swirl in a cylinder, contradicts 2.13 and is wrong. A blowup argument must show where the viscous term is dominated; a scaling check of `νΔω` against `(ω·∇)u` in the proposed profile is mandatory (2.13).
- **T5 (weak-solution test).** A regularity or uniqueness argument that applies to every Leray–Hopf solution using only the energy inequality contradicts 2.11 (ABC forced, Hou et al. unforced) and is wrong; one that applies to every `C^0_t L^2_x` weak solution contradicts Buckmaster–Vicol. A uniqueness argument that does not use smoothness across every time, or the energy inequality, contradicts 2.12. The write-up must state where smoothness of the solution on `[0,T*)` (or of the data) is used.
- **T6 (too-easy test, from the honesty rules).** If a step closes an estimate that would give an a priori bound in any critical or subcritical space (2.15), or proves an LPS-type bound (2.5), or bounds the enstrophy production term by the dissipation with a constant independent of the solution, the default conclusion is that the step is wrong; the step must be re-derived with every constant and every function space made explicit, and then run through T1.

## 4. What this file commits the project to

- Statements attacked: (B) and (D) first (torus; uniqueness in the Clay class is PROVED-HERE in 1.5 so that (D) is exactly "the unique local solution has finite maximal time"), with (A)/(C) as the whole-space translation, always working in Schwartz class in `x` and recording separately what is needed to descend to the weaker Clay class (6)–(7).
- Every attack line will be written with the status legend, and its final section will be a table with one row per no-go item 2.1–2.15 and per transfer test T1–T6, each row saying "passes because ...", "fails at step ...", or "not applicable because ...".
