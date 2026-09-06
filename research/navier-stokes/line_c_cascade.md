# Line C — Tao's cascade program: engineered blowup, and what blocks porting it to true Navier-Stokes

Status legend (one tag per claim): **PROVED-HERE** (complete proof written out, every step checkable), **KNOWN** (cited, not reproved), **COMPUTED** (script plus printed output below), **CONJECTURE**, **FAILED** (attempted, with the exact breaking step). Citation tags: **[checked via WebSearch]** / **[unverified]** as in `00_problem.md`. All no-go items `2.x` and transfer tests `T1`-`T6` are cited by number from `00_problem.md` and are not re-derived here.

## Goal

Tao's 2016 program (no-go item 2.2) has two halves. The first is a *mechanism*: a cascade in which a quadratic, energy-conserving, scale-invariant nonlinearity moves an order-one fraction of the energy from scale `2^n` to scale `2^(n+1)` in a time that shrinks geometrically fast enough to beat viscous damping, so that infinitely many transfers happen before a finite time `T*`. The second is a *carrier*: a bilinear operator `B~` that supports the mechanism while retaining every property of the true Euler nonlinearity `B(u,v) = -P((u.grad)v)` that the known regularity theory uses. This line does the first half concretely on ODE (dyadic / shell) models, measures where it stops, and then states precisely which structural property of `B` blocks the second half for the true equation. The intended output is not "a blowup"; it is a quantitative location of the boundary between the cascade regimes, and a list of the exact properties an engineered nonlinearity must break.

Line A of this project asked whether a coercive Lyapunov functional could rule out blowup and found that none of the natural candidates works; this line asks the complementary question, whether the cascade that would produce blowup can be built, and finds that in the one-mode-per-scale class it provably cannot beat the true dissipation.

## Setting and notation

The dyadic (shell) model studied here is, for `n = 0, 1, 2, ...` with the conventions `a_{-1} = 0` and `a_N = 0` at a truncation to `N` shells,

```text
da_n/dt = -nu * 2^(2 s n) * a_n  +  lam^n * a_{n-1}^2  -  lam^(n+1) * a_n * a_{n+1},     lam = 2^(5/2).
```

`a_n(t)` is the amplitude of the `n`-th dyadic shell in an `L^2`-normalized basis, so that `E = sum_n a_n^2` is the energy analogue, `H^1 = sum_n 4^n a_n^2` is the enstrophy analogue, `k_n = 2^n` is the physical wavenumber, and `s = 1` (dissipation `4^n = k_n^2`) is the true Navier-Stokes Laplacian. The nonlinear coefficient base `lam = 2^(5/2)` is not a free choice: C2 below shows it is the unique base for which the model has the exact Navier-Stokes scaling symmetry of `00_problem.md` section 1.6 together with the dissipation `4^n`.

Two normalizations appear in the literature and they must be kept apart. Writing the model in "flux units", with `K_n := lam^n` as the coefficient of the nonlinearity and the dissipation as `nu * K_n^(2 alpha)`, the exponent `alpha` is the one used by Katz-Pavlovic and Cheskidov; writing it in "physical units", with dissipation `nu * 2^(2 s n) = nu * k_n^(2 s)`, the exponent `s` is the one used here and in the Lions/hyperdissipation literature (`00_problem.md` item 2.10). The dictionary is `2 s n = 2 alpha * (5/2) n`, i.e. `alpha = 2 s / 5`:

| quantity | this file (`s`, dissipation `2^(2 s n)`) | flux units (`alpha`, dissipation `lam^(2 alpha n)`) | ratio to critical |
| --- | --- | --- | --- |
| Cheskidov's proved blowup boundary | `s = 5/6 ~ 0.8333` | `alpha = 1/3` | `2/3` |
| true 3D Navier-Stokes dissipation | `s = 1` | `alpha = 2/5` | `4/5` |
| energy criticality = Lions exponent | `s = 5/4` | `alpha = 1/2` | `1` |

The identification of Cheskidov's `alpha` with `2 s / 5` is **inferred**, not read off the paper: it is forced by two independent consistency checks, namely that his global-regularity threshold `alpha >= 1/2` coincides with the energy-criticality exponent computed in C2 below, and that his remark "`alpha = 1/3` enjoys the same estimates as 4D Navier-Stokes" matches the ratio `alpha / alpha_crit = 2/3`, which is exactly the ratio `1 / ((d+2)/4)` at `d = 4`. Note the last column: the dyadic model with the true `4^n` dissipation sits at `4/5` of critical, exactly where true 3D Navier-Stokes sits, which is the sense in which the model is scale-faithful.

Two profiles will be compared throughout. The **scale-invariant profile** `a_n = c 2^(-n/2)` is the fixed point of the model's scaling map (C2). The **constant-flux (K41) profile** `a_n = c 2^(-5n/6)` is the one that makes the shell-to-shell energy flux `lam^(n+1) a_n^2 a_{n+1}` independent of `n` (C5). The second is strictly steeper than the first, and C10 shows that this one inequality is the whole reason the model does not blow up at `s = 1`.

## Results

### C1. Exact energy conservation (inviscid) and the energy identity — **PROVED-HERE**, confirmed **COMPUTED**

Proof. With `a_{-1} = 0` and `a_N = 0`, and all sums over `0 <= n <= N-1`,

```text
(1/2) d/dt sum_n a_n^2 = sum_n a_n (lam^n a_{n-1}^2 - lam^(n+1) a_n a_{n+1})
                       = sum_n lam^n a_{n-1}^2 a_n  -  sum_n lam^(n+1) a_n^2 a_{n+1}.
```

Re-index the first sum by `m = n-1`: it becomes `sum_{m=-1}^{N-2} lam^(m+1) a_m^2 a_{m+1}`, whose `m = -1` term vanishes because `a_{-1} = 0`, and whose remaining terms are exactly the second sum (its `m = N-1` term vanishes because `a_N = 0`). The two sums cancel term by term. With the dissipation restored, `(1/2) d/dt sum_n a_n^2 = -nu sum_n 2^(2 s n) a_n^2 <= 0`. QED. Script part 1a/1b verifies both identities symbolically for `N = 8` with `lam`, `nu`, `s` free.

This is the property that makes the model a legitimate carrier of Tao's mechanism: it has the exact energy identity of `00_problem.md` item 2.2, so no argument that uses only the energy identity can distinguish it from Navier-Stokes.

### C2. `lam = 2^(5/2)` is the unique base with the Navier-Stokes scaling; the scale-invariant profile — **PROVED-HERE**, confirmed **COMPUTED**

Proposition. Fix the dissipation `2^(2 s n)`. Define `b_n(t) := mu a_{n+1}(sig t)`. Then `b` solves the model with the same `nu` and `s` whenever `a` does, if and only if `sig * 2^(2 s) = 1` and `sig * lam / mu = 1`; i.e. `sig = 2^(-2 s)` and `mu = lam * 2^(-2 s)`.

Proof. Differentiate: `db_n/dt = mu sig (da_{n+1}/dt)(sig t)`. Substituting the equation at shell `n+1` and rewriting `a_{n+1} = b_n/mu`, `a_n = b_{n-1}/mu`, `a_{n+2} = b_{n+1}/mu` gives `db_n/dt = -nu sig 2^(2s) 2^(2 s n) b_n + (sig lam / mu) lam^n b_{n-1}^2 - (sig lam / mu) lam^(n+1) b_n b_{n+1}`, and matching coefficients with the model at shell `n` gives the two stated conditions. QED (script part 1c: residual is `0` after substitution).

Corollary (uniqueness of `lam`). At `s = 1` the first condition forces `sig = 1/4`, which is the parabolic time scaling `t -> t/4` of `u_lam0(x,t) = lam0 u(lam0 x, lam0^2 t)` at `lam0 = 2`. The same map multiplies each shell amplitude by `lam0^(1 - 3/2) = 2^(-1/2)` (because `||u_lam0||_{L^2} = lam0^(-1/2) ||u||_{L^2}` in `d = 3`, `00_problem.md` 1.6), so its inverse — which is the map `b` above, moving content one shell up — must have `mu = 2^(1/2)`. Then `mu = lam/4 = 2^(1/2)` forces `lam = 2^(5/2)`. QED.

Corollary (scale-invariant profile). `a_n = c 2^(-n/2)` is the unique (up to `c`) profile fixed by the map `a -> b` at `s = 1`, since `2^(1/2) c 2^(-(n+1)/2) = c 2^(-n/2)`. It is **not** a stationary solution: substituting it into the nonlinearity gives residual `-2 c^2 2^(3n/2)` (script part 1e, exact). Under this map `E -> 2E` and `H^1 -> H^1/2`, reproducing the supercriticality of the energy and the subcriticality of the enstrophy of `00_problem.md` 1.6.

### C3. The `(nu, amplitude)` collapse — **PROVED-HERE**, confirmed **COMPUTED**

Proposition. If `a` solves the model with viscosity `nu` then `b_n(t) := mu a_n(mu t)` solves it with viscosity `mu nu`, for every `mu > 0`.

Proof. `db_n/dt = mu^2 (da_n/dt)(mu t) = mu^2 [-nu 2^(2 s n) a_n + lam^n a_{n-1}^2 - lam^(n+1) a_n a_{n+1}](mu t)`. The nonlinearity is exactly quadratic, so it equals `-(mu nu) 2^(2 s n) b_n + lam^n b_{n-1}^2 - lam^(n+1) b_n b_{n+1}`. QED.

Consequences used below: only the ratio `amplitude / nu` matters, so the mandatory `(nu, amplitude)` table is one-dimensional; `max H^1` scales by `mu^2` and the stall shell `n_d` (C8) is unchanged; and "blowup for large data at fixed `nu`" and "blowup for fixed data at small `nu`" are the same statement. Verified numerically in part 3 at `mu = 3` over three decades of `nu`.

### C4. The exact self-similar blowup of the inviscid model has infinite energy — **PROVED-HERE**, confirmed **COMPUTED**

Proposition. On the bi-infinite lattice `n in Z`, `a_n(t) = lam^(-n) / ((lam^2 - 1)(T - t))` solves the inviscid model exactly, for every `T`.

Proof. Write `c = 1/(lam^2-1)`. Left side: `c lam^(-n) (T-t)^(-2)`. Right side: `lam^n c^2 lam^(-2(n-1)) (T-t)^(-2) - lam^(n+1) c^2 lam^(-n) lam^(-(n+1)) (T-t)^(-2) = c^2 (T-t)^(-2) (lam^(2-n) - lam^(-n)) = c^2 (lam^2 - 1) lam^(-n) (T-t)^(-2)`. The two agree iff `c = c^2 (lam^2-1)`, which is the definition of `c`. QED (script part 1g: residual `0`).

Its energy `sum_{n in Z} a_n^2` diverges geometrically as `n -> -infinity`, and it violates the boundary condition `a_{-1} = 0` of the half-lattice model. This is the model's analogue of no-go item 2.3 (Necas-Ruzicka-Sverak, Tsai): the exactly self-similar blowup exists but is not in the finite-energy class, so a finite-energy blowup must be a *front*, not an exact self-similar profile. It is also the profile of the far precursor tail seen numerically in C6 (measured slope near `-5/2 = log2(lam^(-1))`).

### C5. The constant-flux (K41) profile is an exact stationary solution and is steeper than scale-invariant — **PROVED-HERE**, confirmed **COMPUTED**

Define the shell-to-shell flux `Pi_n := lam^(n+1) a_n^2 a_{n+1}`, so that the shell energy budget reads `(1/2) d/dt a_n^2 = Pi_{n-1} - Pi_n - nu 2^(2 s n) a_n^2` (immediate from the equation; it is the term-by-term version of C1).

Proposition. `a_n = c 2^(-5n/6)` makes `Pi_n = 2^(5/3) c^3` independent of `n`, and on the bi-infinite lattice it is an exact stationary solution of the inviscid model.

Proof. `Pi_n = 2^(5(n+1)/2) c^3 2^(-5n/3) 2^(-5(n+1)/6) = c^3 2^(5n/2 - 5n/3 - 5n/6) 2^(5/2 - 5/6) = c^3 2^(5/3)`, since `5/2 - 5/3 - 5/6 = 0`. Stationarity is then `Pi_{n-1} - Pi_n = 0`, which holds for all `n in Z`; on the half-lattice it fails only at `n = 0`, where `Pi_{-1} = 0` and energy must be injected. QED (script part 1f, exact).

The comparison that drives everything below: `2^(-5n/6)` decays strictly faster than the scale-invariant `2^(-n/2)` of C2, since `5/6 > 1/2`.

## Breakpoint

(skeleton)

## Transfer check

(skeleton)

## What would be needed next

(skeleton)

## How to reproduce

(skeleton)

## References

(skeleton)
