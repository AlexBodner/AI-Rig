# Line E — Frontier extension: logarithmically supercritical hyperdissipation

Status legend (one tag per claim): **PROVED-HERE** (complete proof written out), **KNOWN** (cited, not reproved), **COMPUTED** (script plus printed output below), **CONJECTURE**, **FAILED** (attempted, with the exact step where it breaks). Citation tags: [checked via WebSearch] / [unverified] as in `00_problem.md`.

## Goal

The Lions exponent `s = 5/4` is the exact place where the 3D Navier-Stokes energy stops being supercritical (`00_problem.md` section 1.6), and the only rigorous frontier that has ever been pushed past a criticality threshold for this equation is the logarithmic one just below it: Tao 2009 and Barbato-Morandin-Romito 2014 prove global regularity for dissipation symbols `|ξ|^{5/2}/g(|ξ|)^2` that are weaker than `|ξ|^{5/2}` by a divergent logarithmic amount. This line does three things: (E1) reproves the Lions result with every inequality visible and with the two independent reasons `5/4` is the threshold made explicit; (E3-E7) redoes the logarithmic estimate from scratch with the weight kept explicit, reproducing Tao's `g^4` condition and then proving that the whole family of weighted-enstrophy Gronwall schemes is *incapable* of reaching Barbato-Morandin-Romito's `g^2` condition, because the weight cancels identically; and (E4, E9) computes the cascade energy budget that says `g^2` is the ceiling of any argument built from the energy identity plus scaling, which is exactly the ceiling that transfer test T1 predicts.

The push is recorded as **FAILED**, with an obstruction that is itself proved. That is the honest outcome: the `g`-scale between `|ξ|^{5/2}/g^2` and `|ξ|^{5/2}` was already closed by Barbato-Morandin-Romito, and the argument below explains why the natural coercive-functional route cannot close it and what structure the successful route uses instead.

## Setting and notation

Throughout, `T^3 = R^3/Z^3`, `u` is smooth, divergence-free and mean-zero, `Λ = (-Δ)^{1/2}`, and `P` is the Leray projector. The generalized (hyperdissipative) system is

```text
∂_t u + (u·∇)u = -ν L u - ∇p,    div u = 0,    u(·,0) = u^0,        (E-1)
```

where `L` is the Fourier multiplier with symbol `ℓ(ξ) ≥ 0`, `ℓ(ξ) = |ξ|^{2s}` in the pure fractional case and `ℓ(ξ) = |ξ|^{5/2}/g(|ξ|)^2` in the logarithmic case. `D` denotes the multiplier with symbol `m(ξ) = ℓ(ξ)^{1/2}`, so `⟨Lu,u⟩ = ‖Du‖_2^2`. Mean-zero is preserved and lets us use homogeneous norms `‖u‖_{Ḣ^σ} = ‖Λ^σ u‖_2` freely.

Littlewood-Paley: `u = Σ_{N} u_N` over dyadic `N = 2^j`, `j ≥ 0`, with the usual annulus projections. Bernstein on `T^3`: `‖u_N‖_{L^q} ≤ C N^{3(1/r - 1/q)}‖u_N‖_{L^r}` for `r ≤ q`, and `‖Λ^σ u_N‖_{L^r} ≍ N^σ ‖u_N‖_{L^r}`. The three bookkeeping quantities used below, all at fixed time, are

```text
E   = ‖u‖_{L^2}                     (energy, bounded by the energy identity)
A^2 = ‖Du‖_{L^2}^2                  (energy-level dissipation, in L^1_t by the energy identity)
Y_h = Σ_N h(N)^2 N^2 ‖u_N‖_{L^2}^2  (weighted enstrophy; h ≡ 1 gives Y = ‖∇u‖_2^2)
```

and the shorthand `a_N = N^{5/4}g(N)^{-1}‖u_N‖_2` (so `A^2 = Σ_N a_N^2` up to a fixed constant) and `c_N = h(N)N‖u_N‖_2` (so `Y_h = Σ_N c_N^2`).

Two standing KNOWN inputs, used but not reproved. **(K1) Local well-posedness and the continuation criterion.** For (E-1) with `ℓ(ξ) ≥ c|ξ|^{2}` and smooth periodic mean-zero data there is a unique smooth solution on a maximal interval `[0,T*)`, and if `T* < ∞` then `‖u(t)‖_{H^1} → ∞` as `t → T*`; hence an a priori `H^1` bound on `[0,T*)` forces `T* = ∞` (standard; for the hyperdissipative case see Lions 1969, Tao 2009 [unverified as to exact statement]). **(K2) The Sobolev product estimate.** For `a, b < 3/2` and `a + b > 0`, `‖fg‖_{Ḣ^{a+b-3/2}(T^3)} ≤ C_{a,b}‖f‖_{Ḣ^a}‖g‖_{Ḣ^b}` (standard paraproduct estimate [unverified as to exact reference]).

## Results

| #   | Statement                                                                                                                                                                          | Status                                                                                                      |
| --- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| E1  | Global regularity for `(-Δ)^s` on `T^3`, `5/4 ≤ s ≤ 5/2`, with the explicit bound (E1-8) and the two independent reasons `5/4` is the threshold                                    | **PROVED-HERE** (uses K1, K2)                                                                               |
| E2  | Tao 2009 (`∫ds/(s g^4)=∞`), the conjecture and Barbato-Morandin-Romito 2014 (`∫ds/(s g^2)=∞`), Colombo-Haffter 2019 (`α ∈ (5/4-ε(M,δ), 5/4]`)                                      | **KNOWN**                                                                                                   |
| E3  | Osgood lemma plus the change of variables that turns a `g^p` loss at frequency `N_0 ≍ Y^κ` into the condition `∫ds/(s g^p)=∞`                                                      | **PROVED-HERE**                                                                                             |
| E4  | The cascade energy budget: `E_j^{1/2}` drops by `(ν/2)g(2^j)^{-2}` per stage, so `g^2` is the natural exponent; and the same computation returns the Lions exponent at general `α` | **PROVED-HERE** (scaling computation); sharpness **CONJECTURE**                                             |
| E5  | Frequency-localized bookkeeping with the weight explicit: `d/dt Y_h ≤ C ν^{-1} g(N_0)^4 Y_h A^2`, recovering Tao's condition from scratch                                          | **PROVED-HERE** under hypothesis (D); (D) itself **PROVED-HERE** on the diagonal, **CONJECTURE** in general |
| E6  | The push to `g^2` by a weighted enstrophy: the weight cancels identically, so every weight gives `g^4`                                                                             | **FAILED**; the obstruction is **PROVED-HERE**                                                              |
| E7  | The push via the commutator gain in the transport term: the absorption threshold is unchanged                                                                                      | **FAILED** (uses K3)                                                                                        |
| E8  | Shell-flux iteration: the structure that can beat the per-frequency threshold                                                                                                      | **CONJECTURE** with proof plan                                                                              |
| E9  | Script output: exponent arithmetic, weight cancellation, Osgood dichotomy, budgets, dyadic shell model                                                                             | **COMPUTED**                                                                                                |

### E1. Lions: global regularity for `(-Δ)^s`, `5/4 ≤ s ≤ 5/2` — **PROVED-HERE**

**Theorem E1.** Let `5/4 ≤ s ≤ 5/2` and `ν > 0`. For every smooth mean-zero divergence-free `u^0` on `T^3`, the solution of `∂_t u + (u·∇)u = -ν(-Δ)^s u - ∇p`, `div u = 0`, `u(·,0)=u^0` is global and smooth, and satisfies the explicit bound (E1-8) below.

Proof, in eight steps; every constant is named and every hypothesis of every cited inequality is checked. Write `E = ‖u^0‖_{L^2}`, `Y(t) = ‖Λu(t)‖_{L^2}^2`.

(1) **Energy identity.** Pair the equation with `u` in `L^2(T^3)`. `∫((u·∇)u)·u = (1/2)∫ u·∇|u|^2 = -(1/2)∫(div u)|u|^2 = 0`; `∫∇p·u = -∫ p div u = 0`; `⟨(-Δ)^su,u⟩ = ‖Λ^su‖_2^2`. Hence

```text
(1/2) d/dt ‖u‖_2^2 + ν ‖Λ^s u‖_2^2 = 0,   so   ‖u(t)‖_2 ≤ E   and   ∫_0^∞ ‖Λ^s u‖_2^2 dt ≤ E^2/(2ν).
```

(2) **Enstrophy identity.** Pair with `Λ^2 u`. The Leray projector `P` commutes with `Λ^2` and is self-adjoint with `Pu = u`, so `⟨Λ^2 u, P F⟩ = ⟨Λ^2 u, F⟩`; and `⟨Λ^2u,∇p⟩ = -⟨div Λ^2u, p⟩ = 0`. Hence, with `N := -⟨Λ^2u,(u·∇)u⟩`,

```text
(1/2) d/dt ‖Λ u‖_2^2 + ν ‖Λ^{1+s} u‖_2^2 = N.
```

(3) **Duality.** `div u = 0` gives `(u·∇)u = div(u⊗u)`, and `∫ div(u⊗u) dx = 0`, so its `Ḣ^σ` norms are defined for negative `σ` on the torus. By Plancherel `Λ^2 = Λ^{1+s}Λ^{1-s}` and `‖div F‖_{Ḣ^{1-s}} ≤ ‖F‖_{Ḣ^{2-s}}`, so

```text
|N| ≤ ‖Λ^{1+s}u‖_2 · ‖u⊗u‖_{Ḣ^{2-s}}.
```

(4) **Product estimate.** For `s > 2` the output index `2 - s` is negative, so first replace `u⊗u` by its mean-zero part `F := u⊗u - ∫_{T^3}u⊗u dx`; this changes nothing in step (3) because `div` annihilates constants, and it makes `‖F‖_{Ḣ^{2-s}}` well defined. All norms of `u⊗u` below are to be read as norms of `F`. Apply (K2) with `a = 1`, `b = 5/2 - s`. Its hypotheses: `a + b - 3/2 = 2 - s` (the required output index, exactly); `a = 1 < 3/2`; `b < 3/2 ⟺ s > 1`, true since `s ≥ 5/4`; `a + b = 7/2 - s > 0` since `s ≤ 5/2`. Hence `‖u⊗u‖_{Ḣ^{2-s}} ≤ C_1‖u‖_{Ḣ^1}‖u‖_{Ḣ^{5/2-s}}`.

(5) **Interpolation.** For `5/4 ≤ s ≤ 5/2` one has `0 ≤ 5/2-s ≤ s`, so `5/2-s = θ s` with `θ := 5/(2s) - 1 ∈ [0,1]`, and `‖u‖_{Ḣ^{θ s}} ≤ ‖u‖_{L^2}^{1-θ}‖Λ^s u‖_2^{θ}` (Hölder on the Fourier side).

(6) **Young.** Combining (3)-(5) and `xy ≤ (ν/2)x^2 + y^2/(2ν)`,

```text
|N| ≤ C_1 E^{1-θ} ‖Λ u‖_2 ‖Λ^s u‖_2^{θ} ‖Λ^{1+s}u‖_2
     ≤ (ν/2)‖Λ^{1+s}u‖_2^2 + (C_1^2/(2ν)) E^{2(1-θ)} ‖Λ^s u‖_2^{2θ} · Y.
```

(7) **Gronwall.** With `φ(t) := (C_1^2/ν) E^{2(1-θ)}‖Λ^s u(t)‖_2^{2θ}`, steps (2) and (6) give `d/dt Y + ν‖Λ^{1+s}u‖_2^2 ≤ φ(t) Y`. By Hölder in time with exponents `1/θ` and `1/(1-θ)`, and then (1),

```text
∫_0^T φ ≤ (C_1^2/ν) E^{2(1-θ)} T^{1-θ} ( ∫_0^T ‖Λ^s u‖_2^2 dt )^{θ} ≤ (C_1^2/ν) E^{2(1-θ)} T^{1-θ} (E^2/(2ν))^{θ} < ∞.
```

(8) **Conclusion.** `Y(t) ≤ Y(0) exp(∫_0^t φ)` on `[0,T*)`, so `‖u(t)‖_{H^1}` is bounded on every bounded subinterval; by (K1) `T* = ∞`. Explicitly, at `s = 5/4` (`θ = 1`) the bound is time-uniform:

```text
(E1-8)    ‖∇u(t)‖_2^2 ≤ ‖∇u^0‖_2^2 · exp( C_1^2 ‖u^0‖_2^2 / (2ν^2) )   for all t ≥ 0.       ∎
```

**Where `s = 5/4` enters — two independent places, both verified exactly in Part 1 of the script.**

- *Time integrability.* `θ = 5/(2s) - 1 ≤ 1` iff `s ≥ 5/4`. The energy identity puts `‖Λ^su‖_2^2` in `L^1_t` and nothing better; `‖Λ^su‖_2^{2θ}` is in `L^1_t(0,T)` for `θ ≤ 1` and is not controlled for `θ > 1`. Step (7) is the only step that fails for `s < 5/4`.
- *Scaling.* Under `u ↦ λ^{2s-1}u(λ·)` the quantity `‖u‖_{Ḣ^σ}^2` has degree `d(σ) = 4s - 5 + 2σ`. The Gronwall coefficient in (6) has degree `d(s) + d(1) = 10s - 8`; the enstrophy dissipation has degree `d(1+s) = 6s - 3`; they agree iff `s = 5/4`. So `5/4` is exactly the exponent at which the Gronwall factor can be the energy dissipation itself, with no room to spare and no loss taken; `s > 5/4` is genuinely subcritical and the `T^{1-θ}` factor in (7) is the price of the slack.

**Not proved here.** The range `s > 5/2` (choose `a = b = (7/2-s)/2 < 3/2` in step (4) and interpolate both factors; strictly easier, **KNOWN**). The whole-space version (the same proof works verbatim for Schwartz data on `R^3`; the only torus-specific point is the mean-zero convention in step (3)). Both (K1) and (K2) are used, not reproved.

### E2. The logarithmic frontier: Tao 2009, Barbato-Morandin-Romito 2014, Colombo-Haffter 2019 — **KNOWN**

**E2a (Tao 2009).** Consider `∂_t u + (u·∇)u = -D^2 u - ∇p`, `div u = 0` in dimension `d ≥ 3` (the snippet did not say whether the paper is stated on `R^d` or on the torus; the estimates involved are the same), where `D` is the Fourier multiplier with symbol `m(ξ)`. If `m(ξ) ≥ |ξ|^{(d+2)/4}/g(|ξ|)` for all sufficiently large `|ξ|`, for some nondecreasing `g:R^+ → R^+` with `∫_1^∞ ds/(s g(s)^4) = +∞`, then Schwartz data give a global smooth solution. The example given in the paper is `m(ξ) = |ξ|^{(d+2)/4}/log(2+|ξ|^2)^{1/4}`. In `d = 3` this is `m(ξ) = |ξ|^{5/4}/g(|ξ|)`, dissipation symbol `|ξ|^{5/2}/g(|ξ|)^2`. Source: T. Tao, Analysis & PDE 2 (2009) 361-366, arXiv:0906.3070. \[checked via WebSearch: the hypothesis `m(ξ) ≥ |ξ|^{(d+2)/4}/g(|ξ|)`, the condition `∫_1^∞ ds/(s g(s)^4) = +∞`, and the example `log^{1/4}(2+|ξ|^2)` were all returned verbatim in a search snippet; the proof was not read.\]

**E2b (the conjecture, and Barbato-Morandin-Romito 2014).** On the basis of a heuristic comparing the speed of propagation of a putative blowup up the frequency axis with the rate of dissipation, Tao conjectured that regularity should hold under the weaker condition `∫_1^∞ ds/(s g(s)^2) = +∞`. Barbato, Morandin and Romito proved global existence of smooth solutions "under the optimal condition on the correction to the dissipation", proving that conjecture. Source: D. Barbato, F. Morandin, M. Romito, Analysis & PDE 7 (2014) 2009-2027, arXiv:1407.6734. \[checked via WebSearch: the venue, volume, pages, the phrase "under the optimal condition on the correction to the dissipation", the attribution to Tao's conjecture, and the form `∫_1^∞ ds/(s G(s)^2) = ∞` of that conjecture together with its heuristic justification were all returned in search snippets; the proof was not read.\] A search snippet also describes the method: "an iterative estimate on the amount of energy stored in all shells larger than `j` at time `t` and the energy dissipated by these shells up to time `t`" [checked via WebSearch; the sentence describes the companion dyadic paper, D. Barbato, F. Morandin, M. Romito, "Global regularity for a logarithmically supercritical hyperdissipative dyadic equation", Dynamics of PDE 11 (2014) 39-52, arXiv:1403.2852 — checked via WebSearch: title, venue, volume, pages].

**E2c (what is between the two).** With `g(s) = log^{β}(2+s)`, `∫ ds/(s g^4) = ∞` iff `β ≤ 1/4` and `∫ ds/(s g^2) = ∞` iff `β ≤ 1/2`. So `g = log^{1/2}`, i.e. dissipation symbol `|ξ|^{5/2}/log(2+|ξ|)`, is the smallest natural example that BMR covers and Tao 2009 does not. Table A of Part 4 of the script tabulates the partial sums and confirms the split. **COMPUTED**.

**E2d (the only known route strictly below the Lions power).** M. Colombo, S. Haffter, "Global regularity for the hyperdissipative Navier-Stokes equation below the critical order", arXiv:1911.02600, J. Differential Equations (2020): for every `M` and every `δ > 0` there is an explicit `ε = ε(M,δ) > 0` such that for divergence-free `u^0` with `‖u^0‖_{H^δ} ≤ M`, the fractional Navier-Stokes system with order `α ∈ (5/4 - ε, 5/4]` has a unique smooth solution; this comes from a stability statement saying that the set of pairs (datum, order) giving smooth solutions is open in `H^{5/4} × (3/4, 5/4]`. [checked via WebSearch: authors, title, arXiv id, and the abstract as paraphrased here.] The crucial point for this line: `ε` depends on the size of the datum, so this is *not* a statement about the equation with a fixed order `α < 5/4`; the class of data covered shrinks as `α` decreases. Nobody has a regularity theorem for a fixed `α < 5/4` with arbitrary data.

**E2e (the negative side of the same frontier).** "Non-uniqueness of weak solutions for a logarithmically supercritical hyperdissipative Navier-Stokes system", arXiv:2406.05853 [checked via WebSearch: title and arXiv id only; authors and contents not verified]. Whatever its exact statement, it is a reminder that the *weak* solution theory at this frontier is not rigid, so E1 and everything below are statements about the smooth solution on `[0,T*)`, not about arbitrary weak solutions (transfer test T5).

### E3. The Osgood lemma that converts a log loss into a divergence condition — **PROVED-HERE**

Every result in this area has the same shape: the enstrophy-type quantity satisfies a differential inequality with a logarithmically bad coefficient, and the divergence condition on `g` is exactly the Osgood condition making that inequality globally solvable. Isolating this makes it visible that the choice of `g^4` versus `g^2` is decided entirely upstream, in the estimate, and not by any cleverness in the Gronwall step.

**Lemma E3.** Let `y_0 > 0`, let `G:[y_0,∞) → (0,∞)` be nondecreasing, let `φ ≥ 0` be in `L^1_{loc}[0,∞)`, and let `y:[0,T*) → (0,∞)` be locally absolutely continuous with `y'(t) ≤ φ(t) y(t) G(y(t))` for a.e. `t` with `y(t) ≥ y_0`. If `∫_{y_0}^∞ dr/(r G(r)) = ∞` then `sup_{[0,T]∩[0,T*)} y < ∞` for every `T < ∞`.

Proof. Replace `y` by `ỹ = max(y, y_0)`, still locally absolutely continuous, still satisfying the inequality a.e. (on `{y < y_0}` the derivative of `ỹ` is `0 ≤ φ ỹ G(ỹ)`). Set `Ψ(z) = ∫_{y_0}^z dr/(rG(r))`. `Ψ` is continuous, strictly increasing (the integrand is positive), `Ψ(y_0) = 0`, and `Ψ(z) → ∞` as `z → ∞` by hypothesis; so `Ψ:[y_0,∞) → [0,∞)` is a bijection with increasing inverse `Ψ^{-1}`. By the chain rule for absolutely continuous functions, `(Ψ∘ỹ)'(t) = ỹ'(t)/(ỹ G(ỹ)) ≤ φ(t)` a.e., and `Ψ∘ỹ` is absolutely continuous, so `Ψ(ỹ(t)) ≤ Ψ(ỹ(0)) + ∫_0^t φ`. Applying `Ψ^{-1}`, `y(t) ≤ ỹ(t) ≤ Ψ^{-1}(Ψ(ỹ(0)) + ∫_0^T φ) < ∞` for `t ≤ T`. ∎

**Corollary E3'** (change of variables; this is why the condition is stated on `g` and not on `G`). Let `κ > 0`, `p > 0`, and suppose `c_1 g(r^κ)^p ≤ G(r) ≤ c_2 g(r^κ)^p` for `r` large, with `g ≥ 1` nondecreasing. Then `∫^∞ dr/(rG(r)) = ∞` if and only if `∫^∞ ds/(s g(s)^p) = ∞`. Proof: `∫^∞ dr/(r g(r^κ)^p) = (1/κ)∫^∞ ds/(s g(s)^p)` under `s = r^κ`, and multiplying `G` by a constant does not change convergence. ∎

So: an estimate that loses `g(N_0)^p` at the frequency `N_0 ≍ Y^{κ}` where the nonlinearity stops being absorbable yields exactly the condition `∫^∞ ds/(s g(s)^p) = ∞`, whatever `κ > 0` is. **The exponent `p` is the only thing that matters, and it is decided by the estimate, not by the Gronwall argument.** Tao's theorem has `p = 4`; the cascade budget of E4 and Barbato-Morandin-Romito's theorem have `p = 2`.

### E4. The cascade energy budget: why `g^2` and not `g^4` is the natural threshold — **PROVED-HERE** (scaling computation) plus **CONJECTURE** (sharpness)

This is the heuristic that Tao's conjecture rests on, written out with every exponent, because it is what tells us where the ceiling of the whole method is. Fix the dissipation symbol `ℓ(ξ) = |ξ|^{5/2}/g(|ξ|)^2` and consider the only blowup scenario that the energy identity plus scaling does not forbid: a *concentrating frequency cascade*, in which at stage `j` essentially all of the remaining energy `E_j` occupies the single dyadic shell `N_j = 2^j`.

**Computation E4 — PROVED-HERE** (as a scaling computation about the four displayed quantities; it is not a theorem about solutions).

```text
energy in the shell            ‖u_{N_j}‖_2^2 = E_j
velocity amplitude             A_j ≍ N_j^{3/2} E_j^{1/2}          (Bernstein, saturated by a single bump of width 1/N_j)
eddy-turnover / transfer time  T_j ≍ (1/N_j)/A_j = N_j^{-5/2} E_j^{-1/2}
energy dissipated over T_j     ΔE_j ≍ ν ℓ(N_j) E_j T_j = ν N_j^{5/2} g(N_j)^{-2} · E_j · N_j^{-5/2} E_j^{-1/2}
                                    = ν E_j^{1/2} / g(N_j)^2
```

Since `E_{j+1} = E_j - ΔE_j`, the last line says `E_{j+1}^{1/2} = E_j^{1/2} - (ν/2) g(N_j)^{-2} + O(ΔE_j^2/E_j^{3/2})`, i.e. **`E_j^{1/2}` decreases by `(ν/2)g(2^j)^{-2}` per stage**. Two consequences.

- The cascade reaches infinite frequency with energy left over — a finite-time blowup, because `Σ_j T_j ≍ Σ_j 2^{-5j/2}E_j^{-1/2} < ∞` as long as `E_j` stays bounded below — if and only if `(ν/2) Σ_{j≥0} g(2^j)^{-2} < E_0^{1/2}`. It is drained at some finite stage for every initial energy and every `ν > 0` if and only if `Σ_j g(2^j)^{-2} = ∞`, which by the integral test is `∫_1^∞ ds/(s g(s)^2) = ∞`.
- The exponent `2` on `g` in that condition comes from one single place: `ΔE_j/E_j^{1/2}` is `ν ℓ(N_j) T_j` and both `ℓ(N_j) = N_j^{5/2}g^{-2}` and `T_j ≍ N_j^{-5/2}` are pinned, the first by the definition of the equation and the second by Bernstein plus the energy normalization. There is no freedom in it.

So `∫ ds/(s g(s)^2) = ∞` is the *natural* threshold, and Tao's 2009 theorem, with `g^4`, is a factor of two off in the exponent — the loss being an artifact of the estimate, which E5 pins down exactly. Barbato-Morandin-Romito reach the natural threshold.

**Sharpness — CONJECTURE.** The converse, that `Σ_j g(2^j)^{-2} < ∞` permits an actual blowup (for the true equation, or even for an averaged or dyadic model), is not proved here and I could not confirm a theorem of that form. What is confirmed: A. Cheskidov, "Blow-up in finite time for the dyadic model of the Navier-Stokes equations", Trans. Amer. Math. Soc. 360 (2008) 5101-5120, proves finite-time blowup for the dyadic model when the dissipation degree is `α < 1/3` in his normalization \[checked via WebSearch: the venue, pages, and the threshold `α < 1/3`; the normalization was not verified, and it is *not* the 3D Navier-Stokes normalization — see the note in Part 4 Table B\]. That is a power-loss blowup, not a logarithmic one. The word "optimal" in the BMR abstract is the strongest evidence available here for logarithmic sharpness, and it is a word in an abstract, not a checked theorem.

**Relation to T1.** Everything in Computation E4 uses only the energy identity, Bernstein's inequality, and dimensional analysis. All three survive Tao's averaging (item 2.2), so the ceiling it identifies is a ceiling for *abstract* arguments — see the transfer check.

### E5. The dyadic bookkeeping with the weight explicit: the scheme delivers `g^4` — **PROVED-HERE** (under hypothesis (D))

Now redo the E1 estimate for `ℓ(ξ) = |ξ|^{5/2}/g(|ξ|)^2` with the logarithm kept explicit at every frequency, and with a general weight in the coercive quantity so that the whole family of enstrophy-type Gronwall schemes is covered at once. The continuous version of E1 cannot be repeated: step (3) needs `‖Λ^{9/4}u‖_2`, and `‖Λ^{9/4}u‖_2^2 = Σ_N g(N)^2 · (N^{9/2}g(N)^{-2}‖u_N‖_2^2)` is not bounded by any multiple of the dissipation because `g` is unbounded. The frequencies must therefore be split, and the split point is the whole story.

**Standing hypotheses (H).** `g:[1,∞) → [1,∞)` nondecreasing with `g(2N) ≤ C_g g(N)`; `h` a positive weight on dyadic frequencies such that `N ↦ N h(N)/g(N)^2` is nondecreasing with limit `+∞` and doubling constant `C_1`, and `N ↦ N g(N)^2 h(N)` is nondecreasing. (Both hold for `h(N) = N^a`, `a > -1`, and `g` of log type.)

**Notation.** `c_N = h(N) N ‖u_N‖_2`, `Y_h = Σ_N c_N^2`; `a_N = N^{5/4}g(N)^{-1}‖u_N‖_2`, `A^2 = Σ_N a_N^2 ≍ ‖Du‖_2^2`; and the weighted dissipation `B_h^2 = Σ_N h(N)^2 N^{9/2}g(N)^{-2}‖u_N‖_2^2 = Σ_N N^{5/2}g(N)^{-2} c_N^2`, which is what `⟨Λ^2 u, Lu⟩` supplies after weighting. The energy identity gives `∫_0^∞ A(t)^2 dt ≤ ‖u^0‖_2^2/(2ν)`, exactly as in E1(1).

**Hypothesis (D) (the diagonal trilinear bound).** Writing `T_h` for the nonlinear contribution to `d/dt Y_h`, assume `|T_h| ≤ C_0 Σ_N h(N)^2 N^{9/2} ‖u_N‖_2^3`.

*(D) holds for the frequency-diagonal part of the trilinear form — **PROVED-HERE**.* For `N_1 ≍ N_2 ≍ N_3 = N`, using `(u·∇)u = div(u⊗u)`, one integration by parts, Hölder with exponents `(2,4,4)`, and Bernstein `‖f_N‖_{L^4(T^3)} ≤ C N^{3/4}‖f_N‖_{L^2}`:

```text
|⟨Λ^2 u_N, div(u_{N_1}⊗u_{N_2})⟩| = |⟨∇Λ^2 u_N, u_{N_1}⊗u_{N_2}⟩|
      ≤ ‖∇Λ^2 u_N‖_2 ‖u_{N_1}‖_4 ‖u_{N_2}‖_4  ≤ C N^3 · N^{3/4} · N^{3/4} ‖u_N‖_2 ‖u_{N_1}‖_2 ‖u_{N_2}‖_2 = C N^{9/2}‖u_N‖_2^3.
```

*(D) is NOT proved here for the off-diagonal parts* (high-high to low, and high-low to high). The high-low-to-high part is treated separately in E7 and, remarkably, gives the same threshold; the high-high-to-low part is better behaved because the output frequency carries the small factor `N_3^3`. Everything below is therefore conditional on (D), and (D) is flagged as **CONJECTURE** as a statement about the full trilinear form.

**Two exact rewritings of the summand** (both verified symbolically in Part 1 of the script, residual `0`):

```text
(i)   h(N)^2 N^{9/2} ‖u_N‖_2^3 = ( N^{5/2} g(N)^{-2} c_N^2 ) · ( g(N)^2 ‖u_N‖_2 )
(ii)  h(N)^2 N^{9/2} ‖u_N‖_2^3 = ( N g(N)^2 h(N) ) · a_N^2 · c_N
```

**Proposition E5.** Assume (H) and (D). Let `N_0 = N_0(t)` be the least dyadic frequency with `N h(N)/g(N)^2 ≥ C_0 Y_h(t)^{1/2}/ν`. Then, at every time,

```text
(E5-1)   d/dt Y_h + ν B_h^2  ≤  C_0 · ( N_0 g(N_0)^2 h(N_0) ) · Y_h^{1/2} · A^2
(E5-2)   d/dt Y_h            ≤  C_0^2 C_1 ν^{-1} · g(N_0)^4 · Y_h · A(t)^2.
```

Proof. `d/dt Y_h + 2ν B_h^2 = T_h`, and `|T_h| ≤ C_0 Σ_N h^2 N^{9/2}‖u_N‖_2^3` by (D). Split the sum at `N_0`.

*High frequencies `N > N_0`, via (i).* `‖u_N‖_2 = c_N/(h(N)N) ≤ Y_h^{1/2}/(h(N)N)`, so `g(N)^2‖u_N‖_2 ≤ g(N)^2 Y_h^{1/2}/(N h(N)) ≤ ν/C_0`, the last step being the definition of `N_0` together with the monotonicity of `N h(N)/g(N)^2` in (H). Hence `C_0 Σ_{N>N_0} (i) ≤ ν Σ_{N>N_0} N^{5/2}g(N)^{-2}c_N^2 ≤ ν B_h^2`, which the dissipation absorbs.

*Low frequencies `N ≤ N_0`, via (ii).* By the monotonicity of `N g(N)^2 h(N)` in (H) and `c_N ≤ Y_h^{1/2}`, `C_0 Σ_{N ≤ N_0} (ii) ≤ C_0 (N_0 g(N_0)^2 h(N_0)) Y_h^{1/2} Σ_N a_N^2`. That is (E5-1).

*From (E5-1) to (E5-2).* By minimality of `N_0` and the doubling constant `C_1` of `N h(N)/g(N)^2`, `N_0 h(N_0)/g(N_0)^2 ≤ C_1 C_0 Y_h^{1/2}/ν`, i.e. `N_0 h(N_0) ≤ C_1 C_0 g(N_0)^2 Y_h^{1/2}/ν`. Substituting into (E5-1) turns `N_0 g(N_0)^2 h(N_0) Y_h^{1/2}` into `C_1 C_0 ν^{-1} g(N_0)^4 Y_h`. ∎

**Corollary E5' (Tao's condition falls out, with the right exponent).** Take `h ≡ 1`. Then `N_0 ≍ g(N_0)^2 Y^{1/2}/ν`, so for `g` of log type `log N_0 = (1/2)(1 + o(1)) log Y` and `g(N_0)^4 ≍ g(Y)^4` up to a constant. By Lemma E3 with `G(Y) = c g(Y^{1/2})^4` and Corollary E3' (with `κ = 1/2`, `p = 4`) and `∫_0^∞ A^2 dt ≤ ‖u^0‖_2^2/(2ν) < ∞`, the enstrophy `Y` does not blow up in finite time provided `∫_1^∞ ds/(s g(s)^4) = ∞`. That is exactly Tao's hypothesis, recovered from scratch. **PROVED-HERE modulo (D) and (K1)-(K2).**

### E6. The push, attempt 1: weighted enstrophy. The weight cancels identically — **FAILED**, obstruction **PROVED-HERE**

The push attempted here was to reach the natural `g^2` condition of E4 by choosing the coercive functional better: replace the enstrophy `‖∇u‖_2^2` by a weighted enstrophy `Y_h`, hoping that a weight growing like `g(N)^2` would soak up one of the two logarithmic losses. It does not work, and the reason is an exact cancellation rather than a near miss.

**Proposition E6 (the obstruction) — PROVED-HERE.** Under (H) and (D), the Gronwall coefficient produced by the scheme of E5 is `C_0^2 C_1 ν^{-1} g(N_0)^4`, in which the weight `h` appears only through `N_0`. Consequently, if `h(N) = N^a` with `a > -1` and `g` is of log type, then `N_0 ≍ (g(N_0)^2 Y_h^{1/2}/ν)^{1/(1+a)}`, so `log N_0 = (1 + o(1)) (log Y_h)/(2(1+a))`, so `g(N_0)^4 ≍_a g(Y_h)^4`, and by Corollary E3' the Osgood condition delivered by the scheme is `∫_1^∞ ds/(s g(s)^4) = ∞` **for every such `h`**.

*The motivating weight itself.* The obvious candidate is `h = g^2`, chosen precisely so that the `g^2` in the threshold cancels. It satisfies (H), and the defining relation becomes `N_0 h(N_0)/g(N_0)^2 = N_0 ≍ Y_h^{1/2}/ν`, with no logarithm left in it at all. But (E5-2) is unchanged: `g(N_0)^4 ≍ g(Y_h)^4` for `g` of log type, so the delivered condition is again `∫ds/(s g^4) = ∞`. The logarithm was removed from the threshold and reappeared, undiminished, in the coefficient — which is the cancellation of Proposition E6 seen from the other side.

Proof. The first sentence is (E5-2), where `h` occurs nowhere except inside `N_0`. For the second, the defining relation `N_0^{1+a}/g(N_0)^2 ≍ Y_h^{1/2}/ν` gives, after taking logarithms and using `log g(N_0) = o(log N_0)` for `g` of log type, `(1+a) log N_0 = (1/2) log Y_h + o(log N_0 + log Y_h)`. Then Corollary E3' applies with `κ = 1/(2(1+a)) > 0` and `p = 4`. ∎

**Where exactly the `g^2` is lost, and why no weight can recover it.** Trace the two factors of `g^2` in `g(N_0)^4`.

- One factor comes from rewriting (ii): converting the two energy-dissipation factors `a_N` back into `‖u_N‖_2` costs `g(N)^2`, because `a_N` measures `Du` and the estimate needs `Λ^{5/4}u = g(Λ)·Du`. This is unavoidable: `A^2` is the energy dissipation, and the energy identity supplies nothing else in `L^1_t`.
- The other comes from the absorption threshold, and this is the informative one. Divide the high-frequency summand by the dissipation at the same frequency: from (i), `(nonlinear at N)/(dissipation at N) = C_0 g(N)^2‖u_N‖_2/ν`. **The weight has already cancelled at this point, before any choice of `h` is made.** The absorption condition at frequency `N` is the `h`-free statement

```text
(E6-1)      ‖u_N‖_{L^2} ≤ ν / (C_0 g(N)^2).
```

That is a condition on the solution at one frequency; a weight in the coercive functional is a change of how frequencies are *summed*, and cannot alter a per-frequency ratio. Raising `h` makes the threshold `N_0` smaller for a given value of `Y_h`, but `Y_h` is correspondingly larger, and Proposition E6 says the two effects cancel to the accuracy that the Osgood condition can see. Part 2 of the script shows this numerically: for `h(N) = N^a` with `a ∈ {-1/2, 0, 1/2, 1}`, the coefficient `g(N_0)^4` varies by a factor that converges, as `Y → ∞`, to exactly `((1+a_max)/(1+a_min))^{4β} = 4^{4β}` for `g = log^β` — `4.00` for `β = 1/4`, `16.06` for `β = 1/2` — i.e. by a bounded factor, and multiplying `g` by a constant does not change `∫ ds/(s g^p)`.

**The exact inequality that will not close, with its exponents.** To reach the `g^2` condition the low-frequency sum would have to obey

```text
(E6-2, target)    C_0 Σ_{N ≤ N_0} h(N)^2 N^{9/2}‖u_N‖_2^3   ≤   C ν^{-1} g(N_0)^2 · Y_h · A^2,
```

whereas the scheme proves only

```text
(E6-3, actual)    C_0 Σ_{N ≤ N_0} h(N)^2 N^{9/2}‖u_N‖_2^3   ≤   C_0 ( N_0 g(N_0)^2 h(N_0) ) Y_h^{1/2} A^2 = C ν^{-1} g(N_0)^4 · Y_h · A^2,
```

and the two differ by exactly `g(N_0)^2 = (N_0 h(N_0) ν)/(C_1 C_0 Y_h^{1/2})`, i.e. by exactly the factor that defines `N_0`. Making the low-frequency sum smaller by `g(N_0)^2` would require the absorption threshold to be `N h(N) ≳ Y_h^{1/2}/ν` instead of `N h(N)/g(N)^2 ≳ Y_h^{1/2}/ν`, which by (E6-1) is false at a single frequency. **FAILED**, and the failure is a proved obstruction rather than an unfinished computation. It is also not repaired by choosing the split point elsewhere: moving it below `N_0` leaves non-absorbable frequencies in the Gronwall term, and moving it above `N_0` enlarges the supremum in the low-frequency sum; `N_0` is the optimum of (E5-1).

### E7. The push, attempt 2: commutator gain in the transport term does not move the threshold — **FAILED**

Hypothesis (D) is crude: it treats the trilinear form as if there were no cancellation. The natural next attempt is that the real trilinear form is smaller than (D) at high frequencies, so that the absorption threshold (E6-1) improves. The dominant high-frequency interaction is transport of a high mode by a low mode, and there the divergence-free condition does produce a cancellation. It does not help, and the computation says why.

**Computation E7 — PROVED-HERE modulo the standard commutator estimate (K3).** Fix `M < N` dyadic and consider `T_{M,N} = -⟨Λ^2 u_N, (u_M·∇)u_N⟩`. Since `div u_M = 0`, `⟨Λ u_N, (u_M·∇)Λ u_N⟩ = 0`, so writing `Λ((u_M·∇)u_N) = (u_M·∇)Λ u_N + [Λ, u_M·∇]u_N`,

```text
T_{M,N} = -⟨Λ u_N, [Λ, u_M·∇] u_N⟩.
```

**(K3) KNOWN, in its frequency-localized form.** For `M ≤ N/8`, `‖[Λ, u_M·∇]u_N‖_{L^2} ≤ C ‖∇u_M‖_{L^∞}‖Λ u_N‖_{L^2}` [unverified as to exact reference]. This is the *localized* bound, not the general Kato-Ponce estimate, and the difference matters: the general one at `σ = 1` reads `‖[Λ,b]∂f‖_2 ≤ C(‖∇b‖_∞‖Λ f‖_2 + ‖Λ b‖_2‖∇f‖_∞)`, whose second term would contribute `M N^{7/2}‖u_M‖_2‖u_N‖_2^2` here — larger than the first by `(N/M)^{3/2}`, which would destroy the argument. That term is absent in the localized case because `û_M` is supported in `|η| ≍ M ≤ N/8` while `û_N` is supported in `|ξ-η| ≍ N`, so on the support of the product the commutator symbol `|ξ| - |ξ-η|` is `O(|η|) = O(M)`: all the growth in `ξ` cancels and only the first term survives. With Bernstein `‖∇u_M‖_∞ ≤ C M^{5/2}‖u_M‖_2`,

```text
|T_{M,N}| ≤ C ‖∇u_M‖_∞ · N^2 ‖u_N‖_2^2 ≤ C M^{5/2} ‖u_M‖_2 · N^2 ‖u_N‖_2^2.
```

This is a genuine gain over (D): the factor `N` that a naive bound would put on the transported field has been moved onto the low mode. Now weight and compare with the dissipation at frequency `N`, which is `ν N^{5/2}g(N)^{-2} c_N^2`; since `h(N)^2 N^2 ‖u_N‖_2^2 = c_N^2`, the weighted term is `C (Σ_{M<N} M^{5/2}‖u_M‖_2) c_N^2`, and

```text
Σ_{M<N} M^{5/2}‖u_M‖_2 = Σ_{M<N} M^{3/2} c_M / h(M) ≤ ( Σ_{M<N} M^{3/2}/h(M) ) · Y_h^{1/2} ≤ C N^{3/2} Y_h^{1/2}/h(N),
```

the last step by the geometric-sum dominance of the top term, valid whenever `M ↦ M^{3/2}/h(M)` is nondecreasing (true for `h = N^a`, `a < 3/2`). Absorption into the dissipation at frequency `N` therefore requires

```text
C N^{3/2} Y_h^{1/2}/h(N) ≤ ν N^{5/2} g(N)^{-2}      ⟺      N h(N)/g(N)^2 ≥ (C/ν) Y_h^{1/2},
```

**which is the same threshold `N_0` as in E5, with the same power of `g`.** The commutator gain of one derivative is exactly compensated by the fact that the low mode is measured in `L^∞` rather than `L^2`, and the `g(N)^2` survives untouched because it comes from the dissipation, not from the nonlinearity.

**Conclusion (this is the sharpened breakpoint).** The absorption threshold `N h(N)/g(N)^2 ≍ Y_h^{1/2}/ν` is produced independently by two quite different trilinear estimates — the crude diagonal bound (D) and the commutator-improved transport bound — so it is not an artifact of a lossy estimate. Together with E6, which says no weight moves it, the conclusion is that **what has to change to reach `g^2` is not the trilinear estimate and not the coercive functional, but the Gronwall structure itself**: the scheme compares a single scalar `Y_h(t)^{1/2}`, which knows only the total size of the solution, against the dissipation at each frequency separately, and therefore pays `g^2` at every scale. **FAILED**.

### E8. The push, attempt 3: shell-flux Gronwall (the direction that works in the literature) — **CONJECTURE** with proof plan

E6 and E7 together say the scalar Gronwall structure is the obstruction. The structure that replaces it is, according to the search snippet quoted in E2b, exactly what Barbato-Morandin-Romito use: an iteration over shells rather than a single differential inequality. I did not carry it out; writing it would exceed the two-page limit, so it is recorded as a plan.

**Conjecture E8 (proof plan, not a proof).** For dyadic `2^j` let

```text
E_{>j}(t) = Σ_{N > 2^j} ‖u_N(t)‖_2^2                        (energy above shell j at time t)
D_{>j}(t) = 2ν ∫_0^t Σ_{N > 2^j} N^{5/2} g(N)^{-2}‖u_N‖_2^2 ds   (energy dissipated above shell j up to t).
```

The frequency-localized energy identity gives `E_{>j}(t) + D_{>j}(t) = E_{>j}(0) + ∫_0^t Π_j(s) ds`, where `Π_j` is the energy flux through shell `j`, a trilinear expression supported on triads with at least one frequency `≍ 2^j`. The plan is to prove an iterative bound of the form `E_{>j+1} + D_{>j+1} ≤ (1 - c ν g(2^j)^{-2}) (E_{>j} + D_{>j}) + (lower order)`, so that after `J` steps the quantity above shell `J` is at most `exp(-cν Σ_{j<J} g(2^j)^{-2}) (E_{>0}+D_{>0})`, which tends to `0` as `J → ∞` precisely when `Σ_j g(2^j)^{-2} = ∞`; the energy above every shell then vanishes fast enough to give a subcritical bound and hence regularity by (K1).

**Why this can beat the scalar scheme, in one sentence.** The absorption threshold (E6-1) compares the dissipation at frequency `N` against `‖u_N‖_2 ≤ Y_h^{1/2}/(h(N)N)`, i.e. against the *total* size of the solution; the shell-flux scheme compares it against `E_{>j}^{1/2}`, the energy that has actually reached that part of the spectrum, which is itself already small because the preceding shells drained it — so the `g^2` from the threshold is paid once against a decreasing budget, rather than once per scale against a fixed one. This matches Computation E4 exactly, where the drain per stage is `ν g(2^j)^{-2}` applied to `E_j^{1/2}`.

**What I did not do.** The flux `Π_j` has a high-high-to-low part that is not obviously controlled by the two quantities above; the `(lower order)` term is not identified; and the passage from `E_{>j} → 0` at a quantified rate to a subcritical a priori bound is not written. Status **CONJECTURE**; the corresponding theorem is KNOWN (E2b) and the reconstruction here is not offered as a substitute for it.

### E9. Shell-model numerics for the cascade budget — **COMPUTED**

All four tables are produced by `line_e_frontier.py`; the numbers below are copied from its printed output.

**E9a — the exponent arithmetic of E1 and E5, exact in sympy (Part 1).** `d(σ) = 4s + 2σ - 5`; `d(s) + d(1) = 10s - 8`; `d(1+s) = 6s - 3`; equality iff `s = 5/4`. The E1 chain has `b = 5/2 - s`, `a + b - 3/2 = 2 - s` identically, `θ = (5/2-s)/s`, `θ ≤ 1` iff `s ≥ 5/4`, `θ ≥ 0` iff `s ≤ 5/2`. The nonlinear-to-dissipation ratio at frequency `N` in the E5 bookkeeping simplifies to `g(N)^2 ‖u_N‖_2` with the weight `h` absent, and the rewriting (ii) is verified with residual exactly `0`.

**E9b — the weight cancellation of E6 (Part 2).** `g(N_0)^4` for `h(N) = N^a`, and the ratio of the largest to the smallest entry in each row:

| `g`         | `log10 Y` | `a = -1/2` | `a = 0` | `a = 1/2` | `a = 1` | spread |
| ----------- | --------- | ---------- | ------- | --------- | ------- | ------ |
| `log^{1/4}` | 10        | 26.34      | 13.17   | 9.011     | 6.242   | 4.22   |
| `log^{1/4}` | 640       | 1481       | 740.3   | 493.5     | 370.1   | 4.002  |
| `log^{1/2}` | 10        | 930.2      | 211.9   | 94.17     | 48.07   | 19.35  |
| `log^{1/2}` | 640       | 2.217e6    | 5.532e5 | 2.456e5   | 1.38e5  | 16.06  |

The spread converges to `((1+a_max)/(1+a_min))^{4β} = 4^{4β}` (`4` and `16`), i.e. the weight changes the Gronwall coefficient only by a bounded factor, which the Osgood condition cannot see. This is Proposition E6, numerically.

**E9c — the Osgood dichotomy of E3 (Part 3).** For `y' = y log(2+y)^p`, `y(0)=1`: `p = 0.5` reaches `t = 20` with `y ≈ 9.9e48`; `p = 1.00` reaches `t = 5` (`y ≈ 6.8e35`) and overflows before `t = 20`, consistent with global existence and a very fast growth; `p = 1.02, 1.2, 2, 4` overflow at or before `t = 5`. The truncated Osgood integral (cutoffs `60` and `600`) grows (13.6 to 47.1 at `p=0.5`; 4.12 to 6.42 at `p=1`) exactly for `p ≤ 1` and is essentially constant (0.895 to 0.910 at `p=2`) for `p > 1`.

**E9d — the budgets, and why the frontier is numerically invisible (Part 4, Table A).** Partial sums `S_p(M) = Σ_{m<M} g(2^m)^{-p}`:

| `g`         | `S_2(10)` | `S_2(10^2)` | `S_2(10^3)` | `S_2(10^4)` | `S_4(10)` | `S_4(10^2)` | `S_4(10^3)` | `S_4(10^4)` |
| ----------- | --------- | ----------- | ----------- | ----------- | --------- | ----------- | ----------- | ----------- |
| `1`         | 10        | 100         | 1000        | 1e4         | 10        | 100         | 1000        | 1e4         |
| `log^{1/4}` | 6.097     | 22.65       | 74.64       | 238.9       | 4.038     | 7.426       | 10.75       | 14.08       |
| `log^{1/2}` | 4.038     | 7.426       | 10.75       | 14.08       | 2.206     | 2.404       | 2.423       | 2.425       |
| `log^{1}`   | 2.206     | 2.404       | 2.423       | 2.425       | 1.117     | 1.119       | 1.119       | 1.119       |
| `x^{0.05}`  | 7.466     | 14.92       | 14.93       | 14.93       | 5.794     | 7.725       | 7.725       | 7.725       |

`g = log^{1/2}` is the smallest natural example inside BMR's hypothesis and outside Tao's: `S_2` still grows at `M = 10^4` while `S_4` has already converged to `2.425`. The divergent sums grow like `log M`, which is the quantitative form of the statement that the Tao/BMR frontier cannot be seen numerically.

**E9e — dyadic shell model (Part 4, Table B).** The model `a_n' = c_{n-1}a_{n-1}^2 - c_n a_n a_{n+1} - ν d_n a_n` with `c_n = λ^{β(n+1)}`, `λ = 2`, `β = 5/2`, `d_n = λ^{(2β/3)n}/g(λ^n)^2`, `ν = 0.05`, `a_0(0) = 1`. **The exponent is `2β/3`, not `β`**: the dyadic model's own criticality is set by its constant-flux (Kolmogorov) cascade `a_n ≍ λ^{-βn/3}`, not by the concentrating cascade of E4, so its critical dissipation order is `2β/3` and it is *not* the Lions exponent. `log10 max_t Σ_n λ^{2n}a_n^2` versus truncation `M`, and the growth rate in decades per shell (a constant-flux inertial range gives `log10 λ^{1/3} = 0.1003`):

| `g`         | `M=10` | `M=14` | `M=18` | `M=22` | dec/shell | vs inertial |
| ----------- | ------ | ------ | ------ | ------ | --------- | ----------- |
| `1`         | 2.361  | 2.732  | 3.098  | 3.462  | 0.0910    | 0.91        |
| `log^{1/4}` | 2.769  | 3.240  | 3.687  | 4.124  | 0.1091    | 1.09        |
| `log^{1/2}` | 3.175  | 3.735  | 4.249  | 4.739  | 0.1225    | 1.22        |
| `log^{1}`   | 3.974  | 4.706  | 5.341  | 5.926  | 0.1462    | 1.46        |
| `x^{0.10}`  | 2.910  | 3.546  | 4.182  | 4.820  | 0.1596    | 1.59        |
| `x^{0.25}`  | 3.733  | 4.743  | 5.747  | 6.751  | 0.2509    | 2.50        |

Reading, including the negative part. The critical case `g = 1` sits at the inertial rate (ratio `0.91`, slightly damped) and the power losses sit well above it, so at 22 shells the model separates a power loss from criticality. It does **not** separate the three logarithmic cases from each other or from a power loss: the ratios interpolate smoothly, exactly as E9d predicts. So these numerics support the bookkeeping of E4-E7 and are **not** evidence for or against the divergence conditions. This is numerical evidence about a shell model, not a proof, and not about Navier-Stokes.

## Breakpoint

**Where this line stops short of the Clay problem, stated as a computation rather than as a mood.** Everything above is about a *modified* equation. Fefferman's (A)-(D) fix `ν Δ`, i.e. `α = 1`, and Lions' threshold is `α = 5/4`; the gap is a full power `|ξ|^{1/2}`, not a logarithm. Run Computation E4 at a general order `α`, with dissipation symbol `|ξ|^{2α}`:

```text
ΔE_j / E_j^{1/2} ≍ ν ℓ(N_j) T_j = ν N_j^{2α} · N_j^{-5/2} = ν N_j^{2α - 5/2} = ν 2^{(2α-5/2) j}.
```

The total drain over the whole cascade is `ν Σ_j 2^{(2α-5/2)j}`, which **diverges iff `α ≥ 5/4` and converges iff `α < 5/4`** — the Lions exponent, recovered as an integral test. At `α = 1` the series is `ν Σ_j 2^{-j/2} = ν(2+√2)`, a *finite number of order `ν`*. So a concentrating cascade at `α = 1` loses only `O(ν)` of `E_0^{1/2}` in total, no matter how many decades it traverses; the energy identity does not drain it at all. The logarithmic frontier of Tao and Barbato-Morandin-Romito lives entirely inside the borderline case `2α = 5/2` where that series is exactly `Σ_j g(2^j)^{-2}`, i.e. marginally divergent, and improving `g^4` to `g^2` moves the boundary from `log^{1/4}` to `log^{1/2}`. Measured against the missing power `|ξ|^{1/2}`, that is not progress towards `α = 1`; it is the sharpening of a boundary case.

**The internal breakpoint of the push.** The push (E6, E7) is **FAILED** at a specific inequality: to obtain the natural exponent `g^2` the low-frequency sum would have to satisfy (E6-2) whereas the scheme proves only (E6-3), the two differing by exactly the factor `g(N_0)^2` that defines the absorption threshold. That threshold is the per-frequency inequality `‖u_N‖_2 ≤ ν/(C_0 g(N)^2)` of (E6-1), in which no weight appears (E6) and which is unchanged by the commutator gain in the transport term (E7). What must change is the Gronwall structure — a scalar `Y_h(t)` compared against the dissipation at every frequency separately pays `g^2` at every scale — and the structure that does change it is the shell-flux iteration of E8, recorded as **CONJECTURE** and known in the literature (E2b).

## Transfer check

Rows are the no-go items of `00_problem.md` section 2 and the transfer tests of section 3. "Object tested" is E1 plus E5-E7 (the regularity scheme) unless stated.

| Item                                  | Verdict                                       | Because                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| ------------------------------------- | --------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 2.1 supercriticality barrier          | passes, and is the point                      | The scheme uses only the energy identity plus scaling, which is exactly why it stops at a modified dissipation; at `α = 1` the Breakpoint computation shows the cascade drain is `O(ν)`, i.e. the barrier bites                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| 2.2 averaged blowup                   | passes                                        | The argument would apply verbatim to `∂_t u = -D^2 u + B̃(u,u)` (see T1); this is consistent, since Tao's blowup is at `α = 1` and the argument is empty there                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| 2.3 no backward self-similar profiles | not applicable                                | No blowup ansatz is proposed                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| 2.4 `L^∞_t L^3_x` (ESS)               | not applicable                                | Statements are for modified dissipation; ESS is an `α = 1` theorem                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| 2.5 LPS and Leray lower bounds        | passes                                        | E1 proves an `H^1` bound only for `s ≥ 5/4`, where `H^1` is strictly subcritical; nothing is claimed at `s = 1`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| 2.6 partial regularity (CKN)          | not applicable                                | No local analysis is used                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| 2.7 Beale-Kato-Majda                  | not applicable                                | No vorticity criterion is used                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| 2.8 2D global regularity              | passes vacuously                              | The line proves regularity, not blowup; in 2D the analogue of E1 holds at `s ≥ 1` by the same chain with `d = 2` (not written out)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| 2.9 axisymmetric swirl-free           | not applicable                                | No symmetry reduction is used                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| 2.10 hyperdissipative regularity      | this line IS this item                        | E1 reproves the `s ≥ 5/4` half; E2 states the logarithmic half; E5 recovers Tao's condition **under hypothesis (D)**, which is proved here only on the frequency diagonal; E6-E7 record the failed push to BMR's                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| 2.11 non-uniqueness of weak solutions | passes                                        | Everything is about the unique smooth solution on `[0,T*)` from (K1); no statement is made about Leray-Hopf solutions                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| 2.12 instantaneous Type I blowup      | passes                                        | Same reason; E1 step (8) uses smoothness at every `t < T*` through the continuation criterion                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| 2.13 Euler blowup                     | passes                                        | Every estimate has a factor `ν^{-1}` (see (E1-8)) and degenerates as `ν → 0`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| 2.14 unstable singularities (PINN)    | not applicable                                | No numerical profile search is attempted; E9e is a shell model, and its negative reading is reported                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| 2.15 critical spaces, dyadic models   | passes                                        | E9e's dyadic model is used only as an analogy and its *different* criticality (`2β/3`, not `β`) is stated explicitly                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| **T1 averaging**                      | **passes, and this is the ceiling statement** | E1, E4, E5, E6 use only: the energy identity, scaling, Bernstein/Littlewood-Paley, Hölder, and the size bound (D) on the trilinear form. Every one of these survives Tao's averaging, so the argument applies verbatim to the averaged equation with the same dissipation. That is not a contradiction — Tao's blowup is at `α = 1`, and at `α = 1` the argument proves nothing (Breakpoint) — but it does fix the ceiling: **no strengthening of this line can ever reach an order at which the averaged equation blows up.** One exception is worth flagging: E7 uses the frequency-localized cancellation `⟨B(u_M,f),f⟩ = 0`, which is strictly stronger than the energy cancellation `⟨B(u,u),u⟩ = 0` that averaging preserves; whether Tao's `B̃` retains it is a question, tagged **CONJECTURE**, and nothing here is built on the answer |
| T2 comparison systems                 | not applicable                                | No blowup mechanism is proposed                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| T3 blowup profile constraints         | not applicable                                | Same                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| T4 inviscid limit                     | passes                                        | See 2.13                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| T5 weak solutions                     | passes                                        | E1 is a statement about the unique smooth solution; the `H^1` bound is propagated by Gronwall along that solution, and the continuation criterion (K1) is the only bridge to global existence. The non-uniqueness result at this very frontier (E2e) confirms that this restriction is not cosmetic                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| T6 too-easy                           | passes                                        | E1 does close an `H^1` bound — but only for `s ≥ 5/4`, where it is Lions' 1969 theorem, KNOWN, not new. The push E6-E7 does **not** close, and is recorded as FAILED with the exact inequality                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |

## What would be needed next

- **Finish E8 or read BMR.** The shell-flux iteration is the only structure identified here that can beat the per-frequency threshold (E6-1). The two missing pieces are a bound on the high-high-to-low part of the flux `Π_j` in terms of `E_{>j}` and `D_{>j}` alone, and the passage from a quantified decay of `E_{>j}` to a subcritical a priori bound. Both are presumably in arXiv:1407.6734, which this environment could not fetch.
- **Settle hypothesis (D) for the full trilinear form.** E5 and E6 are conditional on it. The high-high-to-low part should be strictly better (the output frequency supplies `N_3^3` with `N_3 ≪ N_1 ≍ N_2`), and E7 handles high-low-to-high, but neither is written as a complete paraproduct sum here.
- **Decide the sharpness question of E4.** Is there a blowup theorem, for the averaged equation or for a dyadic model at its own critical order, when `Σ_j g(2^j)^{-2} < ∞`? If yes, then the `g`-scale is provably closed and this whole direction is finished; if no, then "optimal" in the BMR abstract is a statement about the method, not about the equation, and there may be room below `g^2`. I could not resolve this from search snippets.
- **Does Tao's averaged operator retain `⟨B̃(u_M,f),f⟩ = 0`?** Rotation conjugation preserves it; conjugation by a non-unitary order-zero multiplier need not. If it does not, then E7 (and any argument using the localized transport cancellation) is *not* abstract in the T1 sense, and the T1 ceiling for such arguments is higher than stated. **CONJECTURE**, worth one page.
- **The Colombo-Haffter direction (E2d) is the only one that moves the power.** The interesting question it raises is whether the `ε(M,δ)` can be made independent of `M` for data in a fixed *critical* ball rather than a subcritical one; if so the argument would iterate. Nothing here bears on it.

## How to reproduce

```text
cd research/navier-stokes
python3 line_e_frontier.py            # about 3 minutes on 4 CPUs (under the 5-minute cap); prints Parts 1-4 in order
ruff check --line-length 120 line_e_frontier.py
ruff format --check --line-length 120 line_e_frontier.py
```

The script is deterministic (a seed is set and printed but no randomness is used). Part 1 is sympy, exact. Part 2 works entirely in log space so that frequencies up to `2^4000` are handled without overflow. Parts 3 and 4 use `scipy.integrate.solve_ivp` with the stiff `Radau` method; "overflow" in Part 3 means the solver could not reach the final time, which is the numerical signature of finite-time blowup.

## References

- J.-L. Lions, *Quelques méthodes de résolution des problèmes aux limites non linéaires*, Dunod/Gauthier-Villars, 1969 — global regularity for `α ≥ 5/4`. \[checked via WebSearch: book and the `α ≥ 5/4` attribution, via `00_problem.md` item 2.10; the proof in E1 is written independently and is not claimed to be Lions'.\] <!-- codespell:ignore methodes,problemes,resolution -->
- T. Tao, "Global regularity for a logarithmically supercritical hyperdissipative Navier-Stokes equation", Analysis & PDE 2 (2009) 361-366, arXiv:0906.3070. \[checked via WebSearch: venue, volume, pages, arXiv id, and the hypothesis `m(ξ) ≥ |ξ|^{(d+2)/4}/g(|ξ|)` with `∫_1^∞ ds/(s g(s)^4) = +∞` and the example `log^{1/4}(2+|ξ|^2)`.\]
- D. Barbato, F. Morandin, M. Romito, "Global regularity for a slightly supercritical hyperdissipative Navier-Stokes system", Analysis & PDE 7 (2014) 2009-2027, arXiv:1407.6734. \[checked via WebSearch: venue, volume, pages, arXiv id, the abstract phrase "under the optimal condition on the correction to the dissipation", and that it proves a conjecture of Tao whose form `∫_1^∞ ds/(s G(s)^2) = ∞` and heuristic justification were also returned in a snippet.\]
- D. Barbato, F. Morandin, M. Romito, "Global regularity for a logarithmically supercritical hyperdissipative dyadic equation", Dynamics of Partial Differential Equations 11 (2014) 39-52, arXiv:1403.2852. \[checked via WebSearch: title, venue, volume, pages, arXiv id, and the description of the method as an iterative estimate on the energy above shell `j` and the energy dissipated by those shells.\]
- M. Colombo, S. Haffter, "Global regularity for the hyperdissipative Navier-Stokes equation below the critical order", arXiv:1911.02600, J. Differential Equations (2020). \[checked via WebSearch: authors, title, arXiv id, and the abstract, including the `ε = ε(M,δ)` dependence and the openness statement in `H^{5/4} × (3/4, 5/4]`.\]
- A. Cheskidov, "Blow-up in finite time for the dyadic model of the Navier-Stokes equations", Trans. Amer. Math. Soc. 360 (2008) 5101-5120. \[checked via WebSearch: venue, volume, pages, and the threshold `α < 1/3`; the normalization behind that number was not verified.\]
- N. H. Katz, N. Pavlović, "Finite time blow-up for a dyadic model of the Navier-Stokes equations", Trans. Amer. Math. Soc. 357 (2005). [unverified; cited only as the origin of the shell model used in E9e.]
- T. Tao, "Finite time blowup for an averaged three-dimensional Navier-Stokes equation", J. Amer. Math. Soc. 29 (2016) 601-674, arXiv:1402.0290 — the source of transfer test T1; see `00_problem.md` item 2.2. \[checked via WebSearch, in `00_problem.md`.\]
- "Non-uniqueness of weak solutions for a logarithmically supercritical hyperdissipative Navier-Stokes system", arXiv:2406.05853. [checked via WebSearch: title and arXiv id only; authors and content not verified.]
- Kato-Ponce commutator estimate (K3) and the Sobolev product estimate (K2). [unverified as to exact references; both are standard and are used as black boxes.]
