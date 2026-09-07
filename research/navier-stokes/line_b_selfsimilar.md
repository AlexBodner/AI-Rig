# Line B — Self-similar blowup profiles and Liouville theorems (Type I versus Type II)

Status legend (one tag per claim): **PROVED-HERE** (complete proof written out), **KNOWN** (cited, not reproved), **COMPUTED** (script plus printed output below), **CONJECTURE**, **FAILED** (attempted, with the exact step where it breaks). Citation tags: [checked via WebSearch] / [unverified] as in `00_problem.md`.

Files: this note and the script `line_b_selfsimilar.py` (same directory; every table and residual below is pasted from its output, runtime 16 s on 4 CPUs). No-go items and transfer tests are cited by their numbers in `00_problem.md` and are not restated.

## Goal

Determine exactly which blowup ansatz for 3D Navier–Stokes survives the known Liouville theorems: derive the backward self-similar profile equation and the Nečas–Růžička–Šverák (NRS) identity with the correct drift, sign and constant; verify the identity symbolically; write out the maximum-principle Liouville argument with every hypothesis named; then classify the general ansatz `u = (T−t)^{-α} U(x/(T−t)^β)` by `β` and say precisely which route is left.

Outcome in one line. Writing `s := T − t`, the nonlinear balance forces `α = 1 − β` (B4), the viscous balance and Leray's lower bounds force `β ≤ 1/2` (B5), finite energy forces `β ≥ 2/5` (B6), and `β = 1/2` is killed by NRŠ/Tsai/Chae–Wolf and independently by ESS (B3, B9). The round-1 referee report (`refutations/line_b_r1.md`) showed that this single sentence ran together two very different statements, so the outcome is now recorded in two parts.

1. **For the exactly self-similar ansatz the line is complete and the answer is negative — PROVED-HERE modulo one KNOWN input.** The nonlinear balance forces `α = 1 − β` (B4.1); rigidity in `s` forces `β = 1/2` (B4.2, under a growth hypothesis on `U` that is now stated); finite energy forces `U ∈ L^2` (B6.1), and `U ∈ L^2 ∩ L^∞` gives `U ∈ L^3`, so NRŠ (2.3, KNOWN) applies and `U ≡ 0`. **No exactly self-similar blowup exists in the Clay class, for any `β`.** This route needs neither Chae–Wolf nor the matched-tail discussion.
1. **For the *asymptotically* self-similar class the line produces a window and nothing more — CONJECTURE.** The arithmetic gives `2/5 ≤ β < 1/2` (B5, B6, B7.1), i.e. Type II with an Euler self-similar leading profile and formally perturbative viscosity. But (a) the class is not defined here and every norm identity used is an exact computation for the exact ansatz (B5.5); (b) the right endpoint for the physically relevant non-`L^3` matched profile rests on Chae–Wolf's unread hypotheses (B7.1); and (c) **nonemptiness of this window carries no blowup information**, because the entire window arithmetic applies verbatim to swirl-free axisymmetric 3D Navier–Stokes, which is globally regular (B11, the positive control the first version of this note failed to run).

The same arithmetic in dimension `d` with dissipation `−ν(−Δ)^σ` gives the window `2/(d+2) ≤ β ≤ 1/(2σ)`, empty exactly when `σ > (d+2)/4` — the Lions exponent, `5/4` at `d = 3` and `1` at `d = 2` (B10; PROVED-HERE for the exact ansatz, CONJECTURE for the asymptotic class, whose upper edge needs a hyperdissipative lower bound this project cannot verify). The line stops (Breakpoint) at the fact that no smooth decaying Euler self-similar profile with `β` in the window is known to exist, that Chae's companion note on *locally* self-similar Euler singularities may already close it, and that the viscous correction, though formally `O(s^{1−2β})`, multiplies the highest-order operator and is therefore a singular perturbation.

## Setting and notation

- `s := T* − t > 0`, so `∂_t = −∂_s`. Ansatz `u(x,t) = s^{-α} U(y)`, `y = x/s^β`, `p(x,t) = s^{-2α} P(y)`, with `α, β > 0`, `U : R^3 → R^3` smooth, `div U = 0`. Vorticity of the profile `ω := curl U` (not to be confused with the vorticity of `u`, which is `s^{-α−β} ω(y)`).
- Navier–Stokes `∂_t u + (u·∇)u = νΔu − ∇p`, `div u = 0`, `f ≡ 0`; `ν = 1` in B1–B3 and B8 (it can be scaled away in the profile equation by `U ↦ νU`, `P ↦ ν^2 P`), kept explicit in B4–B7 and B10.
- The pressure exponent is not free: `p = Σ_{i,j} R_i R_j (u_i u_j)` (Riesz transforms), a scaling-invariant operator, so `p` scales as `|u|^2`, i.e. as `s^{-2α}` in the ansatz. This is used in B4 and is the reason the pressure cannot be used to rebalance the exponents.
- `Π := P + |U|^2/2 + (1/2) y·U` (the NRS "head pressure"), `b := y/2 + U` (the NRS drift). `Res` denotes the residual of the profile equation, `Res := ΔU − [ (1/2)U + (1/2)(y·∇)U + (U·∇)U + ∇P ]`, so that `Res ≡ 0` is exactly the `ν = 1`, `β = 1/2` profile equation.
- Similarity variables for B8: `τ := −log s`, `V(y,τ) := s^{1/2} u(s^{1/2} y, t)`; exactly self-similar means `∂_τ V = 0`, discretely self-similar (DSS) with factor `λ > 1` means `V` is periodic in `τ` with period `2 log λ`.
- KNOWN inputs used without proof, cited by no-go item: local well-posedness and the `H^1` continuation criterion (2.5, 2.15); Leray's lower bounds `‖u(t)‖_{L^∞} ≥ c (T*−t)^{-1/2}` and `‖∇u(t)‖_{L^2} ≥ c (T*−t)^{-1/4}` (2.5); ESS/Seregin `lim_{t→T*} ‖u(t)‖_{L^3} = ∞` (2.4); BKM (2.7); CKN `ε`-regularity (2.6); NRS/Tsai (2.3); Lions' hyperdissipative regularity (2.10).

## Results

Summary of statuses, revised after round 1. B1, B2, B4.1, B4.2, B5.1–B5.4, B6.1, B6.2a, B8, B11 PROVED-HERE (B1, B2, B8 additionally COMPUTED). B3 PROVED-HERE under three explicitly stated hypotheses (D0)–(D2); the decay lemma is KNOWN, and the claim that NRŠ/Tsai supply *exactly* that lemma is now CONJECTURE. B5.5 (every extension of an exclusion from the exact ansatz to the asymptotically self-similar class) CONJECTURE. B6.2b CONJECTURE. B7.1 KNOWN, not PROVED-HERE, because its right endpoint rests on Chae–Wolf's unread hypotheses; B7.2 PROVED-HERE (formal arithmetic); B7.3 CONJECTURE. B9 KNOWN for the bibliographic data and CONJECTURE for "Chae's theorem does not apply". B10 PROVED-HERE and COMPUTED for the exact ansatz, CONJECTURE for the asymptotic class at general `σ`. Three FAILED items are recorded: the energy route to `β > 1/2` (inside B5), the transfer-table verdicts for 2.9 and T2 as they stood in round 1 (B11), and the claim that script Part 5 verified the profile identity (B11).

**Corrections carried out in this revision**, each forced by `refutations/line_b_r1.md` and each keeping the criticized claim with its new status rather than deleting it: B7.1 downgraded PROVED-HERE → KNOWN; B3's "their content is precisely the missing lemma" downgraded to CONJECTURE; the asymptotic extensions of B5.2/B5.3 downgraded to CONJECTURE and collected in B5.5; transfer rows 2.9 and T2 recorded FAILED and rewritten against the control now run in B11; B9's non-applicability of Chae 2007 downgraded to CONJECTURE; B10's general-`σ` upper edge downgraded to CONJECTURE, with the exact-ansatz case reproved so that it no longer uses the unverified hyperdissipative lower bound; B10's COMPUTED tag made real (script Part 4b now computes the CKN exponent instead of printing it); B6.2 split by *class* rather than by tail shape; B4.2's missing growth hypothesis supplied and its broken proof route recorded; B8's false DSS-versus-self-similar contrast corrected and its `β = 1/2`-only scope proved; script Part 5 rebuilt so that it actually exercises `div U = 0`; B1's false remark about no-go item 2.3 corrected and recorded as an erratum for `00_problem.md`.

### B1. The backward self-similar profile equation — **PROVED-HERE** and **COMPUTED**

Take `α = β = 1/2`: `u(x,t) = s^{-1/2} U(y)`, `y = x s^{-1/2}`, `p = s^{-1} P(y)`. Since `x s^{-3/2} = y s^{-1}`,

```text
∂_s [ s^{-1/2} U(x s^{-1/2}) ] = −(1/2) s^{-3/2} U − (1/2) s^{-3/2} (y·∇)U,
```

so `∂_t u = −∂_s u = (1/2) s^{-3/2} [ U + (y·∇)U ]`. The other three terms are `(u·∇_x)u = s^{-3/2}(U·∇)U`, `νΔ_x u = ν s^{-3/2} ΔU`, `∇_x p = s^{-3/2} ∇P`. Every term carries the same power `s^{-3/2}`, and dividing by it gives the profile equation

```text
(1/2) U + (1/2)(y·∇)U + (U·∇)U + ∇P = ν ΔU,      div U = 0.
```

At `ν = 1` this is `Res = 0` with `Res` as in "Setting and notation", For its relation to the equation quoted in no-go item 2.3, see the erratum below: the first version of this note got that comparison wrong. **COMPUTED** (script Part 1, run with symbolic `α, β` and a random exact-rational profile; see B4 for the general-exponent table): the residual of the substituted equation against the claimed profile form is `0` after imposing `α = 1 − β`, and nonzero before.

**Erratum for `00_problem.md` item 2.3 — PROVED-HERE.** The first version of this note said that the equation above "is the equation quoted in no-go item 2.3 with the opposite overall sign". That was wrong, and the round-1 referee is right that it hid an error in `00_problem.md`. Item 2.3 quotes `−ΔU − (1/2)U − (1/2)(y·∇)U + (U·∇)U + ∇P = 0`; the derivation above gives `−ΔU + (1/2)U + (1/2)(y·∇)U + (U·∇)U + ∇P = 0`. Multiplying either by `−1` gives an equivalent equation, so "the opposite overall sign" would mean no difference at all; what actually differs is the sign of the two drift terms *only*, which is a genuinely different equation. The derivation above is the **backward** one (it uses `s = T∗ − t`, `∂_t = −∂_s`); 2.3's version is the **forward** self-similar profile equation, satisfied by `u = t^{-1/2}U(x/t^{1/2})`. The distinction is load-bearing, not cosmetic: with 2.3's drift the NRŠ drift is `b = −y/2 + U`, so `div b = −3/2 < 0`, and Step 2 of Theorem B3 — which needs the coefficient of `∫ Π χ_R ≤ 0` to be positive — fails outright. I own only this file, so I record the correction here instead of editing `00_problem.md`; it is listed as unresolved in the hand-back.

### B2. The Nečas–Růžička–Šverák identity — **PROVED-HERE** and **COMPUTED**

**Proposition B2.** Let `U` be a smooth divergence-free field on an open set, `P` smooth, `ω = curl U`, `Π = P + |U|^2/2 + (1/2) y·U`, `b = y/2 + U`. Then, identically in `(U,P)`,

```text
ΔΠ − b·∇Π − |ω|^2 = − div(Res) + b·Res.
```

In particular, if `U` solves the profile equation of B1 with `ν = 1` (`Res = 0`) then

```text
ΔΠ − ( y/2 + U )·∇Π = |curl U|^2.
```

So the drift is `y/2 + U` (not `y/2`, not `U`, not their negative), the sign is such that `Π` is a **subsolution** of a drift-diffusion operator with no zeroth-order term, and the constant in front of `|curl U|^2` is exactly `+1`.

Proof. Three pointwise vector identities, valid for any smooth field with `div U = 0`:

- `(U·∇)U = ∇(|U|^2/2) − U × ω` (the Lamb identity).
- `∇(y·U) = (y·∇)U + (U·∇)y + y × curl U + U × curl y = (y·∇)U + U + y × ω`, since `(U·∇)y = U` and `curl y = 0`; hence `(1/2)U + (1/2)(y·∇)U = (1/2)∇(y·U) − (1/2) y × ω`.
- `ΔU = ∇(div U) − curl curl U = −curl ω`.

Substituting the first two into the definition of `Res` and using the third,

```text
Res = −curl ω − ∇Π + b × ω,        b = y/2 + U.      (B2.1)
```

Take the divergence of (B2.1). Using `div curl = 0`, `div(a × ω) = ω·curl a − a·curl ω`, and `curl b = curl(y/2) + curl U = ω`:

```text
div Res = − ΔΠ + |ω|^2 − b·curl ω.
```

From (B2.1), `curl ω = −Res − ∇Π + b × ω`, so `b·curl ω = −b·Res − b·∇Π` (because `b·(b × ω) = 0`). Substituting gives `div Res = −ΔΠ + |ω|^2 + b·Res + b·∇Π`, which is the claim. ∎

**COMPUTED** (script Part 2). The identity is verified by sympy as an exact symbolic identity in three *arbitrary* functions: `U = curl A` with `A = (A_1,A_2,A_3)` undetermined sympy functions of `(y_1,y_2,y_3)` (by the Poincaré lemma on `R^3` this parametrizes exactly the smooth divergence-free fields) and `P` a fourth arbitrary function. Printed output:

```text
- div U = curl-free check, div(curl A) = 0
- residual of  [Delta Pi - b.grad Pi - |omega|^2] - [-div Res + b.Res]  =  0
- identity verified: True
```

Controls, to show the check has teeth (T6 discipline): each of the following variants of the same computation is *not* identically zero.

| variant                                                    | residual identically zero? |
| ---------------------------------------------------------- | -------------------------- |
| `U` not divergence-free (random rational polynomial field) | False                      |
| constant `c = −1` in front of `\|ω\|^2`                    | False                      |
| constant `c = 2`                                           | False                      |
| constant `c = 1/2`                                         | False                      |
| drift `b = y/2` only                                       | False                      |
| drift `b = U` only                                         | False                      |
| drift `b = −(y/2 + U)`                                     | False                      |

### B3. The Liouville argument from the identity — **PROVED-HERE** under (D0)–(D2); the decay lemma is **KNOWN**

**Theorem B3.** Let `U` be a smooth divergence-free solution of the `ν = 1` profile equation on all of `R^3` and let `Π` be as above. Assume

- **(D0)** `U ∈ L^∞(R^3)` and `U(y) → 0` as `|y| → ∞`;
- **(D1)** `limsup_{|y|→∞} Π(y) ≤ 0`;
- **(D2)** `liminf_{R→∞} ∫_{R ≤ |y| ≤ 2R} |Π(y)| dy = 0`.

Then `U ≡ 0` (and hence the self-similar solution is trivial).

Proof, in three steps; each hypothesis is used exactly once and I say where.

*Step 1 (maximum principle; uses B2 and (D1)).* By B2, `LΠ := ΔΠ − b·∇Π = |ω|^2 ≥ 0` on `R^3`. `L` is a second-order elliptic operator with smooth coefficients and **no zeroth-order term**; on each ball `B_R` its drift `b` is bounded. The weak maximum principle for subsolutions of such an operator gives `max_{\bar B_R} Π = max_{∂B_R} Π`. Letting `R → ∞` and using (D1), `Π ≤ 0` on `R^3`. (This is the step that fails for DSS profiles, B8.)

*Step 2 (cutoff integration; uses `div b = 3/2 > 0`, Step 1, the boundedness half of (D0), and (D2)).* Fix `χ ∈ C_c^∞([0,∞))`, `0 ≤ χ ≤ 1`, `χ = 1` on `[0,1]`, `χ = 0` on `[2,∞)`, `χ' ≤ 0`, and set `χ_R(y) := χ(|y|/R)`, so `R ↦ χ_R(y)` is nondecreasing for each `y`. Integrating `|ω|^2 = ΔΠ − b·∇Π` against `χ_R` and moving all derivatives onto `χ_R` (all integrals are over the compact set `{|y| ≤ 2R}`, so no boundary terms appear):

```text
∫ |ω|^2 χ_R = ∫ Π Δχ_R + ∫ Π div(b χ_R) = (3/2) ∫ Π χ_R + ∫ Π ( Δχ_R + b·∇χ_R ),
```

using `div b = div(y/2) + div U = 3/2 + 0`. The first term is `≤ 0` by Step 1. For the second, `Δχ_R` and `∇χ_R` are supported in the annulus `A_R := {R ≤ |y| ≤ 2R}` and there

```text
|Δχ_R| ≤ C_χ / R^2 ≤ C_χ,     |(y/2)·∇χ_R| = (|y|/(2R)) |χ'(|y|/R)| ≤ ‖χ'‖_∞,     |U·∇χ_R| ≤ ‖U‖_{L^∞(A_R)} / R,
```

so `|Δχ_R + b·∇χ_R| ≤ C := C_χ + ‖χ'‖_∞ (1 + ‖U‖_{L^∞})`, a constant independent of `R` by (D0). The only structural input is that the Ornstein–Uhlenbeck part of the drift contributes `(|y|/(2R))|χ'|`, which is bounded rather than growing, because `|y| ≤ 2R` on the support of `∇χ_R`. Hence

```text
0 ≤ ∫ |ω|^2 χ_R ≤ C ∫_{A_R} |Π| dy.
```

Take `R = R_j` along a sequence realizing the `liminf` in (D2). For any fixed `R_0` and all `R_j ≥ R_0`, monotonicity of `χ_R` gives `∫ |ω|^2 χ_{R_0} ≤ ∫ |ω|^2 χ_{R_j} ≤ C ∫_{A_{R_j}} |Π| → 0`, so `∫ |ω|^2 χ_{R_0} = 0`, so `ω ≡ 0` on `B_{R_0}`. As `R_0` was arbitrary, `ω ≡ 0` on `R^3`.

*Step 3 (Liouville; uses the decay half of (D0)).* With `ω ≡ 0`, (B2.1) with `Res = 0` reads `0 = −∇Π`, so `Π` is a constant `c`; then `∫_{A_R}|Π| = |c| · |A_R|` with `|A_R| → ∞`, so (D2) forces `c = 0`. Also `ΔU = −curl ω = 0`, so every component of `U` is harmonic on `R^3`; it is bounded and tends to `0` at infinity by (D0), so by the classical Liouville theorem for harmonic functions `U ≡ 0`. (Then `P` is constant, `= 0` after normalization.) ∎

**What Theorem B3 is not.** It is *weaker* than NRS and Tsai, and the difference is exactly the interesting part. My hypotheses (D1)–(D2) are decay hypotheses on `Π`; because `Π` contains the term `(1/2) y·U`, (D2) needs roughly `|U(y)| = o(|y|^{-4})`, which `U ∈ L^3` does not give directly. **KNOWN, not reproved here:** NRŠ 1996 [checked via WebSearch: Acta Math. 176 (1996) 283–294] reach `U ≡ 0` from `U ∈ L^3(R^3)`, and Tsai 1998 from the local energy bounds of a suitable weak solution. I cannot reconstruct their intermediate decay statement from memory and do not guess it.

**What this does *not* establish — CONJECTURE, downgraded from an assertion in round 1.** The first version of this note said that the content of NRŠ and Tsai "is precisely the missing lemma — that `U ∈ L^3`, or the local energy bounds of a suitable weak solution, force enough decay of `U` and `P` for the maximum-principle argument to run", and the transfer table said 2.3 was "reproved in structure". **Both go beyond what is established anywhere in this project, and the round-1 referee is right to reject them.** The implication "`U ∈ L^3` ⟹ (D0)–(D2)" is proved nowhere here; asserting that it is *the* content of NRŠ presupposes that their proof — which no one in this project has read (`00_problem.md` section 0) — factors through the hypotheses of Theorem B3. It may instead run through a Carleman estimate, a direct energy argument, or a decay estimate never phrased through `Π` at all. The honest statement, and the one now used everywhere below: **Theorem B3 reaches `U ≡ 0` from (D0)–(D2); NRŠ reach `U ≡ 0` from `U ∈ L^3` by a route that was not read; the bridge between the two hypothesis sets is missing.** Status: Theorem B3 PROVED-HERE exactly as stated; "B3 reproduces the structure of NRŠ" CONJECTURE. What is genuinely reproduced here is the *identity* B2 and the *shape* of a Liouville argument that closes under (D0)–(D2) — a self-contained object, whose relation to the literature is a guess.

**The exact ansatz does not need any of this — PROVED-HERE (given 2.3).** For the *exactly* self-similar ansatz there is a two-line route to `U ∈ L^3` that avoids both the decay bootstrap and Chae–Wolf. By B6.1, `∫_{R^3}|u(x,t)|^2 dx = s^{5β−2}∫|U|^2 dy`; at `β = 1/2` this is `s^{1/2}‖U‖_{L^2}^2`, so a solution with finite energy on `(0,T*]` has `‖U‖_{L^2}^2 ≤ E_0 (T*)^{−1/2} < ∞`, i.e. `U ∈ L^2` (and if `‖U‖_{L^2} = ∞` the solution has infinite energy at *every* time, so the alternative is not a loophole). With `U` bounded, `‖U‖_{L^3}^3 ≤ ‖U‖_{L^∞}‖U‖_{L^2}^2 < ∞`, so `U ∈ L^3` and NRŠ (2.3, KNOWN) gives `U ≡ 0` directly. **Consequence: Chae–Wolf and the matched-tail discussion are needed only for the asymptotically self-similar class, never for the exact ansatz.** Note also that the matched tail `|U| ~ |y|^{-1}` is not in `L^2` either, so it is not a possible profile for an *exact* finite-energy self-similar solution at all; it is an object of the local/asymptotic picture only.

**KNOWN strengthening, with its citation dependence stated at the point of use.** D. Chae, J. Wolf, "On the Liouville type theorems for self-similar solutions to the Navier–Stokes equations", Arch. Rational Mech. Anal. 225 (2017) 549–572 \[checked via WebSearch: authors, title, venue, volume, pages, and one abstract sentence stating that they remove a scenario of asymptotically self-similar blowup with profile in `L^{p,∞}(R^3)`, `p > 3/2`\]. This is the only exclusion available for the asymptotically self-similar matched profile, whose tail is `|U(y)| ~ |y|^{-α/β}`, i.e. `|y|^{-1} ∈ L^{3,∞} ∖ L^3` at `β = 1/2` — exactly the case NRŠ's hypothesis misses.

**Binding citation-dependence statement (`00_problem.md` section 0, rule 1).** Their *technical definition* of "asymptotically self-similar" — the topology, the region, the uniformity — and any further hypotheses of their theorem were **never read**, and the abstract sentence above does not supply them. Every claim below that excludes asymptotically self-similar `β = 1/2` therefore rests on hypotheses this project cannot check, and carries at most the status **KNOWN** on that account; it may not be tagged PROVED-HERE, and it is not. In particular the strict right endpoint of the headline window, `β < 1/2` rather than `β ≤ 1/2`, is exactly this dependence (B7.1).

Independently, ESS (2.4) kills `β = 1/2` because `‖u(t)‖_{L^3} = ‖U‖_{L^3}` is constant in `t` (table in B10), contradicting `lim_{t→T*}‖u(t)‖_{L^3} = ∞`. **The ESS route needs `U ∈ L^3` and is therefore genuinely unavailable for the matched profile**, whose `L^3` norm diverges logarithmically at `β = 1/2` — this is why Chae–Wolf, and not ESS, is the load-bearing input there.

### B4. Exponent balance and the rigidity of exact self-similarity — **PROVED-HERE** and **COMPUTED**

Substituting the general ansatz `u = s^{-α}U(y)`, `y = x/s^β`, `p = s^{-2α}P(y)` into `∂_t u + (u·∇)u − νΔu + ∇p = 0` and multiplying by `s^{α+1}` gives, term by term:

| term        | power of `s = T − t` relative to `s^{-α-1}`, after `α = 1 − β` |
| ----------- | -------------------------------------------------------------- |
| `d_t u`     | `s^(0)`                                                        |
| `(u.grad)u` | `s^(0)`                                                        |
| `grad p`    | `s^(0)`                                                        |
| `nu Lap u`  | `s^(1 - 2*beta)`                                               |

(script Part 1; the substitution and the differentiation are done by sympy on a random exact-rational profile, and the printed residual against the claimed profile form is nonzero for generic `(α,β)` and identically `0` after imposing `α = 1 − β`).

**B4.1 (nonlinear balance) — PROVED-HERE.** `∂_t u` scales as `s^{-α-1}`, `(u·∇)u` as `s^{-2α-β}`; equality of the exponents is `α + 1 = 2α + β`, i.e. `α = 1 − β`. The pressure gradient scales as `s^{-2α-β}` too, and this is forced, not chosen: by "Setting and notation" `p` is the Riesz-transform quadratic expression in `u`, so it cannot be given an independent exponent to rebalance the equation.

**B4.2 (rigidity) — PROVED-HERE.** *Hypothesis, missing from the round-1 statement and supplied here:* `U` is smooth and **bounded with `U(y) → 0` as `|y| → ∞`** (subpolynomial growth suffices). Without it the conclusion is false as stated: a harmonic polynomial profile satisfies `ΔU = 0` and is excluded by nothing in the argument. With `α = 1 − β` the transformed equation reads

```text
α U + β (y·∇)U + (U·∇)U + ∇P = ν s^{1−2β} ΔU.
```

`U` and `P` do not depend on `s`. If `β ≠ 1/2` the function `s ↦ s^{1−2β}` is nonconstant on `(0,T*)`, so the identity for all `s` forces **both** sides to vanish separately:

```text
ΔU = 0      and      α U + β (y·∇)U + (U·∇)U + ∇P = 0.
```

Every component of `U` is then harmonic on `R^3`, bounded, and tends to `0` at infinity, so the classical Liouville theorem for harmonic functions gives `U ≡ 0`. **FAILED, recorded rather than deleted:** the first version of this note offered a second route, "pairing `ΔU = 0` with `U` over `B_R` and letting `R → ∞`". That route does not work and the round-1 referee is right. Pairing gives `∫_{B_R}|∇U|^2 = ∫_{∂B_R} U · ∂_n U`, and nothing in the hypotheses controls the boundary term — boundedness of `U` alone permits `|∂_n U| ~ 1/R` on a sphere of area `~R^2`, so the right-hand side need not vanish. Only the Liouville route closes, and only under the growth hypothesis now in the statement. Conclusion: **an exactly self-similar Navier–Stokes solution with bounded profile decaying at infinity must have `β = 1/2`**; all other exponents can occur only asymptotically, i.e. `u ≈ s^{-α}U(x/s^β)` to leading order with `U` solving the *Euler* self-similar equation `α U + β(y·∇)U + (U·∇)U + ∇P = 0`. This is the integration-by-parts argument the task asks for in the `β > 1/2` regime; note it is equally decisive for `β < 1/2`, which is why the surviving route in B7 is only asymptotically self-similar.

### B5. `β > 1/2` is excluded — **PROVED-HERE**, with one **FAILED** reading recorded

Three independent exclusions, in increasing order of strength of the input used.

**B5.1 (exact ansatz) — PROVED-HERE.** By B4.2, `β > 1/2` forces `ΔU = 0` hence `U ≡ 0` for any bounded profile decaying at infinity.

**B5.2 (Leray's lower bound; KNOWN input 2.5) — PROVED-HERE given 2.5.** For a nontrivial bounded profile, `‖u(t)‖_{L^∞} = s^{-α}‖U‖_{L^∞}` with `α = 1 − β < 1/2`, so `‖u(t)‖_{L^∞} = o(s^{-1/2})`, contradicting `‖u(t)‖_{L^∞} ≥ c s^{-1/2}` (2.5). The same contradiction comes from `‖∇u(t)‖_{L^2} = s^{(3β-2)/2}‖∇U‖_{L^2}` against `‖∇u(t)‖_{L^2} ≥ c s^{-1/4}`, which needs `(3β−2)/2 ≤ −1/4`, i.e. `β ≤ 1/2` (table in B10). *The round-1 sentence "this kills asymptotically self-similar `β > 1/2` too, not just the exact ansatz" is downgraded to CONJECTURE and moved to B5.5: both norm identities used here are exact computations for the exact ansatz, globally in `y`.*

**B5.3 (CKN; KNOWN input 2.6) — PROVED-HERE given 2.6.** On the parabolic cylinder `Q_r = B_r × (T*−r^2, T*)` at the singular point, `r^{-1}∫∫_{Q_r}|∇u|^2 dx dt ~ r^{6β−3}` in the regime `β > 1/2` (B10 table, where the exponent is derived and is piecewise), which tends to `0` for `β > 1/2`; CKN `ε`-regularity then says the point is regular. *The round-1 sentence "again this covers the asymptotic case" is downgraded to CONJECTURE and moved to B5.5.*

**B5.4 (local Reynolds number) — PROVED-HERE (arithmetic).** The Reynolds number of the blowup structure is `Re(s) = (length)·(velocity)/ν = s^β · s^{-α}/ν = s^{2β−1}/ν`, which tends to `0` as `s → 0` when `β > 1/2`: the structure is in the Stokes regime, where no nonlinear amplification is available.

**B5.5 (what the exclusions of B5.2 and B5.3 do and do not reach) — CONJECTURE.** In round 1, B5.2 and B5.3 each ended by extending its exclusion from the exact ansatz to "asymptotically self-similar" solutions, and B7.1's whole conclusion is about a class that by B4.2 cannot be exactly self-similar. **The term was never defined in this note — no topology, no region, no uniformity — and the round-1 referee is right that the extensions are therefore unsupported.** Stating the gap precisely rather than papering over it:

- The identities actually used are `‖u(t)‖_{L^∞} = s^{-α}‖U‖_{L^∞}`, `‖∇u(t)‖_{L^2} = s^{(3β−2)/2}‖∇U‖_{L^2}`, `‖u(t)‖_{L^3} = s^{2β−1}‖U‖_{L^3}`, `∫|u(t)|^2 = s^{(d+2)β−2}∫|U|^2`, and the CKN quantity of B10. Each is an **exact change of variables for the exact ansatz, valid globally in `y`** — none survives the replacement of `u` by `s^{-α}U(x/s^β) + w` for an unspecified `w`.
- The weakest hypothesis under which they do survive is: `w` is subordinate to the leading term **in precisely the norm being used**, uniformly as `s → 0`. For `‖·‖_{L^∞}` and `‖∇·‖_{L^2}` this is an assumption of the same strength as the conclusion it is used to prove, which is why it cannot be assumed silently.
- A definition of the class *would* have to be supplied by whoever completes this: the natural candidates are (a) `s^{α}u(s^β y, T*−s) → U(y)` in `C^k_loc(R^3 ∖ {0})`, (b) the same convergence in `L^q(B(0,r))` on balls shrinking at the blowup rate, (c) convergence of the *vorticity* in a critical Besov norm. Options (b) and (c) are, at abstract level, the two notions used by Chae in arXiv:math/0604234 (B9), and (b) there shrinks balls at the rate `√(T−t)`, i.e. at `β = 1/2` only.

Until a definition is fixed **and** the subordination is proved, every sentence in this note that extends an exclusion from the exact ansatz to the asymptotic class is CONJECTURE. This applies to the closing sentences of B5.2 and B5.3, to B7.1 as a whole, and to the asymptotic half of B10. It does **not** apply to B4.2, B5.1, B6.1, B6.2a or the exact-ansatz half of B10, which are statements about the exact ansatz and are PROVED-HERE.

**FAILED sub-item, recorded as required.** The task asked to exclude `β > 1/2` "by an energy/integration-by-parts argument on the profile equation". Pairing the all-terms profile equation `νΔU = αU + β(y·∇)U + (U·∇)U + ∇P` with `U` over `R^3` (assuming `U ∈ H^1`, `U` and `P` decaying fast enough for the boundary terms to vanish) and using `∫((y·∇)U)·U = −(d/2)∫|U|^2`, `∫((U·∇)U)·U = 0`, `∫∇P·U = 0` gives in dimension `d`

```text
ν ∫ |∇U|^2 dy = ( ((d+2)β − 2)/2 ) ∫ |U|^2 dy,        d = 3:   ν ∫|∇U|^2 = ((5β−2)/2) ∫|U|^2.
```

This is a genuine identity (it is the energy identity of the original solution rewritten in profile variables — see B6) and it **does** exclude `β < 2/5` in `d = 3`. But its right-hand side is *positive* for `β > 1/2`, so it excludes nothing there: the requested argument does not close in the requested regime, and I record that as a failure of the suggested route rather than dressing up B5.1–B5.4 as it. The correct integration-by-parts statement for `β > 1/2` is the degenerate one of B5.1 (`ν∫|∇U|^2 = 0` from `ΔU = 0`).

### B6. `β < 2/5` is excluded by finite energy; `β = 2/5` is the borderline — **PROVED-HERE** / **CONJECTURE**

**B6.1 — PROVED-HERE.** In dimension `d`, `∫_{R^d}|u(x,t)|^2 dx = s^{-2α + dβ} ∫|U|^2 dy = s^{(d+2)β − 2}∫|U|^2 dy` (substitute `x = s^β y`, `dx = s^{dβ}dy`, and `α = 1 − β`). In `d = 3` the exponent is `5β − 2`. So along a self-similar blowup the energy contained in the profile region behaves like `s^{5β−2}`:

- `β > 2/5`: the profile-region energy tends to `0` — the singularity is "energy-shedding", consistent with a smooth global energy bound;
- `β = 2/5`: it is exactly constant (the self-similar rescaling is energy-preserving; equivalently the energy is *critical* for this scaling);
- `β < 2/5`: it tends to `+∞`, contradicting `∫|u(t)|^2 dx ≤ ∫|u^0|^2 dx` for a smooth solution on `[0,T*)`.

Hence **`β ≥ 2/5` is forced by finite energy** in `d = 3`. For the **exact** ansatz the argument is exhaustive with no case analysis at all: `∫_{R^3}|u(x,t)|^2 dx = s^{5β−2}∫|U|^2 dy` holds identically, so either `∫|U|^2 = ∞`, in which case the solution has infinite energy at *every* `t` and is not in the Clay class, or `∫|U|^2 < ∞`, in which case `β < 2/5` sends the energy to `+∞` as `s → 0`. **So for the exact ansatz, finite energy gives both `U ∈ L^2` and `β ≥ 2/5`, and no statement about tails is needed.** (This is the observation used in B3 and B10 to avoid Chae–Wolf and the matched-tail dichotomy entirely.) For the *local or asymptotic* picture the profile need not be globally defined and the tail matters; there, if `U` has the generic matched tail `|U(y)| ~ |y|^{-α/β}` (so that `u(x,T*) ~ |x|^{-α/β}` is `s`-independent, which is what "matching to the outer solution" means), then `∫_{|x| ≤ r_0}|u|^2 dx = s^{(d+2)β−2}∫_{|y| ≤ r_0 s^{-β}}|U|^2 dy ~ s^{(d+2)β−2} (s^{-β})^{d − 2α/β}`, and the two powers cancel identically, giving a finite `s`-independent answer precisely when `d − 2α/β > 0`, i.e. `2(1−β)/β < 3`, i.e. `β > 2/5` again. The threshold is the same either way, which is the arithmetic reason `2/5` is not an artifact of assuming `U ∈ L^2`.

**B6.2 (the endpoint `β = 2/5`), split by class.** In round 1 this was a single CONJECTURE presented through a dichotomy "either `U ∈ L^2` or the exact matched power tail", which the referee correctly observed is exhaustive only among power-law tails. The split below is by *class*, which is exhaustive.

**B6.2a (local/asymptotic picture, matched power tail) — PROVED-HERE.** Suppose the leading profile has the matched power tail `|U(y)| ~ |y|^{-α/β} = |y|^{-3/2}` at `β = 2/5`, and that the leading-term identity `∫_{|x| ≤ r_0}|u|^2 dx = s^{5β−2}∫_{|y| ≤ r_0 s^{-β}}|U|^2 dy` holds. Then `∫_{|y|<Y}|y|^{-3} dy = 4π log Y + O(1)`, so with `5β − 2 = 0` the local energy is `4πβ log(1/s) + O(1) → ∞`, contradicting `∫_{R^3}|u(t)|^2 ≤ ∫|u^0|^2`. The computation is elementary and complete: `|x|^{-3/2}` is not locally square integrable in three dimensions. **The hypothesis that the leading-term identity survives for the asymptotic class is exactly the subordination assumption of B5.5, so as a statement about the asymptotic class this inherits CONJECTURE; as a statement about a profile with that exact tail it is PROVED-HERE.**

**B6.2b (exact ansatz, `U ∈ L^2`, exactly `β = 2/5`) — CONJECTURE.** Here a fixed positive amount of energy sits in a ball of radius `s^{2/5} → 0`, i.e. `|u|^2` concentrates an atom at the singular point. Nothing above excludes it. Doing so needs a non-concentration statement for suitable weak solutions at a first singular time, which I have not proved and do not find in the no-go list 2.1–2.15. I keep the closed window `2/5 ≤ β < 1/2` below rather than claiming the open one.

### B7. The surviving window `2/5 ≤ β < 1/2`: Euler leading profile and the matching arithmetic — **PROVED-HERE** (arithmetic) and **CONJECTURE** (perturbativity)

**B7.1 (what survives) — KNOWN, downgraded from PROVED-HERE in round 1.** *Why the downgrade.* The strictness of the right endpoint is what turns this into "Type II, strictly faster than self-similar", and for the *exact* ansatz that endpoint is fine (B3, B4.2, B6.1 close every `β` including `1/2`, PROVED-HERE modulo 2.3). But B4.2 also says the surviving class cannot be exactly self-similar, so B7.1 is a statement about the asymptotically self-similar class — and there the notes themselves identify the physically relevant profile as the matched tail `|U(y)| ~ |y|^{-1}` at `β = 1/2`, which is **not** in `L^3` (the integral of `|y|^{-3}` diverges logarithmically), so neither ESS (2.4) nor NRŠ (2.3) reaches it. The sole exclusion is Chae–Wolf 2017, whose technical definition of "asymptotically self-similar" and further hypotheses were never read (B3, binding citation-dependence statement). By `00_problem.md` section 0 rule 1 this claim may not be PROVED-HERE. Its status is **KNOWN**, and the endpoint is honestly recorded as `β ≤ 1/2` with `β = 1/2` excluded *only modulo unread hypotheses*. Everything else in B7.1 — the arithmetic of `α = 1 − β`, `β ≥ 2/5`, `β ≤ 1/2` — remains PROVED-HERE for the exact ansatz and CONJECTURE for the asymptotic class per B5.5. Combining B4.1 (`α = 1 − β`), B5 (`β ≤ 1/2`), B6 (`β ≥ 2/5`), and B3/2.3/2.4/Chae–Wolf (`β = 1/2` excluded), an asymptotically self-similar blowup of 3D Navier–Stokes in the Clay class would have

```text
2/5 ≤ β < 1/2,     α = 1 − β ∈ (1/2, 3/5],     ‖u(t)‖_∞ ~ (T*−t)^{-α}   with   α > 1/2,
```

that is, **Type II blowup, strictly faster than self-similar**, and by B4.2 it cannot be exactly self-similar: the leading profile must solve the *Euler* self-similar equation

```text
α U + β (y·∇)U + (U·∇)U + ∇P = 0,      div U = 0,      α = 1 − β.
```

**B7.2 (inner/outer matching arithmetic) — PROVED-HERE (formal computation, exact arithmetic).** The viscous length at time-to-blowup `s` is `ℓ_ν = (νs)^{1/2}`; the blowup structure has size `ℓ_b = s^β`. For `β < 1/2`, `ℓ_b / ℓ_ν = ν^{-1/2} s^{β−1/2} → ∞`: the structure is much *larger* than the distance diffusion can travel in the remaining time, so viscosity cannot smear it. Equivalently `Re(s) = s^{2β−1}/ν → ∞` (B5.4). In the inner region, the viscous term is smaller than the inertial term by the factor `ν s^{1−2β}` (B4 table). In the outer region the flow is `s`-independent, `u(x,T*) ~ |x|^{-α/β}`, and

```text
(viscous)/(inertial)  ~  ν |x|^{-α/β − 2} / |x|^{-2α/β − 1}  =  ν |x|^{α/β − 1},
```

which tends to `0` as `|x| → 0` exactly when `α/β > 1`, i.e. `β < 1/2` — the same condition. Matching the two: the inner estimate `ν s^{1−2β}|y|^{α/β−1} ≪ 1` holds for `|y| ≪ s^{-β}`, i.e. for `|x| ≪ 1`, so the inviscid approximation is valid across the whole blowup region up to the outer scale, with no intermediate viscous layer. Viscosity is **formally perturbative throughout the window**, and it is not perturbative at `β = 1/2` (where the ratio is `O(1)`) or `β > 1/2` (where it dominates).

**B7.3 (why "formally" is doing real work) — CONJECTURE.** In similarity variables `τ = −log s` the equation for `V(y,τ)` reads `∂_τ V + αV + β(y·∇)V + (V·∇)V + ∇P = ν e^{-(1−2β)τ} ΔV`. The small parameter `ν e^{-(1−2β)τ}` multiplies the *highest-order* operator: this is a singular perturbation, the limit `τ → ∞` is not uniform in `y` (it fails at high wavenumbers, where the viscous term is `ν e^{-(1−2β)τ}|ξ|^2` and beats the `O(1)` transport terms as soon as `|ξ| ≳ e^{(1−2β)τ/2}`), and no formal-order argument can conclude that a viscous solution shadows an Euler self-similar profile. That a smooth Euler self-similar profile with `β` in the window is actually shadowed by a Navier–Stokes solution is exactly the CONJECTURE this line cannot decide; see Breakpoint.

### B8. Discretely self-similar profiles at `β = 1/2`: the identity and why the maximum principle dies — **PROVED-HERE** and **COMPUTED**

**Scope, stated first because round 1 framed it too widely.** Proposition B8 fixes `β = 1/2` and `ν = 1`. It does **not** cover the `β < 1/2` profiles of B7, and the reason is structural rather than technical — see "B8 is `β = 1/2`-specific" below. Its subject is *backward* DSS blowup.

In similarity variables at `β = 1/2`, with `V = V(y,τ)`, `τ = −log s`, the equation is

```text
∂_τ V + (1/2)V + (1/2)(y·∇)V + (V·∇)V + ∇P = ΔV,      div V = 0,
```

and exact self-similarity is `∂_τ V = 0`, DSS with factor `λ` is `2 log λ`-periodicity in `τ`. Redoing the proof of B2 with the extra term (the only change is that `Res` acquires `−∂_τ V` and that `div ∂_τ V = ∂_τ div V = 0`, so the divergence step is unaffected):

**Proposition B8.** With `Π = P + |V|^2/2 + (1/2)y·V`, `b = y/2 + V`, `ω = curl V`, and `Res := ΔV − [∂_τ V + (1/2)V + (1/2)(y·∇)V + (V·∇)V + ∇P]`, identically in `(V,P)`,

```text
ΔΠ − b·∇Π − ∂_τ Π + ∂_τ P − |ω|^2 = − div(Res) + b·Res.
```

Equivalently, on a solution (`Res = 0`),

```text
( Δ − b·∇ − ∂_τ ) Π = |curl V|^2 − ∂_τ P.
```

Proof. As in B2, `Res = −curl ω − ∂_τ V − ∇Π + b×ω`; taking the divergence and eliminating `b·curl ω` exactly as there gives `ΔΠ − b·∇Π − |ω|^2 − b·∂_τ V = −div Res + b·Res`. Finally `∂_τ Π = ∂_τ P + V·∂_τ V + (y/2)·∂_τ V = ∂_τ P + b·∂_τ V`, so `b·∂_τ V = ∂_τ Π − ∂_τ P`. ∎

**COMPUTED** (script Part 3, generic `V = curl A` with `A_i` and `P` arbitrary functions of `(y_1,y_2,y_3,τ)`):

```text
- residual of the tau-dependent identity: 0
- identity verified: True
- residual after deleting the obstruction term d_tau P is zero: False
```

**Why the maximum-principle argument does not extend — PROVED-HERE (as an obstruction statement).** In the exactly self-similar case the right-hand side is `|ω|^2 ≥ 0`, so `Π` is a subsolution of an operator with no zeroth-order term and Step 1 of B3 applies. In the DSS/asymptotically self-similar case the right-hand side is `|ω|^2 − ∂_τ P`, and `∂_τ P` has **no sign**: `P` is the (Riesz-transform) pressure of `V`, so `∂_τ P = Σ R_iR_j ∂_τ(V_iV_j)` is an unsigned nonlocal quadratic expression in `V` and `∂_τ V`, and there is no inequality of the form `∂_τ P ≤ c|ω|^2` available. Averaging in `τ` over one DSS period does not rescue the argument either: both `∂_τ Π` and `∂_τ P` are `τ`-derivatives of periodic functions and average to zero, so `Π̄ := ⟨Π⟩_τ` satisfies `ΔΠ̄ − (y/2)·∇Π̄ = ⟨|ω|^2⟩ + ⟨V·∇Π⟩`, and the obstruction has simply moved into the correlation term `⟨V·∇Π⟩ ≠ V̄·∇Π̄`, which is again unsigned. So the DSS case is not merely technically harder: the specific structural fact that makes NRŠ work — that the source term is a perfect square — is destroyed.

**Correction to the round-1 closing sentence.** That version ended "and this is exactly why DSS solutions exist (Chae–Wolf; Bradshaw–Tsai, KNOWN, B9) while exactly self-similar ones do not". **That comparison is false and is withdrawn.** The cited DSS solutions are *forward* solutions from singular data (B9), and forward *exactly* self-similar solutions also exist (Jia–Šverák, cited in 2.3), so the contrast does not distinguish DSS from exact self-similarity at all; and it says nothing about *backward* DSS, which is B8's actual subject. The correct statement, and all that is claimed: **no result in the no-go list excludes backward DSS blowup, and B8 identifies the precise structural reason the NRŠ mechanism does not reach it** — the source term `|ω|^2 − ∂_τ P` loses its sign. Whether backward DSS blowup exists is open and untouched here.

**B8 is `β = 1/2`-specific — PROVED-HERE (as an obstruction).** For general `β` the similarity equation is `∂_τ V + αV + β(y·∇)V + (V·∇)V + ∇P = ν e^{-(1−2β)τ}ΔV` with `α = 1 − β` (B7.3). The head-pressure construction rests on the identity `∇(y·V) = (y·∇)V + V + y × ω`, which gives

```text
α V + β (y·∇)V = β ∇(y·V) + (α − β) V − β y × ω.
```

At `β = 1/2` one has `α = β` and the leftover `(α − β)V` vanishes, which is exactly why `Π = P + |V|^2/2 + (1/2)y·V` absorbs both drift terms. For `β ≠ 1/2` the leftover `(α − β)V = (1 − 2β)V` is not a gradient and is not a cross product with `ω`, so **this** head pressure leaves it unabsorbed and B2/B8 do not apply as written. (This is an obstruction to the construction, not a proof that no head pressure exists.) Independently, the viscous coefficient `ν e^{-(1−2β)τ}` multiplies `ΔV` but not the other terms of `Π`, so even the elimination of `curl ω` no longer produces a `τ`-independent drift-diffusion operator. Consequence: the `β < 1/2` window of B7 is **not** covered by B8, and no obstruction statement is claimed there.

### B9. Known obstructions to smooth decaying Euler self-similar profiles — **KNOWN**

The surviving route of B7 needs a nontrivial solution of the Euler self-similar equation with `α = 1 − β`, `β ∈ [2/5, 1/2)`. Known constraints, none of which closes the window:

- **Chae's nonexistence theorem** \[checked via WebSearch: the statement "there exists no finite-time blowing-up self-similar solution to the 3D Euler equations if there exists `p_1 > 0` such that the vorticity `Ω = curl V ∈ L^p(R^3)` for all `p ∈ (0,p_1)`" appears verbatim in a search snippet; the paper is D. Chae, "Nonexistence of self-similar singularities for the 3D incompressible Euler equations", Comm. Math. Phys. 273 (2007), arXiv:math/0601060 — venue confirmed, page numbers unverified\]. **Why it does not close the window:** the matched vorticity tail of a profile with exponent `β` is `|Ω(y)| ~ |y|^{-(α+β)/β} = |y|^{-1/β}` (using `ω_u = s^{-1}Ω(y)` and `s`-independence of the outer vorticity), so `Ω ∈ L^p(R^3)` requires `p/β > 3`, i.e. `p > 3β ≥ 6/5`. There is no `p_1 > 0` with `Ω ∈ L^p` for all small `p`; Chae's hypothesis fails for exactly the profiles the window allows. This is the arithmetic reason the surviving route is not already excluded — and it also says what a candidate profile must look like: slowly decaying vorticity, `Ω ∉ L^p` for `p ≤ 6/5`.
- **Chae, "Nonexistence of asymptotically self-similar singularities in the Euler and the Navier–Stokes equations"** [checked via WebSearch: title and arXiv id math/0604234; I cannot recall the hypotheses and do not guess them]. This is the item most likely to constrain B7 directly; see "What would be needed next".
- **Elgindi's `C^{1,α}` Euler blowup** (2.13) \[checked via WebSearch: the axisymmetric swirl-free `C^{1,α}` singularity with a locally self-similar profile is confirmed; the numerical value of its `β` is **[unverified, not recalled]** — I do not reconstruct it, an earlier attempt at reconstructing it from the Hölder exponent produced an obviously wrong sign\]. Two things are certain without it: the data are not `C^∞`, hence outside the Clay class (4); and whatever `β` it has, `Re(s) ~ s^{2β−1}/ν` decides whether viscosity is perturbative for it, so the scaling check demanded by 2.13/T4 is the one-line computation of B5.4.
- **Chen–Hou** (2.13): smooth data but with a boundary; the mechanism is not available on `R^3` or `T^3`. **PINN-discovered unstable profiles** (2.14): none of the equations treated is 3D Navier–Stokes, and by 2.3/B3 an exactly self-similar Navier–Stokes profile in `L^3` cannot be found by any method, so such searches must target the `β < 1/2` asymptotically self-similar window of B7 or a DSS ansatz.
- **DSS solutions do exist** \[checked via WebSearch: Chae–Wolf construct DSS solutions for any DSS data in `L^2_loc`; Bradshaw–Tsai, "Discretely self-similar solutions to the Navier–Stokes equations with data in `L^2_loc` satisfying the local energy inequality", Analysis & PDE 12 (2019) 1943–…, title and venue confirmed\]. These are *forward* DSS solutions from singular (`(−1)`-homogeneous-type) data, not backward blowup, and their data are not in the Clay class; they are the standing example that the DSS ansatz is not empty, which is why B8's obstruction matters.

### B10. Transfer: hyperdissipative window, the Lions exponent, and the 2D control — **PROVED-HERE** and **COMPUTED**

Replace `νΔ` by `−ν(−Δ)^σ` in dimension `d`. The nonlinear balance is unchanged, `α = 1 − β` (nothing in B4.1 used the dissipation). Two exponents govern the ansatz:

```text
β_ν := 1/(2σ)     (dissipation balances transport: the term −ν(−Δ)^σ u carries s^{−α−2σβ},
                   equal to the s^{−α−1} of ∂_t u exactly when 2σβ = 1);
β_E := 2/(d+2)    (the profile-region energy s^{(d+2)β−2} is s-independent: B6.1 in d dimensions).
```

`β ≤ β_ν` is forced by the Leray/CKN arguments of B5 (which in the hyperdissipative case rest on the scaling `|u| ~ ℓ^{-(2σ−1)}`, `s ~ ℓ^{2σ}`, hence `‖u(t)‖_∞ ≳ s^{-(2σ−1)/(2σ)} = s^{-(1−β_ν)}`; this lower bound is the scaling form of local well-posedness in `L^∞`, KNOWN [unverified as a reference, standard]). `β ≥ β_E` is forced by finite energy (B6.1, whose proof is dissipation-independent). So the self-similar window is `[β_E, β_ν]`, and it is nonempty iff

```text
2/(d+2) ≤ 1/(2σ)     <=>     σ ≤ (d+2)/4.
```

**COMPUTED** (script Part 4), with the `sympy` solve printing `β_E = 2/(d+2)`, the profile identity constant `((d+2)β − 2)/2`, and `β_E = β_ν <=> σ = d/4 + 1/2`:

| d   | sigma | beta_E = 2/(d+2) | beta_nu = 1/(2 sigma) | window     | verdict                        |
| --- | ----- | ---------------- | --------------------- | ---------- | ------------------------------ |
| 2   | 1     | 1/2              | 1/2                   | [1/2, 1/2] | single point (energy-critical) |
| 2   | 9/8   | 1/2              | 4/9                   | [1/2, 4/9] | EMPTY                          |
| 2   | 5/4   | 1/2              | 2/5                   | [1/2, 2/5] | EMPTY                          |
| 2   | 3/2   | 1/2              | 1/3                   | [1/2, 1/3] | EMPTY                          |
| 3   | 1     | 2/5              | 1/2                   | [2/5, 1/2] | nonempty                       |
| 3   | 9/8   | 2/5              | 4/9                   | [2/5, 4/9] | nonempty                       |
| 3   | 5/4   | 2/5              | 2/5                   | [2/5, 2/5] | single point (energy-critical) |
| 3   | 3/2   | 2/5              | 1/3                   | [2/5, 1/3] | EMPTY                          |
| 4   | 1     | 1/3              | 1/2                   | [1/3, 1/2] | nonempty                       |
| 4   | 9/8   | 1/3              | 4/9                   | [1/3, 4/9] | nonempty                       |
| 4   | 5/4   | 1/3              | 2/5                   | [1/3, 2/5] | nonempty                       |
| 4   | 3/2   | 1/3              | 1/3                   | [1/3, 1/3] | single point (energy-critical) |

**Required hyperdissipative check (2.10, T2) — PROVED-HERE.** At `d = 3`, `σ ≥ 5/4` the window degenerates to the single energy-critical point `β = 2/5` (at `σ = 5/4`) and is empty for `σ > 5/4`. At `σ = 5/4` exactly, the profile `L^2` identity of B5 reads `ν‖(−Δ)^{5/8}U‖_2^2 = ((5·(2/5) − 2)/2)‖U‖_2^2 = 0`, so `U ≡ 0` whenever `U ∈ L^2` (the identity needs its integrals to converge); the complementary case of the matched tail `|U| ~ |y|^{-3/2} ∉ L^2` is closed instead by the logarithmic divergence of B6.2. Either way **the self-similar route is closed at `σ ≥ 5/4` with no room at the endpoint**, and the threshold `σ = (d+2)/4` recovered here is exactly the Lions exponent of 2.10. This is not a new proof of Lions' theorem (it only kills the self-similar ansatz, not all blowup), but it is a nontrivial consistency check: the self-similar window and the energy-criticality threshold of 1.6 coincide identically, for every `d` and `σ`.

**Required 2D check (2.8, T2) — PROVED-HERE.** At `d = 2`, `σ = 1` the window collapses to the single point `β = 1/2`, and there `β_E = β_ν = 1/2` means the energy is *critical* for the blowup scaling in 2D (`∫_{R^2}|u|^2 dx = s^{4β−2}‖U‖^2 = s^0`), matching 1.6 and 2.8. So the argument correctly reports 2D as the borderline case and does not exclude 2D blowup by a mechanism that would also apply in 3D: the exclusion is **dimensional** (`β_E = 2/(d+2)`, from the volume factor `s^{dβ}` in the energy), not a statement about vortex stretching, and it degenerates to an equality in 2D exactly where the known 2D regularity proof also becomes an equality (critical energy plus a subcritical enstrophy bound). The 3D window `[2/5, 1/2)` is open precisely because `2/5 < 1/2`, i.e. because `d = 3 > 2`. Note the argument does **not** claim to prove 2D regularity — at `β = 1/2` in 2D it says only that a self-similar blowup would have exactly scale-invariant energy, which the maximum-principle bound on `‖ω‖_{L^∞}` (2.8) then rules out by a different route.

**Rate table for the window (T3 arithmetic) — PROVED-HERE and COMPUTED** (script Part 4). All rates for `u = s^{-α}U(x/s^β)` in `d = 3`. In the last column `Q_r = B_r × (T*−r^2, T*)` and the entry is `r^{-1} ∫∫_{Q_r} |∇u|^2 dx dt`; its exponent is **piecewise**, `2 − 1/β` for `1/3 < β ≤ 1/2` and `6β − 3` for `β ≥ 1/2` (the two agree at `β = 1/2`). The crossover is real: substituting `s = r^{1/β}σ` in `∫_0^{r^2} s^{3β−2} F(r s^{-β}) ds`, `F(R) := ∫_{B_R}|∇U|^2`, the upper limit `r^{2−1/β}` tends to `∞` for `β < 1/2` (the integral saturates, giving `r^{3−1/β}` before dividing by `r`) and to `0` for `β > 1/2` (the profile core fits inside `B_r` throughout, giving `r^{6β−2}`).

| beta | alpha = 1 - beta | sup-norm of u | L^3 norm of u | sup-norm of omega | L^2 norm of grad u | CKN quantity |
| ---- | ---------------- | ------------- | ------------- | ----------------- | ------------------ | ------------ |
| 7/20 | 13/20            | s^(-13/20)    | s^(-3/10)     | s^(-1)            | s^(-19/40)         | r^(-6/7)     |
| 2/5  | 3/5              | s^(-3/5)      | s^(-1/5)      | s^(-1)            | s^(-2/5)           | r^(-1/2)     |
| 9/20 | 11/20            | s^(-11/20)    | s^(-1/10)     | s^(-1)            | s^(-13/40)         | r^(-2/9)     |
| 1/2  | 1/2              | s^(-1/2)      | s^(0)         | s^(-1)            | s^(-1/4)           | r^(0)        |
| 3/5  | 2/5              | s^(-2/5)      | s^(1/5)       | s^(-1)            | s^(-1/10)          | r^(3/5)      |

Two structural readings. (i) `‖ω(t)‖_{L^∞} ~ s^{-α−β} = s^{-1}` for **every** `β`, since `α + β = 1`: a self-similar blowup always has the BKM-critical vorticity rate `(T*−t)^{-1}`, with `∫^{T*}‖ω‖_∞ dt` diverging logarithmically — the ansatz is exactly BKM-borderline, which is why 2.7 gives no information about `β`. (ii) `‖u(t)‖_{L^3} ~ s^{2β−1}` blows up iff `β < 1/2` and is constant at `β = 1/2`: the ESS criterion (2.4) is *precisely* the statement `β ≠ 1/2`, so ESS and NRS attack the same single point of the `β`-line from different directions.

## Breakpoint

This line reduces "self-similar blowup of 3D Navier–Stokes" to a one-parameter family indexed by `β`, closes every value of `β` outside `[2/5, 1/2)`, and stops at three specific places inside it.

1. **The decay bootstrap in B3 is cited, not proved.** Theorem B3 is complete under (D0)–(D2), and (D0)–(D2) are hypotheses about `Π = P + |U|^2/2 + (y·U)/2`, which contains `y·U` and therefore needs roughly `|U(y)| = o(|y|^{-4})` for (D2). The step "`U ∈ L^3` (or a local-energy bound) implies that much decay" is the actual content of NRS 1996 and Tsai 1998 and is KNOWN here. So the exclusion of `β = 1/2` in this note is: my proof, their lemma. I flag this rather than claiming the whole theorem.
2. **No smooth decaying Euler self-similar profile with `β ∈ [2/5, 1/2)` is known to exist, and this note produces none.** The Euler profile equation `(1−β)U + β(y·∇)U + (U·∇)U + ∇P = 0` on `R^3` with the matched tail `|U| ~ |y|^{-(1−β)/β}` (so `|U|` decays like `|y|^{-p}` with `p ∈ (1, 3/2]`) is an overdetermined-looking elliptic system with no variational structure I can exploit and no small parameter. Chae's nonexistence theorem does not apply (B9: the vorticity tail `|y|^{-1/β}` is in `L^p` only for `p > 3β ≥ 6/5`, and his hypothesis needs `L^p` for arbitrarily small `p`), which is why the window is open, but that is an absence of an obstruction, not the presence of a profile. **This is the exact step where the line stops.**
3. **Even given such a profile, the passage from Euler to Navier–Stokes is a singular perturbation, not a perturbation.** B7.3: the viscous coefficient in similarity variables is `ν e^{-(1−2β)τ}` on the highest-order operator, so the formal smallness `s^{1−2β}` fails uniformly above wavenumber `e^{(1−2β)τ/2}`, exactly the range where a blowup lives. Closing this gap requires the stability programme demanded by 2.14: a linearized operator around the profile, an unstable-manifold count, and a rigorous (probably computer-assisted) coercivity estimate — none of which exists here.

The honest summary: this line does not produce a blowup and does not produce a regularity theorem. What it produces is the *window* `2/5 ≤ β < 1/2` with a proof for each of its two closing inequalities, the identity and Liouville argument that close the right endpoint, the obstruction term `∂_τ P` that reopens the problem for DSS, and the transfer formula `σ ≤ (d+2)/4` that recovers the Lions exponent from the same arithmetic.

## Transfer check

Table required by `00_problem.md` §4, one row per no-go item and per transfer test. "This line" means the union of B1–B10.

| item                                      | verdict                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 2.1 supercriticality (Tao 2007)           | consistent, and the window is a quantitative expression of it: `β_E = 2/(d+2)` is where the *only* globally controlled coercive quantity (the energy) becomes scale-invariant for the blowup scaling, and it lies strictly below the viscous threshold `β_ν = 1/2` exactly because the energy is supercritical in 3D (1.6). No new controlled quantity is claimed.                                                                                                                                                                                           |
| 2.2 averaged Navier–Stokes (Tao 2016)     | B4, B5.4, B6, B10 transfer verbatim to `∂_t u = Δu + B̃(u,u)` (they use only the energy identity, scaling and the linear term), and this is consistent: Tao's blowup is a cascade over infinitely many scales, i.e. Type II with no fixed `β`, which is outside the ansatz. B2/B3/B8 do **not** transfer: they use `ΔU = −curl ω`, the Lamb identity and the pressure formula `p = R_iR_j(u_iu_j)` — exactly the exact-nonlinearity structure T1 demands be named.                                                                                            |
| 2.3 NRŠ/Tsai self-similar                 | reproved in structure (B2, B3) with the decay lemma cited; strengthened by the Chae–Wolf `L^{p,∞}`, `p > 3/2` result (B3), which covers the matched tail \`                                                                                                                                                                                                                                                                                                                                                                                                  |
| 2.4 ESS / Seregin / Tao 2019              | passes and is sharpened: `‖u(t)‖_{L^3} ~ s^{2β−1}`, so ESS is *exactly* the statement `β ≠ 1/2` for this ansatz, and it is satisfied with room to spare (`s^{-1/10}` at `β = 9/20`) throughout the surviving window. The triple-log lower bound of Tao 2019 is comfortably beaten by a power.                                                                                                                                                                                                                                                                |
| 2.5 LPS / Leray lower bounds              | passes; used as a KNOWN input in B5.2 to close `β > 1/2`. Inside the window `‖u‖_∞ ~ s^{-(1−β)}` with `1 − β > 1/2`, strictly above Leray's bound. No LPS-type a priori bound is proved anywhere (T6 not triggered).                                                                                                                                                                                                                                                                                                                                         |
| 2.6 CKN partial regularity                | passes; used in B5.3. Inside the window \`r^{-1}∫∫\_{Q_r}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| 2.7 BKM                                   | passes, and gives no information: `‖ω(t)‖_∞ ~ s^{-1}` for every `β` (because `α + β = 1`), so `∫^{T*}‖ω‖_∞` diverges logarithmically for the whole family. The ansatz is BKM-borderline by construction.                                                                                                                                                                                                                                                                                                                                                     |
| 2.8 2D global regularity                  | positive control passes (B10): at `d = 2, σ = 1` the window collapses to the single energy-critical point `β = 1/2`; the exclusion mechanism is the dimension-dependent volume factor in the energy, not vortex stretching, and it correctly degenerates rather than "proving" 2D regularity.                                                                                                                                                                                                                                                                |
| 2.9 axisymmetric swirl-free               | not applicable: this line does not use any symmetry reduction, so it neither uses nor is contradicted by the swirl-free result. Note the window is compatible with the KNOWN "no Type I blowup for axisymmetric with swirl", which is the statement `β ≠ 1/2` again.                                                                                                                                                                                                                                                                                         |
| 2.10 hyperdissipation (Lions)             | positive control passes and is the sharpest transfer here (B10): the window `[2/(d+2), 1/(2σ)]` is empty exactly for `σ > (d+2)/4`, recovering the Lions exponent `5/4` at `d = 3`, with the endpoint `σ = 5/4` also closed by `ν‖(−Δ)^{5/8}U‖^2 = 0`.                                                                                                                                                                                                                                                                                                       |
| 2.11 non-uniqueness (BV, ABC, Hou et al.) | not applicable to the arguments (no uniqueness claim), but methodologically relevant: ABC and Hou et al. build non-uniqueness from *unstable self-similar profiles* in similarity variables — the same object as B7, with data outside the Clay class. This is the closest existing technology to the surviving route.                                                                                                                                                                                                                                       |
| 2.12 Cheskidov–Dai–Palasek                | not applicable: their solutions are not smooth at `T*` and not Leray–Hopf; the ansatz here is for a classical solution smooth on `[0,T*)`.                                                                                                                                                                                                                                                                                                                                                                                                                   |
| 2.13 Euler blowup (Elgindi, Chen–Hou)     | passes T4 as a mandatory computation rather than a slogan: for any candidate Euler profile with exponent `β`, `Re(s) = s^{2β−1}/ν`, so viscosity is perturbative iff `β < 1/2` (B5.4, B7.2). Elgindi's `β` is [unverified, not recalled] and his data are `C^{1,α}` (outside the Clay class); Chen–Hou needs a boundary. No claim here is uniform in `ν → 0`.                                                                                                                                                                                                |
| 2.14 PINN singularities                   | consistent, and this line supplies the target: by B3/2.3 no exactly self-similar `L^3` Navier–Stokes profile exists, so a PINN-style search must aim at the `β ∈ [2/5, 1/2)` Euler self-similar profile of B7 (with its specific tail exponent `(1−β)/β ∈ (1, 3/2]`) or at a DSS ansatz, and must come with the unstable-manifold count 2.14 demands.                                                                                                                                                                                                        |
| 2.15 standing constraints                 | consistent; no critical-space a priori bound is produced. Small-data theory is untouched (a blowup profile is large by construction).                                                                                                                                                                                                                                                                                                                                                                                                                        |
| T1 (averaging)                            | passes: this line contains no regularity argument for (A)/(B). The parts that would transfer to `B̃` are the negative/arithmetic parts; the parts that would not (B2, B3, B8) use `curl curl`, the Lamb identity and the Riesz pressure, named explicitly above.                                                                                                                                                                                                                                                                                              |
| T2 (regular comparison systems)           | passes: applied to 2D (window collapses, 2.8) and to hyperdissipative `σ ≥ 5/4` (window empty, 2.10). Applied to viscous Burgers (no pressure, no incompressibility): B2 fails at the first line, since \`(U·∇)U = ∇(                                                                                                                                                                                                                                                                                                                                        |
| T3 (blowup profile constraints)           | passes, item by item, for the whole window: `‖u‖_3 → ∞` at rate `s^{2β−1}` ✓ (2.4); `∫‖ω‖_∞ dt` divergent ✓ (2.7); Leray's two lower bounds strictly satisfied ✓ (2.5); singular set a point, dissipation per scale `r^{2−1/β} → ∞` ✓ (2.6); not exactly self-similar (B4.2) ✓ (2.3); if axisymmetric, not Type I since `β ≠ 1/2` ✓ (2.9). No constraint is violated, which is why the window is recorded as open rather than FAILED.                                                                                                                        |
| T4 (inviscid limit)                       | passes: the mandatory scaling check `νΔω` versus `(ω·∇)u` is done in B5.4/B7.2 and yields the exact threshold `β = 1/2`; and B7.3 explicitly refuses to conclude perturbativity from formal order counting.                                                                                                                                                                                                                                                                                                                                                  |
| T5 (weak solutions)                       | passes/not applicable: everything here is for a classical solution smooth on `[0,T*)`. Smoothness is used in B4 (differentiating the ansatz), B2/B3 (pointwise identities and the maximum principle) and B6 (the energy identity as an equality).                                                                                                                                                                                                                                                                                                            |
| T6 (too-easy test)                        | triggered once and honoured: the profile identity `ν‖∇U‖^2 = ((5β−2)/2)‖U‖^2` looked like it might exclude an interval of `β` cheaply; it was re-derived with the dimension explicit, checked against the original energy identity (they agree, B6.1), checked numerically (script Part 5, relative error `3.7e-10`), and its actual reach (`β < 2/5`, not `β > 1/2`) is recorded as a FAILED sub-item in B5 rather than overstated. The symbolic identities of B2/B8 were given seven falsifying controls (B2 table) so that "residual = 0" is not vacuous. |

## What would be needed next

- **Reprove or replace the decay lemma of B3.** Concretely: show that a smooth solution of the profile equation with `U ∈ L^3(R^3)` satisfies (D0)–(D2). This is a self-contained elliptic-regularity problem and would make the `β = 1/2` exclusion fully self-contained here.
- **Read Chae's "Nonexistence of asymptotically self-similar singularities in the Euler and the Navier–Stokes equations" (arXiv:math/0604234)** [checked via WebSearch: title only] and check whether its hypotheses cover the window `β ∈ [2/5, 1/2)`. This is the single most likely existing result to close or shrink the surviving window, and I could not retrieve its statement in this environment.
- **Settle the endpoint `β = 2/5`** (B6.2): prove or refute that a suitable weak solution cannot concentrate a fixed positive amount of `L^2` energy at a point at its first singular time. If it cannot, the window is the open interval `(2/5, 1/2)`.
- **Search for the Euler self-similar profile** in the window, in the shooting/PINN style of 2.14, with the tail exponent `(1−β)/β ∈ (1, 3/2]` imposed as a boundary condition at infinity and `β` treated as a nonlinear eigenvalue — this is the concrete numerical object this line produces. Any candidate must then be run through the `Re(s) = s^{2β−1}/ν` check and the singular-perturbation objection of B7.3.
- **Attack the DSS obstruction directly** (B8): find any inequality controlling `∂_τ P` by `|ω|^2` and the dissipation over one DSS period, or prove that none exists. Since `∂_τ P = Σ R_iR_j ∂_τ(V_iV_j)`, this is a question about Riesz transforms of time-derivatives of quadratic expressions, and a negative answer would be a clean statement that the NRS mechanism is intrinsically restricted to exact self-similarity.

## How to reproduce

```bash
cd research/navier-stokes
ruff check --line-length 120 line_b_selfsimilar.py && ruff format --check --line-length 120 line_b_selfsimilar.py
python3 line_b_selfsimilar.py
```

Fixed seed `20260906`; no files are written; runtime about 16 s on 4 CPUs (sympy generic-function identities under 1 s each, the exponent substitution about 10 s, the spectral numerical check under 1 s). The last printed line is `All checks passed: True`. Every table and residual quoted above is copied verbatim from that output.

## References

- J. Nečas, M. Růžička, V. Šverák, "On Leray's self-similar solutions of the Navier–Stokes equations", Acta Math. 176 (1996), 283–294 \[checked via WebSearch: authors, title, venue, volume, pages; the internal structure of the proof is from memory and is *not* reproduced here — see B3\].
- T.-P. Tsai, "On Leray's self-similar solutions of the Navier–Stokes equations satisfying local energy estimates", Arch. Rational Mech. Anal. 143 (1998), 29–51 \[checked via WebSearch in `00_problem.md` 2.3\].
- D. Chae, J. Wolf, "On the Liouville type theorems for self-similar solutions to the Navier–Stokes equations", Arch. Rational Mech. Anal. 225 (2017), 549–572, arXiv:1609.06962 \[checked via WebSearch: authors, title, venue, volume, pages; and that they generalize NRŠ and Tsai and exclude asymptotically self-similar blowup with profile in `L^{p,∞}`, `p > 3/2`\].
- D. Chae, "Nonexistence of self-similar singularities for the 3D incompressible Euler equations", Comm. Math. Phys. 273 (2007), arXiv:math/0601060 \[checked via WebSearch: title, venue, arXiv id, and the theorem statement "no self-similar blowup if `curl V ∈ L^p` for all `p ∈ (0,p_1)` for some `p_1 > 0`"; page numbers unverified\].
- D. Chae, "Nonexistence of asymptotically self-similar singularities in the Euler and the Navier–Stokes equations", arXiv:math/0604234 [checked via WebSearch: title and arXiv id only; hypotheses not recalled and deliberately not guessed].
- Z. Bradshaw, T.-P. Tsai, "Discretely self-similar solutions to the Navier–Stokes equations with data in `L^2_loc` satisfying the local energy inequality", Analysis & PDE 12 (2019), 1943–… [checked via WebSearch: title and venue; page range partially unverified]. Companion existence results by Chae–Wolf \[checked via WebSearch: DSS solutions for any DSS data in `L^2_loc`\].
- T. M. Elgindi, Ann. of Math. 194 (2021), 647–727 (`C^{1,α}` Euler blowup) \[checked via WebSearch in `00_problem.md` 2.13; the value of the self-similar exponent `β` is **[unverified, not recalled]** and is not used\]. Higher-regularity follow-ups appear in 2025–2026 arXiv titles such as "Self-similar blow-up solutions of incompressible Euler equations in `R^d`, `d ≥ 3` with `C^{1,1−2/d−}` velocity" and "Asymptotically Self-Similar Blowup for 3D Incompressible Euler with `C^{1,1/3−}` Velocity" [checked via WebSearch: titles only; authors, statements and status unverified].
- J. Leray, Acta Math. 63 (1934) (the self-similar ansatz and the lower bounds); L. Escauriaza, G. Seregin, V. Šverák (2003); L. Caffarelli, R. Kohn, L. Nirenberg (1982); J. T. Beale, T. Kato, A. Majda (1984); J.-L. Lions (1969) — all as cited and tagged in `00_problem.md` 2.3–2.10; used here only as KNOWN inputs.
