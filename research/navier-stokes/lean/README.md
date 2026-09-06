# Machine-checked certificates for Line A (Lean 4.15.0, core only, no mathlib)

This directory contains one Lean file, `Certificates.lean`, which re-proves inside the Lean kernel the *finite algebraic certificates* that `../line_a_lyapunov.md` currently takes on the word of the Python script `../line_a_lyapunov.py`. Nothing else. The point is narrow and worth stating precisely: for these particular numbers the status tag moves from "a script printed this" to "the Lean kernel evaluated this, and the proof term depends on no axioms". Everything analytic in Line A remains exactly as unformalized as it was.

`Certificates.lean` is checked by running `lean` on it directly. There is no lakefile and no Lake project, on purpose: a single file with no dependencies is the most reproducible artifact available here.

## What is kernel-checked

Fifty-six theorems, all proved by `decide` (kernel evaluation of closed arithmetic). They fall into five groups.

**1. The A2 triad certificate (3D).** The triad is `k = (1,0,0)`, `p = (0,1,1)`, `q = (-1,-1,-1)`; the two Gaussian-integer amplitude choices are copied from A2 of the notes and from the `choices` list of `exact_part` in the script. Lean recomputes the energy transfers from the amplitudes — it does not take them as data:

| theorem                                               | statement                                                                                                                                                                                                                                                                                |
| ----------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `triad_sums_to_zero`, `triad_norms`                   | `k + p + q = 0`; `(abs k)^2, (abs p)^2, (abs q)^2 = 1, 2, 3`                                                                                                                                                                                                                             |
| `triadA_divFree`, `triadB_divFree`                    | every amplitude is orthogonal to its wavevector (the fields are divergence-free)                                                                                                                                                                                                         |
| `triad_nondegenerate`                                 | among the 56 unordered triples with repetition drawn from `k, p, q, -k, -p, -q`, exactly two sum to zero: `(0,1,2)` and `(3,4,5)`. This is the script's `triad_nondegenerate`, and it is what licenses "the nonlinear term at a mode of the support comes only from the other two modes" |
| `triadA_transfers`, `triadB_transfers`                | `(T_k, T_p, T_q) = (-1, 1, 0)` and `(1, -3, 2)`, in units of `2 pi`                                                                                                                                                                                                                      |
| `triadA_transfer_sum`, `triadB_transfer_sum`          | `T_k + T_p + T_q = 0` (detailed energy conservation on the triad, K5 restricted to one triad)                                                                                                                                                                                            |
| `triadA_transfer_conj`, `triadB_transfer_conj`        | `T_{-n} = T_n` (reality), which the script asserts                                                                                                                                                                                                                                       |
| `triadA_proj_agrees`, `triadB_proj_agrees`            | the Leray projector may be dropped: `(abs n)^2 T_n(projected) = (abs n)^2 T_n` at all six modes, in the division-free form                                                                                                                                                               |
| `triadMatrix_value`, `triad_det`, `triad_det_ne_zero` | `det [1 0; -3 2] = 2`, nonzero — the `(T_p, T_q)` vectors span the plane                                                                                                                                                                                                                 |

**2. The 2D negative control.** Without it the 3D determinant means nothing, so it is formalized on the same footing: planar triad `k = (1,0,0)`, `p = (0,1,0)`, `q = (-1,-1,0)`, planar fields, same recomputation. `planarA_transfers` and `planarB_transfers` give `(1,-1,0)` and `(7,-7,0)`; `planarA_transfer_sum`, `planarB_transfer_sum` give the energy sum `0`; `planarA_enstrophy_conserved`, `planarB_enstrophy_conserved` give the *extra* 2D constraint `sum (abs n)^2 T_n = 0`; and `planar_det_zero` gives the vanishing determinant `det [-1 0; -7 0] = 0`. The contrast between `triad_det_ne_zero` and `planar_det_zero` is the mechanism the notes point at: in 3D the enstrophy constraint is lost.

**3. The A3.1 enstrophy budget.** Also recomputed from the amplitudes, in the units of the notes: `triadA_budget` gives `(energy, enstrophy, Euler flux, heat flux) = (14, 28, 4, -128)` (that is, `int abs(u)^2 = 14`, `int abs(omega)^2 = 28 (2 pi)^2`, Euler flux `4 (2 pi)^3`, heat flux `-128 (2 pi)^4`), and `triadB_budget` gives `(26, 58, 4, -292)`. The Euler flux is obtained through the notes' own formula `F_m = 4 sum m(n) T_n` with symbol `m(n) = (abs n)^2`, so it is a consequence of the transfers of group 1 rather than a separate datum; the agreement of the resulting `4` with the script's independently computed `2 <-Lap u, B(u,u)>` is evidence that the transcribed definition is the intended one. `energy_flux_zero` gives the energy Euler flux `0` (K5) for all four fields. With `A = a (2 pi)`, `ensRate` is the initial enstrophy rate in units of `(2 pi)^6`; `triadA_threshold/above/below` check that it vanishes at `a = 32` (that is `A* = 64 pi`) and changes sign across it, and `triadB_threshold/above/below` do the same at `a = 73` (`A* = 146 pi`).

**4. `L^4` thresholds (data, not recomputation).** The `L^4` fluxes need the pressure and are *not* recomputed in Lean; they enter as the script's exact output. What is checked is the rational arithmetic that produces the threshold amplitude the notes quote: `thresholdOK F G a*` is `F a* = -G` by cross-multiplication of integer pairs (core Lean 4.15 has no usable `Rat` numerals). `l4_threshold_choice1` checks `(32/9)(1305/2) = 2320`, that is `A* = 1305 pi` for the field `-u` of choice 1; `l4_threshold_choice2` checks `(1280/9)(12987/160) = 11544`; `ens_threshold_rational` and `ens_threshold_rational2` re-check the two enstrophy thresholds in the same rational format. `l4RateA_below` and `l4RateA_above` exhibit the sign change of the `L^4` rate at integer amplitudes `652` and `653` straddling `a* = 652.5`.

**5. The A6 certificate.** The five matrices `N(w^(m))` need the pressure too, so they are data, taken from the script. `A6rows_value` checks that the coordinate rows `(N11, N22, N12, N13, N23)` are exactly the five rows printed in A6(4) of the notes; `A6_symmetric` checks symmetry; `A6_traceless` checks `tr N(w) = 0` for each; `A6_det` and `A6_det_ne_zero` check `det = 7936 != 0`; and `A6_left_inverse` checks `A6adj * A6rows = 7936 I`, exhibiting an explicit integer matrix that becomes a left inverse after dividing by `7936`.

Two remarks on group 5. The trace check is not a restatement of the coordinate convention: the notes print only the five coordinates, and `N33` (which the trace needs) is *not* in the notes — it was obtained by rerunning the script's exact machinery, so `A6_traceless` genuinely tests that third, independently computed integral against the two printed diagonal entries. And a deliberate mutation test confirms this: changing the `N33` entry of `N3` from `-6` to `-5` breaks `A6_traceless` and nothing else.

The determinant function `det` is pinned down by seven sanity theorems (`det_id2`, `det_id3`, `det_id5`, `det_repeated_row`, `det_swap`, `det_2x2`, `det_3x3`) before it is used, since a kernel-checked proof about a wrongly defined determinant would prove nothing. The same `det` is used for the `2 x 2` triad matrices and the `5 x 5` A6 matrix.

## What is NOT formalized

This is the important half of this README. The Lean file is a certificate checker, not a proof of any theorem of Line A.

- **The analytic core is absent.** A1 and its relatives (A1b, A1c, A1c', A1d, A1e, A1f, A1g) — the lemmas that do the actual work, turning monotonicity of a functional into invariance under Euler — are not formalized in any form. They involve Frechet derivatives on function spaces, local existence for Navier-Stokes, homogeneity limits `lambda -> infinity` and `lambda -> 0`, Danskin's theorem, and integration by parts on the torus. Formalizing them needs mathlib (real and functional analysis, measure theory, Sobolev spaces) plus a substantial development of periodic Navier-Stokes that does not exist there today; it is a multi-person-year project, not something this file approaches. Nothing here should be read as evidence for A1.
- **A6 is not formalized as a theorem.** Only its final linear-algebra certificate is. Steps (1)-(3) and (5) of the A6 proof — the flux formula `DQ(u)[B(u,u)] = int p tr(D^2 phi(u) grad u)`, the reduction to the traceless Hessian, the perturbation `u_eps = v_0 + eps w` with `p(u_eps) = eps^2 p(w)`, and the integration of `D^2 phi = c I` to a quadratic — are pure prose here, as is A2's reduction of the Euler flux to `F_m = 4((m(p) - m(k)) T_p + (m(q) - m(k)) T_q)` and the passage from "some `u` on the triad has nonzero flux" to statements about `H^s` and `Hdot^s` for real `s` (which needs real exponentials: `2^s - 1` and `3^s - 1` nonzero for `s != 0`).
- **The definition of the transfer is transcribed, not derived.** Lean computes `T_n = Im(conj(a_n) . sum_{n1+n2=n} (u_{n1} . n2) u_{n2})`, which is the notes' formula for `Re(conj(u_n) B(u,u)_n)` in units of `2 pi`, with the Leray projector handled by `transferProj`. That this formula *is* the Fourier coefficient of `-P((u.grad)u)` is mathematics done by hand in the notes and by the script's spectral algebra; Lean checks the arithmetic of the formula, not its derivation from the PDE. Three independent things make the transcription trustworthy: it reproduces the script's printed transfers for all four fields, it reproduces the planar control, and the enstrophy Euler flux it derives (`4`, for both choices) matches the value the script computes by two other routes (`2 <-Lap u, B>` and `2 int omega . (omega . grad) u`).
- **Spanning over the rationals is not formalized.** `A6_left_inverse` gives `A6adj * A6rows = 7936 I`; the one-line implication "hence the rows are linearly independent over `Q`, hence the `N(w)` span `Sym_0`, hence `tr(H_0(v_0) N) = 0` for all `N` forces `H_0(v_0) = 0`" is left to the reader. Formalizing it would need either vector spaces over a field (mathlib) or the `omega` tactic, and `omega` costs axioms (see below).
- **The sign conclusions hold only at the integer amplitudes exhibited.** `ensRate 33 > 0` and `ensRate 31 < 0` are two points; "the enstrophy increases for every real `A > 64 pi`" is not a statement this file can even express, since core Lean has no reals.
- **The `L^4` fluxes and the A6 matrices are inputs.** They come from the script's exact `Fraction` arithmetic. If the script's spectral algebra were wrong there, this file would faithfully certify wrong numbers. The A2/A3.1 group does not have this weakness: it is recomputed from the amplitudes.
- **Everything else in the project.** A3.2, A3.3, A3.4, A4, A5, the FFT computations, the sup-norm searches and every floating-point number in the notes are untouched here.

## Why `decide` only, and no axioms

All 56 proofs use `decide` (or would use `rfl`): the Lean kernel evaluates the closed arithmetic term and checks it reduces to `true`. Three alternatives are deliberately avoided.

- `native_decide` would compile the computation and trust the compiled code, introducing the axiom `Lean.ofReduceBool` — precisely the "trust me, it ran" that this exercise is meant to eliminate.
- `omega` (linear integer arithmetic) would let us state the spanning lemma for all integer vectors, but its proof terms depend on `propext` and `Quot.sound`. Those are standard parts of Lean's logic, not classical choice, and relying on them would be perfectly respectable — but the stronger and simpler claim "depends on no axioms at all" is available here, so it is what the file delivers.
- `sorry` is not used anywhere; any use of it would make the file report `sorryAx`.

Mathlib is not available in this environment (its cache hosts are unreachable and a source build is not feasible), which is the reason the file uses `Int` pairs for rationals and lists of `Int` for matrices instead of `Rat`, `Matrix` and `Polynomial`.

## Cross-validation against the notes and the script, and discrepancies

Every number in `Certificates.lean` was taken from `../line_a_lyapunov.md` and checked against a fresh run of Part 2 (`exact_part`) of `../line_a_lyapunov.py`. **No discrepancy was found**: the transfers `(-1,1,0)` and `(1,-3,2)`, the determinant `2`, the planar transfers `(1,-1,0)` and `(7,-7,0)` with determinant `0`, the budget `(14, 28, 4, -128)` and `(26, 58, 4, -292)`, the thresholds `64 pi`, `146 pi`, `1305 pi` and `(12987/160)(2 pi)`, the five A6 rows and `det = 7936` agree between notes, script and Lean.

Three items in the Lean file come from the script's output only and are *not* printed in the notes: the entry `N33` of each `N(w^(m))` (needed for `A6_traceless`), the choice-2 `L^4` fluxes `1280/9` and `-11544` (the notes print only the resulting threshold `(12987/160)(2 pi)`), and the choice-2 energy `26` and `int abs(u)^4 = 1286`. They are labelled as such in the Lean comments. This is the same class of gap the round-1 referee flagged for the A6 amplitudes, which the notes then pasted in.

## Reproducing the check from scratch

Lean 4.15.0, no other dependency, about ten seconds of CPU:

```text
curl -LO https://github.com/leanprover/lean4/releases/download/v4.15.0/lean-4.15.0-linux.tar.zst
tar --zstd -xf lean-4.15.0-linux.tar.zst          # or: unzstd lean-4.15.0-linux.tar.zst && tar -xf lean-4.15.0-linux.tar
export PATH="$PWD/lean-4.15.0-linux/bin:$PATH"    # the extracted toolchain's bin directory
lean --version                                    # Lean (version 4.15.0, x86_64-unknown-linux-gnu, commit 11651562caae, Release)
cd research/navier-stokes/lean
lean Certificates.lean
```

`lean Certificates.lean` must exit `0` and print nothing except the axiom audit below — no errors, no warnings. Any failure of any certificate would surface as `error: tactic 'decide' proved that the proposition ... is false`.

To confirm the theorems are load-bearing rather than vacuous, mutate the file in a scratch copy: changing one Gaussian-integer amplitude of `triadA` (say `a_p = (1, i, -i)` to `(1, i, i)`) breaks eight theorems including `triadA_divFree`, `triadA_transfers` and `triadA_budget`, and changing `N33` of `N3` breaks `A6_traceless` alone.

## Axiom audit (verbatim output)

The complete output of `lean Certificates.lean`, one line per theorem:

```text
'NSCert.det_id2' does not depend on any axioms
'NSCert.det_id3' does not depend on any axioms
'NSCert.det_id5' does not depend on any axioms
'NSCert.det_repeated_row' does not depend on any axioms
'NSCert.det_swap' does not depend on any axioms
'NSCert.det_2x2' does not depend on any axioms
'NSCert.det_3x3' does not depend on any axioms
'NSCert.triad_sums_to_zero' does not depend on any axioms
'NSCert.triad_norms' does not depend on any axioms
'NSCert.triadA_divFree' does not depend on any axioms
'NSCert.triadB_divFree' does not depend on any axioms
'NSCert.triad_nondegenerate' does not depend on any axioms
'NSCert.triadA_transfers' does not depend on any axioms
'NSCert.triadB_transfers' does not depend on any axioms
'NSCert.triadA_transfer_sum' does not depend on any axioms
'NSCert.triadB_transfer_sum' does not depend on any axioms
'NSCert.triadA_transfer_conj' does not depend on any axioms
'NSCert.triadB_transfer_conj' does not depend on any axioms
'NSCert.triadA_proj_agrees' does not depend on any axioms
'NSCert.triadB_proj_agrees' does not depend on any axioms
'NSCert.triadMatrix_value' does not depend on any axioms
'NSCert.triad_det' does not depend on any axioms
'NSCert.triad_det_ne_zero' does not depend on any axioms
'NSCert.planar_triad_sums_to_zero' does not depend on any axioms
'NSCert.planarA_divFree' does not depend on any axioms
'NSCert.planarB_divFree' does not depend on any axioms
'NSCert.planar_nondegenerate' does not depend on any axioms
'NSCert.planarA_transfers' does not depend on any axioms
'NSCert.planarB_transfers' does not depend on any axioms
'NSCert.planarA_transfer_sum' does not depend on any axioms
'NSCert.planarB_transfer_sum' does not depend on any axioms
'NSCert.planarA_enstrophy_conserved' does not depend on any axioms
'NSCert.planarB_enstrophy_conserved' does not depend on any axioms
'NSCert.planarMatrix_value' does not depend on any axioms
'NSCert.planar_det_zero' does not depend on any axioms
'NSCert.triadA_budget' does not depend on any axioms
'NSCert.triadB_budget' does not depend on any axioms
'NSCert.energy_flux_zero' does not depend on any axioms
'NSCert.triadA_threshold' does not depend on any axioms
'NSCert.triadA_above' does not depend on any axioms
'NSCert.triadA_below' does not depend on any axioms
'NSCert.triadB_threshold' does not depend on any axioms
'NSCert.triadB_above' does not depend on any axioms
'NSCert.triadB_below' does not depend on any axioms
'NSCert.ens_threshold_rational' does not depend on any axioms
'NSCert.ens_threshold_rational2' does not depend on any axioms
'NSCert.l4_threshold_choice1' does not depend on any axioms
'NSCert.l4_threshold_choice2' does not depend on any axioms
'NSCert.l4RateA_below' does not depend on any axioms
'NSCert.l4RateA_above' does not depend on any axioms
'NSCert.A6_symmetric' does not depend on any axioms
'NSCert.A6_traceless' does not depend on any axioms
'NSCert.A6rows_value' does not depend on any axioms
'NSCert.A6_det' does not depend on any axioms
'NSCert.A6_det_ne_zero' does not depend on any axioms
'NSCert.A6_left_inverse' does not depend on any axioms
```
