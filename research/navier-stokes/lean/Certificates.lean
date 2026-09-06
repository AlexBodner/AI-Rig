/-
  Certificates.lean -- kernel-checked finite algebraic certificates for Line A
  (coercive Lyapunov functionals) of the Navier-Stokes transparency project.

  Scope.  This file formalizes ONLY the exact finite arithmetic that the notes
  `research/navier-stokes/line_a_lyapunov.md` currently take from the Python
  script `research/navier-stokes/line_a_lyapunov.py`.  Nothing analytic is
  formalized here: no Frechet derivatives, no PDE, no function spaces.  See
  `README.md` in this directory for the precise boundary.

  Conventions (matching the notes and the script).  Unit torus T^3 = R^3/Z^3,
  Fourier modes exp(2 pi i n.x), a real field is u = sum_n 2 Re[a_n e^{2 pi i n.x}]
  over a set of representative modes n, so the full coefficient family is
  u_n = a_n and u_{-n} = conj(a_n).  Derivatives contribute a factor 2 pi, which
  is carried symbolically: every quantity below is an integer coefficient of a
  power of 2 pi, exactly as the script's `pi_pow` bookkeeping does.

  Built with Lean 4.15.0, core only (no mathlib).  Every proof is `decide` or
  `rfl`, so the kernel evaluates the arithmetic; `#print axioms` at the end of
  the file reports that no theorem depends on any axiom.  `native_decide`,
  `omega` and `sorry` are deliberately not used (see README).
-/

namespace NSCert

/-! ## Integer 3-vectors (wavevectors) -/

abbrev V3 := Int × Int × Int

def vadd (a b : V3) : V3 := (a.1 + b.1, a.2.1 + b.2.1, a.2.2 + b.2.2)
def vneg (a : V3) : V3 := (-a.1, -a.2.1, -a.2.2)
def vdotZ (a b : V3) : Int := a.1 * b.1 + a.2.1 * b.2.1 + a.2.2 * b.2.2
def normSqZ (n : V3) : Int := vdotZ n n

/-! ## Gaussian integers and Gaussian-integer 3-vectors -/

structure GI where
  re : Int
  im : Int
deriving DecidableEq, Repr

def gzero : GI := ⟨0, 0⟩
def gadd (a b : GI) : GI := ⟨a.re + b.re, a.im + b.im⟩
def gsub (a b : GI) : GI := ⟨a.re - b.re, a.im - b.im⟩
def gmul (a b : GI) : GI := ⟨a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re⟩
def gconj (a : GI) : GI := ⟨a.re, -a.im⟩
def gsmul (c : Int) (a : GI) : GI := ⟨c * a.re, c * a.im⟩
def gnormSq (a : GI) : Int := a.re * a.re + a.im * a.im

structure GV where
  x : GI
  y : GI
  z : GI
deriving DecidableEq, Repr

def gvzero : GV := ⟨gzero, gzero, gzero⟩
def gvadd (a b : GV) : GV := ⟨gadd a.x b.x, gadd a.y b.y, gadd a.z b.z⟩
def gvsub (a b : GV) : GV := ⟨gsub a.x b.x, gsub a.y b.y, gsub a.z b.z⟩
def gvconj (a : GV) : GV := ⟨gconj a.x, gconj a.y, gconj a.z⟩
def gvsmulG (c : GI) (a : GV) : GV := ⟨gmul c a.x, gmul c a.y, gmul c a.z⟩
def gvsmulZ (c : Int) (a : GV) : GV := ⟨gsmul c a.x, gsmul c a.y, gsmul c a.z⟩
/-- The Gaussian vector `c * n` for an integer wavevector `n`. -/
def gvOfVec (c : GI) (n : V3) : GV := ⟨gsmul n.1 c, gsmul n.2.1 c, gsmul n.2.2 c⟩
/-- Bilinear (unconjugated) pairing of a Gaussian vector with an integer vector. -/
def gvdotZ (a : GV) (n : V3) : GI := gadd (gadd (gsmul n.1 a.x) (gsmul n.2.1 a.y)) (gsmul n.2.2 a.z)
/-- Hermitian pairing `conj(a) . b`. -/
def gvhdot (a b : GV) : GI := gadd (gadd (gmul (gconj a.x) b.x) (gmul (gconj a.y) b.y)) (gmul (gconj a.z) b.z)
def gvNormSq (a : GV) : Int := gnormSq a.x + gnormSq a.y + gnormSq a.z

/-! ## Trigonometric-polynomial fields given by Gaussian-integer amplitudes

A field is a list of `(mode, amplitude)` pairs, one entry per representative
mode; `closure` adds the conjugate modes, giving the full coefficient family. -/

structure Mode where
  n : V3
  a : GV
deriving DecidableEq, Repr

abbrev Field := List Mode

def closure (F : Field) : Field := F ++ F.map (fun m => ⟨vneg m.n, gvconj m.a⟩)

def ampAt (F : Field) (n : V3) : GV :=
  match F.find? (fun m => decide (m.n = n)) with
  | some m => m.a
  | none => gvzero

/-- Divergence-free test: every amplitude is orthogonal to its wavevector. -/
def divFree (F : Field) : Bool := F.all (fun m => decide (gvdotZ m.a m.n = gzero))

/-! ### The nonlinear term on a single triad

For the Euler nonlinearity `B(u,u) = -P((u.grad)u)` the Fourier coefficient at
mode `n` is `B_n = -P_n [ i sum_{n1+n2=n} (u_{n1}.n2) u_{n2} ]` in units of
`2 pi`, with `P_n v = v - n (n.v)/|n|^2` the Leray projector (notes, A2, "Hand
check of the certificate").  `convSum` is the inner sum. -/

def convSum (F : Field) (n : V3) : GV :=
  F.foldl
    (fun acc m1 =>
      F.foldl
        (fun acc2 m2 =>
          if vadd m1.n m2.n = n then gvadd acc2 (gvsmulG (gvdotZ m1.a m2.n) m2.a) else acc2)
        acc)
    gvzero

/-- `|n|^2` times the Leray-projected sum, so that the projector needs no division. -/
def convSumProj (F : Field) (n : V3) : GV :=
  gvsub (gvsmulZ (normSqZ n) (convSum F n)) (gvOfVec (gvdotZ (convSum F n) n) n)

/-- Energy transfer `T_n = Re(conj(u_n) . B_n)` in units of `2 pi`; since
`Re(conj(a) . (-i) v) = Im(conj(a) . v)` this is an imaginary part. -/
def transfer (F : Field) (n : V3) : Int := (gvhdot (ampAt (closure F) n) (convSum (closure F) n)).im

/-- The same with the Leray projector kept, multiplied by `|n|^2`. -/
def transferProj (F : Field) (n : V3) : Int :=
  (gvhdot (ampAt (closure F) n) (convSumProj (closure F) n)).im

/-! ### Spectral functionals of the field, in the units of the notes -/

def sumOver (F : Field) (f : Mode → Int) : Int := F.foldl (fun acc m => acc + f m) 0

/-- `integral |u|^2`. -/
def energy (F : Field) : Int := 2 * sumOver F (fun m => gvNormSq m.a)
/-- `integral |curl u|^2` in units of `(2 pi)^2`. -/
def enstrophy (F : Field) : Int := 2 * sumOver F (fun m => normSqZ m.n * gvNormSq m.a)
/-- Enstrophy heat flux `2 <-Lap u, Lap u>` in units of `(2 pi)^4`. -/
def enstrophyHeatFlux (F : Field) : Int :=
  -4 * sumOver F (fun m => normSqZ m.n * normSqZ m.n * gvNormSq m.a)
/-- Euler flux of `Q_m = sum m(n) |u_n|^2`, i.e. `4 sum_{k,p,q} m(n) T_n` (notes, A2). -/
def eulerFluxSymbol (F : Field) (m : V3 → Int) : Int := 4 * sumOver F (fun md => m md.n * transfer F md.n)
/-- Enstrophy Euler flux in units of `(2 pi)^3` (symbol `m(n) = |n|^2`). -/
def enstrophyEulerFlux (F : Field) : Int := eulerFluxSymbol F normSqZ
/-- Energy Euler flux (symbol `m(n) = 1`); zero by K5. -/
def energyEulerFlux (F : Field) : Int := eulerFluxSymbol F (fun _ => 1)


/-! ## Determinants of small integer matrices (cofactor expansion along the first row) -/

def minorRows (m : List (List Int)) (j : Nat) : List (List Int) := (m.drop 1).map (fun row => row.eraseIdx j)

def det : Nat → List (List Int) → Int
  | 0, _ => 1
  | n + 1, m =>
    (List.range (n + 1)).foldl
      (fun acc j => acc + (if j % 2 = 0 then 1 else -1) * ((m.headD []).getD j 0) * det n (minorRows m j)) 0

/-- Sanity checks fixing the meaning of `det` (a kernel-checked wrong definition would prove nothing). -/
theorem det_id2 : det 2 [[1, 0], [0, 1]] = 1 := by decide
theorem det_id3 : det 3 [[1, 0, 0], [0, 1, 0], [0, 0, 1]] = 1 := by decide
theorem det_id5 : det 5 [[1,0,0,0,0], [0,1,0,0,0], [0,0,1,0,0], [0,0,0,1,0], [0,0,0,0,1]] = 1 := by decide
theorem det_repeated_row : det 3 [[1, 2, 3], [1, 2, 3], [4, 5, 6]] = 0 := by decide
theorem det_swap : det 2 [[3, 4], [1, 2]] = -det 2 [[1, 2], [3, 4]] := by decide
theorem det_2x2 : det 2 [[1, 2], [3, 4]] = -2 := by decide
theorem det_3x3 : det 3 [[1, 2, 3], [4, 5, 6], [7, 8, 10]] = -3 := by decide

/-! ## A2: the 3D triad certificate

Triad `k = (1,0,0)`, `p = (0,1,1)`, `q = (-1,-1,-1)`, with `k + p + q = 0` and
`|k|^2 = 1`, `|p|^2 = 2`, `|q|^2 = 3`.  Two Gaussian-integer amplitude choices,
copied from `line_a_lyapunov.md` (A2, "Certificate (script, Part 2)") and from
`line_a_lyapunov.py` (`exact_part`, list `choices`). -/

def k3 : V3 := (1, 0, 0)
def p3 : V3 := (0, 1, 1)
def q3 : V3 := (-1, -1, -1)

/-- choice 1: `a_k = (0, 1, i)`, `a_p = (1, i, -i)`, `a_q = (1, -1, 0)`. -/
def triadA : Field :=
  [ ⟨k3, ⟨⟨0, 0⟩, ⟨1, 0⟩, ⟨0, 1⟩⟩⟩,
    ⟨p3, ⟨⟨1, 0⟩, ⟨0, 1⟩, ⟨0, -1⟩⟩⟩,
    ⟨q3, ⟨⟨1, 0⟩, ⟨-1, 0⟩, ⟨0, 0⟩⟩⟩ ]

/-- choice 2: `a_k = (0, 1+i, -1)`, `a_p = (1-i, 1, -1)`, `a_q = (1, 1, -2)`. -/
def triadB : Field :=
  [ ⟨k3, ⟨⟨0, 0⟩, ⟨1, 1⟩, ⟨-1, 0⟩⟩⟩,
    ⟨p3, ⟨⟨1, -1⟩, ⟨1, 0⟩, ⟨-1, 0⟩⟩⟩,
    ⟨q3, ⟨⟨1, 0⟩, ⟨1, 0⟩, ⟨-2, 0⟩⟩⟩ ]

theorem triad_sums_to_zero : vadd (vadd k3 p3) q3 = ((0 : Int), (0 : Int), (0 : Int)) := by decide
theorem triad_norms : (normSqZ k3, normSqZ p3, normSqZ q3) = ((1 : Int), (2 : Int), (3 : Int)) := by decide
theorem triadA_divFree : divFree triadA = true := by decide
theorem triadB_divFree : divFree triadB = true := by decide

/-- Non-degeneracy: among the 56 unordered triples (with repetition) drawn from the six
modes `k, p, q, -k, -p, -q` (in this order), exactly two sum to zero, namely `(k,p,q)`
and `(-k,-p,-q)`.  This is the script's `triad_nondegenerate`. -/
def zeroSumTriples (S : List V3) : List (Nat × Nat × Nat) :=
  let len := S.length
  (List.range len).foldl (fun acc i =>
    acc ++ (List.range len).foldl (fun acc2 j =>
      if j < i then acc2 else
        acc2 ++ (List.range len).foldl (fun acc3 l =>
          if l < j then acc3 else
            if vadd (vadd (S.getD i (0, 0, 0)) (S.getD j (0, 0, 0))) (S.getD l (0, 0, 0))
                = ((0 : Int), (0 : Int), (0 : Int))
              then acc3 ++ [(i, j, l)] else acc3) []) []) []

def modesOf (F : Field) : List V3 := (closure F).map (fun m => m.n)

theorem triad_nondegenerate : zeroSumTriples (modesOf triadA) = [(0, 1, 2), (3, 4, 5)] := by decide

/-! ### Energy transfers on the triad (units of `2 pi`) -/

theorem triadA_transfers :
    (transfer triadA k3, transfer triadA p3, transfer triadA q3) = ((-1 : Int), (1 : Int), (0 : Int)) := by
  decide

theorem triadB_transfers :
    (transfer triadB k3, transfer triadB p3, transfer triadB q3) = ((1 : Int), (-3 : Int), (2 : Int)) := by
  decide

/-- Detailed energy conservation on the triad (K5 restricted to one triad): the transfers sum to zero. -/
theorem triadA_transfer_sum : transfer triadA k3 + transfer triadA p3 + transfer triadA q3 = 0 := by decide
theorem triadB_transfer_sum : transfer triadB k3 + transfer triadB p3 + transfer triadB q3 = 0 := by decide

/-- Reality: `T_{-n} = T_n`, as the script asserts. -/
theorem triadA_transfer_conj :
    (transfer triadA (vneg k3), transfer triadA (vneg p3), transfer triadA (vneg q3))
      = (transfer triadA k3, transfer triadA p3, transfer triadA q3) := by decide
theorem triadB_transfer_conj :
    (transfer triadB (vneg k3), transfer triadB (vneg p3), transfer triadB (vneg q3))
      = (transfer triadB k3, transfer triadB p3, transfer triadB q3) := by decide

/-- The Leray projector does not change the transfers (it is killed by `a_n . n = 0`);
checked here in the division-free form `|n|^2 T_n^{proj} = |n|^2 T_n` at all six modes. -/
theorem triadA_proj_agrees :
    (modesOf triadA).all (fun n => decide (transferProj triadA n = normSqZ n * transfer triadA n)) = true := by
  decide
theorem triadB_proj_agrees :
    (modesOf triadB).all (fun n => decide (transferProj triadB n = normSqZ n * transfer triadB n)) = true := by
  decide

/-- The certificate proper: the two vectors `(T_p, T_q)` are linearly independent,
`det [1 0; -3 2] = 2` (notes, A2). -/
def triadMatrix : List (List Int) :=
  [[transfer triadA p3, transfer triadA q3], [transfer triadB p3, transfer triadB q3]]

theorem triadMatrix_value : triadMatrix = [[1, 0], [-3, 2]] := by decide
theorem triadMatrix_shape :
    triadMatrix.length = 2 ∧ triadMatrix.all (fun r => decide (r.length = 2)) = true := by decide
theorem triad_det : det 2 triadMatrix = 2 := by decide
theorem triad_det_ne_zero : det 2 triadMatrix ≠ 0 := by decide

/-! ## A2: the planar (2D) negative control

Planar triad `k = (1,0,0)`, `p = (0,1,0)`, `q = (-1,-1,0)` with planar fields.  Here the
`(T_p, T_q)` vectors are confined to a line and the corresponding determinant vanishes:
this is what makes the 3D certificate meaningful (notes, A2, "2D positive control"). -/

def k2 : V3 := (1, 0, 0)
def p2 : V3 := (0, 1, 0)
def q2 : V3 := (-1, -1, 0)

/-- planar choice 1: `a_k = (0, 1, 0)`, `a_p = (i, 0, 0)`, `a_q = (1, -1, 0)`. -/
def planarA : Field :=
  [ ⟨k2, ⟨⟨0, 0⟩, ⟨1, 0⟩, ⟨0, 0⟩⟩⟩,
    ⟨p2, ⟨⟨0, 1⟩, ⟨0, 0⟩, ⟨0, 0⟩⟩⟩,
    ⟨q2, ⟨⟨1, 0⟩, ⟨-1, 0⟩, ⟨0, 0⟩⟩⟩ ]

/-- planar choice 2: `a_k = (0, 2-i, 0)`, `a_p = (1+i, 0, 0)`, `a_q = (1+2i, -1-2i, 0)`. -/
def planarB : Field :=
  [ ⟨k2, ⟨⟨0, 0⟩, ⟨2, -1⟩, ⟨0, 0⟩⟩⟩,
    ⟨p2, ⟨⟨1, 1⟩, ⟨0, 0⟩, ⟨0, 0⟩⟩⟩,
    ⟨q2, ⟨⟨1, 2⟩, ⟨-1, -2⟩, ⟨0, 0⟩⟩⟩ ]

theorem planar_triad_sums_to_zero : vadd (vadd k2 p2) q2 = ((0 : Int), (0 : Int), (0 : Int)) := by decide
theorem planarA_divFree : divFree planarA = true := by decide
theorem planarB_divFree : divFree planarB = true := by decide
theorem planar_nondegenerate : zeroSumTriples (modesOf planarA) = [(0, 1, 2), (3, 4, 5)] := by decide

theorem planarA_transfers :
    (transfer planarA k2, transfer planarA p2, transfer planarA q2) = ((1 : Int), (-1 : Int), (0 : Int)) := by
  decide
theorem planarB_transfers :
    (transfer planarB k2, transfer planarB p2, transfer planarB q2) = ((7 : Int), (-7 : Int), (0 : Int)) := by
  decide

/-- Energy conservation on the planar triad. -/
theorem planarA_transfer_sum : transfer planarA k2 + transfer planarA p2 + transfer planarA q2 = 0 := by decide
theorem planarB_transfer_sum : transfer planarB k2 + transfer planarB p2 + transfer planarB q2 = 0 := by decide

/-- 2D enstrophy conservation on the triad: `sum |n|^2 T_n = 0` as well (this is the extra
constraint that 3D loses; equivalently, the enstrophy Euler flux vanishes for planar fields). -/
theorem planarA_enstrophy_conserved : enstrophyEulerFlux planarA = 0 := by decide
theorem planarB_enstrophy_conserved : enstrophyEulerFlux planarB = 0 := by decide

/-- The negative control: the planar `(T_p, T_q)` determinant is zero. -/
def planarMatrix : List (List Int) :=
  [[transfer planarA p2, transfer planarA q2], [transfer planarB p2, transfer planarB q2]]

theorem planarMatrix_value : planarMatrix = [[-1, 0], [-7, 0]] := by decide
theorem planarMatrix_shape :
    planarMatrix.length = 2 ∧ planarMatrix.all (fun r => decide (r.length = 2)) = true := by decide
theorem planar_det_zero : det 2 planarMatrix = 0 := by decide

/-! ## A3.1: the exact enstrophy budget of the triad fields

All four quantities are computed here from the Gaussian-integer amplitudes, in the units
of the notes: energy, enstrophy in `(2 pi)^2`, Euler flux in `(2 pi)^3`, heat flux in
`(2 pi)^4`.  The Euler flux is obtained through `eulerFluxSymbol` with the symbol
`m(n) = |n|^2`, i.e. through the formula `F_m = 4 sum m(n) T_n` of A2. -/

theorem triadA_budget :
    (energy triadA, enstrophy triadA, enstrophyEulerFlux triadA, enstrophyHeatFlux triadA)
      = ((14 : Int), (28 : Int), (4 : Int), (-128 : Int)) := by decide

theorem triadB_budget :
    (energy triadB, enstrophy triadB, enstrophyEulerFlux triadB, enstrophyHeatFlux triadB)
      = ((26 : Int), (58 : Int), (4 : Int), (-292 : Int)) := by decide

/-- The energy Euler flux vanishes (K5), for both choices and for the planar fields. -/
theorem energy_flux_zero :
    (energyEulerFlux triadA, energyEulerFlux triadB, energyEulerFlux planarA, energyEulerFlux planarB)
      = ((0 : Int), (0 : Int), (0 : Int), (0 : Int)) := by decide

/-- The `Hdot^{1/2}` Euler flux has symbol `m(n) = 2 pi |n|`, which is not rational; but its
decomposition over the three shells is integral.  `hHalfCoeff F j` is the coefficient of
`(2 pi)^2 sqrt(j)` in that flux, i.e. `4 T_n` summed over the modes with `|n|^2 = j`.  The
notes' value for choice 1 is `(2 pi)^2 (-4 + 4 sqrt 2)`; the irrational combination itself
is out of reach of core Lean (there are no reals), the integer coefficients are not. -/
def hHalfCoeff (F : Field) (j : Int) : Int :=
  eulerFluxSymbol F (fun n => if normSqZ n = j then 1 else 0)

theorem triadA_hhalf_coeffs :
    (hHalfCoeff triadA 1, hHalfCoeff triadA 2, hHalfCoeff triadA 3) = ((-4 : Int), (4 : Int), (0 : Int)) := by
  decide

theorem triadB_hhalf_coeffs :
    (hHalfCoeff triadB 1, hHalfCoeff triadB 2, hHalfCoeff triadB 3) = ((4 : Int), (-12 : Int), (8 : Int)) := by
  decide

/-- With `A = a (2 pi)`, `d/dt integral |omega|^2` at `t = 0` equals `(2 pi)^6` times
`ensRate F a` (notes, A3.1: `-128 (2 pi)^4 A^2 + 4 (2 pi)^3 A^3`). -/
def ensRate (F : Field) (a : Int) : Int :=
  enstrophyEulerFlux F * (a * a * a) + enstrophyHeatFlux F * (a * a)

/-- Choice 1 threshold `a = 32`, i.e. `A* = 64 pi`: the rate vanishes there and changes sign. -/
theorem triadA_threshold : ensRate triadA 32 = 0 := by decide
theorem triadA_above : ensRate triadA 33 > 0 := by decide
theorem triadA_below : ensRate triadA 31 < 0 := by decide

/-- Choice 2 threshold `a = 73`, i.e. `A* = 146 pi`. -/
theorem triadB_threshold : ensRate triadB 73 = 0 := by decide
theorem triadB_above : ensRate triadB 74 > 0 := by decide
theorem triadB_below : ensRate triadB 72 < 0 := by decide

/-! ## A3.1: the `L^4` thresholds (rational cross-multiplication)

The `L^4` fluxes involve the pressure and are NOT recomputed here; they are the script's
exact output, taken as data.  What is checked is the arithmetic that turns them into the
threshold amplitude printed in the notes: `F a* + G = 0` with `F`, `G`, `a*` rational,
in units where `A = a (2 pi)`, `F` counts `(2 pi)`, `G` counts `(2 pi)^2`.
Rationals are pairs of integers compared by cross-multiplication (core Lean 4.15 has no
usable `Rat` numerals). -/

structure Rat0 where
  num : Int
  den : Int
deriving DecidableEq, Repr

def qmul (a b : Rat0) : Rat0 := ⟨a.num * b.num, a.den * b.den⟩
def qneg (a : Rat0) : Rat0 := ⟨-a.num, a.den⟩
def qEqB (a b : Rat0) : Bool := a.num * b.den == b.num * a.den
/-- A pair of integers denotes a rational only when its denominator is positive; without this
guard `qEqB` would accept degenerate pairs such as `(0,0)`. -/
def qWF (a : Rat0) : Bool := decide (a.den > 0)
/-- `F * a* = -G`, i.e. the Navier-Stokes rate `F A^{d+1} + G A^d` vanishes at `A = a* (2 pi)`.
The three inputs are checked to be well formed, so this is a genuine rational identity. -/
def thresholdOK (F G astar : Rat0) : Bool :=
  qWF F && qWF G && qWF astar && qEqB (qmul F astar) (qneg G)

/-- Enstrophy, choice 1: `F = 4`, `G = -128`, `a* = 32` (`A* = 64 pi`), consistent with `ensRate`. -/
theorem ens_threshold_rational : thresholdOK ⟨4, 1⟩ ⟨-128, 1⟩ ⟨32, 1⟩ = true := by decide
/-- Enstrophy, choice 2: `F = 4`, `G = -292`, `a* = 73` (`A* = 146 pi`). -/
theorem ens_threshold_rational2 : thresholdOK ⟨4, 1⟩ ⟨-292, 1⟩ ⟨73, 1⟩ = true := by decide
/-- `L^4` at the field `-u`, choice 1: `F = 32/9`, `G = -2320`, `a* = 1305/2` (`A* = 1305 pi`). -/
theorem l4_threshold_choice1 : thresholdOK ⟨32, 9⟩ ⟨-2320, 1⟩ ⟨1305, 2⟩ = true := by decide
/-- `L^4` at the field `-u`, choice 2: `F = 1280/9`, `G = -11544`, `a* = 12987/160`.
(The two fluxes are script output only; the notes print just the threshold.) -/
theorem l4_threshold_choice2 : thresholdOK ⟨1280, 9⟩ ⟨-11544, 1⟩ ⟨12987, 160⟩ = true := by decide

/-- Sign change of the `L^4` rate of choice 1 (nine times the rate, in units `(2 pi)^6`,
for the field `-u`): negative below `a* = 652.5`, positive above. -/
def l4RateA (a : Int) : Int := 32 * (a * a * a * a * a) - 20880 * (a * a * a * a)

theorem l4RateA_below : l4RateA 652 < 0 := by decide
theorem l4RateA_above : l4RateA 653 > 0 := by decide

/-! ## A6: the five-field certificate for zeroth-order local functionals

`N(w) := sym part of the matrix M(w)_{ij} = integral p_w d_j w_i dx`, in units of `2 pi`,
for the five Gaussian-integer fields `w^(1), ..., w^(5)` of the notes (A6, step (4)).
These matrices are NOT recomputed here (they need the pressure, hence rational spectral
algebra): the entries are data from `line_a_lyapunov.py`.  The notes print the five
coordinate rows `(N11, N22, N12, N13, N23)` and the determinant; the diagonal entry `N33`
is not printed in the notes and was obtained by rerunning the script's exact machinery,
which is why the trace identity below is a genuine check on the data rather than a
restatement of the coordinate convention. -/

def N1 : List (List Int) := [[-6, 5, 6], [5, 4, 5], [6, 5, 2]]
def N2 : List (List Int) := [[-10, -1, 4], [-1, 4, -1], [4, -1, 6]]
def N3 : List (List Int) := [[2, -2, -6], [-2, 4, -10], [-6, -10, -6]]
def N4 : List (List Int) := [[-2, -1, -2], [-1, -2, -7], [-2, -7, 4]]
def N5 : List (List Int) := [[-6, 12, 10], [12, -10, 14], [10, 14, 16]]

def entry (m : List (List Int)) (i j : Nat) : Int := (m.getD i []).getD j 0
def trace3 (m : List (List Int)) : Int := entry m 0 0 + entry m 1 1 + entry m 2 2
def isSym3 (m : List (List Int)) : Bool :=
  decide (entry m 0 1 = entry m 1 0) && decide (entry m 0 2 = entry m 2 0) &&
    decide (entry m 1 2 = entry m 2 1)
/-- Coordinates on the 5-dimensional space `Sym_0` of traceless symmetric matrices. -/
def coords (m : List (List Int)) : List Int :=
  [entry m 0 0, entry m 1 1, entry m 0 1, entry m 0 2, entry m 1 2]

def A6mats : List (List (List Int)) := [N1, N2, N3, N4, N5]
def A6rows : List (List Int) := A6mats.map coords

/-- Each `N(w^(m))` is symmetric. -/
theorem A6_symmetric : A6mats.all isSym3 = true := by decide
/-- Each `N(w^(m))` is traceless: `tr N(w) = integral p_w div w = 0` (notes, A6 step (3)). -/
theorem A6_traceless : A6mats.all (fun m => decide (trace3 m = 0)) = true := by decide

/-- The five coordinate rows are exactly those printed in the notes. -/
theorem A6rows_value :
    A6rows = [[-6, 4, 5, 6, 5], [-10, 4, -1, 4, -1], [2, 4, -2, -6, -10], [-2, -2, -1, -2, -7],
      [-6, -10, 12, 10, 14]] := by decide

/-- Shape checks: `det`, `dotL` and `matMul` read missing entries as `0`, so the data layer
is pinned down explicitly -- five rows of five, and `3 x 3` for the matrices they come from. -/
theorem A6mats_shape :
    A6mats.all (fun m => decide (m.length = 3) && m.all (fun r => decide (r.length = 3))) = true := by decide
theorem A6rows_shape :
    A6rows.length = 5 ∧ A6rows.all (fun r => decide (r.length = 5)) = true := by decide

/-- The certificate proper: the `5 x 5` coordinate matrix has determinant `7936`. -/
theorem A6_det : det 5 A6rows = 7936 := by decide
theorem A6_det_ne_zero : det 5 A6rows ≠ 0 := by decide

/-! ### An explicit integer left inverse (invertibility without linear algebra)

`A6adj` is the adjugate of the coordinate matrix; `A6adj * A6rows = 7936 * I` exhibits
`(1/7936) A6adj` as a left inverse over the rationals, which is the content of "the
`N(w^(m))` span `Sym_0`" used in A6 step (4).  The step from this identity to spanning
over `Q` is a one-line argument that is NOT formalized here. -/

def dotL (u v : List Int) : Int :=
  (List.range u.length).foldl (fun acc i => acc + u.getD i 0 * v.getD i 0) 0
def colOf (m : List (List Int)) (j : Nat) : List Int := m.map (fun row => row.getD j 0)
def matMul (a b : List (List Int)) : List (List Int) :=
  a.map (fun row => (List.range (b.headD []).length).map (fun j => dotL row (colOf b j)))
def scaledId (c : Int) (n : Nat) : List (List Int) :=
  (List.range n).map (fun i => (List.range n).map (fun j => if i = j then c else 0))

def A6adj : List (List Int) :=
  [[2640, -2528, -2448, 2768, -1488],
   [632, -16, 424, -648, -248],
   [280, -208, 1544, -488, 744],
   [5408, -3904, -5664, 6560, -2976],
   [-2520, 1872, 1976, -3544, 1240]]

theorem A6adj_shape :
    A6adj.length = 5 ∧ A6adj.all (fun r => decide (r.length = 5)) = true := by decide

theorem A6_left_inverse : matMul A6adj A6rows = scaledId 7936 5 := by decide

/-! ## Axiom audit

Every theorem above is proved by `decide`, i.e. by kernel evaluation of the arithmetic.
Each line below prints "does not depend on any axioms": no `sorry`, no `native_decide`
(which would introduce `Lean.ofReduceBool`), and no tactic that pulls in `propext` or
`Quot.sound` (`omega` would). -/

#print axioms det_id2
#print axioms det_id3
#print axioms det_id5
#print axioms det_repeated_row
#print axioms det_swap
#print axioms det_2x2
#print axioms det_3x3
#print axioms triad_sums_to_zero
#print axioms triad_norms
#print axioms triadA_divFree
#print axioms triadB_divFree
#print axioms triad_nondegenerate
#print axioms triadA_transfers
#print axioms triadB_transfers
#print axioms triadA_transfer_sum
#print axioms triadB_transfer_sum
#print axioms triadA_transfer_conj
#print axioms triadB_transfer_conj
#print axioms triadA_proj_agrees
#print axioms triadB_proj_agrees
#print axioms triadMatrix_value
#print axioms triadMatrix_shape
#print axioms triad_det
#print axioms triad_det_ne_zero
#print axioms planar_triad_sums_to_zero
#print axioms planarA_divFree
#print axioms planarB_divFree
#print axioms planar_nondegenerate
#print axioms planarA_transfers
#print axioms planarB_transfers
#print axioms planarA_transfer_sum
#print axioms planarB_transfer_sum
#print axioms planarA_enstrophy_conserved
#print axioms planarB_enstrophy_conserved
#print axioms planarMatrix_value
#print axioms planarMatrix_shape
#print axioms planar_det_zero
#print axioms triadA_budget
#print axioms triadB_budget
#print axioms energy_flux_zero
#print axioms triadA_hhalf_coeffs
#print axioms triadB_hhalf_coeffs
#print axioms triadA_threshold
#print axioms triadA_above
#print axioms triadA_below
#print axioms triadB_threshold
#print axioms triadB_above
#print axioms triadB_below
#print axioms ens_threshold_rational
#print axioms ens_threshold_rational2
#print axioms l4_threshold_choice1
#print axioms l4_threshold_choice2
#print axioms l4RateA_below
#print axioms l4RateA_above
#print axioms A6_symmetric
#print axioms A6_traceless
#print axioms A6rows_value
#print axioms A6mats_shape
#print axioms A6rows_shape
#print axioms A6_det
#print axioms A6_det_ne_zero
#print axioms A6adj_shape
#print axioms A6_left_inverse

end NSCert
