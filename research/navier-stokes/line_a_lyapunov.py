"""Line A (coercive Lyapunov functionals) -- computational certificates for line_a_lyapunov.md.

Setting: unit torus T^3 = R^3/Z^3, Fourier modes exp(2*pi*i*n.x), viscosity nu = 1, no force.
Navier-Stokes: d_t u = Lap u + B(u,u),  B(u,v) = -P((u.grad) v),  P = Leray projector.
For a functional Q the "heat flux" is DQ(u)[Lap u], the "Euler flux" is DQ(u)[B(u,u)]; the
Navier-Stokes derivative of Q at t = 0 is their sum (Lemma A1, step 1 of the notes).

Parts:
  1. sympy verification of the pointwise identities behind the enstrophy / L^p / helicity budgets;
  2. exact rational (Fraction) spectral algebra: triad certificate (A2), enstrophy and L^4 fluxes,
     the zeroth-order local certificate (A6);
  3. FFT path (dealiased, exact for polynomial quantities): fluxes of all candidates on explicit
     Gaussian-integer trigonometric polynomials, positive controls (2D, Beltrami), L^inf via Danskin,
     the F(|u|) = |u|^2 log(1+|u|) scans (A5), and a small optimisation of the enstrophy threshold.

Runs in well under five minutes with `python3 line_a_lyapunov.py`; prints Markdown tables.
"""

from __future__ import annotations

import itertools
import time
from fractions import Fraction

import numpy as np
import sympy as sp
from scipy.optimize import minimize

TWO_PI = 2.0 * np.pi
SEED = 20260906

# ----------------------------------------------------------------------------------------------
# Part 1: sympy verification of pointwise identities
# ----------------------------------------------------------------------------------------------


def sympy_identities() -> list[tuple[str, bool]]:
    """Verify, for generic smooth u(x,y,z), the identities used in the flux computations."""
    x, y, z = sp.symbols("x y z", real=True)
    X = (x, y, z)
    u = [sp.Function(f"u{i}")(x, y, z) for i in range(3)]

    def d(f, j):
        return sp.diff(f, X[j])

    def curl(v):
        return [d(v[2], 1) - d(v[1], 2), d(v[0], 2) - d(v[2], 0), d(v[1], 0) - d(v[0], 1)]

    def adv(a, b):  # (a.grad) b
        return [sum(a[j] * d(b[i], j) for j in range(3)) for i in range(3)]

    def dot(a, b):
        return sum(a[i] * b[i] for i in range(3))

    def lap(f):
        return sum(d(d(f, j), j) for j in range(3))

    def div(v):
        return sum(d(v[j], j) for j in range(3))

    w = curl(u)
    w2 = dot(w, w)
    S = [[(d(u[i], j) + d(u[j], i)) / 2 for j in range(3)] for i in range(3)]
    checks = []
    # (I1) w.((u.grad)w) = (1/2) u.grad|w|^2
    checks.append(
        (
            "I1  w.(u.grad)w = (1/2) u.grad|w|^2",
            sp.simplify(dot(w, adv(u, w)) - sum(u[j] * d(w2, j) for j in range(3)) / 2) == 0,
        )
    )
    # (I2) (1/2) u.grad|w|^2 = div(u |w|^2/2) - (|w|^2/2) div u
    lhs = sum(u[j] * d(w2, j) for j in range(3)) / 2
    rhs = div([u[j] * w2 / 2 for j in range(3)]) - (w2 / 2) * div(u)
    checks.append(("I2  (1/2) u.grad|w|^2 = div(u|w|^2/2) - (|w|^2/2) div u", sp.simplify(lhs - rhs) == 0))
    # (I3) curl((u.grad)u) = (u.grad)w - (w.grad)u + w div u   (vorticity equation)
    c = curl(adv(u, u))
    r = [adv(u, w)[i] - adv(w, u)[i] + w[i] * div(u) for i in range(3)]
    checks.append(
        ("I3  curl((u.grad)u) = (u.grad)w - (w.grad)u + w div u", all(sp.simplify(c[i] - r[i]) == 0 for i in range(3)))
    )
    # (I4) w.(w.grad)u = w.S w
    checks.append(
        (
            "I4  w.(w.grad)u = w.Sw",
            sp.simplify(dot(w, adv(w, u)) - sum(w[i] * S[i][j] * w[j] for i in range(3) for j in range(3))) == 0,
        )
    )
    # (I5) w.Lap w = Lap(|w|^2/2) - |grad w|^2
    gw2 = sum(d(w[i], j) ** 2 for i in range(3) for j in range(3))
    checks.append(
        (
            "I5  w.Lap w = Lap(|w|^2/2) - |grad w|^2",
            sp.simplify(dot(w, [lap(w[i]) for i in range(3)]) - (lap(w2 / 2) - gw2)) == 0,
        )
    )
    # (I6) div((u.grad)u) = d_i u_j d_j u_i + u.grad(div u)   (pressure equation)
    lhs = div(adv(u, u))
    rhs = sum(d(u[j], i) * d(u[i], j) for i in range(3) for j in range(3)) + sum(u[j] * d(div(u), j) for j in range(3))
    checks.append(("I6  div((u.grad)u) = d_i u_j d_j u_i + u.grad(div u)", sp.simplify(lhs - rhs) == 0))
    # (I7) (u.grad)u = grad(|u|^2/2) - u x w   (used for helicity and Beltrami fields)
    u2 = dot(u, u)
    uxw = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
    checks.append(
        (
            "I7  (u.grad)u = grad(|u|^2/2) - u x w",
            all(sp.simplify(adv(u, u)[i] - (d(u2, i) / 2 - uxw[i])) == 0 for i in range(3)),
        )
    )
    # (I8) w.(u x w) = 0
    checks.append(("I8  w.(u x w) = 0", sp.simplify(dot(w, uxw)) == 0))
    # (I9) for phi(v) = F(|v|): traceless part of the Hessian is (F'' - F'/r)(vv^T/r^2 - I/3)
    v = sp.symbols("v0:3", real=True)
    rr = sp.symbols("rr", positive=True)
    F = sp.Function("F")
    r = sp.sqrt(sum(vi**2 for vi in v))
    H = sp.Matrix(3, 3, lambda i, j: sp.diff(F(r), v[i], v[j]))
    H0 = (H - (H.trace() / 3) * sp.eye(3)).applyfunc(lambda e: e.doit())
    F1, F2 = sp.diff(F(rr), rr), sp.diff(F(rr), rr, 2)
    vv = sp.Matrix(v)
    target = ((F2 - F1 / rr) * (vv * vv.T / rr**2 - sp.eye(3) / 3)).subs(rr, r)
    diff = (H0 - target).applyfunc(sp.simplify)
    checks.append(("I9  Hess F(|v|) traceless part = (F''-F'/r)(vv^T/r^2 - I/3)", diff == sp.zeros(3, 3)))
    # (I10) flux of a local functional: grad(phi(u)).(u.grad)u = (u.grad)(phi(u)) for phi = F(|u|)
    phi_u = F(sp.sqrt(dot(u, u)))
    gradphi = [sp.diff(phi_u, u[i]) for i in range(3)]
    lhs = sum(gradphi[i] * adv(u, u)[i] for i in range(3))
    rhs = sum(u[j] * d(phi_u, j) for j in range(3))
    checks.append(("I10 grad(F(|u|)).(u.grad)u = (u.grad)F(|u|)", sp.simplify((lhs - rhs).doit()) == 0))
    return checks


# ----------------------------------------------------------------------------------------------
# Part 2: exact rational spectral algebra
# ----------------------------------------------------------------------------------------------


class CQ:
    """Complex number with Fraction real and imaginary parts."""

    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re = Fraction(re)
        self.im = Fraction(im)

    def __add__(self, o):
        return CQ(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        return CQ(self.re - o.re, self.im - o.im)

    def __mul__(self, o):
        if isinstance(o, CQ):
            return CQ(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
        return CQ(self.re * o, self.im * o)

    def __neg__(self):
        return CQ(-self.re, -self.im)

    def conj(self):
        return CQ(self.re, -self.im)

    def is_zero(self):
        return self.re == 0 and self.im == 0

    def __repr__(self):
        return f"({self.re}{'+' if self.im >= 0 else '-'}{abs(self.im)}i)"


class ExactField:
    """Real trigonometric polynomial on T^3 with complex-rational coefficients.

    coef: dict wavevector n (tuple of ints) -> list of CQ (ncomp entries). The physical field is
    sum_n coef[n] exp(2 pi i n.x) times (2 pi)^pi_pow; pi_pow counts derivatives taken.
    """

    def __init__(self, ncomp: int, coef=None, pi_pow: int = 0):
        self.ncomp = ncomp
        self.coef = {} if coef is None else coef
        self.pi_pow = pi_pow

    def copy(self):
        return ExactField(self.ncomp, {n: list(v) for n, v in self.coef.items()}, self.pi_pow)

    @staticmethod
    def from_modes(modes: dict) -> ExactField:
        """modes: n -> list of CQ (coefficient of exp(2 pi i n.x)); conjugate modes are added."""
        coef = {}
        for n, a in modes.items():
            coef.setdefault(n, [CQ() for _ in a])
            coef[n] = [coef[n][i] + a[i] for i in range(len(a))]
            m = tuple(-c for c in n)
            coef.setdefault(m, [CQ() for _ in a])
            coef[m] = [coef[m][i] + a[i].conj() for i in range(len(a))]
        ncomp = len(next(iter(modes.values())))
        f = ExactField(ncomp, coef, 0)
        f.prune()
        return f

    def prune(self):
        self.coef = {n: v for n, v in self.coef.items() if not all(c.is_zero() for c in v)}
        return self

    def add(self, other, scale=1):
        assert self.pi_pow == other.pi_pow and self.ncomp == other.ncomp
        out = self.copy()
        for n, v in other.coef.items():
            cur = out.coef.setdefault(n, [CQ() for _ in range(self.ncomp)])
            out.coef[n] = [cur[i] + v[i] * Fraction(scale) for i in range(self.ncomp)]
        return out.prune()

    def scale(self, s):
        return ExactField(self.ncomp, {n: [c * Fraction(s) for c in v] for n, v in self.coef.items()}, self.pi_pow)

    def deriv(self, j):
        out = {}
        for n, v in self.coef.items():
            f = CQ(0, n[j])  # symbol of d_j / (2 pi) is i n_j
            out[n] = [c * f for c in v]
        return ExactField(self.ncomp, out, self.pi_pow + 1).prune()

    def laplacian(self):
        out = {n: [c * (-(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)) for c in v] for n, v in self.coef.items()}
        return ExactField(self.ncomp, out, self.pi_pow + 2).prune()

    def component(self, i):
        return ExactField(1, {n: [v[i]] for n, v in self.coef.items()}, self.pi_pow).prune()

    def curl(self):
        assert self.ncomp == 3
        d = [[self.component(i).deriv(j) for j in range(3)] for i in range(3)]  # d[i][j] = d_j u_i
        out = {}
        comps = [d[2][1].add(d[1][2], -1), d[0][2].add(d[2][0], -1), d[1][0].add(d[0][1], -1)]
        for i, c in enumerate(comps):
            for n, v in c.coef.items():
                out.setdefault(n, [CQ() for _ in range(3)])[i] = v[0]
        return ExactField(3, out, self.pi_pow + 1).prune()

    def divergence(self):
        out = self.component(0).deriv(0)
        for j in (1, 2):
            out = out.add(self.component(j).deriv(j))
        return out

    @staticmethod
    def product(a: ExactField, b: ExactField) -> ExactField:
        """Pointwise product of two scalar fields (convolution of coefficients)."""
        assert a.ncomp == 1 and b.ncomp == 1
        out = {}
        for n1, v1 in a.coef.items():
            for n2, v2 in b.coef.items():
                n = (n1[0] + n2[0], n1[1] + n2[1], n1[2] + n2[2])
                cur = out.setdefault(n, [CQ()])
                cur[0] = cur[0] + v1[0] * v2[0]
        return ExactField(1, out, a.pi_pow + b.pi_pow).prune()

    @staticmethod
    def dot_field(a: ExactField, b: ExactField) -> ExactField:
        """Scalar field a.b."""
        out = ExactField.product(a.component(0), b.component(0))
        for i in (1, 2):
            out = out.add(ExactField.product(a.component(i), b.component(i)))
        return out

    @staticmethod
    def scalar_times_vector(s: ExactField, v: ExactField) -> ExactField:
        comps = [ExactField.product(s, v.component(i)) for i in range(3)]
        out = {}
        for i, c in enumerate(comps):
            for n, val in c.coef.items():
                out.setdefault(n, [CQ() for _ in range(3)])[i] = val[0]
        return ExactField(3, out, s.pi_pow + v.pi_pow).prune()

    @staticmethod
    def advect(a: ExactField, b: ExactField) -> ExactField:
        """(a.grad) b for vector fields a, b."""
        comps = []
        for i in range(3):
            acc = ExactField.product(a.component(0), b.component(i).deriv(0))
            for j in (1, 2):
                acc = acc.add(ExactField.product(a.component(j), b.component(i).deriv(j)))
            comps.append(acc)
        out = {}
        for i, c in enumerate(comps):
            for n, val in c.coef.items():
                out.setdefault(n, [CQ() for _ in range(3)])[i] = val[0]
        return ExactField(3, out, a.pi_pow + b.pi_pow + 1).prune()

    def leray(self):
        out = {}
        for n, v in self.coef.items():
            n2 = n[0] ** 2 + n[1] ** 2 + n[2] ** 2
            if n2 == 0:
                out[n] = list(v)
                continue
            nv = CQ()
            for i in range(3):
                nv = nv + v[i] * n[i]
            out[n] = [v[i] - nv * Fraction(n[i], n2) for i in range(3)]
        return ExactField(3, out, self.pi_pow).prune()

    @staticmethod
    def pressure(u: ExactField) -> ExactField:
        """p with Lap p = -d_i d_j (u_i u_j), mean zero; p_hat_n = -(n_i n_j/|n|^2)(u_i u_j)^_n."""
        out = {}
        for i in range(3):
            for j in range(3):
                pij = ExactField.product(u.component(i), u.component(j))
                for n, v in pij.coef.items():
                    n2 = n[0] ** 2 + n[1] ** 2 + n[2] ** 2
                    if n2 == 0:
                        continue
                    cur = out.setdefault(n, [CQ()])
                    cur[0] = cur[0] - v[0] * Fraction(n[i] * n[j], n2)
        return ExactField(1, out, u.pi_pow * 2).prune()

    @staticmethod
    def euler_B(u: ExactField) -> ExactField:
        """B(u,u) = -P((u.grad)u)."""
        return ExactField.advect(u, u).leray().scale(-1)

    @staticmethod
    def inner(a: ExactField, b: ExactField) -> tuple[Fraction, int]:
        """Integral over T^3 of a.b: returns (rational, power of 2 pi). Imaginary part must vanish."""
        tot = CQ()
        for n, va in a.coef.items():
            vb = b.coef.get(n)
            if vb is None:
                continue
            for i in range(a.ncomp):
                tot = tot + va[i].conj() * vb[i]
        assert tot.im == 0, tot
        return tot.re, a.pi_pow + b.pi_pow

    def energy_transfer(self, B: ExactField) -> dict:
        """T_n = Re(conj(u_n).B_n) for each mode n of u (rational part; power of 2 pi = B.pi_pow)."""
        out = {}
        for n, v in self.coef.items():
            b = B.coef.get(n)
            if b is None:
                out[n] = Fraction(0)
                continue
            s = CQ()
            for i in range(3):
                s = s + v[i].conj() * b[i]
            out[n] = s.re
        return out


def perp_basis(n):
    """Two integer vectors spanning n^perp (n a nonzero integer vector)."""
    n = np.array(n)
    cands = []
    for e in np.eye(3, dtype=int):
        c = np.cross(n, e)
        if np.any(c != 0):
            cands.append(c // np.gcd.reduce(c[c != 0]))
    b1 = cands[0]
    for c in cands[1:]:
        if np.any(np.cross(b1, c) != 0):
            b2 = c
            break
    return tuple(int(t) for t in b1), tuple(int(t) for t in b2)


def gaussian_integer_field(modes, coeffs) -> ExactField:
    """Real field sum_n 2 Re[(alpha_n b_n^1 + beta_n b_n^2) e^{2 pi i n.x}] with Gaussian-integer alpha, beta."""
    d = {}
    for n, (alpha, beta) in zip(modes, coeffs):
        b1, b2 = perp_basis(n)
        a = [CQ(alpha[0], alpha[1]) * b1[i] + CQ(beta[0], beta[1]) * b2[i] for i in range(3)]
        d[n] = a
    return ExactField.from_modes(d)


def triad_nondegenerate(k, p, q) -> bool:
    """True iff the only zero-sum unordered triples from {+-k,+-p,+-q} are +-(k,p,q)."""
    vecs = [np.array(v) for v in (k, p, q)]
    S = vecs + [-v for v in vecs]
    zero_sum = set()
    for i, j, m in itertools.combinations_with_replacement(range(6), 3):
        if np.all(S[i] + S[j] + S[m] == 0):
            zero_sum.add(tuple(sorted((i, j, m))))
    return zero_sum == {(0, 1, 2), (3, 4, 5)}


def exact_part():
    lines = []
    k, p, q = (1, 0, 0), (0, 1, 1), (-1, -1, -1)
    lines.append(
        f"Triad k={k}, p={p}, q={q}: non-degenerate (only +-(k,p,q) sum to zero): {triad_nondegenerate(k, p, q)}"
    )
    kp = (1, 0, 0), (0, 1, 0), (-1, -1, 0)
    lines.append(f"Planar triad {kp}: non-degenerate: {triad_nondegenerate(*kp)}")

    # two amplitude choices on the 3D triad
    choices = [
        {k: [CQ(0), CQ(1), CQ(0, 1)], p: [CQ(1), CQ(0, 1), CQ(0, -1)], q: [CQ(1), CQ(-1), CQ(0)]},
        {k: [CQ(0), CQ(1, 1), CQ(-1)], p: [CQ(1, -1), CQ(1), CQ(-1)], q: [CQ(1), CQ(1), CQ(-2)]},
    ]
    rows = []
    for c in choices:
        for n, a in c.items():  # divergence-free check
            assert all((a[0] * n[0] + a[1] * n[1] + a[2] * n[2]).is_zero() for _ in [0])
        u = ExactField.from_modes(c)
        assert not u.divergence().coef, "not divergence free"
        B = ExactField.euler_B(u)
        T = u.energy_transfer(B)
        Tk, Tp, Tq = T[k], T[p], T[q]
        Tmk, Tmp, Tmq = T[tuple(-t for t in k)], T[tuple(-t for t in p)], T[tuple(-t for t in q)]
        assert (Tk, Tp, Tq) == (Tmk, Tmp, Tmq)
        rows.append((Tk, Tp, Tq))
    det = rows[0][1] * rows[1][2] - rows[0][2] * rows[1][1]
    lines.append(
        "Energy transfers T_n = Re(conj(u_n).B(u,u)_n) on the 3D triad, units of 2*pi (each mode and its conjugate):"
    )
    for i, (Tk, Tp, Tq) in enumerate(rows):
        lines.append(f"  choice {i + 1}: T_k = {Tk}, T_p = {Tp}, T_q = {Tq}, sum = {Tk + Tp + Tq}")
    lines.append(f"  det[(T_p,T_q) choice 1; (T_p,T_q) choice 2] = {det}  (nonzero => (T_p,T_q) span R^2)")

    # enstrophy Euler flux two ways, heat flux, threshold, for choice 1 and 2
    for i, c in enumerate(choices):
        u = ExactField.from_modes(c)
        B = ExactField.euler_B(u)
        w = u.curl()
        # Euler flux of enstrophy: 2 <-Lap u, B> = 2 int w.(w.grad)u
        mlap = u.laplacian().scale(-1)
        f1, pw1 = ExactField.inner(mlap, B)
        f2, pw2 = ExactField.inner(w, ExactField.advect(w, u))
        assert pw1 == pw2 == 3 and f1 == f2, (f1, f2, pw1, pw2)
        g, pwg = ExactField.inner(mlap, u.laplacian())
        assert pwg == 4
        ens, pwe = ExactField.inner(w, w)
        en, _ = ExactField.inner(u, u)
        F = 2 * f1
        G = 2 * g
        sgn = "u" if F > 0 else "-u"
        thr = -Fraction(G) / abs(Fraction(F)) if F != 0 else None
        lines.append(
            f"choice {i + 1}: energy = {en}, enstrophy = {ens}*(2pi)^2; enstrophy Euler flux F = {F}*(2pi)^3 (both formulas agree), "
            f"heat flux G = {G}*(2pi)^4; for the field {sgn}: threshold amplitude A* = -G/|F| = {thr}*(2pi) = {float(thr) * TWO_PI:.6f}"
            if thr is not None
            else f"choice {i + 1}: enstrophy Euler flux vanishes"
        )
        # L^4: Q = int |u|^4, DQ[h] = 4 int |u|^2 u.h ; Euler flux two ways
        u2 = ExactField.dot_field(u, u)
        u2u = ExactField.scalar_times_vector(u2, u)
        f4a, pw4a = ExactField.inner(u2u, B)
        pr = ExactField.pressure(u)
        gradp = ExactField(3, {}, 0)
        gp = {}
        for j in range(3):
            dj = pr.deriv(j)
            for n, v in dj.coef.items():
                gp.setdefault(n, [CQ() for _ in range(3)])[j] = v[0]
        gradp = ExactField(3, gp, pr.pi_pow + 1)
        f4b, pw4b = ExactField.inner(u2u, gradp)
        assert pw4a == pw4b == 1 and f4a == -f4b, (f4a, f4b)
        g4, pwg4 = ExactField.inner(u2u, u.laplacian())
        F4, G4 = 4 * f4a, 4 * g4
        sgn4 = "u" if F4 > 0 else "-u"
        thr4 = -Fraction(G4) / abs(Fraction(F4)) if F4 != 0 else None
        l4, _ = ExactField.inner(u2, u2)
        lines.append(
            f"choice {i + 1}: int|u|^4 = {l4}; L^4 Euler flux F = {F4}*(2pi) (= -4 int |u|^2 u.grad p, agrees), heat flux G = {G4}*(2pi)^2; "
            f"for the field {sgn4}: A* = -G/|F| = {thr4}*(2pi) = {float(thr4) * TWO_PI:.6f}"
            if thr4 is not None
            else f"choice {i + 1}: L^4 Euler flux vanishes"
        )
        # H^{1/2}: F = 2 sum_n (2 pi |n|) T_n = (2 pi)^2 * 2 * sum |n| T_n  (T_n carries one factor 2 pi)
        T = u.energy_transfer(B)
        coeff = {}
        for n, t in T.items():
            n2 = n[0] ** 2 + n[1] ** 2 + n[2] ** 2
            coeff[n2] = coeff.get(n2, Fraction(0)) + 2 * t
        val = sum(float(cf) * np.sqrt(n2) for n2, cf in coeff.items()) * TWO_PI**2
        lines.append(
            f"choice {i + 1}: Hdot^(1/2) Euler flux = (2pi)^2 * ["
            + " + ".join(f"({cf})*sqrt({n2})" for n2, cf in sorted(coeff.items()))
            + f"] = {val:.6f}"
        )

    # planar control: the 2D triad; (T_p,T_q) must be confined to a line
    kp, pp, qp = kp
    choices2d = [
        {kp: [CQ(0), CQ(1), CQ(0)], pp: [CQ(0, 1), CQ(0), CQ(0)], qp: [CQ(1), CQ(-1), CQ(0)]},
        {kp: [CQ(0), CQ(2, -1), CQ(0)], pp: [CQ(1, 1), CQ(0), CQ(0)], qp: [CQ(1, 2), CQ(-1, -2), CQ(0)]},
    ]
    rows2 = []
    for c in choices2d:
        u = ExactField.from_modes(c)
        assert not u.divergence().coef
        B = ExactField.euler_B(u)
        T = u.energy_transfer(B)
        rows2.append((T[kp], T[pp], T[qp]))
    det2 = rows2[0][1] * rows2[1][2] - rows2[0][2] * rows2[1][1]
    lines.append("Planar triad (2D control): transfers T_k, T_p, T_q, energy sum, enstrophy sum |n|^2 T_n:")
    for i, (Tk, Tp, Tq) in enumerate(rows2):
        lines.append(
            f"  choice {i + 1}: {Tk}, {Tp}, {Tq}; sum = {Tk + Tp + Tq}; enstrophy-weighted sum = {1 * Tk + 1 * Tp + 2 * Tq}"
        )
    lines.append(
        f"  det[(T_p,T_q) choice 1; choice 2] = {det2}  (zero => confined to a line, as 2D enstrophy conservation forces)"
    )

    # A6 certificate: N(w) = sym part of M_ij = int p_w d_j w_i, for five Gaussian-integer fields
    rng = np.random.default_rng(SEED)
    modes = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, -1, 0), (1, 0, -1), (0, 1, -1)]
    best = None
    for _trial in range(40):
        fields = []
        for _m in range(5):
            coeffs = [
                (
                    (int(rng.integers(-1, 2)), int(rng.integers(-1, 2))),
                    (int(rng.integers(-1, 2)), int(rng.integers(-1, 2))),
                )
                for _ in modes
            ]
            fields.append((coeffs, gaussian_integer_field(modes, coeffs)))
        mat = []
        for _coeffs, w in fields:
            assert not w.divergence().coef
            pr = ExactField.pressure(w)
            M = [[ExactField.inner(pr, w.component(i).deriv(j))[0] for j in range(3)] for i in range(3)]
            N = [[(M[i][j] + M[j][i]) / 2 for j in range(3)] for i in range(3)]
            assert N[0][0] + N[1][1] + N[2][2] == 0
            mat.append([N[0][0], N[1][1], N[0][1], N[0][2], N[1][2]])
        detN = sp.Matrix([[sp.Rational(v.numerator, v.denominator) for v in row] for row in mat]).det()
        if detN != 0:
            best = (fields, mat, detN)
            break
    assert best is not None
    fields, mat, detN = best
    lines.append(
        "A6 certificate: five Gaussian-integer fields w^(m) on modes "
        + str(modes)
        + " (coefficients (alpha,beta) per mode, alpha,beta in Z[i] written as (re,im)):"
    )
    for m, (coeffs, _w) in enumerate(fields):
        lines.append(f"  w^({m + 1}): " + " ".join(f"{c}" for c in coeffs))
    lines.append("  rows (N11, N22, N12, N13, N23) of N(w^(m)) = sym part of int p_w d_j w_i dx, units of 2*pi:")
    for m, row in enumerate(mat):
        lines.append(f"  N(w^({m + 1})) = ({', '.join(str(v) for v in row)})")
    lines.append(
        f"  det of the 5x5 coordinate matrix = {detN} (nonzero => the N(w) span all traceless symmetric matrices)"
    )
    return lines


# ----------------------------------------------------------------------------------------------
# Part 3: FFT machinery (numpy), exact for polynomial quantities when N is large enough
# ----------------------------------------------------------------------------------------------


class Spectral:
    def __init__(self, N: int):
        self.N = N
        n = np.fft.fftfreq(N, 1.0 / N).astype(int)
        self.n = np.stack(np.meshgrid(n, n, n, indexing="ij"))  # (3, N, N, N) integer wavevectors
        self.n2 = np.sum(self.n**2, axis=0)
        self.absn = np.sqrt(self.n2)
        self.ik = 2j * np.pi * self.n  # symbol of grad

    def coef_from_modes(self, modes: dict) -> np.ndarray:
        """modes: n -> complex 3-vector (coefficient of e^{2 pi i n.x}); returns (3,N,N,N) coefficients of a real field."""
        c = np.zeros((3, self.N, self.N, self.N), dtype=complex)
        for n, a in modes.items():
            a = np.asarray(a, dtype=complex)
            c[(slice(None),) + tuple(int(t) % self.N for t in n)] += a
            c[(slice(None),) + tuple((-int(t)) % self.N for t in n)] += np.conj(a)
        return c

    def phys(self, c):  # coefficients -> grid samples (real field)
        return np.real(np.fft.ifftn(c, axes=(-3, -2, -1)) * self.N**3)

    def spec(self, f):  # grid samples -> coefficients
        return np.fft.fftn(f, axes=(-3, -2, -1)) / self.N**3

    def grad(self, c_scalar):
        return self.ik * c_scalar[None]

    def laplacian(self, c):
        return -(TWO_PI**2) * self.n2 * c

    def curl(self, c):
        ik = self.ik
        return np.stack([ik[1] * c[2] - ik[2] * c[1], ik[2] * c[0] - ik[0] * c[2], ik[0] * c[1] - ik[1] * c[0]])

    def leray(self, c):
        n2 = np.where(self.n2 == 0, 1, self.n2)
        nv = np.sum(self.n * c, axis=0)
        out = c - self.n * (nv / n2)[None]
        return out

    def advect(self, ca, cb):
        """coefficients of (a.grad) b (exact if N > 2*K_max + ... i.e. products resolved)."""
        a = self.phys(ca)
        out = np.zeros_like(a)
        for i in range(3):
            for j in range(3):
                out[i] += a[j] * self.phys(self.ik[j] * cb[i])
        return self.spec(out)

    def B(self, c):
        return -self.leray(self.advect(c, c))

    def pressure(self, c):
        u = self.phys(c)
        n2 = np.where(self.n2 == 0, 1, self.n2)
        p = np.zeros((self.N,) * 3, dtype=complex)
        for i in range(3):
            for j in range(3):
                p -= self.n[i] * self.n[j] / n2 * self.spec(u[i] * u[j])
        p[0, 0, 0] = 0
        return p

    def integral(self, f):
        return float(np.mean(f))

    def inner(self, ca, cb):
        return float(np.real(np.sum(np.conj(ca) * cb)))


def field_modes_from_gaussian(modes, coeffs) -> dict:
    d = {}
    for n, (alpha, beta) in zip(modes, coeffs):
        b1, b2 = perp_basis(n)
        a = complex(alpha[0], alpha[1]) * np.array(b1) + complex(beta[0], beta[1]) * np.array(b2)
        d[n] = a
    return d


def fluxes(sp_: Spectral, c: np.ndarray, which: str, N_quad: int | None = None) -> tuple[float, float, float]:
    """Return (Q(u), heat flux DQ(u)[Lap u], Euler flux DQ(u)[B(u,u)]) for candidate `which`.

    Polynomial candidates use the exact-N grid `sp_`; non-polynomial ones re-sample on N_quad.
    """
    lap = sp_.laplacian(c)
    B = sp_.B(c)
    if which == "energy":
        return sp_.inner(c, c), 2 * sp_.inner(c, lap), 2 * sp_.inner(c, B)
    if which == "enstrophy":
        m = TWO_PI**2 * sp_.n2
        return sp_.inner(m * c, c), 2 * sp_.inner(m * c, lap), 2 * sp_.inner(m * c, B)
    if which == "H12":
        m = TWO_PI * sp_.absn
        return sp_.inner(m * c, c), 2 * sp_.inner(m * c, lap), 2 * sp_.inner(m * c, B)
    if which == "helicity":
        w = sp_.curl(c)
        return sp_.inner(c, w), 2 * sp_.inner(w, lap), 2 * sp_.inner(w, B)
    if which == "local_energy":  # int |u|^2 chi, chi = 1 + cos(2 pi x_1)
        xg = np.arange(sp_.N) / sp_.N
        chi = 1.0 + np.cos(TWO_PI * xg)[:, None, None]
        u = sp_.phys(c)
        return (
            sp_.integral(chi * np.sum(u**2, 0)),
            2 * sp_.integral(chi * np.sum(u * sp_.phys(lap), 0)),
            2 * sp_.integral(chi * np.sum(u * sp_.phys(B), 0)),
        )
    if which == "L4":
        u = sp_.phys(c)
        u2 = np.sum(u**2, axis=0)
        wgt = u2[None] * u
        return (
            sp_.integral(u2**2),
            4 * sp_.integral(np.sum(wgt * sp_.phys(lap), 0)),
            4 * sp_.integral(np.sum(wgt * sp_.phys(B), 0)),
        )
    # non-polynomial candidates: resample everything on a finer grid
    sq = Spectral(N_quad) if N_quad else sp_

    def upsample(cc):
        out = np.zeros((cc.shape[0], sq.N, sq.N, sq.N), dtype=complex)
        K = sp_.N // 2
        for idx in itertools.product(range(-K, K + 1), repeat=3):
            src = tuple(t % sp_.N for t in idx)
            dst = tuple(t % sq.N for t in idx)
            out[(slice(None),) + dst] = cc[(slice(None),) + src]
        return out

    cU, lapU, Bfield = upsample(c), upsample(lap), upsample(B)
    u, lu, bphys = sq.phys(cU), sq.phys(lapU), sq.phys(Bfield)
    if which == "L3":
        au = np.sqrt(np.sum(u**2, 0))
        return sq.integral(au**3), 3 * sq.integral(au * np.sum(u * lu, 0)), 3 * sq.integral(au * np.sum(u * bphys, 0))
    if which == "Flog":  # int |u|^2 log(1+|u|)
        au = np.sqrt(np.sum(u**2, 0))
        g = 2 * np.log1p(au) + au / (1 + au)  # F'(r)/r
        return (
            sq.integral(au**2 * np.log1p(au)),
            sq.integral(g * np.sum(u * lu, 0)),
            sq.integral(g * np.sum(u * bphys, 0)),
        )
    if which == "omega32":  # int |w|^{3/2}
        w, lw, bw = sq.phys(sq.curl(cU)), sq.phys(sq.curl(lapU)), sq.phys(sq.curl(Bfield))
        aw = np.sqrt(np.sum(w**2, 0))
        awm = np.where(aw > 0, aw, 1.0) ** (-0.5) * (aw > 0)
        return (
            sq.integral(aw**1.5),
            1.5 * sq.integral(awm * np.sum(w * lw, 0)),
            1.5 * sq.integral(awm * np.sum(w * bw, 0)),
        )
    raise ValueError(which)


def eval_modes(modes: dict, x: np.ndarray, deriv: str = "") -> np.ndarray:
    """Evaluate a real field given by modes at point x; deriv = "" (value) or "lap" (Laplacian)."""
    val = np.zeros(3, dtype=complex)
    for n, a in modes.items():
        n = np.asarray(n, dtype=float)
        ph = np.exp(2j * np.pi * np.dot(n, x))
        fac = 1.0
        if deriv == "lap":
            fac = -(TWO_PI**2) * np.dot(n, n)
        val += fac * (np.asarray(a) * ph + np.conj(np.asarray(a)) * np.conj(ph))
    return np.real(val)


def grad_modes(modes: dict, x: np.ndarray) -> np.ndarray:
    """Matrix g[j, i] = d_j f_i (x) for the real field f given by modes."""
    g = np.zeros((3, 3), dtype=complex)
    for n, a in modes.items():
        nn = np.asarray(n, dtype=float)
        ph = np.exp(2j * np.pi * np.dot(nn, x))
        g += np.outer(2j * np.pi * nn, np.asarray(a) * ph) + np.outer(
            -2j * np.pi * nn, np.conj(np.asarray(a)) * np.conj(ph)
        )
    return np.real(g)


def sup_norm_fluxes(sp_: Spectral, modes: dict, kind: str) -> dict:
    """Danskin fluxes for Q = ||u||_inf (kind='u') or ||omega||_inf (kind='w')."""
    c = sp_.coef_from_modes(modes)
    if kind == "w":
        cw = sp_.curl(c)
        target_modes = {n: cw[(slice(None),) + tuple(int(t) % sp_.N for t in n)] for n in modes}
    else:
        target_modes = modes
    fine = Spectral(48)
    f = fine.phys(fine.coef_from_modes(target_modes))
    mag = np.sqrt(np.sum(f**2, 0))
    idx = np.unravel_index(np.argmax(mag), mag.shape)
    x0 = np.array(idx) / fine.N

    def negsq(x):
        return -np.sum(eval_modes(target_modes, x) ** 2)

    def negsq_grad(x):
        return -2.0 * grad_modes(target_modes, x) @ eval_modes(target_modes, x)

    res = minimize(negsq, x0, jac=negsq_grad, method="BFGS", options={"gtol": 1e-10})
    xs = res.x % 1.0
    gradnorm = float(np.linalg.norm(negsq_grad(xs)))
    # second-best local maximum on the grid (uniqueness check)
    loc = np.ones_like(mag, dtype=bool)
    for ax in range(3):
        for sh in (1, -1):
            loc &= mag >= np.roll(mag, sh, axis=ax)
    vals = np.sort(mag[loc])[::-1]
    second = vals[1] if len(vals) > 1 else float("nan")
    fx = eval_modes(target_modes, xs)
    mag_s = np.linalg.norm(fx)
    lap_s = eval_modes(target_modes, xs, "lap")
    G = np.dot(fx, lap_s) / mag_s
    if kind == "u":
        p = sp_.pressure(c)
        nz = np.argwhere(np.abs(p) > 1e-14)
        gradp = np.zeros(3, dtype=complex)
        for ii in nz:
            n = np.array([sp_.n[a][tuple(ii)] for a in range(3)], dtype=float)
            gradp += 2j * np.pi * n * p[tuple(ii)] * np.exp(2j * np.pi * np.dot(n, xs))
        F = -np.dot(fx, np.real(gradp)) / mag_s
    else:
        gu = grad_modes(modes, xs)  # gu[j, i] = d_j u_i
        S = (gu + gu.T) / 2
        F = fx @ S @ fx / mag_s
    return {"x*": xs, "max": mag_s, "second_local_max": second, "G": G, "F": F, "gradnorm": gradnorm}


def fft_part():
    lines = []
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    modes = [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 1, 0),
        (1, 0, 1),
        (0, 1, 1),
        (1, -1, 0),
        (1, 0, -1),
        (0, 1, -1),
        (1, 1, 1),
        (1, 1, -1),
        (1, -1, 1),
        (-1, 1, 1),
    ]
    N_exact = 16  # K = 1 in each component here (|n|_inf <= 1): quartic products have degree 5 < 16
    S16 = Spectral(N_exact)

    # --- consistency checks against the exact part (choice 1 of the triad) ---
    k, p, q = (1, 0, 0), (0, 1, 1), (-1, -1, -1)
    tri = {k: [0, 1, 1j], p: [1, 1j, -1j], q: [1, -1, 0]}
    c = S16.coef_from_modes(tri)
    _, G, F = fluxes(S16, c, "enstrophy")
    lines.append(f"FFT cross-check (triad choice 1): enstrophy heat flux G = {G:.10f}, Euler flux F = {F:.10f}")
    _, G4, F4 = fluxes(S16, c, "L4")
    lines.append(f"FFT cross-check (triad choice 1): L^4 heat flux G = {G4:.10f}, Euler flux F = {F4:.10f}")
    _, Gh, Fh = fluxes(S16, c, "H12")
    lines.append(f"FFT cross-check (triad choice 1): Hdot^(1/2) heat flux G = {Gh:.10f}, Euler flux F = {Fh:.10f}")

    # --- helicity: ABC flow (curl u = 2 pi u) and its mirror image ---
    abc = {(0, 0, 1): [-0.5j, 0.5, 0], (1, 0, 0): [0, -0.5j, 0.5], (0, 1, 0): [0.5, 0, -0.5j]}
    # u = (sin 2pi z + cos 2pi y, sin 2pi x + cos 2pi z, sin 2pi y + cos 2pi x): sin t = (e^{it}-e^{-it})/(2i) -> coef -i/2 at +n
    c = S16.coef_from_modes(abc)
    w = S16.curl(c)
    lines.append(f"ABC field: max |curl u - 2 pi u| over coefficients = {np.max(np.abs(w - TWO_PI * c)):.2e}")
    Qh, Gh, Fh = fluxes(S16, c, "helicity")
    lines.append(
        f"ABC: helicity = {Qh:.6f} (= 2 pi * 3 = {TWO_PI * 3:.6f}), heat flux = {Gh:.6f} (= -48 pi^3 = {-48 * np.pi**3:.6f}), Euler flux = {Fh:.2e}"
    )
    mirror = {n: np.conj(a) for n, a in abc.items()}  # u(-x) has coefficients conj(a_n) at n
    cm = S16.coef_from_modes(mirror)
    Qm, Gm, Fm = fluxes(S16, cm, "helicity")
    lines.append(
        f"mirrored ABC: helicity = {Qm:.6f}, heat flux = {Gm:.6f} (= +48 pi^3), Euler flux = {Fm:.2e}  -> helicity increases"
    )
    _, Ge, Fe = fluxes(S16, c, "enstrophy")
    lines.append(f"ABC: enstrophy Euler flux = {Fe:.2e} (Beltrami: steady Euler solution), heat flux = {Ge:.4f}")

    # --- 2D positive control: planar fields have zero enstrophy Euler flux ---
    planar_modes = [(1, 0, 0), (0, 1, 0), (1, 1, 0), (1, -1, 0), (2, 1, 0), (1, 2, 0)]
    worst = 0.0
    for _ in range(5):
        d = {}
        for n in planar_modes:
            amp = complex(rng.normal(), rng.normal())
            d[n] = amp * np.array([-n[1], n[0], 0.0])
        c = Spectral(32).coef_from_modes(d)
        _, _, F2 = fluxes(Spectral(32), c, "enstrophy")
        _, _, F2b = fluxes(Spectral(32), c, "L3", N_quad=32)
        worst = max(worst, abs(F2))
        _, _, Fh2 = fluxes(Spectral(32), c, "H12")
        worst = max(worst, abs(F2))
    lines.append(
        f"2D control: max |enstrophy Euler flux| over 5 random planar fields = {worst:.2e} (zero: no stretching in 2D); L^3 and H^(1/2) fluxes are NOT zero in 2D (e.g. {F2b:.4f}, {Fh2:.4f})"
    )

    # --- candidate table on Gaussian-integer fields ---
    cands = ["enstrophy", "L3", "H12", "omega32", "L4", "local_energy", "Flog"]
    n_fields = 60
    fields = []
    for _ in range(n_fields):
        coeffs = [
            ((int(rng.integers(-1, 2)), int(rng.integers(-1, 2))), (int(rng.integers(-1, 2)), int(rng.integers(-1, 2))))
            for _ in modes
        ]
        fields.append(coeffs)
    results = {w: [] for w in cands}
    for fi, coeffs in enumerate(fields):
        d = field_modes_from_gaussian(modes, coeffs)
        c = S16.coef_from_modes(d)
        en = S16.inner(c, c)
        if en == 0:
            continue
        for w in cands:
            Q, G, F = fluxes(S16, c, w, N_quad=32 if w in ("L3", "omega32", "Flog") else None)
            results[w].append((fi, Q, G, F, en))
    lines.append("")
    lines.append(
        "| candidate | fields with nonzero Euler flux | best field | Q(u) | heat flux G | Euler flux F | ‖u‖_2 | A* = -G/F | A* ‖u‖_2 |"
    )
    lines.append("|---|---|---|---|---|---|---|---|---|")
    best_fields = {}
    for w in cands:
        rows = results[w]
        nz = sum(1 for r in rows if abs(r[3]) > 1e-9 * max(1.0, abs(r[2])))
        # choose sign so that F > 0 (Euler flux is odd under u -> -u), then minimise A*||u||_2
        scored = []
        for fi, Q, G, F, en in rows:
            if abs(F) < 1e-9 * max(1.0, abs(G)):
                continue
            sign = 1 if F > 0 else -1
            Astar = -G / (sign * F)
            scored.append((Astar * np.sqrt(en), fi, sign, Q, G, sign * F, np.sqrt(en), Astar))
        scored.sort()
        s = scored[0]
        best_fields[w] = (s[1], s[2])
        lines.append(
            f"| {w} | {nz}/{len(rows)} | #{s[1]} (sign {s[2]:+d}) | {s[3]:.6g} | {s[4]:.6g} | {s[5]:.6g} | {s[6]:.6g} | {s[7]:.6g} | {s[0]:.6g} |"
        )

    lines.append("")
    lines.append(
        "Coefficients (alpha, beta) per mode (Gaussian integers as (re,im)) of the best fields; u = sum_n 2 Re[(alpha_n b_n^1 + beta_n b_n^2) e^{2 pi i n.x}]:"
    )
    lines.append("modes: " + str(modes))
    lines.append("perp bases b^1, b^2: " + str([perp_basis(n) for n in modes]))
    for w in cands:
        fi, sign = best_fields[w]
        lines.append(f"  {w}: field #{fi}, sign {sign:+d}: " + " ".join(str(cc) for cc in fields[fi]))

    # --- verification at A = 2 A* for the polynomial candidates, and convergence in N for the others ---
    lines.append("")
    lines.append("| candidate | field | N | G | F | flux at A = 2A* (should be > 0) |")
    lines.append("|---|---|---|---|---|---|")
    for w in cands:
        if w == "Flog":
            continue  # not homogeneous: the polynomial-in-A structure does not apply (see A5)
        fi, sign = best_fields[w]
        d = field_modes_from_gaussian(modes, fields[fi])
        d = {n: sign * a for n, a in d.items()}
        for Nq in [16] if w in ("enstrophy", "H12", "L4", "local_energy") else [32, 64, 96]:
            c = S16.coef_from_modes(d)
            Q, G, F = fluxes(S16, c, w, N_quad=Nq)
            Astar = -G / F
            A = 2 * Astar
            cA = S16.coef_from_modes({n: A * a for n, a in d.items()})
            QA, GA, FA = fluxes(S16, cA, w, N_quad=Nq)
            lines.append(f"| {w} | #{fi} | {Nq} | {G:.8g} | {F:.8g} | {GA + FA:.8g} |")

    # --- sup norms via Danskin ---
    lines.append("")
    lines.append(
        "| candidate | field | x* (maximiser) | max | second grid local max | ‖grad of the squared magnitude at x*‖ | G = u.Lap u/‖u‖ at x* | F = -u.grad p/‖u‖ (or w.Sw/‖w‖) at x* | A* |"
    )
    lines.append("|---|---|---|---|---|---|---|---|---|")
    sup_fields = []
    for kind, name in (("u", "Linf_u"), ("w", "Linf_omega")):
        bestk = None
        for fi in range(min(30, len(fields))):
            d = field_modes_from_gaussian(modes, fields[fi])
            if S16.inner(S16.coef_from_modes(d), S16.coef_from_modes(d)) == 0:
                continue
            r = sup_norm_fluxes(S16, d, kind)
            if r["gradnorm"] > 1e-8 * r["max"] ** 2:
                continue
            sign = 1 if r["F"] > 0 else -1
            gap = r["max"] - r["second_local_max"]
            if abs(r["F"]) < 1e-8 or gap < 0.05 * r["max"]:
                continue
            Astar = -r["G"] / (sign * r["F"])
            if bestk is None or Astar < bestk[0]:
                bestk = (Astar, fi, sign, r)
        Astar, fi, sign, r = bestk
        xs = ", ".join(f"{t:.6f}" for t in r["x*"])
        lines.append(
            f"| {name} | #{fi} (sign {sign:+d}) | ({xs}) | {r['max']:.6g} | {r['second_local_max']:.6g} | {r['gradnorm']:.1e} | {r['G']:.6g} | {sign * r['F']:.6g} | {Astar:.6g} |"
        )
        sup_fields.append(f"  {name}: field #{fi}, sign {sign:+d}: " + " ".join(str(cc) for cc in fields[fi]))

    lines.append("Coefficients of the sup-norm fields (same modes and bases as above):")
    lines.extend(sup_fields)

    # --- A5: F(|u|) = |u|^2 log(1+|u|): Pi(w), the amplitude scan on the torus, and the asymptotics ---
    fi, sign = best_fields["Flog"]
    d = field_modes_from_gaussian(modes, fields[fi])
    d = {n: sign * a for n, a in d.items()}
    lines.append("")
    lines.append("A5: Q = int |u|^2 log(1+|u|) on the best Flog field w above (field #%d, sign %+d)." % (fi, sign))
    lines.append(
        "Pi(w) = int p_w w.grad log|w| dx = int p_w w_j w_i d_j w_i / |w|^2 dx; large-A asymptotics: Euler flux = 2 Pi A^3 + o(A^3), heat flux = -(2 log A + 1) A^2 ‖grad w‖_2^2 + A^2 L1(w)[Lap w] + o(A^2)."
    )
    lines.append(
        "| N | ‖grad w‖_2^2 | Pi(w) | Euler flux / A^3 at A = 1e3 | at A = 1e5 | heat flux / (A^2 log A) at A = 1e5 |"
    )
    lines.append("|---|---|---|---|---|---|")
    for Nq in (32, 64, 96):
        Sq2 = Spectral(Nq)
        wc2 = Sq2.coef_from_modes(d)
        w2 = Sq2.phys(wc2)
        pw = Sq2.phys(Sq2.pressure(wc2))
        gw = np.stack([[Sq2.phys(Sq2.ik[j] * wc2[i]) for i in range(3)] for j in range(3)])  # gw[j, i] = d_j w_i
        aw2 = np.sum(w2**2, 0)
        num = np.einsum("jxyz,ixyz,jixyz->xyz", w2, w2, gw)
        Pi = Sq2.integral(pw * np.where(aw2 > 0, num / np.where(aw2 > 0, aw2, 1.0), 0.0))
        enst = -Sq2.inner(wc2, Sq2.laplacian(wc2))
        Bw2 = Sq2.phys(Sq2.B(wc2))
        lw2 = Sq2.phys(Sq2.laplacian(wc2))
        vals = []
        for A in (1e3, 1e5):
            au = A * np.sqrt(aw2)
            g = 2 * np.log1p(au) + au / (1 + au)
            vals.append(Sq2.integral(g * np.sum(w2 * Bw2, 0)) * A**3 / A**3)  # = int g(A|w|) w.B(w,w)
        A = 1e5
        au = A * np.sqrt(aw2)
        g = 2 * np.log1p(au) + au / (1 + au)
        heat_norm = Sq2.integral(g * np.sum(w2 * lw2, 0)) / np.log(A)
        lines.append(f"| {Nq} | {enst:.6g} | {Pi:.6g} | {vals[0]:.6g} | {vals[1]:.6g} | {heat_norm:.6g} |")
    Sq2 = Spectral(96)
    wc2 = Sq2.coef_from_modes(d)
    w2 = Sq2.phys(wc2)
    aw = np.sqrt(np.sum(w2**2, 0))
    Bw2 = Sq2.phys(Sq2.B(wc2))
    lw2 = Sq2.phys(Sq2.laplacian(wc2))
    enst = -Sq2.inner(wc2, Sq2.laplacian(wc2))
    gw = np.stack([[Sq2.phys(Sq2.ik[j] * wc2[i]) for i in range(3)] for j in range(3)])
    gradw2 = np.sum(gw**2, axis=(0, 1))
    aw2 = aw**2
    C1 = Sq2.integral(np.where(aw2 > 0, np.log(np.where(aw2 > 0, aw, 1.0)), 0.0) * gradw2)  # int log|w| |grad w|^2
    grad_abs2 = np.sum(np.einsum("ixyz,jixyz->jxyz", w2, gw) ** 2, 0)  # |w_i d_j w_i|^2 summed over j
    C2 = Sq2.integral(np.where(aw2 > 0, grad_abs2 / np.where(aw2 > 0, aw2, 1.0), 0.0))  # int |grad |w||^2
    lines.append(
        f"Asymptotic constants at N = 96: ‖grad w‖_2^2 = {enst:.6g}, C1 = int log|w| |grad w|^2 = {C1:.6g}, C2 = int |grad|w||^2 = {C2:.6g}, Pi = {Pi:.6g}"
    )
    lines.append(
        "Amplitude scan (N_quad = 96), flux = DQ(Aw)[Lap(Aw) + B(Aw,Aw)]; prediction = 2 Pi A^3 - A^2 [(2 log A + 1) ‖grad w‖_2^2 + 2 C1 + 2 C2] (error o(A^2)):"
    )
    lines.append("| A | heat flux | Euler flux | total flux | prediction |")
    lines.append("|---|---|---|---|---|")
    first_pos = None
    for A in [1.0, 10.0, 100.0, 200.0, 300.0, 400.0, 500.0, 700.0, 1000.0, 3000.0]:
        au = A * aw
        g = 2 * np.log1p(au) + au / (1 + au)
        heat = Sq2.integral(g * np.sum(w2 * lw2, 0)) * A**2
        eul = Sq2.integral(g * np.sum(w2 * Bw2, 0)) * A**3
        pred = 2 * Pi * A**3 - A**2 * ((2 * np.log(A) + 1) * enst + 2 * C1 + 2 * C2)
        if eul + heat > 0 and first_pos is None:
            first_pos = A
        lines.append(f"| {A:g} | {heat:.6g} | {eul:.6g} | {heat + eul:.6g} | {pred:.6g} |")
    lines.append(f"A5: first scanned amplitude with positive Navier-Stokes derivative of Q: A = {first_pos:g}")
    # convergence in N of the total flux at A = 1000
    conv = []
    for Nq in (48, 64, 96, 128):
        Sq3 = Spectral(Nq)
        wc3 = Sq3.coef_from_modes(d)
        w3 = Sq3.phys(wc3)
        A = 1000.0
        au = A * np.sqrt(np.sum(w3**2, 0))
        g = 2 * np.log1p(au) + au / (1 + au)
        tot = (
            Sq3.integral(g * np.sum(w3 * Sq3.phys(Sq3.laplacian(wc3)), 0)) * A**2
            + Sq3.integral(g * np.sum(w3 * Sq3.phys(Sq3.B(wc3)), 0)) * A**3
        )
        conv.append(f"N = {Nq}: {tot:.8g}")
    lines.append("A5: total flux at A = 1000 versus quadrature resolution: " + "; ".join(conv))

    # helicity: Euler flux always zero, heat flux of either sign
    hel_pos, hel_tot = 0, 0
    hel_maxF = 0.0
    for coeffs in fields[:30]:
        c = S16.coef_from_modes(field_modes_from_gaussian(modes, coeffs))
        if S16.inner(c, c) == 0:
            continue
        _, G, F = fluxes(S16, c, "helicity")
        hel_tot += 1
        hel_pos += G > 0
        hel_maxF = max(hel_maxF, abs(F))
    lines.append(
        f"Helicity on {hel_tot} Gaussian-integer fields: max |Euler flux| = {hel_maxF:.1e} (Euler invariant), heat flux positive for {hel_pos}/{hel_tot} fields (no sign)"
    )

    # --- small optimisation of the enstrophy threshold A*||u||_2 over real coefficients (|n|_inf <= 1 modes) ---
    def unpack(x):
        d = {}
        for i, n in enumerate(modes):
            b1, b2 = perp_basis(n)
            a = complex(x[4 * i], x[4 * i + 1]) * np.array(b1) + complex(x[4 * i + 2], x[4 * i + 3]) * np.array(b2)
            d[n] = a
        return d

    def objective(x):
        c = S16.coef_from_modes(unpack(x))
        en = S16.inner(c, c)
        _, G, F = fluxes(S16, c, "enstrophy")
        if F <= 1e-12 or en < 1e-12:
            return 1e6
        return -G / F * np.sqrt(en)

    fi, sign = best_fields["enstrophy"]
    x0 = []
    for alpha, beta in fields[fi]:
        x0 += [sign * alpha[0], sign * alpha[1], sign * beta[0], sign * beta[1]]
    x0 = np.array(x0, dtype=float)
    res = minimize(
        objective, x0, method="Nelder-Mead", options={"maxiter": 4000, "maxfev": 4000, "xatol": 1e-6, "fatol": 1e-9}
    )
    lines.append("")
    c_opt = S16.coef_from_modes(unpack(res.x))
    _, G_opt, F_opt = fluxes(S16, c_opt, "enstrophy")
    en_opt = np.sqrt(S16.inner(c_opt, c_opt))
    lines.append(
        f"Enstrophy threshold optimisation (Nelder-Mead over 52 real coefficients, modes |n|_inf <= 1): start A*‖u‖_2 = {objective(x0):.6g}, "
        f"end {res.fun:.6g} ({res.nfev} evaluations); optimised field: ‖u‖_2 = {en_opt:.6g}, G = {G_opt:.6g}, F = {F_opt:.6g}, A* = {-G_opt / F_opt:.6g}"
    )
    lines.append(f"FFT part wall-clock: {time.time() - t0:.1f} s")
    return lines


def main():
    t0 = time.time()
    print("# line_a_lyapunov.py output\n")
    print("## Part 1: sympy pointwise identities\n")
    for name, ok in sympy_identities():
        print(f"- {name}: {'OK' if ok else 'FAILED'}")
    print(f"\n(sympy part: {time.time() - t0:.1f} s)\n")
    t1 = time.time()
    print("## Part 2: exact rational certificates\n")
    for line in exact_part():
        print(line)
    print(f"\n(exact part: {time.time() - t1:.1f} s)\n")
    print("## Part 3: FFT computations\n")
    for line in fft_part():
        print(line)
    print(f"\nTotal wall-clock: {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
