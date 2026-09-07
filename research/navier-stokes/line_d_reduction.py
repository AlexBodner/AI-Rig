"""Line D: the L^3 wall for 3D Navier-Stokes, and numerics bracketing its threshold.

Part 1: explicit constants (Talenti's sharp Sobolev constant; the sharp Gagliardo-Nirenberg
constant for ||w||_4 <= C ||grad w||_2^{3/4} ||w||_2^{1/4}, via the ground state of
-Delta Q + Q = Q^3 on R^3) and the resulting threshold kappa of D2.
Part 2: verification on T^3 of the exact identity D1
    d/dt int|u|^3 = -3 nu int|u||grad u|^2 - (4nu/3) int|grad(|u|^{3/2})|^2 + 3 int p (u.grad|u|).
Part 2b: the antiperiodic class of D3b, where the pressure work vanishes identically.
Part 3: partial descent on the T^3 functionals Theta and K, and the T^3 pressure ratio.
Part 4: an explicit T^3 field with d/dt ||u||_{L^3}^3 > 0 at nu = 1.
Part 5: the pressure work is not identically zero on T^3 (D3c).  Exact symbolic evaluation of
    J(v) = int q d_1 v_1 with q = R_i R_j (v_i v_j), for the three-mode divergence-free field
    v = A e_a cos x_1 + B e_b cos x_2 + C e_c sin(x_1+x_2), plus a grid check of the expansion
    P(c e_1 + eps v) = 3 c eps^3 J(v) + O(eps^4).

Run: python3 line_d_reduction.py   (about 55 s on 4 CPUs, no files written)
"""

import numpy as np
import sympy as sp
from mpmath import gamma as mp_gamma
from mpmath import mp, mpf, pi
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

SEED = 20260906
NU = 1.0
mp.dps = 30


# Part 1: constants
def talenti_sobolev_constant(n=3, p=2):
    """Sharp C in ||u||_{L^{p*}(R^n)} <= C ||grad u||_{L^p(R^n)}, Talenti 1976."""
    n_, p_ = mpf(n), mpf(p)
    pref = pi ** (-mpf(1) / 2) * n_ ** (-1 / p_) * ((p_ - 1) / (n_ - p_)) ** (1 - 1 / p_)
    br = mp_gamma(1 + n_ / 2) * mp_gamma(n_) / (mp_gamma(n_ / p_) * mp_gamma(1 + n_ - n_ / p_))
    return pref * br ** (1 / n_)


def _shoot(alpha, rmax=30.0):
    """Radial shooting for -Q'' - (2/r)Q' + Q = Q^3 from r = eps with Q(0) = alpha."""
    eps = 1e-8
    hit_zero, blow_up = (lambda r, y: y[0]), (lambda r, y: y[0] - 10.0 * alpha)
    hit_zero.terminal = blow_up.terminal = True
    sol = solve_ivp(
        lambda r, y: [y[1], -2.0 * y[1] / r + y[0] - y[0] ** 3],
        (eps, rmax),
        [alpha + (alpha - alpha**3) * eps**2 / 6.0, (alpha - alpha**3) * eps / 3.0],
        events=(hit_zero, blow_up),
        rtol=1e-12,
        atol=1e-14,
        dense_output=True,
        max_step=0.05,
    )
    return len(sol.t_events[0]) > 0, sol


def ground_state(n_bisect=200):
    """Bisect on Q(0) for the positive decaying ground state of -Delta Q + Q = Q^3 on R^3."""
    lo, hi = 1.0, 20.0  # lo: blows up to +infinity; hi: crosses zero
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        if _shoot(mid)[0]:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-15:
            break
    alpha = 0.5 * (lo + hi)
    return alpha, _shoot(alpha)[1]


def ground_state_norms(sol, rmax):
    """Return (M, G, F4) = (||Q||_2^2, ||grad Q||_2^2, ||Q||_4^4) with the 4 pi r^2 measure."""
    r = np.linspace(1e-8, rmax, 400001)
    y = sol.sol(r)
    q, dq = y[0], y[1]
    w = 4.0 * np.pi * r**2
    trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return trap(w * q**2, r), trap(w * dq**2, r), trap(w * q**4, r)


def part1():
    print("=" * 78)
    print("PART 1 - explicit constants and the threshold kappa")
    print("=" * 78)
    c_sob = float(talenti_sobolev_constant())
    print(f"Talenti sharp Sobolev C_S (||w||_6 <= C_S ||grad w||_2, R^3) = {c_sob:.10f}")
    print(f"  S = C_S^-2 = {c_sob**-2:.10f} (= 3(pi/2)^(4/3) = {3 * (np.pi / 2) ** (4 / 3):.10f})")
    c_gn_rig = c_sob**0.75
    print(f"rigorous GN constant via Hoelder+Sobolev: C_GN <= C_S^(3/4) = {c_gn_rig:.10f}")

    alpha, sol = ground_state()
    rend = sol.t[-1]
    m2, g2, f4 = ground_state_norms(sol, min(rend, 30.0))
    print(f"\nground state of -Delta Q + Q = Q^3 on R^3: Q(0) = {alpha:.10f}, integrated to r = {rend:.3f}")
    print(f"  M = ||Q||_2^2 = {m2:.8f};  G = ||grad Q||_2^2 = {g2:.8f}  (Pohozaev G = 3M: {abs(g2 - 3 * m2) / m2:.1e})")
    print(f"  F = ||Q||_4^4 = {f4:.8f}  (Pohozaev F = 4M: {abs(f4 - 4 * m2) / m2:.1e})")
    c_gn_sharp = f4**0.25 / (g2**0.375 * m2**0.125)
    print(f"  sharp C_GN = ||Q||_4/(||grad Q||_2^(3/4) ||Q||_2^(1/4)) = {c_gn_sharp:.10f}")
    print(f"  closed form 4^(1/4) 3^(-3/8) M^(-1/4)                  = {4**0.25 * 3**-0.375 * m2**-0.25:.10f}")

    c_p_rig = 3.0 ** (8.0 / 3.0)
    print(f"\nrigorous pressure constant (Iwaniec-Martin ||R_j||_{{L^3}} = sqrt3): C_p <= 3^(8/3) = {c_p_rig:.6f}")
    print("\nkappa / nu  =  4 / (3 C_p C_GN^{4/3})   [D2]")
    print(f"| {'C_GN':<22} | {'C_p':<19} | {'C_GN^{4/3}':>10} | {'kappa/nu':>9} |")
    print(f"| {'-' * 22} | {'-' * 19} | {'-' * 10} | {'-' * 9} |")
    for name, gn in (("C_S^(3/4) (rigorous)", c_gn_rig), ("sharp (ground state)", c_gn_sharp)):
        g43 = gn ** (4.0 / 3.0)
        print(f"| {name:<22} | {'3^(8/3) (rigorous)':<19} | {g43:>10.6f} | {4.0 / (3.0 * c_p_rig * g43):>9.6f} |")
    return {"C_S": c_sob, "C_GN_rig": c_gn_rig, "C_GN_sharp": c_gn_sharp, "C_p_rig": c_p_rig, "M": m2}


# Spectral machinery on T^3 = [0,2pi)^3
class Grid:
    def __init__(self, n=32):
        self.n = n
        k = np.fft.fftfreq(n, d=1.0 / n)
        self.kx, self.ky, self.kz = np.meshgrid(k, k, k, indexing="ij")
        self.k2 = self.kx**2 + self.ky**2 + self.kz**2
        self.k2inv = np.where(self.k2 == 0, 1.0, self.k2)
        self.kvec = (self.kx, self.ky, self.kz)
        self.dV = (2.0 * np.pi / n) ** 3

    def integral(self, f):
        return float(np.sum(f) * self.dV)

    def grad(self, f):
        fh = np.fft.fftn(f)
        return [np.real(np.fft.ifftn(1j * kj * fh)) for kj in self.kvec]

    def project(self, v):
        """Leray projection of a real vector field given in physical space."""
        vh = [np.fft.fftn(c) for c in v]
        div = sum(kj * ch for kj, ch in zip(self.kvec, vh))
        out = []
        for kj, ch in zip(self.kvec, vh):
            ph = ch - kj * div / self.k2inv
            ph[0, 0, 0] = ch[0, 0, 0]
            out.append(np.real(np.fft.ifftn(ph)))
        return out

    def pressure(self, u):
        """p with -Delta p = d_i d_j (u_i u_j), zero mean."""
        num = np.zeros_like(self.k2, dtype=complex)
        for i, ki in enumerate(self.kvec):
            for j, kj in enumerate(self.kvec):
                num += ki * kj * np.fft.fftn(u[i] * u[j])
        # -Delta p = d_i d_j (u_i u_j)  =>  p-hat = -k_i k_j (u_i u_j)-hat / |k|^2
        ph = -num / self.k2inv
        ph[0, 0, 0] = 0.0
        return np.real(np.fft.ifftn(ph))


def functionals(g, u, nu=NU, full=True):
    """Pieces of d/dt int|u|^3 for a solenoidal field u on T^3.  D = 3 int|u||grad u|^2 +
    3 int|u||grad|u||^2 = 3 int|u||grad u|^2 + (4/3) X^2 with X^2 = int|grad(|u|^{3/2})|^2;
    P = 3 int p (u.grad|u|).  With full=True also the direct d/dt and the transport term."""
    absu = np.sqrt(sum(c**2 for c in u))
    safe = np.maximum(absu, 1e-300)
    du = [g.grad(c) for c in u]  # du[i][j] = d_j u_i
    gradsq = sum(du[i][j] ** 2 for i in range(3) for j in range(3))
    grad_absu = [sum(u[i] * du[i][j] for i in range(3)) / safe for j in range(3)]
    p = g.pressure(u)
    xsq = (9.0 / 4.0) * g.integral(absu * sum(c**2 for c in grad_absu))
    diss = 3.0 * g.integral(absu * gradsq) + (4.0 / 3.0) * xsq
    press = 3.0 * g.integral(p * sum(u[j] * grad_absu[j] for j in range(3)))
    out = {"diss": diss, "press": press, "X2": xsq, "L3": g.integral(absu**3) ** (1.0 / 3.0)}
    out["rhs"] = -nu * diss + press
    if not full:
        return out
    conv = [sum(u[j] * du[i][j] for j in range(3)) for i in range(3)]
    lap = [np.real(np.fft.ifftn(-g.k2 * np.fft.fftn(c))) for c in u]
    ut = g.project([nu * lap[i] - conv[i] for i in range(3)])
    out["lhs"] = 3.0 * g.integral(absu * sum(u[i] * ut[i] for i in range(3)))
    out["transport"] = 3.0 * g.integral(absu * sum(u[i] * conv[i] for i in range(3)))
    out["L6sq"] = g.integral(absu**6) ** (1.0 / 3.0)
    out["p3"] = g.integral(np.abs(p) ** 3) ** (1.0 / 3.0)
    return out


def part2(rng):
    """Verify D1 by refining the grid: the exact identity is independent of N, the residual
    is the quadrature error of the non-band-limited |u| and it goes to zero with N."""
    print()
    print("=" * 78)
    print("PART 2 - spectral verification of the exact identity D1 on T^3, refining N")
    print("=" * 78)
    g0 = Grid(16)
    c = rng.standard_normal(len(basis_fields(g0, kmax2=2)))
    print("| field                  |  N |  d/dt int|u|^3 (direct) |  -nu D + P (identity D1) | rel. err | transport |")
    print("| ---------------------- | -- | ----------------------- | ------------------------ | -------- | --------- |")
    out = []
    for mean, label in ((0.0, "zero mean, |k|^2 <= 2"), (3.0, "mean (3,0,0) added   ")):
        for n in (24, 32, 48, 64):
            g = Grid(n)
            u = assemble(basis_fields(g, kmax2=2), c)
            u[0] = u[0] + mean
            f = functionals(g, u)
            rel = abs(f["lhs"] - f["rhs"]) / abs(f["lhs"])
            out.append(rel)
            print(
                f"| {label} | {n:>2d} | {f['lhs']:>23.12e} | {f['rhs']:>24.12e} | {rel:>8.1e} | {f['transport']:>9.1e} |"
            )
    print("\nBoth the residual and the transport term (exactly zero in the continuum) decrease with N:")
    print("the identity D1 is confirmed; the residual is grid quadrature of |u|, not a defect of D1.")
    return out


def part2b(rng):
    """D3b: if u(x + a) = -u(x) for a half-period a, the pressure work P(u) vanishes."""
    print()
    print("=" * 78)
    print("PART 2b - the antiperiodic class: P(u) = 0 whenever u(x + pi*(1,1,1)) = -u(x)")
    print("=" * 78)
    g = Grid(32)
    sh = g.n // 2
    print("| field                                     | antiperiodicity defect |     D(u) |        P(u) |")
    print("| ----------------------------------------- | ---------------------- | -------- | ----------- |")
    cases = [("shell |k| = 1 (12 modes), random", 1), ("shells |k|^2 <= 2 (36 modes), random", 2)]
    vals = []
    for label, kmax2 in cases:
        fields = basis_fields(g, kmax2=kmax2)
        u = assemble(fields, rng.standard_normal(len(fields)))
        nrm = np.sqrt(sum(g.integral(ci**2) for ci in u))
        u = [ci / nrm for ci in u]
        us = [np.roll(np.roll(np.roll(ci, sh, 0), sh, 1), sh, 2) for ci in u]
        defect = np.sqrt(sum(float(((us[i] + u[i]) ** 2).sum()) for i in range(3)) / g.n**3)
        f = functionals(g, u, full=False)
        vals.append((defect, f["press"]))
        print(f"| {label:<41} | {defect:>22.3e} | {f['diss']:>8.5f} | {f['press']:>11.3e} |")
    print("\nThe |k| = 1 shell is antiperiodic (defect at round-off) and has P = 0 at round-off;")
    print("mixing |k|^2 = 1 with |k|^2 = 2 breaks the parity and P becomes nonzero. See D3b.")
    return vals


# Parts 3-4: optimizing the scale-invariant ratios over a divergence-free basis
def basis_fields(g, kmax2=2):
    """Real divergence-free cos/sin building blocks with |k|^2 <= kmax2, on T^3."""
    x = np.arange(g.n) * (2.0 * np.pi / g.n)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    seen, fields = set(), []
    rng_k = range(-2, 3)
    for k1 in rng_k:
        for k2 in rng_k:
            for k3 in rng_k:
                k = np.array([k1, k2, k3], dtype=float)
                n2 = int(k1 * k1 + k2 * k2 + k3 * k3)
                if n2 == 0 or n2 > kmax2:
                    continue
                key = tuple(int(v) for v in k)
                if tuple(-v for v in key) in seen:
                    continue
                seen.add(key)
                a = np.array([1.0, 0.0, 0.0])
                if abs(np.dot(a, k)) > 0.9 * np.linalg.norm(k):
                    a = np.array([0.0, 1.0, 0.0])
                e1 = np.cross(k, a)
                e1 /= np.linalg.norm(e1)
                e2 = np.cross(k, e1)
                e2 /= np.linalg.norm(e2)
                phase = k1 * X + k2 * Y + k3 * Z
                for e in (e1, e2):
                    for trig in (np.cos(phase), np.sin(phase)):
                        fields.append([e[i] * trig for i in range(3)])
    return fields


def assemble(fields, c):
    return [sum(ci * f[i] for ci, f in zip(c, fields)) for i in range(3)]


def theta(g, fields, c):
    """Theta(u) = nu D(u) ||u||_3 / |P(u)|: the L^3 norm at which the amplitude family a*u
    first has d/dt ||a u||_3^3 > 0. Scale-invariant in a (D is 3- and P is 4-homogeneous);
    P is odd in u and D is even, so the sign of P costs nothing."""
    u = assemble(fields, c)
    nrm = np.sqrt(sum(g.integral(ci**2) for ci in u))
    if nrm < 1e-8:
        return 1e12
    u = [ci / nrm for ci in u]
    f = functionals(g, u, full=False)
    if abs(f["press"]) < 1e-14:
        return 1e12
    return NU * f["diss"] * f["L3"] / abs(f["press"])


def optimize_theta(gopt, kmax2, rng, nrestart, maxfev, method="Nelder-Mead"):
    fields = basis_fields(gopt, kmax2=kmax2)
    best, bestc = np.inf, None
    for _ in range(nrestart):
        c0 = rng.standard_normal(len(fields))
        res = minimize(
            lambda c: theta(gopt, fields, c),
            c0,
            method=method,
            options={"maxiter": maxfev, "maxfev": maxfev},
        )
        if res.fun < best:
            best, bestc = res.fun, res.x
    return len(fields), best, bestc


def evaluate(gfine, kmax2, c):
    fields = basis_fields(gfine, kmax2=kmax2)
    u = assemble(fields, c)
    nrm = np.sqrt(sum(gfine.integral(ci**2) for ci in u))
    u = [ci / nrm for ci in u]
    f = functionals(gfine, u)
    if f["press"] < 0:
        u = [-ci for ci in u]
        f = functionals(gfine, u)
    return u, f


def part3(gopt, gfine, rng, consts):
    print()
    print("=" * 78)
    print(f"PART 3 - partial descent on Theta over T^3 fields (search N={gopt.n}, evaluation N={gfine.n})")
    print("=" * 78)
    rows, bestc = [], {}
    for kmax2, method, nrestart, maxfev in ((2, "Nelder-Mead", 2, 1200), (3, "Powell", 1, 2200)):
        m, best_coarse, c = optimize_theta(gopt, kmax2, rng, nrestart, maxfev, method)
        _, f = evaluate(gfine, kmax2, c)
        thr = NU * f["diss"] * f["L3"] / f["press"]
        rows.append((kmax2, m, method, best_coarse, thr, f["press"] / (f["X2"] * f["L3"]), f["p3"] / f["L6sq"]))
        bestc[kmax2] = c
    print(
        "| basis    | modes | method      | Theta at N=%2d | Theta at N=%2d | K = P/(X^2||u||_3) | ||p||_3/|| |u|^2 ||_3 |"
        % (gopt.n, gfine.n)
    )
    print(
        "| -------- | ----- | ----------- | ------------- | ------------- | ------------------ | --------------------- |"
    )
    for kmax2, m, meth, bc, thr, kval, cp in rows:
        print(f"| |k|^2<={kmax2} | {m:>5d} | {meth:<11} | {bc:>13.5f} | {thr:>13.5f} | {kval:>18.6f} | {cp:>21.6f} |")
    thr_ub = min(r[4] for r in rows)
    k_lb = max(r[5] for r in rows)
    cp_lb = max(r[6] for r in rows)
    kub = 2.0 * consts["C_p_rig"] * consts["C_GN_rig"] ** (4.0 / 3.0)
    kub_s = 2.0 * consts["C_p_rig"] * consts["C_GN_sharp"] ** (4.0 / 3.0)
    kap, kap_s = (8.0 / 3.0) / kub, (8.0 / 3.0) / kub_s
    print()
    print("EVERY NUMBER BELOW THAT COMES FROM THIS SEARCH IS A T^3 NUMBER.  The proved bounds are")
    print("R^3 statements (D2).  No transference between the two is proved anywhere in this line,")
    print("so the two columns are NOT two ends of one bracket; they are bounds on two functionals.")
    print()
    print(f"T^3:  K*_T3 >= {k_lb:.6f} (partial descent)      kappa_sharp_T3 <= {thr_ub:.5f} nu (partial descent)")
    print(f"      ||p||_3/|| |u|^2 ||_3 >= {cp_lb:.6f}  (a lower bound for the T^3 multiplier norm)")
    print()
    print(f"R^3:  K*_R3 <= {kub:.5f} (proved, C_GN = C_S^(3/4)) / {kub_s:.5f} (proved, sharp C_GN)")
    print(f"      kappa_sharp_R3 >= {kap:.6f} nu (rigorous C_GN) / {kap_s:.6f} nu (sharp C_GN)")
    print(f"      C_p <= {consts['C_p_rig']:.6f} (proved, given ||R_j||_{{L^3}} = sqrt3, an unread constant)")
    print("      C_p >= the T^3 ratio above ONLY via de Leeuw transference, whose continuity")
    print("      hypothesis the symbol -k_i k_j/|k|^2 fails at the lattice point k = 0 (see D4d).")
    print()
    print(f"If K* were transferable the constant chain would be loose by {kub / k_lb:.0f}x (rigorous C_GN),")
    print(f"{kub_s / k_lb:.0f}x (sharp C_GN); it is not known to be, so this is a comparison, not a bound.")
    return {"thr_ub": thr_ub, "K_lb": k_lb, "C_p_lb": cp_lb, "kap": kap, "kap_s": kap_s, "c": bestc}


def part4(gfine, c, kmax2):
    print()
    print("=" * 78)
    print("PART 4 - an explicit T^3 field with d/dt ||u||_{L^3}^3 > 0 at nu = 1")
    print("=" * 78)
    u, f = evaluate(gfine, kmax2, c)
    a_crit = NU * f["diss"] / f["press"]
    cc = c / np.linalg.norm(c)
    print(f"v: the unit-L^2 field on the |k|^2 <= {kmax2} divergence-free modes found in Part 3,")
    print("signed so P(v) > 0. Coefficients in the basis order of basis_fields():")
    print("  " + ", ".join(f"{v:+.4f}" for v in cc))
    print(f"  D(v) = {f['diss']:.6f}   P(v) = {f['press']:.6f}   ||v||_3 = {f['L3']:.6f}")
    print(f"  d/dt ||a v||_3^3 = -a^3 nu D(v) + a^4 P(v) > 0  iff  a > nu D(v)/P(v) = {a_crit:.6f}")
    print()
    print("| a / a_crit |  ||a v||_3 |     -nu D + P (identity) |  direct d/dt int|u|^3 |  err / nu a^3 D |")
    print("| ---------- | ---------- | ------------------------ | --------------------- | --------------- |")
    for r in (0.5, 1.0, 1.5, 3.0):
        a = r * a_crit
        ff = functionals(gfine, [a * ci for ci in u])
        rel = abs(ff["lhs"] - ff["rhs"]) / (NU * a**3 * f["diss"])
        print(f"| {r:>10.2f} | {ff['L3']:>10.5f} | {ff['rhs']:>24.10e} | {ff['lhs']:>21.10e} | {rel:>15.1e} |")
    print()
    print("=> ||u||_{L^3} is NOT a Lyapunov functional for 3D Navier-Stokes on T^3 at nu = 1;")
    print(f"   for this family the L^3 norm starts to increase at ||u||_3 = {a_crit * f['L3']:.5f}.")
    return a_crit * f["L3"]


# Part 5: the pressure work does not vanish identically on T^3 (D3c)
# v = A e_a cos(x_1) + B e_b cos(x_2) + C e_c sin(x_1 + x_2) with e_a = (0,1,0), e_b = (1,0,0),
# e_c = (1,-1,0); each polarization is orthogonal to its wavevector, so div v = 0, and the three
# wavevectors form a triad (1,0,0) + (0,1,0) = (1,1,0), which is what makes the pressure work
# nonzero: on a single shell there is no triad and P vanishes (D3b).
E_A, E_B, E_C = (0, 1, 0), (1, 0, 0), (1, -1, 0)
K_A, K_B, K_C = (1, 0, 0), (0, 1, 0), (1, 1, 0)


def _vhat():
    """Fourier coefficients {k: [vhat_1, vhat_2, vhat_3]} of v, exact in sympy."""
    a, b, c = sp.symbols("A B C", real=True)
    half = sp.Rational(1, 2)
    modes = {}

    def add(k, vec):
        cur = modes.setdefault(tuple(k), [sp.Integer(0)] * 3)
        for i in range(3):
            cur[i] += vec[i]

    for k, e, amp in ((K_A, E_A, a), (K_B, E_B, b)):
        add(k, [amp * half * e[i] for i in range(3)])  # cos(k.x) -> e/2 at +-k
        add([-ki for ki in k], [amp * half * e[i] for i in range(3)])
    add(K_C, [-sp.I * c * half * E_C[i] for i in range(3)])  # sin(k.x) -> -i e/2 at +k
    add([-ki for ki in K_C], [sp.I * c * half * E_C[i] for i in range(3)])
    return (a, b, c), modes


def exact_J():
    """J(v) = int_{T^3} q d_1 v_1 dx with -Delta q = d_i d_j (v_i v_j), zero mean, computed
    exactly from the Fourier coefficients:  J = (2pi)^3 sum_k qhat(k) conj(i k_1 vhat_1(k))."""
    syms, modes = _vhat()
    ks = list(modes)
    for k in ks:
        assert sp.simplify(sum(k[i] * modes[k][i] for i in range(3))) == 0, "v is not solenoidal"
    total = sp.Integer(0)
    for k in ks:
        k2 = sum(ki * ki for ki in k)
        f = [[sp.Integer(0)] * 3 for _ in range(3)]
        for kp in ks:
            kpp = tuple(k[m] - kp[m] for m in range(3))
            if kpp not in modes:
                continue
            for i in range(3):
                for j in range(3):
                    f[i][j] += modes[kp][i] * modes[kpp][j]
        qh = -sum(k[i] * k[j] * f[i][j] for i in range(3) for j in range(3)) / sp.Integer(k2)
        total += qh * sp.conjugate(sp.I * k[0] * modes[k][0])
    return syms, sp.simplify(sp.expand(total)) * (2 * sp.pi) ** 3


def _v_on_grid(g, amp=(1.0, 1.0, 1.0)):
    a, b, c = amp
    x = np.arange(g.n) * (2.0 * np.pi / g.n)
    xx, yy, _ = np.meshgrid(x, x, x, indexing="ij")
    s = np.sin(xx + yy)
    return [b * np.cos(yy) + c * s, a * np.cos(xx) - c * s, np.zeros_like(xx)]


def part5(gfine):
    print()
    print("=" * 78)
    print("PART 5 - the pressure work is not identically zero on T^3 (D3c)")
    print("=" * 78)
    syms, jsym = exact_J()
    a, b, c = syms
    jnum = float(jsym.subs({a: 1, b: 1, c: 1}))
    assert sp.simplify(jsym + 2 * sp.pi**3 * a * b * c) == 0, jsym
    print(f"exact (sympy, closed form):  J(v) = int q d_1 v_1 dx = {jsym}")
    print(f"                             at A = B = C = 1:  J = -2 pi^3 = {jnum:.10f}")
    print("expansion (D3c):  P(c e_1 + eps v) = 3 c eps^3 J + O(eps^4), so P < 0 for small eps > 0")
    print("                  and P > 0 for small eps < 0; hence P is not identically zero on T^3.")
    print()
    v = _v_on_grid(gfine)
    vmax = float(np.sqrt(sum(ci**2 for ci in v)).max())
    print(f"grid check at N = {gfine.n}, c = 1 (|u| >= 1 - eps*{vmax:.4f} > 0, so |u| is analytic):")
    print("|    eps |            P(c e_1 + eps v) |  P / (3 c eps^3) |   rel. dev. from J |")
    print("| ------ | --------------------------- | ---------------- | ------------------ |")
    for eps in (0.2, 0.1, 0.05, 0.025):
        u = [1.0 + eps * v[0], eps * v[1], eps * v[2]]
        press = functionals(gfine, u, full=False)["press"]
        ratio = press / (3.0 * eps**3)
        print(f"| {eps:>6.3f} | {press:>27.15e} | {ratio:>16.8f} | {abs(ratio - jnum) / abs(jnum):>18.2e} |")
    g2 = Grid(64)
    v2 = _v_on_grid(g2)
    worst = 0.0
    for eps in (0.2, 0.1, 0.05, 0.025):
        r1 = functionals(gfine, [1.0 + eps * v[0], eps * v[1], eps * v[2]], full=False)["press"]
        r2 = functionals(g2, [1.0 + eps * v2[0], eps * v2[1], eps * v2[2]], full=False)["press"]
        worst = max(worst, abs(r1 - r2) / abs(r1))
    print()
    print("The deviation falls like eps: it is the O(eps^4) term of the expansion, not a quadrature")
    print(f"error.  Largest relative change in P between N = {gfine.n} and N = 64: {worst:.1e} (|u| never")
    print("vanishes here, so |u| is analytic and the trapezoidal rule is spectrally accurate).")
    return jnum


def main():
    rng = np.random.default_rng(SEED)
    consts = part1()
    part2(rng)
    part2b(rng)
    gfine = Grid(48)
    res = part3(Grid(24), gfine, rng, consts)
    kmax2 = 3 if 3 in res["c"] else 2
    part4(gfine, res["c"][kmax2], kmax2)
    part5(gfine)
    print()
    print("done.")


if __name__ == "__main__":
    main()
