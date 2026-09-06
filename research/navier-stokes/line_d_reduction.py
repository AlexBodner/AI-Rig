"""Line D: the L^3 wall for 3D Navier-Stokes, and numerics bracketing its threshold.

Parts:
  1. Explicit constants: Talenti's sharp Sobolev constant on R^3, the sharp
     Gagliardo-Nirenberg constant for ||w||_4 <= C ||grad w||_2^{3/4} ||w||_2^{1/4}
     (via the ground state of -Q'' - (2/r)Q' + Q = Q^3), and the resulting
     small-data threshold kappa in the theorem of D2.
  2. Spectral verification on T^3 of the exact identity
        d/dt int|u|^3 = -3 nu int |u||grad u|^2 - (4 nu/3) int |grad(|u|^{3/2})|^2
                        + 3 int p (u . grad|u|)
     including the vanishing of the transport term.
  3. Sampling of random solenoidal fields on T^3 to bracket the sharp threshold
     from above and to lower-bound the pressure multiplier norm (de Leeuw).
  4. An explicit field at which d/dt ||u||_{L^3}^3 > 0 at nu = 1.

Run: python3 line_d_reduction.py   (about 40 s on 4 CPUs, no files written)
"""

import numpy as np
from mpmath import gamma as mp_gamma
from mpmath import mp, mpf, pi
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

SEED = 20260906
NGRID = 32
NU = 1.0
mp.dps = 30


# ----------------------------------------------------------------------------
# Part 1: constants
# ----------------------------------------------------------------------------
def talenti_sobolev_constant(n=3, p=2):
    """Sharp C in ||u||_{L^{p*}(R^n)} <= C ||grad u||_{L^p(R^n)}, Talenti 1976."""
    n_, p_ = mpf(n), mpf(p)
    pref = pi ** (-mpf(1) / 2) * n_ ** (-1 / p_) * ((p_ - 1) / (n_ - p_)) ** (1 - 1 / p_)
    br = mp_gamma(1 + n_ / 2) * mp_gamma(n_) / (mp_gamma(n_ / p_) * mp_gamma(1 + n_ - n_ / p_))
    return pref * br ** (1 / n_)


def _ground_state_rhs(r, y):
    q, dq = y
    return [dq, -2.0 * dq / r + q - q**3]


def _shoot(alpha, rmax=30.0):
    """Integrate the radial ODE from r=eps; return (event, r_end, sol)."""
    eps = 1e-8
    q0 = alpha + (alpha - alpha**3) * eps**2 / 6.0
    dq0 = (alpha - alpha**3) * eps / 3.0

    def hit_zero(r, y):
        return y[0]

    def blow_up(r, y):
        return y[0] - 10.0 * alpha

    hit_zero.terminal = True
    blow_up.terminal = True
    sol = solve_ivp(
        _ground_state_rhs,
        (eps, rmax),
        [q0, dq0],
        events=(hit_zero, blow_up),
        rtol=1e-12,
        atol=1e-14,
        dense_output=True,
        max_step=0.05,
    )
    crossed = len(sol.t_events[0]) > 0
    blew = len(sol.t_events[1]) > 0
    return crossed, blew, sol


def ground_state(n_bisect=200):
    """Bisect on Q(0) for the positive decaying ground state of -Delta Q + Q = Q^3 on R^3."""
    lo, hi = 1.0, 20.0  # lo: blows up to +infinity; hi: crosses zero
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        crossed, blew, _ = _shoot(mid)
        if crossed:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-15:
            break
    alpha = 0.5 * (lo + hi)
    _, _, sol = _shoot(alpha)
    return alpha, sol


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
    print(f"Talenti sharp Sobolev  C_S  (||w||_6 <= C_S ||grad w||_2, R^3) = {c_sob:.10f}")
    print(f"  equivalently S = C_S^-2 = {c_sob**-2:.10f}   (= 3(pi/2)^(4/3) = {3 * (np.pi / 2) ** (4 / 3):.10f})")
    c_gn_rig = c_sob**0.75
    print(f"rigorous GN constant via Hoelder+Sobolev: C_GN <= C_S^(3/4) = {c_gn_rig:.10f}")

    alpha, sol = ground_state()
    rend = sol.t[-1]
    m2, g2, f4 = ground_state_norms(sol, min(rend, 30.0))
    print(f"\nground state of -Delta Q + Q = Q^3 on R^3: Q(0) = {alpha:.10f}, integrated to r = {rend:.3f}")
    print(f"  M = ||Q||_2^2      = {m2:.8f}")
    print(f"  G = ||grad Q||_2^2 = {g2:.8f}    Pohozaev predicts G = 3M = {3 * m2:.8f}  (rel. err {abs(g2 - 3 * m2) / m2:.2e})")
    print(f"  F = ||Q||_4^4      = {f4:.8f}    Pohozaev predicts F = 4M = {4 * m2:.8f}  (rel. err {abs(f4 - 4 * m2) / m2:.2e})")
    c_gn_sharp = f4**0.25 / (g2**0.375 * m2**0.125)
    c_gn_closed = 4**0.25 * 3**-0.375 * m2**-0.25
    print(f"  sharp C_GN = ||Q||_4 / (||grad Q||_2^{{3/4}} ||Q||_2^{{1/4}}) = {c_gn_sharp:.10f}")
    print(f"  closed form 4^{{1/4}} 3^{{-3/8}} M^{{-1/4}}                     = {c_gn_closed:.10f}")

    c_p_rig = 3.0 ** (8.0 / 3.0)
    print(f"\nrigorous pressure constant (Iwaniec-Martin ||R_j||_{{L^3}} = sqrt3): C_p <= 3^(8/3) = {c_p_rig:.6f}")
    print("\nkappa / nu  =  4 / (3 C_p C_GN^{4/3})   [D2]")
    rows = []
    for gn_name, gn in (("C_S^(3/4) (rigorous)", c_gn_rig), ("sharp (ground state)", c_gn_sharp)):
        for cp_name, cp in (("3^(8/3) (rigorous)", c_p_rig),):
            kap = 4.0 / (3.0 * cp * gn ** (4.0 / 3.0))
            rows.append((gn_name, cp_name, gn ** (4.0 / 3.0), kap))
    print(f"| {'C_GN':<22} | {'C_p':<19} | {'C_GN^{4/3}':>10} | {'kappa/nu':>9} |")
    print(f"| {'-' * 22} | {'-' * 19} | {'-' * 10} | {'-' * 9} |")
    for a, b, c, d in rows:
        print(f"| {a:<22} | {b:<19} | {c:>10.6f} | {d:>9.6f} |")
    return {"C_S": c_sob, "C_GN_rig": c_gn_rig, "C_GN_sharp": c_gn_sharp, "C_p_rig": c_p_rig, "M": m2}


# ----------------------------------------------------------------------------
# Spectral machinery on T^3 = [0,2pi)^3
# ----------------------------------------------------------------------------
class Grid:
    def __init__(self, n=NGRID):
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
        ph = num / self.k2inv
        ph[0, 0, 0] = 0.0
        return np.real(np.fft.ifftn(ph))

    def random_solenoidal(self, rng, slope=1.5, kmax=6, aniso=(1.0, 1.0, 1.0)):
        """Random real mean-zero divergence-free field with a decaying spectrum."""
        kmag = np.sqrt(self.k2)
        env = np.where((kmag > 0) & (kmag <= kmax), (kmag + 1e-12) ** (-slope), 0.0)
        v = []
        for a in aniso:
            re = rng.standard_normal(self.k2.shape)
            im = rng.standard_normal(self.k2.shape)
            fh = a * env * (re + 1j * im)
            v.append(np.real(np.fft.ifftn(fh)))
        v = self.project(v)
        for c in v:
            c -= c.mean()
        nrm = np.sqrt(sum(self.integral(c**2) for c in v))
        return [c / nrm for c in v]


def functionals(g, u, nu=NU):
    """Return the pieces of d/dt int |u|^3 for a solenoidal field u on T^3."""
    absu = np.sqrt(sum(c**2 for c in u))
    safe = np.maximum(absu, 1e-300)
    du = [g.grad(c) for c in u]  # du[i][j] = d_j u_i
    gradsq = sum(du[i][j] ** 2 for i in range(3) for j in range(3))
    # grad|u| = (u_i grad u_i)/|u|, computed from spectral derivatives of the smooth u
    grad_absu = [sum(u[i] * du[i][j] for i in range(3)) / safe for j in range(3)]
    p = g.pressure(u)
    # dissipation D = 3 int |u| |grad u|^2 + 3 int |u| |grad|u||^2 = 3 int|u||grad u|^2 + (4/3) X^2
    xsq = (9.0 / 4.0) * g.integral(absu * sum(c**2 for c in grad_absu))
    diss = 3.0 * g.integral(absu * gradsq) + (4.0 / 3.0) * xsq
    press = 3.0 * g.integral(p * sum(u[j] * grad_absu[j] for j in range(3)))
    transport = 3.0 * g.integral(absu * sum(u[i] * sum(u[j] * du[i][j] for j in range(3)) for i in range(3)))
    lap = [np.real(np.fft.ifftn(-g.k2 * np.fft.fftn(c))) for c in u]
    conv = [sum(u[j] * du[i][j] for j in range(3)) for i in range(3)]
    ut = g.project([nu * lap[i] - conv[i] for i in range(3)])
    lhs = 3.0 * g.integral(absu * sum(u[i] * ut[i] for i in range(3)))
    l3 = g.integral(absu**3) ** (1.0 / 3.0)
    l6sq = g.integral(absu**6) ** (1.0 / 3.0)
    p3 = g.integral(np.abs(p) ** 3) ** (1.0 / 3.0)
    return {
        "lhs": lhs,
        "diss": diss,
        "press": press,
        "transport": transport,
        "X2": xsq,
        "L3": l3,
        "L6sq": l6sq,
        "p3": p3,
        "rhs": -nu * diss + press,
    }




def part2(g, rng):
    """Verify D1. With a nonzero mean, |u| is bounded away from 0, hence smooth, and the
    grid quadrature is spectrally accurate; with zero mean, |u| has conical zeros and the
    quadrature degrades - both are shown."""
    print()
    print("=" * 78)
    print("PART 2 - spectral verification of the exact identity D1 on T^3 (N = %d)" % g.n)
    print("=" * 78)
    print("| field                        |  min |u|/max|u| | d/dt int|u|^3 (direct) |  rel. err vs -nu D + P |  transport |")
    print("| ---------------------------- | --------------- | ---------------------- | ---------------------- | ---------- |")
    cases = [
        ("mean (2,0,0), slope 1.5, k<=3", (1.5, 3, (1, 1, 1)), 2.0),
        ("mean (2,1,0), slope 0.5, k<=2", (0.5, 2, (1, 1, 1)), 2.0),
        ("mean (5,0,0), slope 1.0, k<=4", (1.0, 4, (3, 1, 0.3)), 5.0),
        ("zero mean,    slope 1.5, k<=3", (1.5, 3, (1, 1, 1)), 0.0),
    ]
    out = []
    for name, (slope, kmax, aniso), mean in cases:
        u = g.random_solenoidal(rng, slope=slope, kmax=kmax, aniso=aniso)
        u = [2.0 * c for c in u]
        if mean:
            u[0] = u[0] + mean
            if "2,1,0" in name:
                u[1] = u[1] + 1.0
        f = functionals(g, u)
        absu = np.sqrt(sum(c**2 for c in u))
        ratio = absu.min() / absu.max()
        rel = abs(f["lhs"] - f["rhs"]) / max(abs(f["lhs"]), 1e-30)
        out.append(rel)
        print(
            f"| {name:<28} | {ratio:>15.4f} | {f['lhs']:>22.12e} | {rel:>22.3e} | {f['transport']:>10.2e} |"
        )
    print("\nThe three fields with |u| bounded away from zero confirm D1 to spectral accuracy;")
    print("the zero-mean field loses accuracy exactly because |u| is only Lipschitz at its zeros.")
    return out


# ----------------------------------------------------------------------------
# Parts 3-4: optimizing the scale-invariant ratios over a divergence-free basis
# ----------------------------------------------------------------------------
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


def lean_functionals(g, u, nu=NU):
    """Only the pieces the optimizer needs: D, P, X^2, ||u||_3."""
    absu = np.sqrt(sum(c**2 for c in u))
    safe = np.maximum(absu, 1e-300)
    du = [g.grad(c) for c in u]
    gradsq = sum(du[i][j] ** 2 for i in range(3) for j in range(3))
    grad_absu = [sum(u[i] * du[i][j] for i in range(3)) / safe for j in range(3)]
    p = g.pressure(u)
    xsq = (9.0 / 4.0) * g.integral(absu * sum(c**2 for c in grad_absu))
    diss = 3.0 * g.integral(absu * gradsq) + (4.0 / 3.0) * xsq
    press = 3.0 * g.integral(p * sum(u[j] * grad_absu[j] for j in range(3)))
    return diss, press, xsq, g.integral(absu**3) ** (1.0 / 3.0)


def theta(g, fields, c):
    """Theta(u) = nu D(u) ||u||_3 / |P(u)|: the L^3 norm at which the amplitude family a*u
    first has d/dt ||a u||_3^3 > 0. Scale-invariant in a (D is 3- and P is 4-homogeneous);
    P is odd in u and D is even, so the sign of P costs nothing."""
    u = assemble(fields, c)
    nrm = np.sqrt(sum(g.integral(ci**2) for ci in u))
    if nrm < 1e-8:
        return 1e12
    u = [ci / nrm for ci in u]
    diss, press, _, l3 = lean_functionals(g, u)
    if abs(press) < 1e-14:
        return 1e12
    return NU * diss * l3 / abs(press)


def optimize_theta(gopt, kmax2, rng, nrestart, maxfev):
    fields = basis_fields(gopt, kmax2=kmax2)
    best, bestc = np.inf, None
    for _ in range(nrestart):
        c0 = rng.standard_normal(len(fields))
        res = minimize(
            lambda c: theta(gopt, fields, c),
            c0,
            method="Nelder-Mead",
            options={"maxiter": maxfev, "maxfev": maxfev, "xatol": 1e-7, "fatol": 1e-9},
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
    print(f"PART 3 - bracketing the sharp threshold by optimization (N={gopt.n} search, N={gfine.n} evaluation)")
    print("=" * 78)
    rows, bestc = [], {}
    for kmax2, nrestart, maxfev in ((1, 3, 700), (2, 2, 1400)):
        m, best_coarse, c = optimize_theta(gopt, kmax2, rng, nrestart, maxfev)
        _, f = evaluate(gfine, kmax2, c)
        thr = NU * f["diss"] * f["L3"] / f["press"]
        kval = f["press"] / (f["X2"] * f["L3"])
        rows.append((kmax2, m, best_coarse, thr, kval, f["p3"] / f["L6sq"]))
        bestc[kmax2] = c
    print("| basis     | modes | Theta at N=%2d | Theta at N=%2d | K = P/(X^2 ||u||_3) | ||p||_3/|| |u|^2 ||_3 |"
          % (gopt.n, gfine.n))
    print("| --------- | ----- | ------------- | ------------- | ------------------- | --------------------- |")
    for kmax2, m, bc, thr, kval, cp in rows:
        print(f"| |k|^2<={kmax2}   | {m:>5d} | {bc:>13.6f} | {thr:>13.6f} | {kval:>19.6f} | {cp:>21.6f} |")
    thr_ub = min(r[3] for r in rows)
    k_lb = max(r[4] for r in rows)
    cp_lb = max(r[5] for r in rows)
    kub = 2.0 * consts["C_p_rig"] * consts["C_GN_rig"] ** (4.0 / 3.0)
    kub_s = 2.0 * consts["C_p_rig"] * consts["C_GN_sharp"] ** (4.0 / 3.0)
    kap, kap_s = (8.0 / 3.0) / kub, (8.0 / 3.0) / kub_s
    print()
    print(f"K* = sup_u P/(X^2 ||u||_3):   {k_lb:.6f} (optimized lower bd)  <=  K*  <=  {kub:.6f} (proved, C_GN=C_S^(3/4))")
    print(f"                                                                        <=  {kub_s:.6f} (proved, sharp C_GN)")
    print(f"   looseness of the constant chain: {kub / k_lb:.0f}x (rigorous C_GN), {kub_s / k_lb:.0f}x (sharp C_GN)")
    print(f"C_p = ||R_iR_j(u_iu_j)||_3 / || |u|^2 ||_3:  {cp_lb:.6f} (de Leeuw lower bd)  <=  C_p  <=  {consts['C_p_rig']:.6f} (proved)")
    print()
    print("BRACKET on the sharp threshold kappa_sharp = inf_u Theta(u):")
    print(f"   proved   {kap:.6f} nu (rigorous C_GN) / {kap_s:.6f} nu (sharp C_GN)  <=  kappa_sharp  <=  {thr_ub:.6f} nu (optimized, torus)")
    return {"thr_ub": thr_ub, "K_lb": k_lb, "C_p_lb": cp_lb, "kap": kap, "kap_s": kap_s, "c": bestc}


def part4(gfine, c):
    print()
    print("=" * 78)
    print("PART 4 - an explicit field with d/dt ||u||_{L^3}^3 > 0 at nu = 1")
    print("=" * 78)
    u, f = evaluate(gfine, 1, c)
    a_crit = NU * f["diss"] / f["press"]
    cc = c / np.linalg.norm(c)
    print("v: unit-L^2 combination of the twelve |k| = 1 divergence-free modes found in Part 3,")
    print("signed so that P(v) > 0. Coefficients (k = e1,e2,e3; two polarizations; cos then sin):")
    print("  " + ", ".join(f"{v:+.5f}" for v in cc))
    print(f"  D(v) = {f['diss']:.6f}   P(v) = {f['press']:.6f}   ||v||_3 = {f['L3']:.6f}")
    print(f"  d/dt ||a v||_3^3 = -a^3 nu D(v) + a^4 P(v) > 0  iff  a > nu D(v)/P(v) = {a_crit:.6f}")
    print()
    print("| a / a_crit |  ||a v||_3 |      -nu D + P (identity) |   direct d/dt int|u|^3 |  rel. err |")
    print("| ---------- | ---------- | ------------------------- | ---------------------- | --------- |")
    for r in (0.5, 1.0, 1.5, 3.0):
        a = r * a_crit
        ff = functionals(gfine, [a * ci for ci in u])
        rel = abs(ff["lhs"] - ff["rhs"]) / max(abs(ff["lhs"]), abs(ff["rhs"]), 1e-30)
        print(f"| {r:>10.2f} | {ff['L3']:>10.5f} | {ff['rhs']:>25.10e} | {ff['lhs']:>22.10e} | {rel:>9.2e} |")
    print()
    print("=> ||u||_{L^3} is NOT a Lyapunov functional for 3D Navier-Stokes at nu = 1;")
    print(f"   for this family it starts to increase at ||u||_3 = {a_crit * f['L3']:.6f}.")
    return a_crit * f["L3"]


def main():
    rng = np.random.default_rng(SEED)
    consts = part1()
    gfine = Grid(48)
    part2(gfine, rng)
    res = part3(Grid(24), gfine, rng, consts)
    part4(gfine, res["c"][1])
    print()
    print("done.")


if __name__ == "__main__":
    main()
