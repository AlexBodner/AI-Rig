"""Line C: Tao's cascade program on dyadic (shell) models with exact 3D Navier-Stokes scaling.

Model (n = 0..N-1, a_{-1} = a_N = 0), lambda = 2^{5/2}, dissipation exponent s (s = 1 is the true
Navier-Stokes Laplacian):

    da_n/dt = -nu * 2^(2*s*n) * a_n + lam^n * a_{n-1}^2 - lam^(n+1) * a_n * a_{n+1}

Parts:
  1  exact/symbolic checks (energy telescoping, scaling symmetry, uniqueness of lam, self-similar sol)
  2  inviscid front diagnostics: escape time, event-state H^1, front amplitude slope, boundary layer
  3  viscous runs with the true dissipation 4^n: (nu, amplitude) table + the (nu, A) collapse check
  4  dissipation-threshold scan in s, against the predicted stall shell n_d = log2(1/nu)/(2s - 5/3)
  5  engineered multi-mode ("Tao gadget") cascade: exact energy conservation + front-regime transit

Run: python3 line_c_cascade.py   (numpy/scipy/sympy; a few minutes, fixed seeds, no output files)
"""

import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

LAM = 2.0**2.5
SLOPE_K41 = -5.0 / 6.0  # log2 amplitude slope of the constant-flux (K41) profile, C5
SLOPE_SCALEINV = -0.5  # log2 amplitude slope of the scale-invariant profile, C2
SLOPE_FRONT = 0.0  # log2 amplitude slope of the maximal ("front", Tao) regime, C10


def rhs(_t, a, nu, s, n):
    lam_n = LAM ** np.arange(n)
    am = np.concatenate(([0.0], a[:-1]))
    ap = np.concatenate((a[1:], [0.0]))
    return -nu * 2.0 ** (2.0 * s * np.arange(n)) * a + lam_n * am**2 - LAM * lam_n * a * ap


def jac(_t, a, nu, s, n):
    lam_n = LAM ** np.arange(n)
    ap = np.concatenate((a[1:], [0.0]))
    j = np.zeros((n, n))
    d = -nu * 2.0 ** (2.0 * s * np.arange(n)) - LAM * lam_n * ap
    np.fill_diagonal(j, d)
    for i in range(1, n):
        j[i, i - 1] = 2.0 * lam_n[i] * a[i - 1]
    for i in range(n - 1):
        j[i, i + 1] = -LAM * lam_n[i] * a[i]
    return j


def run(nu, s, n, amp=1.0, t_end=200.0, rtol=1e-9, atol=1e-16, esc_tol=1e-9):
    """Integrate from a delta at shell 0. Terminal event: last shell exceeds esc_tol*amp."""
    a0 = np.zeros(n)
    a0[0] = amp

    def escaped(t, a, *args):
        return a[-1] - esc_tol * amp

    escaped.terminal = True
    escaped.direction = 1
    t_grid = np.unique(np.clip(np.concatenate(([0.0], np.logspace(-14.0, np.log10(t_end), 1500))), 0.0, t_end))
    sol = solve_ivp(
        rhs,
        (0.0, t_end),
        a0,
        method="Radau",
        jac=jac,
        args=(nu, s, n),
        rtol=rtol,
        atol=atol,
        t_eval=t_grid,
        events=escaped,
    )
    a, t = sol.y, sol.t
    if len(sol.t_events[0]):  # append the exact event state: the t_eval grid stops just short of it
        a = np.concatenate((a, sol.y_events[0][0][:, None]), axis=1)
        t = np.concatenate((t, sol.t_events[0][:1]))
    kn = np.arange(n)
    h1 = (4.0 ** kn[:, None] * a**2).sum(axis=0)
    dsp = 2.0 ** (2.0 * s * kn)[:, None] * a**2
    diss = 2.0 * nu * dsp.sum(axis=0)
    it = int(np.argmax(diss))
    ener = (a**2).sum(axis=0)
    return {
        "t": t, "a": a, "E": ener, "H1": h1, "n": n, "nu": nu, "s": s, "amp": amp,
        "escaped": len(sol.t_events[0]) > 0, "n_d": int(np.argmax(dsp[:, it])),
        "n_dmax": int(np.argmax(dsp, axis=0).max()), "eps": float(diss.max()),
        "maxH1": float(h1.max()), "amin": float(a.min()), "E0": float(ener[0]),
        "Eend": float(ener[-1]), "tend": float(t[-1]), "nfev": int(sol.nfev),
        "t_esc": float(sol.t_events[0][0]) if len(sol.t_events[0]) else np.inf,
    }  # fmt: skip


def slope_fit(a_end, lo, hi):
    """Least-squares log2 slope of the shell profile over shells lo..hi (inclusive)."""
    idx = np.arange(max(lo, 0), min(hi, len(a_end) - 1) + 1)
    idx = idx[a_end[idx] > 1e-18]
    if len(idx) < 4:
        return np.nan, 0
    return float(np.polyfit(idx, np.log2(a_end[idx]), 1)[0]), len(idx)


def inertial_slope(res):
    """log2 slope per shell of the profile at the instant of maximum dissipation, over 1..n_d-1."""
    a, s, n = res["a"], res["s"], res["n"]
    dsp = 2.0 ** (2.0 * s * np.arange(n))[:, None] * a**2
    it = int(np.argmax(dsp.sum(axis=0)))
    return slope_fit(a[:, it], 1, res["n_d"] - 1)


def part1_exact():
    print("=" * 78)
    print("PART 1 -- exact symbolic checks (sympy)")
    print("=" * 78)
    n = 8
    lam, nu, s = sp.symbols("lam nu s", positive=True)
    a = sp.symbols(f"a0:{n}")
    f = []
    for i in range(n):
        am = a[i - 1] if i >= 1 else 0
        ap = a[i + 1] if i + 1 < n else 0
        f.append(lam**i * am**2 - lam ** (i + 1) * a[i] * ap)
    dE = sp.expand(2 * sum(a[i] * f[i] for i in range(n)))
    print(f"  (1a) inviscid energy: d/dt sum a_n^2 = {dE}   -> telescoping exact: {dE == 0}")
    dEv = sp.expand(2 * sum(a[i] * (f[i] - nu * 2 ** (2 * s * i) * a[i]) for i in range(n)))
    dis = 2 * nu * sum(2 ** (2 * s * i) * a[i] ** 2 for i in range(n))
    print(f"  (1b) viscous: d/dt sum a_n^2 + 2 nu sum 2^(2 s n) a_n^2 = {sp.simplify(dEv + dis)}")
    print("       ((1a) is inviscid; (1b) holds with lam, nu and s all symbolic. The dissipation is")
    print("        diagonal, so it plays no part in the telescoping.)")

    # (1c) discrete scaling symmetry b_n(t) = mu a_{n+1}(sig t) with dissipation 2^{2 s n}
    mu, sig = sp.symbols("mu sig", positive=True)
    t = sp.Symbol("t")
    A = [sp.Function(f"A{i}")(sig * t) for i in range(n)]
    lam_s = sp.Rational(5, 2)  # lam = 2^{5/2}
    two = sp.Integer(2)

    def eq_resid(idx, arr, nu_):  # residual of the model at shell idx for array arr
        am = arr[idx - 1] if idx >= 1 else 0
        ap = arr[idx + 1] if idx + 1 < len(arr) else 0
        return (
            -nu_ * two ** (2 * s * idx) * arr[idx]
            + two ** (lam_s * idx) * am**2
            - two ** (lam_s * (idx + 1)) * arr[idx] * ap
        )

    i0 = 3
    rhs_shift = mu * sig * eq_resid(i0 + 1, A, nu)  # what dA_{n+1}/dt equals, times mu*sig
    B = [mu * A[j + 1] for j in range(n - 1)]
    want = eq_resid(i0, B, nu)
    cond = sp.simplify(sp.expand(rhs_shift - want))
    sol_ms = sp.solve([sp.Eq(sig * two ** (2 * s), 1), sp.Eq(sig * two**lam_s / mu, 1)], [mu, sig], dict=True)[0]
    print(f"  (1c) shift symmetry needs sig = {sol_ms[sig]}, mu = {sol_ms[mu]}")
    chk = sp.simplify(cond.subs(sol_ms))
    print(f"       residual after substitution at the INTERIOR shell n = {i0}: {sp.simplify(chk.subs(s, 1))}  (at s=1)")
    b0 = mu * sig * eq_resid(1, A, nu) - eq_resid(0, B, nu)  # boundary shell n = 0 of the half-lattice
    print(f"       residual at the BOUNDARY shell n = 0 (a_{{-1}} = 0): {sp.simplify(sp.expand(b0).subs(sol_ms))}")
    print("       -> = mu^2 a_0(sig t)^2, nonzero unless a_0 == 0: the map solves the model only on the")
    print("          bi-infinite lattice (or at interior shells).  (lam = 2^(5/2) is hard-coded in (1c).)")
    print(
        f"       at s = 1: sig = {sol_ms[sig].subs(s, 1)} = 1/4, mu = {sp.nsimplify(sol_ms[mu].subs(s, 1))} = 2^(1/2)"
    )
    print("       -> t -> t/4 and a -> 2^(1/2) a per shell: exactly the L^2-normalized NS scaling.")

    # (1d) uniqueness of lam: sig=1/4 forced by 4^n, mu = sig*lam must equal 2^{1/2}
    lam_sol = sp.solve(sp.Eq(sp.Rational(1, 4) * lam, sp.sqrt(2)), lam)[0]
    print(f"  (1d) mu = sig*lam = lam/4 = 2^(1/2)  =>  lam = {lam_sol} = 2^(5/2) = {float(lam_sol):.6f}  (unique)")

    # (1e) scale-invariant profile and constant-flux profile
    c, x = sp.symbols("c x", positive=True)
    prof = c * 2 ** (-sp.Rational(1, 2) * sp.Symbol("nn"))
    print(
        f"  (1e) a_n = c 2^(-n/2) is the fixed point of the shift map: 2^(1/2)*c*2^(-(n+1)/2) = c*2^(-n/2): "
        f"{sp.simplify(sp.sqrt(2) * prof.subs(sp.Symbol('nn'), sp.Symbol('nn') + 1) - prof) == 0}"
    )
    nn = sp.Symbol("nn")
    resid = 2 ** (lam_s * nn) * (c * 2 ** (-(nn - 1) / 2)) ** 2 - 2 ** (lam_s * (nn + 1)) * (c * 2 ** (-nn / 2)) * (
        c * 2 ** (-(nn + 1) / 2)
    )
    print(f"       but it is NOT stationary: residual = {sp.simplify(sp.powsimp(sp.expand(resid), force=True))}")
    flux = sp.simplify(
        2 ** (lam_s * (nn + 1)) * (c * 2 ** (-sp.Rational(5, 6) * nn)) ** 2 * (c * 2 ** (-sp.Rational(5, 6) * (nn + 1)))
    )
    print(
        f"  (1f) constant-flux profile a_n = c 2^(-5n/6): flux lam^(n+1) a_n^2 a_(n+1) = {sp.powsimp(flux, force=True)} (n-independent)"
    )

    # (1g) exact self-similar blowup on the bi-infinite lattice: a_n = lam^{-n}/((lam^2-1)(T-t))
    T = sp.Symbol("T")
    cc = 1 / (lam**2 - 1)

    def an(i):
        return cc * lam ** (-i) / (T - t)

    resid_ss = sp.simplify(sp.diff(an(nn), t) - (lam**nn * an(nn - 1) ** 2 - lam ** (nn + 1) * an(nn) * an(nn + 1)))
    print(f"  (1g) a_n(t) = lam^(-n)/((lam^2-1)(T-t)) on n in Z solves the inviscid model: residual = {resid_ss}")
    print("       its energy sum_{n in Z} a_n^2 = +infinity (geometric divergence as n -> -infinity).")

    # (1h) the power-law comparison lemma of C10: the only part of the regime dichotomy that is a theorem
    b, dd = sp.symbols("b d", positive=True)
    lam_d = 1 + dd / 2  # log2 of lam in dimension d; = 5/2 at d = 3
    print("  (1h) power-law lemma (C10 core).  For A(n) = c 2^(-b n) with c > 0, in dimension d with")
    print("       log2(lam) = 1 + d/2, both statements are geometric series/sequences:")
    ratio_sum = sp.Rational(1, 1) * 2 ** (b - lam_d)  # ratio of sum_n 1/(lam^n A(n))
    ratio_beat = 2 ** (lam_d - b - 2 * s)  # ratio of lam^n A(n) / (nu 2^(2 s n))
    print(f"         sum_n 1/(lam^n A(n)) converges iff ratio {ratio_sum} < 1 iff b < {lam_d} (= 5/2 at d = 3)")
    print(f"         lam^n A(n)/(nu 2^(2 s n)) -> oo iff ratio {ratio_beat} > 1 iff b < {sp.expand(lam_d - 2 * s)}")
    for name, bval in (("constant flux (C5)", (2 + dd) / 6), ("maximal front (C10)", sp.Integer(0))):
        thr = sp.solve(sp.Eq(bval, lam_d - 2 * s), s)[0]
        print(
            f"         b = {sp.simplify(bval)}  [{name}]: beats dissipation iff s < {sp.simplify(thr)}"
            f"   -> d=3: {sp.nsimplify(thr.subs(dd, 3))}, d=2: {sp.nsimplify(thr.subs(dd, 2))},"
            f" d=1: {sp.nsimplify(thr.subs(dd, 1))}"
        )
    print("       This lemma is PROVED-HERE; the identification of A(n) with a cascade front is not.")


def part2_inviscid():
    print("\n" + "=" * 78)
    print("PART 2 -- inviscid model (nu = 0): finite-time front escape")
    print("=" * 78)
    print("  H^1(grid) is the value at the last t_eval point, a FIXED time for every N (the grid does not")
    print("  depend on N), so its constancy across N is automatic; H^1(event) is the value at the escape")
    print("  event itself, and 'last shell' is the single term 4^(N-1) a_(N-1)^2 inside it.")
    print("  N   T_esc(N)   t(grid)    H^1(grid)  H^1(event)  last shell   E/E_0    slope 2..8  slope 9..N-6")
    prof32 = None
    for n in (16, 20, 24, 28, 32, 36):
        r = run(0.0, 1.0, n, t_end=5.0)
        a_ev = r["a"][:, -1]
        si, wi = slope_fit(a_ev, 2, 8)
        so, wo = slope_fit(a_ev, 9, n - 6)
        h1_ev = float((4.0 ** np.arange(n) * a_ev**2).sum())
        last = float(4.0 ** (n - 1) * a_ev[n - 1] ** 2)
        if n == 32:
            prof32 = a_ev
        print(
            f" {n:3d}  {r['t_esc']:.7f} {r['t'][-2]:.7f} {r['H1'][-2]:10.5g} {h1_ev:11.4f} {last:11.4f} "
            f"{r['Eend'] / r['E0']:.12f} {si:8.4f} ({wi:2d}) {so:8.4f} ({wo:2d})"
        )
    print(f"  reference slopes: constant-flux (K41) {SLOPE_K41:.4f}; scale-invariant (C2) {SLOPE_SCALEINV:.4f};")
    print(f"                    maximal front regime (C10) {SLOPE_FRONT:.4f};")
    print(f"                    infinite-energy self-similar tail a_n ~ lam^(-n) (C4): {-2.5:.4f}")
    print("\n  sensitivity of T_esc to the escape threshold esc_tol, at N = 32 (the event is a threshold")
    print("  crossing at the last retained shell, so only the digits stable in this column are real):")
    for et in (1e-4, 1e-6, 1e-9, 1e-12, 1e-14):
        r = run(0.0, 1.0, 32, t_end=5.0, esc_tol=et)
        print(f"    esc_tol = {et:7.1g}:  T_esc = {r['t_esc']:.7f}")
    d32 = np.diff(np.log2(np.maximum(prof32, 1e-300)))
    print("\n  per-shell log2 differences of the N = 32 profile at the event (shells 0->1, 1->2, ...):")
    print("   " + " ".join(f"{x:.2f}" for x in d32[:20]))
    print("   " + " ".join(f"{x:.2f}" for x in d32[20:]))
    print("  -> a flat inner range near -1.55, then a non-monotone shoulder over the last ~10 shells:")
    print("     the truncation boundary layer (shell N-1 has no drain), not a power-law tail.")
    print("  verdict: T_esc(N) converges to 6 digits as N grows -> the front reaches n = infinity in a")
    print("           finite, truncation-independent time (inviscid front escape). The profile at that")
    print("           instant is NOT a clean power law: the inner range sits near -3/2 while the outer")
    print("           range drifts shallower as the truncation is extended, so whether H^1 stays finite")
    print("           at T_esc is NOT decided by these runs (recorded as inconclusive in C6).")


def part3_viscous():
    print("\n" + "=" * 78)
    print("PART 3 -- true Navier-Stokes dissipation 4^n (s = 1): does the front escape?")
    print("=" * 78)
    print("   nu      A    N    T_end     max H^1   n_d  min_n a_n   E_end/E_0   verdict")
    rows = []
    for nu in (0.3, 0.1, 0.03, 0.01):
        for amp in (1.0, 3.0):
            for n in (24, 32):
                r = run(nu, 1.0, n, amp=amp)
                v = "front escapes (check N)" if r["escaped"] else f"dissipation wins (stall at n_d={r['n_d']})"
                print(
                    f" {nu:6.3g} {amp:5.1f} {n:3d} {r['tend']:9.4g} {r['maxH1']:10.4g} {r['n_d']:4d} "
                    f"{r['amin']:10.2e} {r['Eend'] / r['E0']:10.3e}   {v}"
                )
                rows.append(r)
    print(f"  (most negative amplitude over all runs above: min_n a_n = {min(r['amin'] for r in rows):.2e};")
    print("   positivity is preserved up to roundoff.  Whether this is the solution class of the")
    print("   Barbato-Morandin-Romito well-posedness theorem is NOT verified here: that paper was")
    print("   never read, only its title and one abstract snippet.  See C7 in the notes.)")
    print("  (max H^1 below and above is a maximum over the 1500-point log time grid, not over the")
    print("   trajectory; at a stall the peak is broad and the two agree to plotting accuracy.)")
    print("\n  (nu, A) collapse check [C3]: b_n(t) = mu a_n(mu t) solves the model with (nu, A) -> (mu nu, mu A),")
    print("  so max H^1 scales by mu^2 and the stall shell n_d is unchanged.  mu = 3:")
    for nu in (0.1, 0.03, 0.01):
        r1 = run(nu, 1.0, 32, amp=1.0)
        r2 = run(3.0 * nu, 1.0, 32, amp=3.0)
        print(
            f"    (nu={nu:6.3g}, A=1): H^1max={r1['maxH1']:10.5g}, n_d={r1['n_d']:3d}   "
            f"(nu={3 * nu:6.3g}, A=3): H^1max={r2['maxH1']:10.5g}, n_d={r2['n_d']:3d}   "
            f"ratio={r2['maxH1'] / r1['maxH1']:7.4f} (predicted 9)"
        )
    print("\n  inertial-range profile slope at the instant of peak dissipation (predicted -5/6 = -0.8333;")
    print("  the scale-invariant profile of C2 is -1/2, and a front cascade would give slope 0):")
    for nu in (0.03, 0.01, 0.003):
        r = run(nu, 1.0, 34)
        sl, w = inertial_slope(r)
        print(f"    nu = {nu:7.4g}: n_d = {r['n_d']:3d}, slope = {sl:8.4f} over {w} shells")


def part4_threshold():
    print("\n" + "=" * 78)
    print("PART 4 -- dissipation-threshold scan: replace 4^n by 2^(2 s n)")
    print("=" * 78)
    print("  4a. s-scan at nu(s) = 2^(-14 max(2s-5/3, 2/15)), N = 24 and 32 (truncation control)")
    print("    s      nu      n_d(24) n_d(32)  T_esc(24)  T_esc(32)  esc  verdict")
    for s in (0.70, 0.80, 0.8333, 0.90, 1.00, 1.10, 1.25, 1.35):
        nu = 2.0 ** (-14.0 * max(2.0 * s - 5.0 / 3.0, 2.0 / 15.0))
        ra = run(nu, s, 24)
        rb = run(nu, s, 32)
        if ra["escaped"] and rb["escaped"]:
            drift = 100.0 * (rb["t_esc"] - ra["t_esc"]) / ra["t_esc"]
            v = f"escapes at both N; T_esc drift 24->32 = {drift:+.2f}%"
        elif not (ra["escaped"] or rb["escaped"]):
            v = "dissipation wins (no escape at either N)"
        else:
            v = "borderline (escape at one N only)"
        print(
            f" {s:6.4f} {nu:8.3g} {ra['n_d']:6d} {rb['n_d']:7d} {ra['t_esc']:10.5f} {rb['t_esc']:10.5f}"
            f"  {str(ra['escaped'])[0]}/{str(rb['escaped'])[0]}  {v}"
        )
    print("\n  4c. truncation dependence of the escape time at the two sub-threshold rows, nu = 0.274.")
    print("      (Compare the inviscid model in part 2, where T_esc is fixed to 6 digits from N = 24 on.)")
    for s in (0.70, 0.80):
        seq = []
        for n in (18, 22, 26, 30, 34):
            seq.append(run(0.274, s, n)["t_esc"])
        drift = 100.0 * (seq[-1] - seq[0]) / seq[0]
        print(f"    s = {s:4.2f}, N = 18,22,26,30,34: T_esc = " + ", ".join(f"{x:.6f}" for x in seq))
        print(
            f"      -> monotone {'increasing' if all(np.diff(seq) > 0) else 'non-monotone'},"
            f" total drift {drift:+.2f}% and still rising at N = 34."
        )

    print("\n  4b. stall-shell law   n_d = log2(1/nu)/(2s - 5/3) + const   (K41 balance; diverges at s = 5/6)")
    print("     s     nu-range          n_d values           fitted slope   predicted 1/(2s-5/3)")
    grid = (
        (1.00, (0.3, 0.06, 0.012, 0.0024)),
        (1.10, (0.05, 5e-3, 5e-4, 5e-5)),
        (1.25, (1e-2, 2e-4, 4e-6, 1e-7)),
        (1.40, (1e-3, 1e-5, 1e-7, 1e-9)),
    )
    svals, slopes = [], []
    for s, nus in grid:
        xs, ys, flag = [], [], ""
        for nu in nus:
            r = run(nu, s, 34)
            xs.append(np.log2(1.0 / nu))
            ys.append(r["n_d"])
            if r["n_d"] > 30 or r["escaped"]:
                flag = " (*)"
        sl = float(np.polyfit(xs, ys, 1)[0])
        svals.append(s)
        slopes.append(sl)
        print(
            f" {s:5.2f}  {nus[0]:8.3g}..{nus[-1]:8.3g}  {str(ys):22s} {sl:9.3f}      {1.0 / (2 * s - 5.0 / 3.0):9.3f}{flag}"
        )
    sv, iv = np.array(svals), 1.0 / np.array(slopes)
    aa, bb = np.polyfit(sv, iv, 1)
    print(f"\n  linear fit of 1/slope against s:  1/slope = {aa:.4f} s + {bb:.4f}   (predicted 2 s - 5/3)")
    print(
        f"  coefficient errors against the prediction: slope {100 * (aa - 2) / 2:+.1f}%,"
        f" intercept {100 * (bb + 5 / 3) / (-5 / 3):+.1f}%"
    )
    print(
        f"  numerically located blowup threshold s* = {-bb / aa:.4f}   (predicted 5/6 = {5 / 6:.4f},"
        f" error {100 * (-bb / aa - 5 / 6) / (5 / 6):+.1f}%)"
    )
    loo = []
    for k in range(len(sv)):
        m = np.arange(len(sv)) != k
        a2, b2 = np.polyfit(sv[m], iv[m], 1)
        loo.append(-b2 / a2)
    print(f"  leave-one-out range of s* over the four s values: {min(loo):.4f} .. {max(loo):.4f}")
    print("  n_d is an integer (argmax over shells), so each of the 16 inputs carries at least half a")
    print("  shell of quantization, and s* = 0.8446 is an extrapolation below the fitted range")
    print("  s in [1.00, 1.40]: the measurement is consistent with 5/6 at the few-percent level and")
    print("  should not be read as a one-percent determination.")
    print("  Lions exponent for this model (energy criticality) s = 5/4 = 1.2500; s* is well below it.")


def chain_rhs(_t, a, nu, s, m, cpat, nmod):
    """Chain of pairwise transfers, m modes per dyadic scale; c_{j+m} = lam * c_j (scale invariance)."""
    scale = np.arange(nmod) // m
    c = LAM**scale * np.asarray(cpat)[np.arange(nmod) % m]
    cp = np.concatenate((c[1:], [0.0]))
    am = np.concatenate(([0.0], a[:-1]))
    ap = np.concatenate((a[1:], [0.0]))
    return -nu * 2.0 ** (2.0 * s * scale) * a + c * am**2 - cp * a * ap


def run_chain(nu, s, m, cpat, nscales, amp=1.0, t_end=200.0):
    nmod = m * nscales
    scale = np.arange(nmod) // m
    a0 = np.zeros(nmod)
    a0[0] = amp

    def escaped(t, a, *args):
        return a[-1] - 1e-9 * amp

    escaped.terminal = True
    escaped.direction = 1
    t_grid = np.unique(np.clip(np.concatenate(([0.0], np.logspace(-14.0, np.log10(t_end), 1200))), 0.0, t_end))
    sol = solve_ivp(
        chain_rhs,
        (0.0, t_end),
        a0,
        method="Radau",
        args=(nu, s, m, cpat, nmod),
        rtol=1e-9,
        atol=1e-16,
        t_eval=t_grid,
        events=escaped,
    )
    dsp = 2.0 ** (2.0 * s * scale)[:, None] * sol.y**2
    it = int(np.argmax(dsp.sum(axis=0)))
    ener = (sol.y**2).sum(axis=0)
    return {
        "n_d": int(scale[int(np.argmax(dsp[:, it]))]),
        "escaped": len(sol.t_events[0]) > 0,
        "E0": float(ener[0]),
        "Eend": float(ener[-1]),
        "a": sol.y,
        "scale": scale,
        "it": it,
    }


def part5_gadget():
    print("\n" + "=" * 78)
    print("PART 5 -- engineered multi-mode cascade: does 'several modes per scale' change the exponent?")
    print("=" * 78)
    nmod = 7
    cs = sp.symbols(f"c1:{nmod + 1}", positive=True)
    a = sp.symbols(f"b0:{nmod}")
    f = [
        (cs[j] * (a[j - 1] if j >= 1 else 0) ** 2 - (cs[j + 1] * a[j] * a[j + 1] if j + 1 < nmod else 0))
        for j in range(nmod)
    ]
    dE = sp.expand(2 * sum(a[j] * f[j] for j in range(nmod)))
    print(f"  (5a) trilinear form of ANY pairwise chain: d/dt sum_j a_j^2 = {dE}  (vanishes identically: {dE == 0})")
    print("       -> every such chain, with any coefficients c_j, has the exact energy identity of C1.")
    print("\n  (5b) inertial-range slope per SCALE for m modes per scale, s = 1; nu chosen per pattern so")
    print("       that the stall scale lands inside the truncation (the O(1) offset is pattern-dependent)")
    print("     m  pattern C_i           nu        stall scale n_d   log2 slope per scale (pred -0.8333)")
    for m, cpat in ((1, (1.0,)), (2, (1.0, 1.0)), (2, (8.0, 0.125)), (2, (0.125, 8.0)),
                    (3, (1.0, 1.0, 1.0)), (3, (4.0, 1.0, 0.25))):  # fmt: skip
        best = None
        for nu in (3e-3, 1e-4, 3e-6, 1e-7, 3e-9):
            r = run_chain(nu, 1.0, m, cpat, 26)
            if best is None or abs(r["n_d"] - 13) < abs(best[1]["n_d"] - 13):
                best = (nu, r)
            if 8 <= r["n_d"] <= 20:
                break
        nu, r = best
        en = np.array([(r["a"][r["scale"] == n, r["it"]] ** 2).sum() for n in range(26)])
        sl, w = slope_fit(np.sqrt(np.maximum(en, 1e-300)), 1, max(r["n_d"] - 1, 4))
        print(f"     {m}  {str(tuple(cpat)):20s} {nu:9.3g} {r['n_d']:10d}          {sl:8.4f}  ({w} scales)")
    print("       -> the intra-scale pattern moves the O(1) offset only; the slope stays at the")
    print("          constant-flux value -5/6, never reaching the front value 0. See C10-C11.")


if __name__ == "__main__":
    t_start = time.time()
    part1_exact()
    part2_inviscid()
    part3_viscous()
    part4_threshold()
    part5_gadget()
    print(f"\n[total wall time {time.time() - t_start:.1f} s]")
