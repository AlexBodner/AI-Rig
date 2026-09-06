"""Line E: logarithmically supercritical hyperdissipation — checkable computations.

Four independent parts, all printed to stdout:

Part 1  exact exponent arithmetic behind E1 (Lions) and E5 (the dyadic log budget), in sympy.
Part 2  the weight-cancellation identity of E6: for every weight h the Gronwall coefficient is g(N_0)^4.
Part 3  the Osgood ODE y' = phi(t) y g(y)^p: blowup iff int dy/(y g(y)^p) converges.
Part 4  a Katz-Pavlovic / Cheskidov dyadic shell model with dissipation lambda^{beta n}/g(lambda^n)^2,
        which realizes the cascade energy budget of E4 (sum_n g(2^n)^{-2}) directly.

Run:  python3 line_e_frontier.py    (from research/navier-stokes; < 5 minutes, 4 CPUs, no GPU)
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

SEED = 20260906
rng = np.random.default_rng(SEED)


def rule(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# ----------------------------------------------------------------------------------------------
# Part 1: exact exponent arithmetic
# ----------------------------------------------------------------------------------------------


def part1_exponents() -> None:
    rule("PART 1  exponent arithmetic for E1 (Lions) and E5 (dyadic budget), exact in sympy")
    s, sigma, lam = sp.symbols("s sigma lambda", positive=True)

    # Scaling exponent of ||u||_{Hdot^sigma}^2 under u -> lambda^{2s-1} u(lambda x), in 3D.
    d = 2 * (2 * s - 1) + 2 * sigma - 3
    print("d(sigma) = deg of ||u||_{Hdot^sigma}^2 under the (-Delta)^s scaling :", sp.expand(d))
    lhs = d.subs(sigma, s) + d.subs(sigma, 1)  # ||Lambda^s u||^2 * ||Lambda u||^2
    rhs = d.subs(sigma, 1 + s)  # ||Lambda^{1+s} u||^2  (= the enstrophy dissipation)
    print("  deg[ ||Lambda^s u||^2 * ||grad u||^2 ] =", sp.expand(lhs))
    print("  deg[ ||Lambda^{1+s} u||^2 ]           =", sp.expand(rhs))
    print("  equality holds iff s =", sp.solve(sp.Eq(lhs, rhs), s))

    # The product-estimate / interpolation chain of E1.
    a = sp.Integer(1)
    b = sp.Rational(5, 2) - s
    theta = sp.simplify(b / s)  # Hdot^b between L^2 and Hdot^s
    print("\nE1 chain: ||u tensor u||_{Hdot^{2-s}} <= C ||u||_{Hdot^a} ||u||_{Hdot^b},  a =", a, ", b =", b)
    print(
        "  product-estimate index identity  a + b - 3/2 =", sp.simplify(a + b - sp.Rational(3, 2)), " (must equal 2-s)"
    )
    print("  need a < 3/2 : always;   need b < 3/2  <=>  s > 1")
    print("  interpolation exponent theta = b/s =", theta)
    print("  theta <= 1  <=>  s >=", sp.solve(sp.Eq(theta, 1), s))
    print("  theta >= 0  <=>  s <=", sp.solve(sp.Eq(b, 0), s), " (so E1 is stated for 5/4 <= s <= 5/2)")
    for sv in [sp.Rational(5, 4), sp.Rational(13, 10), sp.Rational(3, 2), sp.Rational(2), sp.Rational(5, 2)]:
        print(
            f"    s = {sv}:  b = {b.subs(s, sv)},  theta = {theta.subs(s, sv)},  b<3/2: {b.subs(s, sv) < sp.Rational(3, 2)}"
        )

    # E5 bookkeeping: a_N, b_N, c_N in terms of N, g, h, ||u_N||.
    n, gN, hN, uN = sp.symbols("N g h u", positive=True)
    aN = n ** sp.Rational(5, 4) / gN * uN  # energy-dissipation coefficient
    cN = hN * n * uN  # weighted-enstrophy coefficient
    diss = hN**2 * n ** sp.Rational(9, 2) / gN**2 * uN**2  # weighted enstrophy dissipation at shell N
    nonl = hN**2 * n ** sp.Rational(9, 2) * uN**3  # diagonal trilinear bound at shell N
    print("\nE5: nonlinear/dissipation ratio at frequency N =", sp.simplify(nonl / diss))
    print("  -> absorption at shell N iff  ||u_N||_2 <= nu / g(N)^2   (no h: the weight cancels)")
    print(
        "  low-frequency form  nonl =",
        sp.simplify(sp.factor(nonl.subs(uN, sp.solve(sp.Eq(aN, sp.Symbol("a")), uN)[0]))),
    )
    print("  i.e. nonl = N * g(N)^2 * h(N) * a_N^2 * c_N   with c_N = h N ||u_N||; check:")
    expr = n * gN**2 * hN * sp.Symbol("a") ** 2 * sp.Symbol("c")
    subs = {sp.Symbol("a"): aN, sp.Symbol("c"): cN}
    print("   identity residual =", sp.simplify(expr.subs(subs) - nonl))


# ----------------------------------------------------------------------------------------------
# Part 2: the weight-cancellation identity of E6
# ----------------------------------------------------------------------------------------------


def _g_family() -> dict:
    return {
        "g=1": lambda x: np.ones_like(np.asarray(x, dtype=float)),
        "g=log^{1/4}": lambda x: np.log(2.0 + np.asarray(x, dtype=float)) ** 0.25,
        "g=log^{1/2}": lambda x: np.log(2.0 + np.asarray(x, dtype=float)) ** 0.5,
        "g=log^{1}": lambda x: np.log(2.0 + np.asarray(x, dtype=float)),
        "g=x^{0.05}": lambda x: np.asarray(x, dtype=float) ** 0.05,
        "g=x^{0.10}": lambda x: np.asarray(x, dtype=float) ** 0.10,
        "g=x^{0.25}": lambda x: np.asarray(x, dtype=float) ** 0.25,
    }


def _log_g(name: str, log_n: np.ndarray) -> np.ndarray:
    """log g(N) as a function of log N, evaluated without ever forming N (which reaches 2^4000)."""
    small = log_n < 100.0
    log_log = np.where(small, np.log(np.log(2.0 + np.exp(np.minimum(log_n, 700.0)))), np.log(np.maximum(log_n, 1e-9)))
    if name == "g=1":
        return np.zeros_like(log_n)
    if name == "g=log^{1/4}":
        return 0.25 * log_log
    if name == "g=log^{1/2}":
        return 0.5 * log_log
    if name == "g=log^{1}":
        return log_log
    if name == "g=x^{0.05}":
        return 0.05 * log_n
    if name == "g=x^{0.10}":
        return 0.10 * log_n
    if name == "g=x^{0.25}":
        return 0.25 * log_n
    raise ValueError(name)


def part2_weight_cancellation() -> None:
    rule("PART 2  E6: the absorption threshold N_0 and the Gronwall coefficient, for several weights h")
    print("Scheme: Y_h = sum_N h(N)^2 N^2 ||u_N||^2 ;  N_0 = smallest dyadic N with N h(N)/g(N)^2 >= Y^{1/2}/nu.")
    print("E6 says the Gronwall coefficient is nu^{-1} g(N_0)^4 Y_h A^2 for EVERY weight h, so the only way a")
    print("weight could help is by moving N_0. The table shows g(N_0)^4 for h(N) = N^a; 'spread' is the ratio")
    print("of the largest to the smallest entry in the row. A spread bounded uniformly in Y means no weight")
    print("changes the Osgood condition, since int ds/(s g(s)^4) = infinity is invariant under g -> c g.\n")
    nu = 1.0
    gs = _g_family()
    aexps = {"h=N^{-1/2}": -0.5, "h=1": 0.0, "h=N^{1/2}": 0.5, "h=N^{1}": 1.0}
    lg2 = np.arange(0.0, 4000.0)  # N = 2^lg2, i.e. frequencies up to 2^4000, handled in log space
    logN = lg2 * np.log(2.0)
    for gname in gs:
        log_g = _log_g(gname, logN)
        print(f"  {gname}")
        print("    " + f"{'log10 Y':>9}" + "".join(f"{w:>15}" for w in aexps) + f"{'spread':>12}")
        for log10y in [10.0, 40.0, 160.0, 640.0]:
            target = 0.5 * log10y * np.log(10.0) - np.log(nu)  # log(Y^{1/2}/nu)
            vals = []
            for aexp in aexps.values():
                thr = (1.0 + aexp) * logN - 2.0 * log_g
                hit = np.nonzero(np.maximum.accumulate(thr) >= target)[0]
                vals.append(np.exp(4.0 * log_g[hit[0]]) if len(hit) else np.nan)
            spread = max(vals) / min(vals)
            print("    " + f"{log10y:>9.0f}" + "".join(f"{v:>15.4g}" for v in vals) + f"{spread:>12.4g}")
        print()
    print("  For every log-type g the spread settles to a constant as Y grows (it converges to")
    print("  ((1+a_max)/(1+a_min))^{4 beta} for g = log^beta), so no weight in the family improves the")
    print("  divergence condition. For the power loss g = x^{0.05} the spread grows without bound, but that")
    print("  g fails both Tao's and BMR's hypotheses anyway (Part 4), so it is outside the frontier.")


# ----------------------------------------------------------------------------------------------
# Part 3: the Osgood ODE
# ----------------------------------------------------------------------------------------------


def part3_osgood() -> None:
    rule("PART 3  E3: y' = phi(t) y g(y)^p with phi in L^1; blowup iff int^inf dy/(y g(y)^p) converges")
    print("Test case g(y) = log(2+y), phi(t) = 1 on [0,T].  int^inf dy/(y log(2+y)^p) diverges iff p <= 1,")
    print("so the ODE is global for p <= 1 and blows up in finite time for p > 1.")
    print("The last two columns are the truncated Osgood integral int_{log 2}^{C} du / log(2+e^u)^p at two")
    print("cutoffs C; it keeps growing iff the condition diverges.\n")
    print(f"  {'p':>6}{'y(5)':>16}{'y(20)':>16}{'Osgood C=60':>14}{'Osgood C=600':>14}{'predicted':>12}")
    for p in [0.5, 1.0, 1.02, 1.2, 2.0, 4.0]:

        def f(_t, y, p=p):
            return y * np.log(2.0 + y) ** p

        out = []
        for tf in (5.0, 20.0):
            try:
                sol = solve_ivp(f, (0.0, tf), [1.0], method="Radau", rtol=1e-10, atol=1e-12)
                ok = bool(sol.success) and sol.t[-1] >= tf - 1e-9 and np.isfinite(sol.y[0, -1])
                out.append(sol.y[0, -1] if ok else np.inf)
            except (OverflowError, FloatingPointError, ValueError):
                out.append(np.inf)
        integ = []
        for cutoff in (60.0, 600.0):
            u = np.linspace(np.log(2.0), cutoff, 400001)
            integ.append(float(np.trapezoid(1.0 / np.log(2.0 + np.exp(np.minimum(u, 700.0))) ** p, u)))
        pred = "global" if p <= 1.0 else "blowup"
        sa, sb = (f"{v:.4g}" if np.isfinite(v) else "overflow" for v in out)
        print(f"  {p:>6.2f}{sa:>16}{sb:>16}{integ[0]:>14.5g}{integ[1]:>14.5g}{pred:>12}")
    print("\n  'overflow' = the stiff solver could not reach the final time (y left double range): finite-time blowup.")


# ----------------------------------------------------------------------------------------------
# Part 4: dyadic shell model realizing the cascade budget of E4
# ----------------------------------------------------------------------------------------------

BETA = 2.5  # lambda^{beta n} coupling: the 3D Navier-Stokes dyadic normalization (Lions exponent 2*alpha = beta)
LAM = 2.0


def shell_rhs(_t: float, w: np.ndarray, coup: np.ndarray, diss: np.ndarray) -> np.ndarray:
    """State w = (a, z): a = shell amplitudes, z_n = cumulative energy dissipated by shell n."""
    m = coup.size
    a = w[:m]
    am1 = np.concatenate(([0.0], a[:-1]))
    ap1 = np.concatenate((a[1:], [0.0]))
    cm1 = np.concatenate(([0.0], coup[:-1]))
    return np.concatenate((cm1 * am1**2 - coup * a * ap1 - diss * a, 2.0 * diss * a**2))


ALPHA2_DYADIC = 2.0 * BETA / 3.0  # the dyadic model's own critical dissipation exponent (Kolmogorov balance)


def part4_shell_model(nu: float = 0.05, t_final: float = 60.0) -> None:
    rule("PART 4  E4/E9: the cascade budget, and a dyadic shell model at its own critical exponent")
    print("Table A: the exact budgets. Tao 2009 assumes sum_m g(2^m)^{-4} = infinity; BMR 2014 and the E4")
    print("cascade budget assume sum_m g(2^m)^{-2} = infinity. Partial sums S2, S4 over m < M:\n")
    names = list(_g_family())
    ks = (1, 2, 3, 4)
    print(f"  {'g':>14}" + "".join(f"{f'S2(M=1e{k})':>13}" for k in ks) + "".join(f"{f'S4(M=1e{k})':>13}" for k in ks))
    for gname in names:
        row = []
        for power in (2.0, 4.0):
            for k in ks:
                mm = np.arange(0.0, 10.0**k) * np.log(LAM)
                row.append(float(np.sum(np.exp(-power * _log_g(gname, mm)))))
        print(f"  {gname:>14}" + "".join(f"{v:>13.4g}" for v in row))
    print("\n  Reading: S2 diverges for g = 1, log^{1/4}, log^{1/2} and converges for g = log^{1}; S4 diverges")
    print("  only for g = 1 and log^{1/4}. So g = log^{1/2} lies strictly between Tao's theorem and BMR's:")
    print("  it is the smallest natural example the 2014 improvement covers and the 2009 result does not.")
    print("  The columns also show why the frontier is numerically invisible: the divergent sums grow like")
    print("  log M, so separating log^{1/2} from log^{1} needs dyadic ranges of thousands of shells.\n")

    print(
        f"Table B: shell model a_n' = c_{{n-1}} a_{{n-1}}^2 - c_n a_n a_{{n+1}} - nu d_n a_n, c_n = {LAM}^{{{BETA}(n+1)}},"
    )
    print(f"  d_n = {LAM}^{{{ALPHA2_DYADIC:.4f} n}}/g({LAM}^n)^2, nu = {nu}, a_0(0) = 1, T = {t_final}.")
    print("  The exponent 2*beta/3 is where the dyadic model's OWN constant-flux (Kolmogorov) cascade")
    print("  balances its dissipation; it is the dyadic analogue of the Lions exponent, and it is NOT 2*alpha")
    print("  = 5/2 (the dyadic model is a different system with a different criticality). Reported: log10 of")
    print("  max_t of the enstrophy sum_n lambda^{2n} a_n^2 as the truncation M is refined, and the growth")
    print("  rate in decades per shell. A constant-flux inertial range gives exactly log10(lambda^{1/3}) =")
    print(
        f"  {np.log10(LAM ** (1.0 / 3.0)):.4f} decades per shell with no blowup; a larger rate signals an undamped cascade.\n"
    )
    ms = [10, 14, 18, 22]
    print(f"  {'g':>14}" + "".join(f"{f'M={m}':>10}" for m in ms) + f"{'dec/shell':>11}{'vs inertial':>13}")
    for gname, g in _g_family().items():
        vals = []
        for shells in ms:
            n = np.arange(shells, dtype=float)
            coup = LAM ** (BETA * (n + 1.0))
            diss = nu * LAM ** (ALPHA2_DYADIC * n) / g(LAM**n) ** 2
            w0 = np.zeros(2 * shells)
            w0[0] = 1.0
            sol = solve_ivp(shell_rhs, (0.0, t_final), w0, args=(coup, diss), method="Radau", rtol=1e-9, atol=1e-16)
            ens = np.max(np.sum((LAM ** (2.0 * n))[:, None] * sol.y[:shells] ** 2, axis=0))
            vals.append(np.log10(max(float(ens), 1e-300)))
        rate = (vals[-1] - vals[-2]) / (ms[-1] - ms[-2])
        print(
            f"  {gname:>14}"
            + "".join(f"{v:>10.3f}" for v in vals)
            + f"{rate:>11.4f}{rate / np.log10(LAM ** (1 / 3)):>13.2f}"
        )
    print("\n  Reading, and the honest negative part of it: the critical case g = 1 sits at the inertial rate")
    print("  (ratio ~ 0.9, slightly damped), and the power losses g = x^{0.10}, x^{0.25} sit well above it,")
    print("  so the model does separate a power loss from criticality at 22 shells. It does NOT separate the")
    print("  three logarithmic cases from each other or from a power loss: their ratios interpolate smoothly.")
    print("  That is Table A's resolution statement again. No numerical experiment at any feasible truncation")
    print("  can see the Tao/BMR frontier, so the numerics here support the bookkeeping of E4 and E5 and")
    print("  cannot be evidence for or against the divergence conditions themselves.")
    print("  This is numerical evidence about a shell model, NOT a proof and NOT about Navier-Stokes.")


def main() -> None:
    print("Line E — logarithmically supercritical hyperdissipation. Seed =", SEED, "(no randomness is used).")
    _ = rng.random()
    part1_exponents()
    part2_weight_cancellation()
    part3_osgood()
    part4_shell_model()


if __name__ == "__main__":
    main()
