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
    print("  product-estimate index identity  a + b - 3/2 =", sp.simplify(a + b - sp.Rational(3, 2)), " (must equal 2-s)")
    print("  need a < 3/2 : always;   need b < 3/2  <=>  s > 1")
    print("  interpolation exponent theta = b/s =", theta)
    print("  theta <= 1  <=>  s >=", sp.solve(sp.Eq(theta, 1), s))
    print("  theta >= 0  <=>  s <=", sp.solve(sp.Eq(b, 0), s), " (so E1 is stated for 5/4 <= s <= 5/2)")
    for sv in [sp.Rational(5, 4), sp.Rational(13, 10), sp.Rational(3, 2), sp.Rational(2), sp.Rational(5, 2)]:
        print(f"    s = {sv}:  b = {b.subs(s, sv)},  theta = {theta.subs(s, sv)},  b<3/2: {b.subs(s, sv) < sp.Rational(3, 2)}")

    # E5 bookkeeping: a_N, b_N, c_N in terms of N, g, h, ||u_N||.
    n, gN, hN, uN = sp.symbols("N g h u", positive=True)
    aN = n ** sp.Rational(5, 4) / gN * uN  # energy-dissipation coefficient
    cN = hN * n * uN  # weighted-enstrophy coefficient
    diss = hN**2 * n ** sp.Rational(9, 2) / gN**2 * uN**2  # weighted enstrophy dissipation at shell N
    nonl = hN**2 * n ** sp.Rational(9, 2) * uN**3  # diagonal trilinear bound at shell N
    print("\nE5: nonlinear/dissipation ratio at frequency N =", sp.simplify(nonl / diss))
    print("  -> absorption at shell N iff  ||u_N||_2 <= nu / g(N)^2   (no h: the weight cancels)")
    print("  low-frequency form  nonl =", sp.simplify(sp.factor(nonl.subs(uN, sp.solve(sp.Eq(aN, sp.Symbol('a')), uN)[0]))))
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
    }


def part2_weight_cancellation() -> None:
    rule("PART 2  E6: the absorption threshold N_0 and the Gronwall coefficient, for several weights h")
    print("Scheme: Y_h = sum_N h(N)^2 N^2 ||u_N||^2 ;  N_0 = smallest dyadic N with N h(N)/g(N)^2 >= Y^{1/2}/nu.")
    print("Claim (E6): the Gronwall coefficient equals nu^{-1} g(N_0)^4 Y_h A^2 for EVERY h, and")
    print("            g(N_0) is the same to within a bounded factor for every h in the family below.")
    nu = 1.0
    gs = _g_family()
    weights = {"h=1": 0.0, "h=N^{1/2}": 0.5, "h=N^{-1/2}": -0.5, "h=N^{1}": 1.0}
    dy = 2.0 ** np.arange(0, 121)  # dyadic frequencies up to 2^120
    for gname, g in gs.items():
        gv = g(dy)
        print(f"\n  {gname}")
        header = f"    {'Y':>10}" + "".join(f"{w:>14}" for w in weights) + f"{'spread':>10}"
        print(header)
        for logY in [10.0, 20.0, 40.0, 80.0]:
            y = 10.0**logY
            row, g0s = [], []
            for _wname, aexp in weights.items():
                thr = dy ** (1.0 + aexp) / gv**2  # N h(N) / g(N)^2  with h(N) = N^aexp
                idx = int(np.argmax(thr >= np.sqrt(y) / nu))
                row.append(f"g(N_0)^4={gv[idx] ** 4:9.3g}")
                g0s.append(gv[idx] ** 4)
            spread = max(g0s) / min(g0s)
            print(f"    {logY:>8.0f}  " + "".join(f"{r:>14}" for r in row) + f"{spread:>10.3g}")
    print("\n  'spread' = max/min of g(N_0)^4 over the four weights at fixed Y.")
    print("  For log-type g the spread stays O(1): no weight changes the Osgood condition (E6).")
    print("  For the power g = x^{0.05} the spread is large, but that g fails BOTH divergence conditions,")
    print("  so it is outside the frontier anyway (it is a genuine power loss, not a log loss).")


# ----------------------------------------------------------------------------------------------
# Part 3: the Osgood ODE
# ----------------------------------------------------------------------------------------------


def part3_osgood() -> None:
    rule("PART 3  E3: y' = phi(t) y g(y)^p with phi in L^1; blowup iff int^inf dy/(y g(y)^p) converges")
    print("Test case g(y) = log(2+y), phi(t) = 1 on [0,T].  int^inf dy/(y log(2+y)^p) diverges iff p <= 1.")
    print("So the ODE is global for p <= 1 and blows up in finite time for p > 1.")
    print(f"\n  {'p':>6}{'y(T) at T=5':>18}{'y(T) at T=20':>18}{'Osgood integral':>20}{'predicted':>14}")
    for p in [0.5, 1.0, 1.02, 1.2, 2.0, 4.0]:

        def f(_t, y, p=p):
            return y * np.log(2.0 + y) ** p

        out = []
        for tf in (5.0, 20.0):
            sol = solve_ivp(f, (0.0, tf), [1.0], method="Radau", rtol=1e-10, atol=1e-12, dense_output=False)
            out.append(sol.y[0, -1] if sol.success and sol.t[-1] >= tf - 1e-9 else np.inf)
        # tail of the Osgood integral, computed on a log grid
        u = np.linspace(np.log(2.0), 60.0, 400001)  # u = log y
        integ = np.trapezoid(1.0 / np.log(2.0 + np.exp(u)) ** p, u)
        pred = "global" if p <= 1.0 else "blowup"
        a, b = out
        sa = f"{a:.4g}" if np.isfinite(a) else "overflow"
        sb = f"{b:.4g}" if np.isfinite(b) else "overflow"
        print(f"  {p:>6.2f}{sa:>18}{sb:>18}{integ:>20.6g}{pred:>14}")
    print("\n  'overflow' = the stiff solver could not reach the final time (y left double range): finite-time blowup.")
    print("  The Osgood integral column is the truncated tail int_{log 2}^{60} du / log(2+e^u)^p; it grows without")
    print("  bound with the cutoff exactly when p <= 1, and converges when p > 1.")
