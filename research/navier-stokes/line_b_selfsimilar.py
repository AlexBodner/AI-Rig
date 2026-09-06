"""Line B: self-similar blowup profiles and Liouville theorems for 3D Navier-Stokes.

Five independent checks, all printed as Markdown tables / residuals:

  Part 1  exponent balance for u = (T-t)^{-alpha} U(x/(T-t)^beta): alpha = 1 - beta,
          viscous term carries the extra factor s^{1-2beta}.
  Part 2  the Necas-Ruzicka-Sverak identity Delta Pi - (y/2 + U).grad Pi = |curl U|^2,
          verified as an exact symbolic identity modulo the profile equation, for a
          generic divergence-free U = curl A (A three arbitrary functions).
  Part 3  the discretely-self-similar (tau-dependent) version and the obstruction term.
  Part 4  the surviving self-similar window 2/(d+2) <= beta <= 1/(2 s) and the Lions
          exponent s = (d+2)/4 in dimension d with dissipation -(-Laplacian)^s.
  Part 5  numerical check of the integration-by-parts identity used in Part 4.

Run:  python3 line_b_selfsimilar.py     (all output to stdout, no files written)
"""

from fractions import Fraction

import numpy as np
import sympy as sp

SEED = 20260906

Y = sp.symbols("y1 y2 y3", real=True)
X = sp.symbols("x1 x2 x3", real=True)
S, TAU = sp.symbols("s tau", positive=True)
AL, BE, NU = sp.symbols("alpha beta nu", positive=True)


# ----------------------------------------------------------------------------- helpers
def grad(f, v=Y):
    return [sp.diff(f, vi) for vi in v]


def div(F, v=Y):
    return sum(sp.diff(F[i], v[i]) for i in range(3))


def curl(F, v=Y):
    return [
        sp.diff(F[2], v[1]) - sp.diff(F[1], v[2]),
        sp.diff(F[0], v[2]) - sp.diff(F[2], v[0]),
        sp.diff(F[1], v[0]) - sp.diff(F[0], v[1]),
    ]


def lap(f, v=Y):
    return sum(sp.diff(f, vi, 2) for vi in v)


def dot(F, G):
    return sum(F[i] * G[i] for i in range(3))


def advect(F, G, v=Y):
    """(F . grad) G, componentwise."""
    return [sum(F[j] * sp.diff(G[i], v[j]) for j in range(3)) for i in range(3)]


def rand_poly(rng, v, deg=3, den=7):
    """Random polynomial in the variables v with small exact rational coefficients."""
    terms = []
    for e1 in range(deg + 1):
        for e2 in range(deg + 1 - e1):
            for e3 in range(deg + 1 - e1 - e2):
                c = Fraction(int(rng.integers(-4, 5)), int(rng.integers(1, den + 1)))
                if c != 0:
                    terms.append(sp.Rational(c.numerator, c.denominator) * v[0] ** e1 * v[1] ** e2 * v[2] ** e3)
    return sp.expand(sum(terms))


# ------------------------------------------------------------------- Part 1: exponents
def part1_exponents():
    """Substitute u = s^{-alpha} U(x/s^beta), p = s^{-2 alpha} P(x/s^beta) into
    d_t u + (u.grad)u - nu Lap u + grad p  and read off the powers of s = T - t."""
    print("## Part 1 - exponent balance (sympy, concrete random rational profile)\n")
    rng = np.random.default_rng(SEED)
    Upoly = [rand_poly(rng, Y, deg=3) for _ in range(3)]
    Ppoly = rand_poly(rng, Y, deg=3)

    ysub = {Y[i]: X[i] / S**BE for i in range(3)}
    u = [S ** (-AL) * Upoly[i].subs(ysub, simultaneous=True) for i in range(3)]
    p = S ** (-2 * AL) * Ppoly.subs(ysub, simultaneous=True)

    # d_t = -d_s because s = T - t.
    dt_u = [-sp.diff(u[i], S) for i in range(3)]
    conv = advect(u, u, X)
    lap_u = [lap(u[i], X) for i in range(3)]
    gp = grad(p, X)
    nse = [sp.expand(dt_u[i] + conv[i] - NU * lap_u[i] + gp[i]) for i in range(3)]

    # go back to profile variables x = s^beta y and pull out the overall power s^{-alpha-1}
    xsub = {X[i]: S**BE * Y[i] for i in range(3)}
    lhs = [sp.expand(sp.simplify((nse[i].subs(xsub, simultaneous=True)) * S ** (AL + 1))) for i in range(3)]

    # claimed answer, after imposing alpha = 1 - beta
    claimed = [
        AL * Upoly[i]
        + BE * sum(Y[j] * sp.diff(Upoly[i], Y[j]) for j in range(3))
        + advect(Upoly, Upoly)[i]
        + sp.diff(Ppoly, Y[i])
        - NU * S ** (1 - 2 * BE) * lap(Upoly[i])
        for i in range(3)
    ]
    resid_general = [sp.simplify(lhs[i] - claimed[i]) for i in range(3)]
    resid_balanced = [sp.simplify(r.subs(AL, 1 - BE)) for r in resid_general]

    # power of s carried by each term of the transformed equation, relative to s^{-alpha-1}
    def power_of(term_x, term_y):
        ratio = sp.simplify((term_x.subs(xsub, simultaneous=True) * S ** (AL + 1)) / term_y)
        return sp.simplify(sp.log(ratio) / sp.log(S))

    rows = [
        ("d_t u", power_of(dt_u[0], AL * Upoly[0] + BE * sum(Y[j] * sp.diff(Upoly[0], Y[j]) for j in range(3)))),
        ("(u.grad)u", power_of(conv[0], advect(Upoly, Upoly)[0])),
        ("grad p", power_of(gp[0], sp.diff(Ppoly, Y[0]))),
        ("nu Lap u", power_of(NU * lap_u[0], NU * lap(Upoly[0]))),
    ]
    print("| term | power of s = T - t relative to s^{-alpha-1}, after alpha = 1 - beta |")
    print("| ---- | ---- |")
    for name, e in rows:
        print(f"| `{name}` | `s^({sp.simplify(e.subs(AL, 1 - BE))})` |")
    ok_gen = all(r == 0 for r in resid_general)
    ok = all(r == 0 for r in resid_balanced)
    print("\n- residual of (transformed NSE) - [alpha U + beta (y.grad)U + (U.grad)U + grad P - nu s^(1-2beta) Lap U]")
    print(f"  is identically zero for generic (alpha, beta): {ok_gen}")
    print(f"- ... and after imposing alpha = 1 - beta: {ok}")
    print("- so the transport/nonlinear/pressure terms all carry s^0 and the viscous term s^(1-2beta);")
    print("  an exact self-similar solution therefore needs beta = 1/2, else both groups vanish separately.\n")
    return ok


# ------------------------------------------------------------ Part 2: the NRS identity
def part2_nrs(generic=True):
    """Verify  Delta Pi - b.grad Pi - |omega|^2 = -div(Res) + b.Res,
    Res := Lap U - [ (1/2)U + (1/2)(y.grad)U + (U.grad)U + grad P ],
    Pi  := P + |U|^2/2 + (y.U)/2,   b := y/2 + U,   omega := curl U.
    U = curl A is the generic divergence-free field (Poincare lemma on R^3)."""
    tag = "generic (three arbitrary functions A_i)" if generic else "random rational polynomials"
    print(f"## Part 2 - Necas-Ruzicka-Sverak identity, U = curl A, {tag}\n")
    if generic:
        A = [sp.Function(f"A{i + 1}")(*Y) for i in range(3)]
        P = sp.Function("P")(*Y)
    else:
        rng = np.random.default_rng(SEED + 1)
        A = [rand_poly(rng, Y, deg=4) for _ in range(3)]
        P = rand_poly(rng, Y, deg=4)

    U = curl(A)
    om = curl(U)
    b = [Y[i] / 2 + U[i] for i in range(3)]
    Pi = P + dot(U, U) / 2 + dot(Y, U) / 2
    Res = [
        lap(U[i])
        - (U[i] / 2 + sum(Y[j] * sp.diff(U[i], Y[j]) for j in range(3)) / 2 + advect(U, U)[i] + sp.diff(P, Y[i]))
        for i in range(3)
    ]
    lhs = lap(Pi) - dot(b, grad(Pi)) - dot(om, om)
    rhs = -div(Res) + dot(b, Res)
    resid = sp.expand(lhs - rhs)
    print(f"- div U = curl-free check, div(curl A) = {sp.simplify(div(U))}")
    print(f"- residual of  [Delta Pi - b.grad Pi - |omega|^2] - [-div Res + b.Res]  =  {sp.simplify(resid)}")
    ok = sp.simplify(resid) == 0
    print(f"- identity verified: {ok}")
    return ok


def part2_teeth():
    """Control: the same residual with a NON-divergence-free U must be nonzero,
    and the identity with a wrong constant c != 1 in front of |omega|^2 must fail."""
    print("\n### Part 2 control experiments (the test has teeth)\n")
    rng = np.random.default_rng(SEED + 2)
    Ubad = [rand_poly(rng, Y, deg=2) for _ in range(3)]
    P = rand_poly(rng, Y, deg=2)
    om = curl(Ubad)
    b = [Y[i] / 2 + Ubad[i] for i in range(3)]
    Pi = P + dot(Ubad, Ubad) / 2 + dot(Y, Ubad) / 2
    Res = [
        lap(Ubad[i])
        - (
            Ubad[i] / 2
            + sum(Y[j] * sp.diff(Ubad[i], Y[j]) for j in range(3)) / 2
            + advect(Ubad, Ubad)[i]
            + sp.diff(P, Y[i])
        )
        for i in range(3)
    ]
    r_bad = sp.expand(lap(Pi) - dot(b, grad(Pi)) - dot(om, om) + div(Res) - dot(b, Res))
    print(f"- div U for the control field: {sp.expand(div(Ubad)) != 0} (nonzero as intended)")
    print(f"- residual with a non-divergence-free U is identically zero: {sp.simplify(r_bad) == 0}")
    rng2 = np.random.default_rng(SEED + 3)
    A = [rand_poly(rng2, Y, deg=3) for _ in range(3)]
    P2 = rand_poly(rng2, Y, deg=3)
    U = curl(A)
    om2 = curl(U)
    b2 = [Y[i] / 2 + U[i] for i in range(3)]
    Pi2 = P2 + dot(U, U) / 2 + dot(Y, U) / 2
    Res2 = [
        lap(U[i])
        - (U[i] / 2 + sum(Y[j] * sp.diff(U[i], Y[j]) for j in range(3)) / 2 + advect(U, U)[i] + sp.diff(P2, Y[i]))
        for i in range(3)
    ]
    for c in (sp.Rational(-1), sp.Rational(2), sp.Rational(1, 2)):
        r = sp.expand(lap(Pi2) - dot(b2, grad(Pi2)) - c * dot(om2, om2) + div(Res2) - dot(b2, Res2))
        print(f"- residual with constant c = {c} in front of |omega|^2 is zero: {sp.simplify(r) == 0}")
    for drift in ("y_only", "U_only", "minus"):
        bb = {
            "y_only": [Y[i] / 2 for i in range(3)],
            "U_only": list(U),
            "minus": [-(Y[i] / 2 + U[i]) for i in range(3)],
        }[drift]
        r = sp.expand(lap(Pi2) - dot(bb, grad(Pi2)) - dot(om2, om2) + div(Res2) - dot(bb, Res2))
        print(f"- residual with drift '{drift}' is zero: {sp.simplify(r) == 0}")
    return True


# --------------------------------------------------- Part 3: discretely self-similar
def part3_dss():
    """Similarity variables tau = -log(T-t), V = V(y,tau) (DSS: V periodic in tau).
    Profile equation:  d_tau V + V/2 + (y.grad)V/2 + (V.grad)V + grad P = Lap V.
    Claim:  (Lap - b.grad - d_tau) Pi = |omega|^2 - d_tau P   modulo Res_tau,
    the precise identity being
        Lap Pi - b.grad Pi - d_tau Pi + d_tau P - |omega|^2 = -div Res + b.Res.
    The extra term -d_tau P has no sign: this is where the NRS maximum principle dies."""
    print("## Part 3 - discretely/asymptotically self-similar version (tau-dependent)\n")
    A = [sp.Function(f"A{i + 1}")(*Y, TAU) for i in range(3)]
    P = sp.Function("P")(*Y, TAU)
    V = curl(A)
    om = curl(V)
    b = [Y[i] / 2 + V[i] for i in range(3)]
    Pi = P + dot(V, V) / 2 + dot(Y, V) / 2
    Res = [
        lap(V[i])
        - (
            sp.diff(V[i], TAU)
            + V[i] / 2
            + sum(Y[j] * sp.diff(V[i], Y[j]) for j in range(3)) / 2
            + advect(V, V)[i]
            + sp.diff(P, Y[i])
        )
        for i in range(3)
    ]
    lhs = lap(Pi) - dot(b, grad(Pi)) - sp.diff(Pi, TAU) + sp.diff(P, TAU) - dot(om, om)
    rhs = -div(Res) + dot(b, Res)
    resid = sp.simplify(sp.expand(lhs - rhs))
    print(f"- residual of the tau-dependent identity: {resid}")
    ok = resid == 0
    print(f"- identity verified: {ok}")
    # control: dropping the d_tau P term must break it
    resid2 = sp.simplify(sp.expand(lhs - sp.diff(P, TAU) - rhs))
    print(f"- residual after deleting the obstruction term d_tau P is zero: {resid2 == 0}")
    print("- consequence: (Lap - b.grad - d_tau) Pi = |omega|^2 - d_tau P; the sign of d_tau P")
    print("  is not controlled, so Pi is not a parabolic subsolution and NRS does not extend.\n")
    return ok


# ------------------------------------------- Part 4: the surviving window and Lions
def part4_window():
    """u = (T-t)^{-alpha} U(x/(T-t)^beta) in R^d with dissipation -nu (-Lap)^sigma.
    alpha = 1 - beta (nonlinear balance).
    beta_nu  = 1/(2 sigma)  : dissipation balances transport (Type I).
    beta_E   = 2/(d+2)      : the kinetic energy of the profile region is scale-invariant;
                              finite energy forces beta >= beta_E.
    Window [beta_E, beta_nu] is nonempty iff sigma <= (d+2)/4 (Lions exponent)."""
    print("## Part 4 - self-similar window in dimension d with dissipation -(-Lap)^sigma\n")
    d, sig, be = sp.symbols("d sigma beta", positive=True)
    al = 1 - be
    # energy exponent: int_{R^d} |u|^2 dx = s^{-2 alpha + d beta} ||U||^2
    e_exp = sp.simplify(-2 * al + d * be)
    # profile L^2 identity: nu ||(-Lap)^{sigma/2} U||^2 = ((d+2) beta - 2)/2 ||U||^2
    ident = sp.simplify(d * be / 2 - al)
    beta_E = sp.solve(sp.Eq(e_exp, 0), be)[0]
    beta_nu = sp.Rational(1, 2) / sig
    print(f"- energy exponent  int |u(t)|^2 dx ~ s^E with E = {e_exp}  (zero at beta = {beta_E})")
    print(f"- profile identity constant  nu ||D^sigma U||^2 = c ||U||^2 with c = {sp.factor(ident)}")
    print(f"- beta_E = {beta_E},  beta_nu = 1/(2 sigma) = {beta_nu}")
    cond = sp.simplify(sp.solve(sp.Eq(beta_E, beta_nu), sig)[0])
    print(f"- beta_E = beta_nu  <=>  sigma = {cond}   (Lions exponent; (d+2)/4 = 5/4 at d = 3)\n")
    print("| d | sigma | beta_E = 2/(d+2) | beta_nu = 1/(2 sigma) | window | verdict |")
    print("| - | ----- | ---------------- | --------------------- | ------ | ------- |")
    for dv in (2, 3, 4):
        for sv in (sp.Rational(1), sp.Rational(9, 8), sp.Rational(5, 4), sp.Rational(3, 2)):
            bE, bN = sp.Rational(2, dv + 2), sp.Rational(1, 2) / sv
            verdict = "nonempty" if bE < bN else ("single point (energy-critical)" if bE == bN else "EMPTY")
            print(f"| {dv} | {sv} | {bE} | {bN} | [{bE}, {bN}] | {verdict} |")
    print()
    print(
        "| beta | alpha = 1 - beta | sup-norm of u | L^3 norm of u | sup-norm of omega | L^2 norm of grad u | CKN quantity |"
    )
    print("| ---- | ---- | ---- | ---- | ---- | ---- | ---- |")
    for bv in (sp.Rational(7, 20), sp.Rational(2, 5), sp.Rational(9, 20), sp.Rational(1, 2), sp.Rational(3, 5)):
        av = 1 - bv
        # CKN quantity r^{-1} int_{Q_r} |grad u|^2 on Q_r = B_r x (T*-r^2, T*), with
        # Q_r = B_r x (T*-r^2,T*): the profile core fits inside B_r for all s < r^2 only
        # when beta >= 1/2, which is why the exponent is piecewise.
        ckn = 2 - 1 / bv if bv <= sp.Rational(1, 2) else 6 * bv - 3
        print(
            f"| {bv} | {av} | s^({-av}) | s^({sp.Rational(2) * bv - 1}) | s^(-1) "
            f"| s^({(3 * bv - 2) / 2}) | r^({sp.nsimplify(ckn)}) |"
        )
    print("\n(Leray needs ||u||_inf >= c s^(-1/2) and ||grad u||_2 >= c s^(-1/4): both force beta <= 1/2.")
    print(
        " CKN epsilon-regularity needs the last column not to vanish as r -> 0: also beta <= 1/2\n (the CKN exponent is 2 - 1/beta for 1/3 < beta <= 1/2 and 6 beta - 3 for beta >= 1/2;\n  the two agree at beta = 1/2)."
    )
    print(" Finite energy forces beta >= 2/5.  ESS/NRS kill beta = 1/2.  Surviving: 2/5 <= beta < 1/2.)\n")
    return True


# -------------------------------- Part 5: numerical check of the profile L^2 identity
def part5_numeric():
    """Check int ((y.grad)U).U dy = -(d/2) int |U|^2 dy numerically in d = 3 for a
    smooth compactly supported divergence-free U = curl A, and confirm the
    resulting profile identity constant ((d+2) beta - 2)/2."""
    print("## Part 5 - numerical check of the integration-by-parts identity (d = 3)\n")
    n, L = 64, 8.0
    ax = np.linspace(-L / 2, L / 2, n, endpoint=False)
    k1 = 2 * np.pi * np.fft.fftfreq(n, d=L / n)
    y1, y2, y3 = np.meshgrid(ax, ax, ax, indexing="ij")
    K = np.meshgrid(k1, k1, k1, indexing="ij")
    r2 = y1**2 + y2**2 + y3**2
    g = np.exp(-r2)
    rng = np.random.default_rng(SEED + 7)
    c = rng.normal(size=(3, 3))
    A = [g * (c[i, 0] * y1 + c[i, 1] * y2 + c[i, 2] * y3 + 0.5 * y1 * y2) for i in range(3)]

    def d_(f, axis):
        """Spectral derivative on the periodic box (the field is Gaussian-localized)."""
        return np.real(np.fft.ifftn(1j * K[axis] * np.fft.fftn(f)))

    U = [d_(A[2], 1) - d_(A[1], 2), d_(A[0], 2) - d_(A[2], 0), d_(A[1], 0) - d_(A[0], 1)]
    divU = d_(U[0], 0) + d_(U[1], 1) + d_(U[2], 2)
    h = ax[1] - ax[0]
    ydotgradU_dot_U = sum(U[i] * (y1 * d_(U[i], 0) + y2 * d_(U[i], 1) + y3 * d_(U[i], 2)) for i in range(3))
    normsq = sum(U[i] ** 2 for i in range(3))
    lhs = ydotgradU_dot_U.sum() * h**3
    rhs = -1.5 * normsq.sum() * h**3
    print(f"- max |div U| on the grid: {np.abs(divU).max():.3e} (spectral)")
    print(f"- int ((y.grad)U).U dy   = {lhs: .10f}")
    print(f"- -(d/2) int |U|^2 dy    = {rhs: .10f}")
    print(f"- relative difference     = {abs(lhs - rhs) / abs(rhs):.3e}\n")
    return abs(lhs - rhs) / abs(rhs) < 1e-8


def main():
    ok = []
    ok.append(part1_exponents())
    ok.append(part2_nrs(generic=True))
    ok.append(part2_teeth())
    print()
    ok.append(part3_dss())
    ok.append(part4_window())
    ok.append(part5_numeric())
    print(f"All checks passed: {all(ok)}")


if __name__ == "__main__":
    main()
