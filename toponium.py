"""Everything the toponium study computes, in one file.  The stage scripts 0_ .. 5_ import it.

  [JKT]  M. Jezabek, J.H. Kuhn, T. Teubner, Z. Phys. C 56 (1992) 653  (reference/BF01474740.pdf)
  [Fuks] B. Fuks, K. Hagiwara, K. Ma, Y.-J. Zheng, arXiv:2411.18962v2  (reference/2411.18962v2.pdf)

=====================================================================================
CONTENTS  (six classes, in this order)
=====================================================================================
  Constants     m_t, Gamma_t, C_F, n_F and the free Green's function G0     JKT (6) = Fuks (11)
  --- ours --------------------------------------------------------------------------
  JKTPotential  the potential of JKT Sec. 2 with the cut of Sec. 3           JKT (13)-(15), (19)-(20)
                (JKTPotential.coulomb(alpha): fixed-coupling Coulomb instead)  Fuks (12), alpha fixed
  LSSolver      Lippmann-Schwinger equation solved as in JKT Sec. 3 -> G(E,p) JKT (16)-(18), (21)-(23)
  CoulombExact  exact analytic G for fixed-coupling Coulomb (ground truth)  Fuks below (12)
  --- Fuks's ------------------------------------------------------------------------
  FuksTables    their published tables (swData) and their CALCGREEN routine  Fuks p.3, Sec. 3
  --- both --------------------------------------------------------------------------
  Events        MadGraph events and the per-event weight |G/G0|^2          Fuks (14), (15), (21)

Every G here has the sign of JKT (6) / Fuks (11) (G -> G0 at zero coupling). Fuks's swData
files store -G; FuksTables flips it on loading.

Events.weights takes G from any source (FuksTables, CoulombExact, LSSolver): the events, cuts
and weight are the same for all of them, so differences in the final plots come from G alone.

Where we fill in what JKT leave open (writeup.tex, last section):
  - grid: JKT's subintervals split once more, breaks 0, 4, 16, 60, m_t GeV (0_solver_checks.py (f))
  - A(p) of (22) on a fine auxiliary Gauss-Legendre rule graded towards q = p
  - K between the grid nodes by cubic-spline interpolation
  - alpha_s given at mZ and run with (14)
"""

import gzip
import os

import numpy as np
from scipy import integrate
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq


# =====================================================================================
class Constants:
    """Numbers used everywhere, and the free Green's function.

    EQUATION MAP
      JKT (4), (6)  G0(p, calE) = 1/(E + i Gamma_t - p^2/m_t)  (= Fuks (11))      G0
    """

    M_T = 173.0       # GeV, Fuks
    GAMMA_T = 1.49    # GeV, Fuks
    C_F = 4.0 / 3.0
    N_F = 5

    @staticmethod
    def G0(E, p, m_t=M_T, Gamma_t=GAMMA_T):
        """JKT (6) with JKT (4):  G0(p, calE) = 1 / (E + i Gamma_t - p^2/m_t).  Identical to Fuks (11)."""

        return 1.0 / (np.asarray(E) + 1j * Gamma_t - np.asarray(p) ** 2 / m_t)


C_F, N_F, M_T, GAMMA_T = Constants.C_F, Constants.N_F, Constants.M_T, Constants.GAMMA_T
TWO_PI2 = 2.0 * np.pi ** 2
EIGHT_PI2 = 8.0 * np.pi ** 2


# =====================================================================================
class JKTPotential:
    """The QCD potential in momentum space, JKT Sec. 2, cut at small momentum as in JKT Sec. 3.

    EQUATION MAP  (JKT equation  ->  method), read alongside the paper
      (14)  alpha_MS(Q^2) = 4pi / [b0 ln(Q^2/L^2) + (b1/b0) ln ln(Q^2/L^2)]    alpha_2loop
            Lambda such that (14) gives alpha_s(mZ)                           lambda_msbar_2loop
      (13)  V_pert(Q^2) = -(16pi/3) alpha/Q^2 [1 + (31/3 - 10nF/9) alpha/4pi] V_pert  (= Fuks (12))
      (15)  V = V_pert above 5 GeV; -(4/3)(48pi^2/27)(1/p^2)[1/ln(1+p^2/L_R^2) + C]
            below, C from continuity                                       V_ir, continuity_constant_C, V_jkt
      (19)  V_cut = V above q_cut; -C1/(p^2+q_cut^2) + V0 d3(p) below       V_cut
            C1 from continuity at p = q_cut                                 cut_constant_C1
      (20)  V0 from V_cut(r = 1 GeV^-1) = -1/4 GeV; the d3 term is an
            energy shift E -> E - V0/(2pi)^3                                V_cut_position, energy_shift

    Why the cut (JKT below (18)): the kernel (18) is singular at q = p, like 1/(p-q)^2 for (15),
    whose 1/p^4 small-momentum behaviour is the linear confining potential. Below q_cut = 50 MeV
    JKT replace it by a finite form (19), and fix the d3 term of (19) in position space (20).

    What the solver reads from an instance:
      V(k2)      the potential in momentum space ((19) for JKT, -B/k^2 for Coulomb)
      shift      energy shift from the d3 term of (19)-(20); G0 is taken at E - shift
      coulomb_B  None, or B if V = -B/k^2 exactly (closed-form kernel (18))
      kinks      k where V has a kink (q_cut, 5 GeV), used as breakpoints
      q_cut, label; for JKT also Lam (Lambda_MS), C (continuity constant of (15)),
      C1 (continuity constant of (19)), alpha_mz, Lambda_R, p_match

    JKTPotential(alpha_mz=0.12)       the JKT potential, JKT's q_cut = 50 MeV
    JKTPotential.coulomb(alpha=0.15)  fixed-coupling Coulomb, same interface (no cut, no shift)
    """

    KAPPA = C_F * 48 * np.pi ** 2 / 27   # (4/3)(48 pi^2/27) in JKT (15)

    def __init__(self, alpha_mz=0.12, q_cut=0.05, Lambda_R=0.4, p_match=5.0, nF=N_F, nlo_bracket=True):
        """The full JKT potential: (13)-(15), cut at q_cut by (19) (JKT: 50 MeV), with the d3 term
        of (19) fixed by (20) and applied as an energy shift."""

        Lam = self.lambda_msbar_2loop(alpha_mz, nF=nF)
        C = self.continuity_constant_C(Lam, Lambda_R, p_match, nF, nlo_bracket)
        V = lambda p2: self.V_jkt(p2, Lam, Lambda_R, C, p_match, nF, nlo_bracket)
        self.V = lambda p2: self.V_cut(V, p2, q_cut)
        self.kinks = (q_cut, p_match)
        self.q_cut, self.Lam, self.C, self.C1 = q_cut, Lam, C, self.cut_constant_C1(V, q_cut)
        self.alpha_mz, self.Lambda_R, self.p_match = alpha_mz, Lambda_R, p_match
        self.coulomb_B = None
        self.label = f"JKT alpha_mZ={alpha_mz} q_cut={q_cut * 1e3:.0f} MeV"
        self.shift = self.energy_shift(self.V, q_cut, self.kinks)

    @classmethod
    def coulomb(cls, alpha, nlo_bracket=False, nF=N_F):
        """Fuks (12) / JKT (13) with alpha held fixed:
            V(k^2) = -4 pi C_F alpha_eff / k^2,  alpha_eff = alpha [1 + (31/3 - 10nF/9) alpha/4pi] (bracket optional).
        No cut and no shift: the Coulomb kernel has only a logarithmic singularity, which (21) handles."""

        a_eff = alpha * (1 + (31 / 3 - 10 * nF / 9) * alpha / (4 * np.pi)) if nlo_bracket else alpha
        B = 4 * np.pi * C_F * a_eff
        pot = cls.__new__(cls)
        pot.V = lambda k2: -B / np.asarray(k2, dtype=float)
        pot.coulomb_B, pot.shift, pot.kinks, pot.q_cut = B, 0.0, (), None
        pot.label = f"Coulomb alpha={alpha}" + (" +NLO" if nlo_bracket else "")
        return pot

    # ------------------------------------------------------------------ JKT Sec. 2
    @staticmethod
    def beta_coefficients(nF=N_F):
        """JKT below (14):  b0 = 11 - 2 nF/3,  b1 = 102 - 38 nF/3."""

        return 11 - 2 * nF / 3, 102 - 38 * nF / 3

    @staticmethod
    def alpha_2loop(Q2, Lam, nF=N_F):
        """JKT (14):  alpha_MS(Q^2) = 4 pi / [ b0 ln(Q^2/L^2) + (b1/b0) ln ln(Q^2/L^2) ]."""

        b0, b1 = JKTPotential.beta_coefficients(nF)
        with np.errstate(invalid="ignore", divide="ignore"):
            L = np.log(Q2 / Lam ** 2)
            return 4 * np.pi / (b0 * L + b1 / b0 * np.log(L))

    @staticmethod
    def lambda_msbar_2loop(alpha_mz, mZ=91.1876, nF=N_F):
        """Lambda_MS^(nF) such that JKT (14) at Q = mZ equals alpha_mz (JKT quote alpha_s
        directly; anchoring at mZ connects it to the value Fuks quote).
        Solves b0 x + (b1/b0) ln x = 4 pi/alpha_mz for x = ln(mZ^2/L^2)."""

        b0, b1 = JKTPotential.beta_coefficients(nF)
        x = brentq(lambda x: b0 * x + b1 / b0 * np.log(x) - 4 * np.pi / alpha_mz, 1e-6, 1e4)
        return mZ * np.exp(-x / 2)

    @staticmethod
    def V_pert(Q2, Lam, nF=N_F, nlo_bracket=True):
        """JKT (13):  V_pert(Q^2) = -(16 pi/3) alpha(Q^2)/Q^2 [1 + (31/3 - 10 nF/9) alpha(Q^2)/(4 pi)],
        alpha from JKT (14). nlo_bracket=False drops the square bracket (tree level)."""

        a = JKTPotential.alpha_2loop(Q2, Lam, nF)
        bracket = 1 + (31 / 3 - 10 * nF / 9) * a / (4 * np.pi) if nlo_bracket else 1.0
        return -4 * np.pi * C_F * a / Q2 * bracket

    @staticmethod
    def V_ir(p2, Lambda_R, C):
        """JKT (15), |p| < 5 GeV:  V(p) = -(4/3)(48 pi^2/27) (1/p^2) [ 1/ln(1 + p^2/L_R^2) + C ]."""

        with np.errstate(divide="ignore", invalid="ignore"):
            return -JKTPotential.KAPPA / p2 * (1 / np.log1p(p2 / Lambda_R ** 2) + C)

    @staticmethod
    def continuity_constant_C(Lam, Lambda_R, p_match=5.0, nF=N_F, nlo_bracket=True):
        """JKT (15): C from continuity at p = p_match,  V_ir(p_match^2) = V_pert(p_match^2)."""

        Q2 = p_match ** 2
        return (JKTPotential.V_pert(np.array(Q2), Lam, nF, nlo_bracket) * Q2 / (-JKTPotential.KAPPA)
                - 1 / np.log1p(Q2 / Lambda_R ** 2))

    @staticmethod
    def V_jkt(p2, Lam, Lambda_R, C, p_match=5.0, nF=N_F, nlo_bracket=True):
        """JKT (15): V_pert above p_match = 5 GeV, V_ir below."""

        p2 = np.asarray(p2, dtype=float)
        Q2m = p_match ** 2
        return np.where(p2 >= Q2m, JKTPotential.V_pert(np.maximum(p2, Q2m), Lam, nF, nlo_bracket),
                        JKTPotential.V_ir(np.minimum(p2, Q2m), Lambda_R, C))

    # ------------------------------------------------------------------ JKT Sec. 3: the cut
    @staticmethod
    def cut_constant_C1(V, q_cut):
        """JKT (19): C1 from continuity at p = q_cut,  -C1/(2 q_cut^2) = V(q_cut^2)."""

        qc2 = q_cut ** 2
        return -2 * qc2 * V(np.array([qc2]))[0]

    @staticmethod
    def V_cut(V, p2, q_cut):
        """JKT (19) without the d3 term (that term is the energy shift, see energy_shift):
            V_cut(p) = V(p)                  for |p| > q_cut
                     = -C1/(p^2 + q_cut^2)   for |p| < q_cut,
        C1 from continuity at p = q_cut:  -C1/(2 q_cut^2) = V(q_cut^2)."""

        p2 = np.asarray(p2, dtype=float)
        qc2 = q_cut ** 2
        C1 = JKTPotential.cut_constant_C1(V, q_cut)
        return np.where(p2 > qc2, V(np.maximum(p2, qc2)), -C1 / (p2 + qc2))

    @staticmethod
    def V_cut_position(Vc, r, q_cut, kinks=()):
        """JKT (20) without the d3 term:
            V_cut(r) = Int d^3q/(2pi)^3 e^{-iq.r} V_cut(q) = (1/(2 pi^2 r)) Int_0^inf dq q sin(qr) V_cut(q^2)."""

        f = lambda q: q * Vc(np.array([q * q]))[0]
        edges = [0.0, q_cut] + sorted(k for k in kinks if k > q_cut) + [50.0]
        I = sum(integrate.quad(f, a, b, weight="sin", wvar=r, limit=4000)[0] for a, b in zip(edges[:-1], edges[1:]))
        I += integrate.quad(f, 50.0, np.inf, weight="sin", wvar=r, limlst=200)[0]
        return I / (TWO_PI2 * r)

    @staticmethod
    def energy_shift(Vc, q_cut, kinks=(), r_norm=1.0, V_norm=-0.25):
        """JKT (19)-(20) with JKT (5): the d3 term V0 d3(p) of (19) is the position-space constant
        V0/(2pi)^3, chosen so that V_cut(r = 1 GeV^-1) = -1/4 GeV. Inside the LS equation it gives
        (V0/(2pi)^3) G(p), i.e. E -> E - V0/(2pi)^3. Returns V0/(2pi)^3 in GeV."""

        return V_norm - JKTPotential.V_cut_position(Vc, r_norm, q_cut, kinks)


# =====================================================================================
class LSSolver:
    """The Lippmann-Schwinger equation solved as in JKT Sec. 3:  LSSolver(potential).G(E, p).

    EQUATION MAP  (JKT equation  ->  method), read alongside the paper
      (16)  G(p) = G0(p) K(p)                                                 G
      (17)  K(p) = 1 + Int_0^inf dq Vhat(p,q) G0(q) K(q)                      (solved via (21))
      (18)  Vhat(p,q) = (1/4pi^2)(q/p) Int_{|p-q|}^{p+q} du u V(u)              kernel
      (21)  [1 - A(p)] K(p) + Int dq Vhat G0 [K(p) - K(q)] = 1                solve_K
      (22)  A(p) = Int dq Vhat(p,q) G0(q)                                     script_A (on aux_rule)
      (23)  Gauss-Legendre on (0,a1), (a1,a2), (0,1/a2) in t = 1/q            jkt_grid

    (16)-(17): write G = G0 K. K is smooth in p; the sharp structure of G (the on-shell peak
    for E > 0) is in G0, which is known exactly.
    (21) is (17) with K(p) A(p) added and subtracted -- an identity, not an approximation.
     - The peak of Vhat at q = p now multiplies K(p) - K(q) ~ K'(p)(p - q), which vanishes at
       q = p: what is left is ~1000x smaller and odd about q = p, so the grid nodes on either
       side largely cancel (figure: results/0_solver_checks/kernel_subtraction.png).
     - The peak itself is in A(p), (22), which contains no unknown: "a known function which can
       be evaluated by a numerical integration" (JKT). We use the fine rule aux_rule.

    Built once (independent of E): the node kernel W_ij = w_j Vhat(q_i, q_j), and Vhat on the
    auxiliary rule of A(q_i). Each energy then costs A(q_i) (one dot product per row) and one
    N x N linear solve. Default n = 48 per subinterval (N = 240; JKT used N = 100-400).
    The d3 term of JKT (19) enters as G0 evaluated at E - pot.shift.
    """

    def __init__(self, pot, n=48, m_t=M_T, Gamma_t=GAMMA_T, breaks=None):
        self.pot, self.m, self.Gam = pot, m_t, Gamma_t
        self._build_F()
        self.q, self.w = self.jkt_grid(n, m_t, breaks)
        with np.errstate(divide="ignore", invalid="ignore"):
            self.W = self.w * self.kernel(self.q[:, None], self.q[None, :])
        # q_j = q_i: in (21) this term multiplies K(q_i) - K(q_i) = 0 (and Vhat is infinite there for Coulomb)
        np.fill_diagonal(self.W, 0.0)
        self.aux = [self._aux(p) for p in self.q]

    def G0(self, E, p):
        """JKT (6) at E - shift (the d3 term of JKT (19))."""

        return Constants.G0(E - self.pot.shift, p, self.m, self.Gam)

    # ------------------------------------------------------------------ JKT (18)
    def _build_F(self, s_lo=1e-8, s_hi=1e16, n_cells=6000):
        """F(s) = Int_0^s dk^2 V(k^2) for the kernel (18): tabulated once (10-point Gauss-Legendre
        per cell, 6000 cells in ln s, the kinks of V as cell edges) and splined in ln s."""

        pot = self.pot
        if pot.coulomb_B is not None:
            return
        s = np.unique(np.concatenate([np.geomspace(s_lo, s_hi, n_cells), np.square(pot.kinks)]))
        a, b = s[:-1], s[1:]
        x, w = np.polynomial.legendre.leggauss(10)
        nodes = 0.5 * (b - a)[:, None] * (x + 1) + a[:, None]
        cells = (pot.V(nodes) * w).sum(1) * 0.5 * (b - a)
        F0 = pot.V(np.array([0.5 * s_lo]))[0] * s_lo           # V is finite at k -> 0 after the cut (19)
        self._s, self._F0 = s, F0
        self._Fspl = CubicSpline(np.log(s), np.concatenate([[F0], F0 + np.cumsum(cells)]))

    def F(self, s):
        s = np.asarray(s, dtype=float)
        return np.where(s < self._s[0], self._F0 * s / self._s[0], self._Fspl(np.log(np.maximum(s, self._s[0]))))

    def kernel(self, p, q):
        """JKT (18):  Vhat(p,q) = (1/4pi^2)(q/p) Int_{|p-q|}^{p+q} du u V(u)
                                 = (q/(8 pi^2 p)) [F((p+q)^2) - F((p-q)^2)].
        Coulomb V = -B/k^2 has the closed form  Vhat = -(B q/(8 pi^2 p)) ln[(p+q)^2/(p-q)^2].
        If the k^2 interval is short compared with its position, F would be differenced at nearly
        equal values, so that integral is done directly (16-point Gauss-Legendre)."""

        p, q = np.broadcast_arrays(np.asarray(p, dtype=float), np.asarray(q, dtype=float))
        lo, hi = (p - q) ** 2, (p + q) ** 2
        if self.pot.coulomb_B is not None:
            with np.errstate(divide="ignore"):
                return -self.pot.coulomb_B * q / (EIGHT_PI2 * p) * np.log(hi / lo)
        integral = self.F(hi) - self.F(lo)
        short = (hi - lo) < 0.5 * lo
        if np.any(short):
            x, w = np.polynomial.legendre.leggauss(16)
            l, h = lo[short], hi[short]
            nodes = 0.5 * (h - l)[:, None] * (x + 1) + l[:, None]
            integral[short] = (self.pot.V(nodes) * w).sum(1) * 0.5 * (h - l)
        return q / (EIGHT_PI2 * p) * integral

    # ------------------------------------------------------------------ JKT (23) and the rule for (22)
    @staticmethod
    def gauss_legendre(a, b, n):
        """n-point Gauss-Legendre nodes and weights on (a, b)."""

        x, w = np.polynomial.legendre.leggauss(n)
        return a + 0.5 * (b - a) * (x + 1), 0.5 * (b - a) * w

    @staticmethod
    def gauss_legendre_tail(a, n):
        """Int_a^inf dq f(q) = Int_0^{1/a} dt f(1/t)/t^2, Gauss-Legendre in t = 1/q (JKT (23), last subinterval)."""

        t, wt = LSSolver.gauss_legendre(0.0, 1.0 / a, n)
        return 1 / t, wt / t ** 2

    @staticmethod
    def jkt_grid(n, m_t=M_T, breaks=None):
        """JKT (23): n Gauss-Legendre points on each subinterval between `breaks` and on (breaks[-1], inf)
        as t = 1/q. JKT: breaks = (0, a1 = sqrt(m_t Gamma_t), a2 = m_t). Default here (0, 4, 16, 60, m_t):
        with JKT's three subintervals the on-shell peak of G0 above threshold (width m_t Gamma_t/2p0,
        1.2 GeV at m_t = 120, E = 2) falls between the nodes (0_solver_checks.py (f)).
        Nodes are returned in ascending q."""

        if breaks is None:
            breaks = (0.0, 4.0, 16.0, 60.0, m_t)
        q, w = zip(*[LSSolver.gauss_legendre(a, b, n) for a, b in zip(breaks[:-1], breaks[1:])],
                   LSSolver.gauss_legendre_tail(breaks[-1], n))
        q, w = np.concatenate(q), np.concatenate(w)
        order = np.argsort(q)
        return q[order], w[order]

    @staticmethod
    def aux_rule(p, n=16, grading=24, q_uniform=200.0, dq_uniform=2.0):
        """Fine quadrature rule on (0, inf) for A(p), JKT (22): n-point Gauss-Legendre panels with edges
          - graded geometrically towards q = p (p(1 - 2^-k), p(1 + 2^-k)), resolving the peak / log
            singularity of the kernel at q = p,
          - p (1 + 2^k) further out, for the kernel's fall-off away from q = p,
          - every 2 GeV up to 200 GeV, resolving G0's on-shell peak at q = sqrt(m_t E) for E > 0,
        and a t = 1/q panel beyond the last edge. Returns nodes and weights."""

        k = np.arange(1, grading + 1)
        edges = np.concatenate([[0.0], p * (1 - 2.0 ** -k), [p], p * (1 + 2.0 ** -k), p * (1 + 2.0 ** np.arange(0, 12)),
                                np.arange(dq_uniform, q_uniform + 1e-9, dq_uniform)])
        edges = np.unique(edges[edges <= max(q_uniform, 4097 * p)])
        q, w = zip(*[LSSolver.gauss_legendre(a, b, n) for a, b in zip(edges[:-1], edges[1:])],
                   LSSolver.gauss_legendre_tail(edges[-1], n))
        return np.concatenate(q), np.concatenate(w)

    def _aux(self, p):
        Q, Wq = self.aux_rule(p)
        with np.errstate(divide="ignore", invalid="ignore"):
            Vq = self.kernel(np.full_like(Q, p), Q)
        return Q, Wq * Vq

    # ------------------------------------------------------------------ JKT (22), (21), (16)
    def script_A(self, E, aux=None):
        """JKT (22):  A(p) = Int_0^inf dq Vhat(p,q) G0(q), for every node p = q_i (auxiliary rule)."""

        aux = self.aux if aux is None else aux
        return np.array([np.dot(WV, self.G0(E, Q)) for Q, WV in aux])

    def solve_K(self, E):
        """JKT (21) at the nodes p = q_i, integral -> sum over nodes (23):
            [1 - A_i + sum_j W_ij g_j] K_i - sum_j W_ij g_j K_j = 1,   g_j = G0(q_j)."""

        g = self.G0(E, self.q)
        Wg = self.W * g[None, :]
        M = np.diag(1 - self.script_A(E) + Wg.sum(1)) - Wg
        return np.linalg.solve(M, np.ones(len(self.q)))

    def G(self, E, p_out):
        """JKT (16): G(p) = G0(p) K(p), with K between the nodes from a cubic spline in ln q
        through K(q_i). K is smooth; the sharp structure of G is carried by G0 exactly."""

        K = self.solve_K(E)
        x, lp = np.log(self.q), np.log(np.asarray(p_out, dtype=float))
        K_out = CubicSpline(x, K.real)(lp) + 1j * CubicSpline(x, K.imag)(lp)
        return self.G0(E, p_out) * K_out

    def table(self, E_vals, p_vals):
        """G(E, p) on a grid: shape (len(E_vals), len(p_vals))."""

        return np.array([self.G(E, p_vals) for E in E_vals])


# =====================================================================================
class CoulombExact:
    """Exact S-wave G for a fixed-coupling Coulomb potential: no numerical solver, so it is the
    ground truth. Fuks write below their (12) that for constant alpha_s the LS equation "can be
    solved analytically ... which then yields a Whittaker function". This carries that out and
    ends in a one-dimensional integral (Step 5).

      CoulombExact.G(E, p, alpha)           fast, E and p arrays (160-point Gauss-Legendre)
      CoulombExact.G_adaptive(E, p, alpha)  slow reference, scalar (adaptive quadrature)

    Used to test LSSolver (0_solver_checks.py), to identify swData (2_identify_swdata.py), and
    as the "ours: Coulomb exact" curve (3_fuks_fig1.py, 4_fuks_figs2to5.py).

    DERIVATION (whiteboard order)                                   implemented in
    Step 1  Position-space equation, Fuks (7) with V = -C_F alpha / r, calE = E + i Gamma_t
            (JKT (4)):
                [ -Laplacian/m_t - C_F alpha/r - calE ] G(r) = -delta^3(r)
            Sign: Fuks (7) as printed has +delta^3 on the right. That G is MINUS the G of
            Fuks (11) and JKT (6) (free case: G(p) = 1/(p^2/m_t - calE)). We use the sign of
            Fuks (11)/JKT (6), which needs -delta^3 here. Fuks's swData tables follow (7):
            they store -G.
    Step 2  Solution regular at r -> infinity, with G -> -m_t/(4 pi r) at r -> 0:
                G(r) = -(m_t/(4 pi r)) Gamma(1-lam) W_{lam,1/2}(2 kap r),
                kap = sqrt(-m_t calE)  (Re kap > 0),   lam = C_F alpha m_t / (2 kap)
            (reduced mass m_t/2; lam = 1 is the 1S pole, E_1 = -m_t (C_F alpha)^2/4)      kap_lam
    Step 3  Laplace representation of the Whittaker function (DLMF Sec. 13.16), Re(1-lam) > 0:
                W_{lam,1/2}(z) = z e^{-z/2}/Gamma(1-lam) Int_0^inf dt e^{-zt} t^{-lam} (1+t)^{lam}
            so
                G(r) = -(m_t kap/(2 pi)) Int_0^inf dt t^{-lam} (1+t)^{lam} e^{-kap(1+2t) r}
    Step 4  S-wave Fourier transform, Gt(p) = (4 pi/p) Int_0^inf dr r sin(pr) G(r), with
                Int_0^inf dr r sin(pr) e^{-a r} = 2 p a/(a^2 + p^2)^2,   a = kap (1+2t):
    Step 5  Result:
                Gt(p) = -4 m_t kap^2 Int_0^inf dt  t^{-lam} (1+t)^{lam} (1+2t) / (kap^2 (1+2t)^2 + p^2)^2
                                                                     f_integrand, G, G_adaptive
    Step 6  Check alpha -> 0 (lam = 0): Gt(p) = -m_t/(p^2 + kap^2) = 1/(calE - p^2/m_t),
            exactly JKT (6) = Fuks (11).                            checked in 0_solver_checks.py (e)
    Step 7  Continuation to 1 <= Re lam < 2 (needed below the 1S pole): on [0,1] subtract
            the t -> 0 value,
                Int_0^1 t^{-lam} f(t) dt = Int_0^1 t^{-lam} [f(t) - f(0)] dt + f(0)/(1 - lam)
                                                                     G, G_adaptive
    """

    _u, _wu = np.polynomial.legendre.leggauss(160)
    _u, _wu = 0.5 * (_u + 1), 0.5 * _wu

    @staticmethod
    def kap_lam(E, alpha, m, Gam):
        """Step 2:  kap = sqrt(-m_t (E + i Gamma_t)) with Re kap > 0;  lam = C_F alpha m_t / (2 kap)."""

        kap = np.sqrt(-m * (np.asarray(E) + 1j * Gam))
        kap = np.where(kap.real < 0, -kap, kap)
        return kap, C_F * alpha * m / (2 * kap)

    @staticmethod
    def f_integrand(t, p, kap, lam):
        """Step 5 integrand without t^{-lam}:  f(t) = (1+t)^lam (1+2t) / (kap^2 (1+2t)^2 + p^2)^2."""

        return (1 + t) ** lam * (1 + 2 * t) / (kap ** 2 * (1 + 2 * t) ** 2 + p ** 2) ** 2

    @staticmethod
    def G(E, p, alpha, m=M_T, Gam=GAMMA_T):
        """Steps 5 and 7 with fixed Gauss-Legendre rules; E and p broadcastable arrays.
        [0,1]: t = u^4 (dt = 4u^3 du) tames t^{-lam}; [1,inf): t = 1/v (dt = dv/v^2).
        Agrees with G_adaptive to < 1e-4."""

        u, wu = CoulombExact._u, CoulombExact._wu
        E = np.asarray(E, float)
        p = np.asarray(p, float)
        kap, lam = CoulombExact.kap_lam(E, alpha, m, Gam)
        kap, lam = kap[..., None], lam[..., None]
        pp = np.broadcast_to(p, E.shape)[..., None]
        f = lambda t: CoulombExact.f_integrand(t, pp, kap, lam)
        t = u ** 4
        f0 = f(0.0 * t)
        I_01 = np.sum(wu * 4 * u ** 3 * t ** (-lam) * (f(t) - f0), -1) + f0[..., 0] / (1 - lam[..., 0])
        t2 = 1 / u
        I_1inf = np.sum(wu * t2 ** 2 * t2 ** (-lam) * f(t2), -1)
        return -4 * m * kap[..., 0] ** 2 * (I_01 + I_1inf)

    @staticmethod
    def G_adaptive(E, p, alpha, m=M_T, Gam=GAMMA_T):
        """Steps 5 and 7 by adaptive quadrature (scalar E, p). Sign as JKT (6) / Fuks (11)."""

        kap, lam = CoulombExact.kap_lam(E, alpha, m, Gam)
        kap, lam = complex(kap), complex(lam)
        f = lambda t: CoulombExact.f_integrand(t, p, kap, lam)

        def cquad(fun, a, b):
            re = integrate.quad(lambda t: fun(t).real, a, b, limit=400)[0]
            im = integrate.quad(lambda t: fun(t).imag, a, b, limit=400)[0]
            return re + 1j * im

        f0 = f(0.0)
        I_01 = cquad(lambda t: t ** (-lam) * (f(t) - f0) if t > 0 else 0j, 0, 1) + f0 / (1 - lam)   # Step 7
        I_1inf = cquad(lambda t: t ** (-lam) * f(t), 1, np.inf)
        return -4 * m * kap ** 2 * (I_01 + I_1inf)


# =====================================================================================
class FuksTables:
    """What Fuks et al. published, used exactly as published: their Green's-function tables
    (swData) and their interpolation routine CALCGREEN. Nothing of ours is in this class.

      reference/fuks2024/swData.M20.M5.txt  E = -20 .. -5 GeV (dE = 0.1),  p = 0.01 .. 99.76 GeV (dp = 0.25)
      reference/fuks2024/swData.M5.P20.txt  E =  -5 .. 20 GeV (dE = 0.02), same p grid
      columns: E, p, Re, Im.  The row E = -5 is in both files, with identical values.

    EQUATION MAP  (Fuks paper / code  ->  method)
      Fuks p.3   "the files hold the ratio G~/G~0": they do not; they hold G~ with the sign
                 of Fuks (7), i.e. -G in the convention of Fuks (11)              load_file
      Fuks (11)  convention used everywhere else; blocks(), points() return it    blocks, points, lookup
      Fuks Sec. 3, subroutine CALCGREEN: bilinear interpolation of swData at (E, p*), times
                 (E - p*^2/173 + 1.49 i) = 1/G~0 -> GREEN; GREEN = 1 outside
                 -20 < E < 20, p* < 80 GeV                                        calcgreen_ratio
      Fuks Sec. 3, MATRIX1: |M|^2 times |GREEN|^2 (Fuks (15), Events.weights), zero outside
                 W <= 350 GeV (the samples' generation window) and p* <= 50 GeV
                 (cut in 4_fuks_figs2to5.py)

    CALCGREEN is transcribed as written, including an off-by-one in the table index: the cell
    index is max(int((E - Emin)/dE), 1) (1-based), so it interpolates in the cell below and
    extrapolates (the fraction BB lies in [1, 2)). Effect on event weights: median 0.06%,
    max 1.6%. Fortran 1-based -> Python 0-based.
    """

    FILES = ("reference/fuks2024/swData.M20.M5.txt", "reference/fuks2024/swData.M5.P20.txt")
    E_SPLIT = -5.0   # take E < -5 from the first file, E >= -5 from the second

    @staticmethod
    def load_file(path):
        """One file on its regular grid:  E (n_E,), p (n_p,), stored values (n_E, n_p) = -G."""

        d = np.loadtxt(path)
        E, p = np.unique(d[:, 0]), np.unique(d[:, 1])
        return E, p, (d[:, 2] + 1j * d[:, 3]).reshape(len(E), len(p))

    @staticmethod
    def blocks():
        """The two regular blocks without the duplicated E = -5 row, G in the sign of Fuks (11):
        [(E, p, G) for E < -5, (E, p, G) for E >= -5]."""

        (Ea, pa, Sa), (Eb, pb, Sb) = (FuksTables.load_file(f) for f in FuksTables.FILES)
        keep = Ea < FuksTables.E_SPLIT - 1e-9
        return [(Ea[keep], pa, -Sa[keep]), (Eb, pb, -Sb)]

    @staticmethod
    def points():
        """All table points as flat arrays E, p, G (sign of Fuks (11)), E ascending."""

        E, p, G = [], [], []
        for Eb, pb, Gb in FuksTables.blocks():
            EE, PP = np.meshgrid(Eb, pb, indexing="ij")
            E.append(EE.ravel()); p.append(PP.ravel()); G.append(Gb.ravel())
        return np.concatenate(E), np.concatenate(p), np.concatenate(G)

    @staticmethod
    def lookup():
        """{(round(E, 2), round(p, 2)): G} for exact table points (sign of Fuks (11))."""

        E, p, G = FuksTables.points()
        return {(round(e, 2), round(q, 2)): g for e, q, g in zip(E, p, G)}

    @staticmethod
    def calcgreen_ratio():
        """CALCGREEN as a vectorised function (E, p*) -> GREEN. It multiplies the stored table
        values (-G~) by 1/G~0, so GREEN = -G~/G~0; only |GREEN|^2 enters the event weight."""

        tabs = [FuksTables.load_file(f) for f in FuksTables.FILES]   # stored values = -G~ as in the files
        steps = ((-20.0, 0.1), (-5.0, 0.02))                         # (Emin, dE) hard-coded in CALCGREEN

        def ratio(E, p):
            out = np.ones(len(E), dtype=complex)
            ok = (E > -20) & (E < 20) & (p < 80)
            for (Ed, Pd, G), (emin, de), sel in zip(tabs, steps, (E <= -5, E > -5)):
                m = ok & sel
                II = np.maximum(((E[m] - emin) / de).astype(int), 1) - 1
                JJ = np.maximum((p[m] / 0.25).astype(int), 1) - 1
                AA, BB = (E[m] - Ed[II]) / de, (p[m] - Pd[JJ]) / 0.25
                g = ((1 - AA) * (1 - BB) * G[II, JJ] + AA * (1 - BB) * G[II + 1, JJ]
                     + (1 - AA) * BB * G[II, JJ + 1] + AA * BB * G[II + 1, JJ + 1])
                out[m] = g * (E[m] - p[m] ** 2 / 173.0 + 1j * 1.49)
            return out
        return ratio


# =====================================================================================
class Events:
    """MadGraph events and the re-weighting of Fuks (15). Shared by Fuks's side and ours.

    EQUATION MAP
      Fuks (14)  g g -> t tbar -> b l+ nu b~ l'- nu~; t = b e+ ve, tbar = b~ mu- vm~    read_lhe
      Fuks (21)  p* = |p_t| in the t-tbar rest frame; Lorentz-invariant form
                 p* = sqrt(lambda(W^2, m_t1^2, m_t2^2)) / (2W); W, E = W - 2 m_t, m_tL/m_tH  kinematics
      Fuks (15)  |M|^2 -> |M|^2 |G(E,p*)/G0(E,p*)|^2 per event                       weights

    The MadGraph samples were deleted; kinematics() reads results/4_fuks_figs2to5/_cache_*.npz.
    THOSE FILES ARE THE ONLY COPY OF THE EVENTS: they are only ever read, and a cache is
    written only if none exists.
    """

    @staticmethod
    def weights(G_source, E, pstar):
        """Fuks (15): per-event factor |G(E, p*) / G0(E, p*)|^2 for a source G(E, p) -> complex array."""

        return np.abs(G_source(E, pstar) / Constants.G0(E, pstar)) ** 2

    @staticmethod
    def minv2(p):
        return p[..., 0] ** 2 - p[..., 1] ** 2 - p[..., 2] ** 2 - p[..., 3] ** 2

    @staticmethod
    def read_lhe(path):
        """(weights [pb], 4-momenta of the t and tbar systems) from the final-state partons."""

        wts, pt, ptb = [], [], []
        top_ids, anti_ids = {5, -11, 12}, {-5, 13, -14}
        with gzip.open(path, "rt") as f:
            in_ev, first = False, False
            for line in f:
                if line.startswith("<event"):
                    in_ev, first = True, True
                    a, b = np.zeros(4), np.zeros(4)
                    continue
                if not in_ev:
                    continue
                if line.startswith("</event"):
                    in_ev = False
                    pt.append(a); ptb.append(b)
                    continue
                if first:
                    wts.append(float(line.split()[2]))
                    first = False
                    continue
                c = line.split()
                if len(c) < 10 or c[0].startswith("<") or c[0].startswith("#") or c[1] != "1":
                    continue
                pid = int(c[0])
                mom = np.array([float(c[9]), float(c[6]), float(c[7]), float(c[8])])
                if pid in top_ids:
                    a += mom
                elif pid in anti_ids:
                    b += mom
        return np.array(wts), np.array(pt), np.array(ptb)

    @staticmethod
    def kinematics(lhe_path, cache):
        """Per-event dict w, W, E, pstar, mtL, mtH (see EQUATION MAP). Read from `cache` if it
        exists; otherwise computed from the LHE file and saved there."""

        if os.path.exists(cache):
            d = np.load(cache)
            return {k: d[k] for k in d.files}
        w, pt, ptb = Events.read_lhe(lhe_path)
        m1s, m2s, W2 = Events.minv2(pt), Events.minv2(ptb), Events.minv2(pt + ptb)
        W = np.sqrt(W2)
        lam = (W2 - m1s - m2s) ** 2 - 4 * m1s * m2s
        m1, m2 = np.sqrt(np.maximum(m1s, 0)), np.sqrt(np.maximum(m2s, 0))
        ev = dict(w=w, W=W, E=W - 2 * M_T, pstar=np.sqrt(np.maximum(lam, 0)) / (2 * W),
                  mtL=np.minimum(m1, m2), mtH=np.maximum(m1, m2))
        np.savez(cache, **ev)
        return ev
