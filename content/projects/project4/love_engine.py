"""
love_engine.py -- Love-wave dispersion and attenuation in a stack of N isotropic
viscoelastic layers over a half-space, via the Thomson-Haskell propagator matrix
and complex phase-velocity root finding.

Single source of truth for `Complex Velocity Method.ipynb`: every driver cell
does `from love_engine import *`.

======================================================================
CONVENTIONS  (identical to the write-up; every sign follows from these)

  time factor        e^{+i w t}
  horizontal factor  e^{-i k x},   k = w / c,   Re(k) > 0
  dissipative        Im(mu) > 0,   Q^-1 = Im(mu) / Re(mu)
  trapped mode       Im(nu_halfspace) < 0      (field decays as z -> +inf)
  attenuation        alpha = -Im(k) = w*Im(c)/|c|^2 > 0   =>   Im(c) > 0

  Layer propagator maps  f_bottom -> f_top  (UPWARD):
       T_j = [[      cos(nu h),  -sin(nu h)/(mu nu) ],
              [ mu nu sin(nu h),       cos(nu h)    ]],      det T_j = 1

  WARNING.  The version with the off-diagonal signs swapped is the DOWNWARD
  propagator exp(+A h).  Using it with the chain G = T_1...T_N, together with
  the (also wrong) branch Im(nu_halfspace) > 0, maps D -> -D and therefore
  leaves the dispersion curves correct while making every eigenfunction and
  every intermediate G entry wrong.  Each error alone shifts the roots.
======================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from functools import lru_cache
from scipy.optimize import brentq, minimize_scalar


@dataclass
class Layer:
    h: float       # thickness (m); np.inf for the half-space
    rho: float     # density (kg/m^3)
    beta: float    # shear velocity (m/s) -- see `f_ref` below for its meaning
    Qs: float      # target (frequency-independent) quality factor
    name: str = ""


# ----------------------------------------------------------------------
# 1.  GMB fitted to a constant Q.
#
#     Imposing Im(mu) = Q0^-1 Re(mu) and dividing by mu_u gives a system that is
#     EXACTLY LINEAR in the dimensionless weights
#         ahat_l = dmu * a_l / mu_u ,      mu_0/mu_u = 1 - sum_l ahat_l
#
#         sum_l  (w_i w_l + Q0^-1 w_l^2) / (w_i^2 + w_l^2) * ahat_l = Q0^-1
#
#     No linearisation is needed.  Dropping the Q0^-1 w_l^2 term (the usual
#     "denominator ~ mu_u" shortcut) biases the realised Q LOW by a relative
#     ~ dmu/mu_u (the truncated system enforces Im(mu)/mu_u = 1/Q0, but the
#     realised Q^-1 = Im(mu)/Re(mu) with Re(mu) < mu_u).
#     Prescribing dmu = (2/pi) mu_u / Q with ad-hoc equal weights does not
#     produce Q = Q0 at all -- it can be off by an order of magnitude.
# ----------------------------------------------------------------------
@dataclass
class GMB:
    Q0: float
    omega_l: np.ndarray
    ahat: np.ndarray

    @property
    def a_l(self):           return self.ahat / self.ahat.sum()   # sums to 1
    @property
    def dmu_over_muu(self):  return self.ahat.sum()


@lru_cache(maxsize=None)
def _fit_gmb(Q0, fmin, fmax, L, ncol):
    w_l = 2*np.pi*np.logspace(np.log10(fmin), np.log10(fmax), L)
    w_i = 2*np.pi*np.logspace(np.log10(fmin), np.log10(fmax), ncol)
    A = (w_i[:, None]*w_l[None, :] + (1.0/Q0)*w_l[None, :]**2) / (w_i[:, None]**2 + w_l[None, :]**2)
    ahat, *_ = np.linalg.lstsq(A, np.full(ncol, 1.0/Q0), rcond=None)
    return w_l, ahat


def _mu_unit(g: GMB, omega: float) -> complex:
    """mu(w)/mu_u for the fitted GMB."""
    return (1.0 - g.ahat.sum()) + np.sum(g.ahat*(1j*omega)/(g.omega_l + 1j*omega))


class GMBank:
    """One fitted GMB per distinct Q (each layer carries its own Q).

    `f_ref` pins the model to a measured velocity:
      * f_ref = None  -> `Layer.beta` is the UNRELAXED (infinite-frequency)
        velocity, so mu_u = rho*beta^2.  Then Re(beta(w)) < beta everywhere,
        by as much as 15% at low frequency when Q is small.
      * f_ref = f     -> mu_u is rescaled so that Re(mu(2 pi f)) == rho*beta^2,
        i.e. Re(beta(2 pi f)) == Layer.beta to O(Q^-2).
        This is what you want when the tabulated beta is a measured phase
        velocity at frequency f.  Omitting it offsets the dispersion curves.
    """
    def __init__(self, f_band=(0.05, 200.0), per_decade=4, ncol_per_mech=5, f_ref=1.0):
        self.f_band, self.per_decade, self.ncol_per_mech, self.f_ref = \
            f_band, per_decade, ncol_per_mech, f_ref

    def get(self, Q0) -> GMB:
        fmin, fmax = self.f_band
        L = max(3, int(round(self.per_decade*np.log10(fmax/fmin))))
        w_l, ahat = _fit_gmb(float(Q0), float(fmin), float(fmax), int(L), int(L*self.ncol_per_mech))
        return GMB(Q0, w_l, ahat)

    def mu_u(self, layer) -> float:
        g = self.get(layer.Qs)
        mu = layer.rho*layer.beta**2
        if self.f_ref is None:
            return mu
        return mu/_mu_unit(g, 2*np.pi*self.f_ref).real       # pin Re(mu(f_ref)) = rho beta^2

    def report(self, layers, freqs):
        print(f"GMB fit: band={self.f_band} Hz, f_ref={self.f_ref}")
        print("  target Q |  realised Q at " + " ".join(f"{f:.3g}Hz" for f in freqs))
        for Q0 in sorted({L_.Qs for L_ in layers}):
            g = self.get(Q0)
            vals = [_mu_unit(g, 2*np.pi*f) for f in freqs]
            flag = "   [negative weight present -- passivity checked below]" if (g.ahat < 0).any() else ""
            print(f"   {Q0:7.1f} |  " + "  ".join(f"{v.real/v.imag:7.2f}" for v in vals) + flag)


def complex_modulus(layer: Layer, omega: float, bank: GMBank) -> complex:
    return bank.mu_u(layer)*_mu_unit(bank.get(layer.Qs), omega)


def beta_complex(layer: Layer, omega: float, bank: GMBank) -> complex:
    b = np.sqrt(complex_modulus(layer, omega, bank)/layer.rho)
    return b if b.real > 0 else -b                       # branch: Re(beta) > 0


# ----------------------------------------------------------------------
# 2.  Vertical wavenumber and layer propagator.
#
#     T_j depends on nu_j only through cos(nu h) and sin(nu h)/nu, both EVEN in
#     nu.  The branch of the square root is therefore irrelevant inside the
#     layers; it matters only for the half-space, where the radiation condition
#     selects one exponential.
# ----------------------------------------------------------------------
def vertical_wavenumber(layer, omega, c, bank, halfspace=False):
    mu = complex_modulus(layer, omega, bank)
    nu = (omega/c)*np.sqrt((c/beta_complex(layer, omega, bank))**2 - 1.0 + 0j)
    if halfspace and nu.imag > 0:
        nu = -nu                                          # Im(nu_{N+1}) < 0
    return nu, mu


def _sin_over(a: complex) -> complex:
    """sin(a)/a with the removable singularity at a = 0 (i.e. at c = beta_j)."""
    return 1.0 - a*a/6.0 if abs(a) < 1e-6 else np.sin(a)/a


def layer_propagator(layer, omega, c, bank, normalize=True):
    nu, mu = vertical_wavenumber(layer, omega, c, bank)
    a = nu*layer.h
    C = np.cos(a)
    T = np.array([[C,                  -layer.h*_sin_over(a)/mu],
                  [mu*nu*np.sin(a),     C                      ]], dtype=complex)
    if normalize:
        # A POSITIVE scalar factor leaves the zeros of D unchanged but tames the
        # e^{|Im(nu h)|} growth that destroys Thomson-Haskell at high frequency.
        s = np.exp(abs(a.imag))
        if np.isfinite(s) and s > 1.0:
            T = T/s
    if not np.all(np.isfinite(T)):
        raise FloatingPointError(f"non-finite propagator: nu={nu}, mu={mu}, h={layer.h}")
    return T


def global_propagator(layers, omega, c, bank, normalize=True):
    G = np.eye(2, dtype=complex)
    for L_ in layers:                     # ordered product, j increasing to the right
        G = G @ layer_propagator(L_, omega, c, bank, normalize)
    return G


def dispersion(c, omega, layers, halfspace, bank, normalize=True) -> complex:
    """D(w,c) = G_21 - i mu_{N+1} nu_{N+1} G_22 ;  D = 0 at a Love mode."""
    G = global_propagator(layers, omega, c, bank, normalize)
    nu6, mu6 = vertical_wavenumber(halfspace, omega, c, bank, halfspace=True)
    return G[1, 0] - 1j*mu6*nu6*G[1, 1]


# ----------------------------------------------------------------------
# 3.  Roots.
#
#   * The trapped-mode window is  min_j Re(beta_j(w)) < Re(c) < Re(beta_{N+1}(w)),
#     using the FREQUENCY-DEPENDENT Re(beta_j(w)), not the table value.  A window
#     built from the table velocities admits leaky roots at the top (where
#     nu_{N+1} is propagating, not evanescent) and hides near-cutoff overtones at
#     the bottom.
#   * Sign changes of Re(D) on the real axis are only SEEDS.  Each is promoted to
#     the true complex root by Newton on the complex D.  At small Q the real-axis
#     scan silently loses modes, so roots from the previous frequency are carried
#     forward as additional seeds (continuation).
# ----------------------------------------------------------------------
def mode_window(layers, halfspace, omega, bank):
    lo = min(beta_complex(L_, omega, bank).real for L_ in layers)
    hi = beta_complex(halfspace, omega, bank).real
    return lo, hi


def _newton(f, c0, box, tol=1e-12, itmax=200):
    """Newton on the complex D.  `box` = (c_min, c_max) confines the iterate to the
    physical strip; a step outside it aborts rather than driving nu -> inf.
    Returns nan on failure -- never a silently wrong root."""
    c_min, c_max = box
    c = complex(c0)
    for _ in range(itmax):
        try:
            d = 1e-7*abs(c)
            f0 = f(c)
            fp = (f(c + d) - f(c - d))/(2*d)
        except (FloatingPointError, ZeroDivisionError, OverflowError, ValueError):
            return complex("nan")
        if fp == 0 or not np.isfinite(fp):
            return complex("nan")
        step = f0/fp
        if not np.isfinite(step):
            return complex("nan")
        c -= step
        if not np.isfinite(c) or not (c_min < c.real < c_max) or c.imag < 0:
            return complex("nan")
        if abs(step) < tol*abs(c):
            break
    return c


def find_modes(omega, layers, halfspace, bank, n_modes=5, n_scan=2000, seeds=()):
    """Complex phase velocities of the trapped Love modes, fundamental (slowest) first."""
    lo, hi = mode_window(layers, halfspace, omega, bank)
    if not hi > lo:
        return []

    def D(c):
        return dispersion(c, omega, layers, halfspace, bank)

    def D_real_safe(x):
        try:
            v = D(x).real
            return v if np.isfinite(v) else np.nan
        except (FloatingPointError, ZeroDivisionError, OverflowError, ValueError):
            return np.nan

    box = (0.5*lo, 1.5*hi)                                  # generous but finite
    cands = [_newton(D, s, box) for s in seeds]             # continuation seeds

    grid = np.linspace(lo*(1 + 1e-6), hi*(1 - 1e-6), n_scan)
    v = np.array([D_real_safe(c) for c in grid])
    Q_scale = 2.0*max(L_.Qs for L_ in layers)
    for i in range(len(grid) - 1):
        if not (np.isfinite(v[i]) and np.isfinite(v[i+1])):
            continue
        if v[i] == 0.0 or np.sign(v[i]) == np.sign(v[i+1]):
            continue
        try:
            cr = brentq(D_real_safe, grid[i], grid[i+1], xtol=1e-11, rtol=1e-14)
        except Exception:
            continue
        cands.append(_newton(D, complex(cr, cr/Q_scale), box))

    out = []
    for c in cands:
        if not np.isfinite(c) or c.imag <= 0:
            continue
        if not lo < c.real < hi:                            # reject leaky roots
            continue
        try:
            resid, scale = abs(D(c)), max(abs(D(c*1.01)), 1.0)
        except (FloatingPointError, ZeroDivisionError, OverflowError, ValueError):
            continue
        if resid > 1e-6*scale:                              # scale-free residual test
            continue
        if all(abs(c - z) > 1e-6*abs(c) for z in out):
            out.append(c)

    out.sort(key=lambda z: z.real)      # fundamental is the slowest at a given w
    return out[:n_modes]


def observables(c: complex, omega: float):
    """(phase velocity, attenuation alpha [1/m], Love-mode Q^-1).

    NOTE  c_phase = w/Re(k) = |c|^2/Re(c)  is NOT Re(c); they agree only to
    first order in Im(c)/Re(c).
    """
    k = omega/c
    return omega/k.real, -k.imag, 2.0*c.imag/c.real


# ----------------------------------------------------------------------
# 4.  Self-checks.  Cheap, and they catch every sign error in this file.
# ----------------------------------------------------------------------
def self_test(layers, halfspace, bank, omega, c):
    d = np.linalg.det(global_propagator(layers, omega, c, bank, normalize=False))
    assert abs(d - 1) < 1e-8, f"det(G) = {d}, expected 1"

    nu6, mu6 = vertical_wavenumber(halfspace, omega, c, bank, halfspace=True)
    assert nu6.imag < 0, f"half-space branch wrong: Im(nu) = {nu6.imag}"
    assert abs(np.exp(-1j*nu6*1e3)) < 1.0, "half-space field grows with depth"

    for L_ in list(layers) + [halfspace]:
        mu = complex_modulus(L_, omega, bank)
        assert mu.imag > 0, f"non-dissipative modulus, beta={L_.beta}"
        q = mu.real/mu.imag
        assert abs(q - L_.Qs)/L_.Qs < 0.05, f"realised Q={q:.2f} vs target {L_.Qs}"

    # single-layer reduction: the assembled D must equal the closed form
    #   D = mu1 nu1 sin(nu1 h1) - i mu_hs nu_hs cos(nu1 h1)
    # (the write-up's N=1 check; it pins every sign in T and V at once)
    L1 = layers[0]
    nu1, mu1 = vertical_wavenumber(L1, omega, c, bank)
    D_code = dispersion(c, omega, [L1], halfspace, bank, normalize=False)
    D_ref = mu1*nu1*np.sin(nu1*L1.h) - 1j*mu6*nu6*np.cos(nu1*L1.h)
    assert abs(D_code - D_ref) < 1e-9*abs(D_ref), "single-layer D mismatch"

    print("self_test passed: det(G)=1, half-space decays, Im(mu)>0, Q on target, "
          "single-layer D = closed form.")


def build_model(model_id=1):
    if model_id == 1:      # Table 4
        layers = [Layer(5, 2000, 180, 18),
                  Layer(5, 2000, 300, 30),
                  Layer(5, 2000, 420, 42)]
        halfspace = Layer(np.inf, 2000, 500, 50)
    elif model_id == 2:  
        layers = [Layer(10, 1800,  200, 20),
                  Layer(20, 2200,  600, 40),
                  Layer(25, 2400,  800, 50),
                  Layer(30, 2600, 1000, 60)]
        halfspace = Layer(np.inf, 3000, 1500, 100)
    elif model_id == 3:    # Table S1
        layers = [Layer(10e3, 2500, 2500, 250),
                  Layer(10e3, 2600, 3300, 330),
                  Layer(10e3, 2700, 3600, 360),
                  Layer(10e3, 2800, 4000, 400),
                  Layer(10e3, 2900, 4300, 430)]
        halfspace = Layer(np.inf, 3000, 4500, 450)
    else:  # Table S2
        layers = [Layer(3, 1820, 190, 19),
                  Layer(4, 1860, 320, 32),
                  Layer(5, 1910, 280, 28),
                  Layer(5, 1950, 460, 46),
                  Layer(6, 2000, 630, 63)]
        halfspace = Layer(np.inf, 2100, 750, 75)
    return layers, halfspace


# ----------------------------------------------------------------------
# 5.  Dispersion sweep and plots.
# ----------------------------------------------------------------------
def run_dispersion(model_id=1, n_modes=7, freq_min=0.1, freq_max=80.0, n_freq=200,
                   f_ref=1.0, per_decade=4):
    layers, halfspace = build_model(model_id)
    f_band = (min(freq_min, 0.05)/2.0, freq_max*2.0)     # fit Q beyond the swept band
    bank = GMBank(f_band=f_band, per_decade=per_decade, f_ref=f_ref)

    bank.report(layers + [halfspace], [freq_min, np.sqrt(freq_min*freq_max), freq_max])
    self_test(layers, halfspace, bank, omega=2*np.pi*np.sqrt(freq_min*freq_max),
              c=0.5*(min(L_.beta for L_ in layers) + halfspace.beta) + 1j)

    freq  = np.linspace(freq_min, freq_max, n_freq)
    c_ph  = np.full((n_modes, n_freq), np.nan)
    alpha = np.full((n_modes, n_freq), np.nan)
    Q_L   = np.full((n_modes, n_freq), np.nan)

    prev = []                                    # continuation seeds
    print("=" * 78)
    for i, f in enumerate(freq):
        omega = 2*np.pi*f
        modes = find_modes(omega, layers, halfspace, bank, n_modes, seeds=prev)
        prev = modes
        for m, c in enumerate(modes):
            cp, al, qinv = observables(c, omega)
            c_ph[m, i], alpha[m, i] = cp, al
            Q_L[m, i] = 1.0/qinv if qinv else np.nan
        if i % max(1, n_freq//12) == 0:
            ok = np.isfinite(c_ph[:, i])
            cs  = ", ".join(f"{v:.1f}" for v in c_ph[ok, i])
            als = ", ".join(f"{v*1e6:.2f}" for v in alpha[ok, i])
            print(f"f = {f:7.3f} Hz | {ok.sum()} mode(s) | c = [{cs}] m/s"
                  f" | alpha = [{als}] x 1e-6 /m")
    print("=" * 78)
    return freq, c_ph, alpha, Q_L, layers, halfspace, bank


def plot_results(freq, c_ph, alpha, Q_L, layers, halfspace, bank, model_id):
    """Two-panel figure: phase velocity and attenuation vs frequency."""
    n_modes = c_ph.shape[0]
    colors  = plt.cm.plasma(np.linspace(0, 0.88, n_modes))
    w = 2*np.pi*freq
    b_lo = np.array([min(beta_complex(L_, wi, bank).real for L_ in layers) for wi in w])
    b_hi = np.array([beta_complex(halfspace, wi, bank).real for wi in w])

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), dpi=600)
    for m in range(n_modes):
        ax[0].plot(freq, c_ph[m],  color=colors[m], lw=1.3, label=f"M{m}")
        ax[1].plot(freq, alpha[m], color=colors[m], lw=1.3, label=f"M{m}")
    ax[0].plot(freq, b_lo, "k:",  lw=1, label=r"$\min_j \Re\,\beta_j(\omega)$")
    ax[0].plot(freq, b_hi, "k--", lw=1, label=r"$\Re\,\beta_{N+1}(\omega)$")
    ax[0].set_xlabel("Frequency (Hz)"); ax[0].set_ylabel("Phase velocity (m/s)")
    ax[1].set_xlabel("Frequency (Hz)"); ax[1].set_ylabel(r"Attenuation $\alpha$ (1/m)")
    for a in ax:
        a.grid(True, ls=":", alpha=0.35)
        a.legend(fontsize="x-small", ncol=2)
        a.set_xlim(freq.min(), freq.max())
    fig.suptitle(f"Model {model_id}: Love waves in an isotropic viscoelastic stack",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"Model_{model_id}_PhaseVelocity_Attenuation.png", dpi=600, bbox_inches="tight")
    plt.show()


def plot_results_full(freq, c_ph, alpha, Q_L, layers, halfspace, bank, model_id):
    """Separate figures: (a) phase velocity, (b) attenuation, (c) 3D alpha-c-f."""
    n_modes = c_ph.shape[0]
    n_used  = int(np.isfinite(c_ph).any(axis=1).sum())
    if n_used == 0:
        print("No trapped modes in this band.")
        return
    colors = plt.cm.viridis(np.linspace(0, 0.9, n_used))
    w = 2*np.pi*freq
    b_lo = np.array([min(beta_complex(L_, wi, bank).real for L_ in layers) for wi in w])
    b_hi = np.array([beta_complex(halfspace, wi, bank).real for wi in w])

    # (a) phase velocity ------------------------------------------------
    plt.figure(figsize=(7, 5), dpi=600)
    for m in range(n_used):
        plt.plot(freq, c_ph[m], color=colors[m], lw=1.3, label=f"M{m}")
    plt.plot(freq, b_lo, "k:",  lw=1, label=r"$\min_j \Re\,\beta_j(\omega)$")
    plt.plot(freq, b_hi, "k--", lw=1, label=r"$\Re\,\beta_{N+1}(\omega)$")
    plt.xlabel("Frequency (Hz)"); plt.ylabel("Phase velocity (m/s)")
    plt.xlim(freq.min(), freq.max())
    plt.title(f"Model {model_id}: Love wave dispersion", fontsize=11, fontweight="bold")
    plt.grid(True, ls=":", alpha=0.3)
    plt.legend(fontsize="x-small", ncol=3, loc="upper right", edgecolor="black")
    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_PhaseVelocity.png", dpi=600, bbox_inches="tight")
    plt.show()

    # (b) attenuation ---------------------------------------------------
    plt.figure(figsize=(7, 5), dpi=600)
    for m in range(n_used):
        plt.plot(freq, alpha[m], color=colors[m], lw=1.3, label=f"M{m}")
    plt.xlabel("Frequency (Hz)"); plt.ylabel(r"Attenuation $\alpha$ (1/m)")
    plt.xlim(freq.min(), freq.max())
    plt.title(f"Model {model_id}: Love wave attenuation", fontsize=11, fontweight="bold")
    plt.grid(True, ls=":", alpha=0.3)
    plt.legend(fontsize="x-small", ncol=3, edgecolor="black")
    plt.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))
    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_Attenuation.png", dpi=600, bbox_inches="tight")
    plt.show()

    # (c) 3D alpha-c-f ---------------------------------------------------
    fig = plt.figure(figsize=(9, 6), dpi=600)
    ax = fig.add_subplot(111, projection="3d")
    for m in range(n_used):
        ok = np.isfinite(c_ph[m]) & np.isfinite(alpha[m])
        ax.plot(c_ph[m, ok], freq[ok], alpha[m, ok], color=colors[m], lw=1.3, label=f"M{m}")
    # ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Phase velocity (m/s)")
    ax.set_xlabel("Phase velocity (m/s)"); ax.set_ylabel("Frequency (Hz)") 
    ax.set_zlabel(r"$\alpha$ (1/m)")
    ax.set_title(rf"Model {model_id}: 3D $\alpha$-$c$-$f$", fontsize=11, fontweight="bold")
    ax.legend(fontsize="x-small")
    fig.tight_layout()
    fig.savefig(f"Model_{model_id}_3D_alpha_c_f.png", dpi=600, bbox_inches="tight")
    plt.show()


# ======================================================================
# 6.  Attenuation, three ways.
#
#   1. "Simplified"   alpha = omega / (2 c Q_eff),  Q_eff a thickness-weighted
#      mean of the layer Q's.  Ignores the mode shape entirely: a mode that
#      lives mostly in one layer should feel THAT layer's Q, not the average.
#
#   2. "Perturbation" fix Re(c) at the real-axis root of Re(D) and minimise
#      |D| over Im(c).  This is what the original notebook called the complex
#      velocity method.  It is accurate only while Im(c)/Re(c) is small.
#
#   3. "Exact"        solve D(omega, c) = 0 in the complex c-plane by Newton.
#      This is the actual eigenvalue; 1 and 2 are approximations to it.
#
# In every case  alpha = -Im(k) = omega*Im(c)/|c|^2,  NOT omega*Im(c)/Re(c)^2.
# ======================================================================
def alpha_simplified(omega, c_real, layers):
    total_h = sum(L_.h for L_ in layers if np.isfinite(L_.h))
    Q_eff = sum(L_.Qs*L_.h for L_ in layers if np.isfinite(L_.h))/total_h
    return omega/(2.0*c_real*Q_eff)


def alpha_perturbation(omega, c_real, layers, halfspace, bank):
    def obj(ci):
        try:
            v = abs(dispersion(complex(c_real, ci), omega, layers, halfspace, bank))
            return v if np.isfinite(v) else np.inf
        except (FloatingPointError, ZeroDivisionError, OverflowError, ValueError):
            return np.inf
    Q_min = min(L_.Qs for L_ in layers)
    res = minimize_scalar(obj, bounds=(0.0, c_real/Q_min), method="bounded",
                          options={"xatol": 1e-12})
    c_imag = max(res.x, 0.0)
    return omega*c_imag/(c_real**2 + c_imag**2), c_imag


# ----------------------------------------------------------------------
# 7.  Eigenfunction.
#
#   Start from the half-space state vector V = (1, -i mu_{N+1} nu_{N+1}) with
#   Im(nu_{N+1}) < 0, and walk UP through the layers with the field matrices,
#   solving a_j = E_j(h_j)^{-1} f_bottom at each step.
#
#   Evaluated at the TRUE COMPLEX root, this satisfies sigma_yz(0) = 0 to
#   machine precision.  With the wrong half-space branch the residual is O(1),
#   i.e. the eigenfunction does not satisfy the free surface at all.  Note the
#   branch error does NOT cancel here, because no propagator T_j is involved.
# ----------------------------------------------------------------------
def compute_eigenfunction(omega, c, layers, halfspace, bank, n_points=400, n_hs_depths=3.0):
    nu6, mu6 = vertical_wavenumber(halfspace, omega, c, bank, halfspace=True)
    f = np.array([1.0 + 0j, -1j*mu6*nu6], dtype=complex)     # V * B_{N+1},  B=1

    amps = []
    for j in range(len(layers) - 1, -1, -1):
        L_ = layers[j]
        nu, mu = vertical_wavenumber(L_, omega, c, bank)
        ep, em = np.exp(1j*nu*L_.h), np.exp(-1j*nu*L_.h)
        E_bot = np.array([[ep, em], [1j*mu*nu*ep, -1j*mu*nu*em]], dtype=complex)
        a_j = np.linalg.solve(E_bot, f)
        amps.insert(0, (a_j, nu, mu))
        E_top = np.array([[1.0, 1.0], [1j*mu*nu, -1j*mu*nu]], dtype=complex)
        f = E_top @ a_j

    u0, s0 = f
    nu1, mu1 = vertical_wavenumber(layers[0], omega, c, bank)
    residual = abs(s0)/(abs(mu1*nu1)*abs(u0))          # must be ~ 0

    H = sum(L_.h for L_ in layers)
    z_hs = n_hs_depths/abs(nu6.imag) if nu6.imag != 0 else H
    z = np.linspace(0.0, H + z_hs, n_points)
    u = np.zeros_like(z, dtype=complex)

    z0 = 0.0
    for (a_j, nu, mu), L_ in zip(amps, layers):
        m = (z >= z0) & (z <= z0 + L_.h)
        zl = z[m] - z0
        u[m] = a_j[0]*np.exp(1j*nu*zl) + a_j[1]*np.exp(-1j*nu*zl)
        z0 += L_.h
    # Half-space: B_{N+1} = 1 by construction, so u(z) = exp(-i nu6 (z-H)),
    # which equals 1 at z = H and therefore matches the layer stack continuously.
    m = z > H
    u[m] = np.exp(-1j*nu6*(z[m] - H))

    peak = np.max(np.abs(u))
    if peak > 0:
        u /= peak
    return z, u, residual


# ----------------------------------------------------------------------
# 8.  Three-methods study and its plots.
# ----------------------------------------------------------------------
def run_all_methods(model_id=1, n_modes=6, freq_min=2.0, freq_max=60.0, n_freq=80,
                    f_ref=1.0, per_decade=4):
    layers, halfspace = build_model(model_id)
    bank = GMBank(f_band=(min(freq_min, 0.05)/2.0, freq_max*2.0),
                  per_decade=per_decade, f_ref=f_ref)
    bank.report(layers + [halfspace], [freq_min, np.sqrt(freq_min*freq_max), freq_max])
    self_test(layers, halfspace, bank, omega=2*np.pi*np.sqrt(freq_min*freq_max),
              c=0.5*(min(L_.beta for L_ in layers) + halfspace.beta) + 1j)

    freq = np.linspace(freq_min, freq_max, n_freq)
    shape = (n_modes, n_freq)
    c_ph   = np.full(shape, np.nan)
    a_simp = np.full(shape, np.nan)
    a_pert = np.full(shape, np.nan)
    a_exact= np.full(shape, np.nan)
    Q_L    = np.full(shape, np.nan)

    prev = []
    print("=" * 78)
    for i, f in enumerate(freq):
        omega = 2*np.pi*f
        modes = find_modes(omega, layers, halfspace, bank, n_modes, seeds=prev)
        prev = modes
        for m, c in enumerate(modes):
            cp, al, qinv = observables(c, omega)
            c_ph[m, i]   = cp
            a_exact[m, i] = al
            Q_L[m, i]    = 1.0/qinv if qinv else np.nan
            a_simp[m, i] = alpha_simplified(omega, c.real, layers)
            a_pert[m, i] = alpha_perturbation(omega, c.real, layers, halfspace, bank)[0]
        if i % max(1, n_freq//10) == 0:
            print(f"f = {f:6.2f} Hz -> {len(modes)} mode(s)")
    print("=" * 78)
    return freq, c_ph, a_simp, a_pert, a_exact, Q_L, layers, halfspace, bank


def plot_method_comparison(freq, a_simp, a_pert, a_exact, model_id):
    n_used = int(np.isfinite(a_exact).any(axis=1).sum())
    colors = plt.cm.tab10(np.arange(max(n_used, 1)))
    fig, ax = plt.subplots(2, 2, figsize=(12, 9), dpi=600)

    for data, a, title in [(a_simp, ax[0, 0], "(a) Simplified:  $\\omega/(2cQ_{eff})$"),
                           (a_pert, ax[0, 1], "(b) Perturbation: fixed $\\Re c$"),
                           (a_exact, ax[1, 0], "(c) Exact: complex root of $D$")]:
        for m in range(n_used):
            ok = np.isfinite(data[m])
            if ok.any():
                a.plot(freq[ok], data[m, ok], color=colors[m], lw=1.4, label=f"M{m}")
        a.set_xlabel("Frequency (Hz)"); a.set_ylabel(r"$\alpha$ (1/m)")
        a.set_title(title, fontsize=10, fontweight="bold")
        a.grid(True, ls=":", alpha=0.3); a.legend(fontsize="x-small")
        a.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    a = ax[1, 1]
    for m in range(n_used):
        ok = np.isfinite(a_exact[m]) & np.isfinite(a_simp[m])
        if ok.any():
            a.plot(freq[ok], a_simp[m, ok]/a_exact[m, ok], color=colors[m], lw=1.4,
                   ls="--", label=f"M{m} simplified/exact")
        ok = np.isfinite(a_exact[m]) & np.isfinite(a_pert[m])
        if ok.any():
            a.plot(freq[ok], a_pert[m, ok]/a_exact[m, ok], color=colors[m], lw=1.4,
                   label=f"M{m} perturbation/exact")
    a.axhline(1.0, color="k", ls=":", lw=1)
    a.set_xlabel("Frequency (Hz)"); a.set_ylabel(r"ratio to exact $\alpha$")
    a.set_title("(d) Accuracy of the approximations", fontsize=10, fontweight="bold")
    a.grid(True, ls=":", alpha=0.3); a.legend(fontsize="xx-small", ncol=2)

    fig.suptitle(f"Model {model_id}: attenuation, three methods", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(f"Model_{model_id}_Attenuation_Comparison.png", dpi=600, bbox_inches="tight")
    plt.show()


def plot_eigenfunctions(freq, c_ph, layers, halfspace, bank, model_id, freq_indices=None):
    n_used = int(np.isfinite(c_ph).any(axis=1).sum())
    if n_used == 0:
        print("No modes to plot.")
        return
    n_show = min(n_used, 3)
    if freq_indices is None:
        n = len(freq)
        freq_indices = [n//4, n//2, 3*n//4]
    colors = plt.cm.plasma(np.linspace(0, 0.85, len(freq_indices)))

    fig, axes = plt.subplots(1, n_show, figsize=(4*n_show, 5), dpi=600, squeeze=False)
    axes = axes[0]
    H = sum(L_.h for L_ in layers)
    worst = 0.0

    for m in range(n_show):
        ax = axes[m]
        for k, idx in enumerate(freq_indices):
            cp = c_ph[m, idx]
            if not np.isfinite(cp):
                continue
            omega = 2*np.pi*freq[idx]
            modes = find_modes(omega, layers, halfspace, bank, n_used)
            if m >= len(modes):
                continue
            z, u, res = compute_eigenfunction(omega, modes[m], layers, halfspace, bank)
            worst = max(worst, res)
            ax.plot(np.abs(u), z, color=colors[k], lw=1.8, label=f"f = {freq[idx]:.2f} Hz")
        zc = 0.0
        for L_ in layers:
            zc += L_.h
            ax.axhline(zc, color="gray", ls="--", lw=0.8, alpha=0.6)
        ax.axhline(H, color="k", lw=1.0, alpha=0.8)
        ax.set_xlabel(r"$|u_y(z)|$ (normalised)")
        ax.set_ylabel("Depth (m)")
        ax.set_title(f"Mode {m}", fontsize=11, fontweight="bold")
        ax.invert_yaxis(); ax.grid(True, ls=":", alpha=0.3)
        ax.legend(fontsize="x-small")

    print(f"worst free-surface residual |sigma_yz(0)| / (|mu1 nu1||u(0)|) = {worst:.2e}")
    assert worst < 1e-8, "eigenfunction violates the free-surface condition"
    fig.suptitle(f"Model {model_id}: Love wave eigenfunctions "
                 f"(solid line = top of half-space)", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(f"Model_{model_id}_Eigenfunctions.png", dpi=600, bbox_inches="tight")
    plt.show()


def plot_quality_factor(freq, Q_L, layers, halfspace, model_id):
    """Q(omega) of the MODE, Q_L^-1 = 2 Im(c)/Re(c).

    NOT bounded by the layer Q's: first-order perturbation theory gives
    Q_L^-1 = (c/U) * sum_j w_j Q_j^-1 with w_j >= 0, sum w_j = 1, so where the
    dispersion is steep (group velocity U < c) the mode Q drops BELOW min Q_j
    by the factor U/c.  The min/max Q_j lines below are references, not bounds.
    """
    n_used = int(np.isfinite(Q_L).any(axis=1).sum())
    if n_used == 0:
        return
    colors = plt.cm.viridis(np.linspace(0, 0.9, n_used))
    plt.figure(figsize=(7, 5), dpi=600)
    for m in range(n_used):
        ok = np.isfinite(Q_L[m])
        plt.plot(freq[ok], Q_L[m, ok], color=colors[m], lw=1.4, label=f"M{m}")
    Qs = [L_.Qs for L_ in layers] + [halfspace.Qs]
    plt.axhline(min(Qs), color="k", ls=":",  lw=1, label=f"$\\min Q_j$ = {min(Qs):.0f}")
    plt.axhline(max(Qs), color="k", ls="--", lw=1, label=f"$\\max Q_j$ = {max(Qs):.0f}")
    plt.xlabel("Frequency (Hz)"); plt.ylabel(r"Love-mode quality factor $Q_L$")
    plt.title(f"Model {model_id}: $Q_L(\\omega) = \\Re(c)/(2\\,\\Im c)$",
              fontsize=11, fontweight="bold")
    plt.grid(True, ls=":", alpha=0.3); plt.legend(fontsize="x-small")
    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_Q_of_omega.png", dpi=600, bbox_inches="tight")
    plt.show()
