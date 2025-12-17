#!/usr/bin/env python
# coding: utf-8


# ## Attenuation Coefficient in Wave Propagation
# 
# The attenuation coefficient in wave propagation is directly related to the __imaginary part of the complex wavenumber__, which can be derived from the complex velocity.
# 
# ### 1. The Complex Wavenumber and Complex Velocity
# 
# In a lossy medium, a one-dimensional propagating wave can be described by:
# 
# $$
# p(x, t) = p_0 e^{i(\omega t - kx)}
# $$
# 
# where:
# - $ p $ is the wave field (e.g., pressure for sound, electric field for EM waves)
# - $ \omega $ is the angular frequency
# - $ k $ is the **wavenumber**
# 
# In a dissipative medium, the wavenumber becomes a **complex number**:
# 
# $$
# k = k_r - i \alpha
# $$
# 
# where:
# - $ k_r $ is the **real part** of the wavenumber, related to the phase velocity
# - $ \alpha $ is the **attenuation coefficient** (a positive real number), responsible for the exponential decay of the wave's amplitude
# 
# Substituting this into the wave equation gives:
# 
# $$
# p(x, t) = p_0 e^{i(\omega t - (k_r - i\alpha)x)} = p_0 e^{i(\omega t - k_r x)} \cdot e^{-\alpha x}
# $$
# 
# The term $ e^{-\alpha x} $ clearly shows the exponential amplitude decay.
# 
# The **complex velocity**, $ c^* $, is defined as:
# 
# $$
# c^* = \frac{\omega}{k}
# $$
# 
# ### 2. The Formula for the Attenuation Coefficient
# 
# We can rearrange the definition of the complex velocity to solve for $ k $:
# 
# $$
# k = \frac{\omega}{c^*}
# $$
# 
# Since $ k $ is complex ($ k = k_r - i\alpha $), $ c^* $ must also be complex. Let's express the complex velocity in terms of its real and imaginary parts. The standard convention is:
# 
# $$
# c^* = c_r + i c_i
# $$
# 
# where $ c_i $ is typically negative, representing loss.
# 
# Now, let's find the formula for $ \alpha $:
# 
# $$
# k = k_r - i\alpha = \frac{\omega}{c^*} = \frac{\omega}{c_r + i c_i}
# $$
# 
# To separate the real and imaginary parts, we multiply the numerator and denominator by the complex conjugate of the denominator:
# 
# $$
# k_r - i\alpha = \frac{\omega}{(c_r + i c_i)} \cdot \frac{(c_r - i c_i)}{(c_r - i c_i)} = \frac{\omega (c_r - i c_i)}{c_r^2 + c_i^2}
# $$
# 
# Now, equate the real and imaginary parts from both sides:
# 
# - **Real Part:** $ k_r = \dfrac{\omega c_r}{c_r^2 + c_i^2} $
# - **Imaginary Part:** $ -\alpha = \dfrac{-\omega c_i}{c_r^2 + c_i^2} $
# 
# From the imaginary part, we get the primary formula:
# 
# $$
# \boxed{\alpha = \frac{\omega \, (c_i)}{c_r^2 + c_i^2}}
# $$
# 
# $$
# \alpha = \frac{\omega , \text{Im}(c)}{|c|^2} \approx \frac{\omega c_i}{c_r^2} \quad \text{for small } c_i
# $$
# 
# For most geophysical problems, $c_i \ll c_r$, so this is valid.
# 
# This is the general formula for the attenuation coefficient in terms of the complex velocity $ c^* = c_r + i c_i $
# 
# ---
# 
# ### 3. The "Quality Factor" Convention
# 
# In many fields like geophysics and elastodynamics, a different representation is used. The complex velocity is expressed in terms of a **quality factor**, $ Q $, which is a dimensionless measure of attenuation (low $ Q $ means high attenuation).
# 
# For small attenuation ($ Q \gg 1 $), the standard formulation is:
# 
# $$
# \frac{1}{c^*} \approx \frac{1}{c} \left( 1 - \frac{i}{2Q} \right)
# $$
# 
# where $ c $ is the real-valued phase velocity at low loss.
# 
# Starting from the wavenumber:
# 
# $$
# k = \frac{\omega}{c^*} \approx \frac{\omega}{c} \left( 1 - \frac{i}{2Q} \right) = \frac{\omega}{c} - i \frac{\omega}{2cQ}
# $$
# 
# Comparing this to $ k = k_r - i\alpha $, we immediately identify:
# 
# $$
# \boxed{\alpha \approx \frac{\omega}{2cQ}}
# $$
# 
# This is a very common and useful approximation that clearly shows attenuation increases linearly with frequency for a constant $ Q $, which is observed in many real materials like the Earth's crust.
# 
# This is the standard high-Q attenuation model:
# 
# * attenuation ∝ frequency
# * phase velocity assumed constant
# * causal dispersion ignored (valid only for large (Q))
# * weak-attenuation models (high Q, no dispersion).
# 
# ### Summary
# 
# | Context | Formula for Attenuation Coefficient $$ \alpha $$ |
# | :--- | :--- |
# | **General Complex Velocity**<br>($ c^* = c_r + i c_i $) | $ \alpha = \dfrac{\omega \, (-c_i)}{c_r^2 + c_i^2} $ |
# | **Low-Loss "Q" Approximation**<br>($ Q \gg 1 $) | $ \alpha \approx \dfrac{\omega}{2cQ} $ |
# 
# The choice of formula depends on how the properties of your lossy medium are characterized—whether by a direct complex velocity or by a quality factor ($ Q $).
# 



# ### Model 1

# In[6]:

import numpy as np
import matplotlib.pyplot as plt
import cmath
from scipy.optimize import brentq, minimize_scalar
from dataclasses import dataclass
from typing import List, Tuple

# ================================================================
# 1. Layer and viscoelastic model definitions
# ================================================================

@dataclass
class Layer:
    """Represents a single layer in the model"""
    h: float      # thickness (m)
    rho: float    # density (kg/m³)
    beta: float   # shear-wave velocity (m/s)
    Qs: float     # quality factor
    name: str = ""

@dataclass
class GMBParameters:
    """Generalized Maxwell Body parameters"""
    L: int                   # Number of Maxwell bodies
    weights: List[float]     # Weighting coefficients a_l (must sum to 1)
    relax_freqs: List[float] # Relaxation frequencies (Hz)

# ================================================================
# 2. Complex shear modulus using Generalized Maxwell Body (GMB)
# ================================================================
def complex_modulus_gmb(beta, rho, Qs, omega, gmb):
    """
    Compute complex shear modulus using Generalized Maxwell Body model
    
    μ(ω) = μ_R + δμ * Σ(iω*a_l / (ω_l + iω))
    
    Args:
        beta: Shear wave velocity (m/s)
        rho: Density (kg/m³)
        Qs: Quality factor
        omega: Angular frequency (rad/s)
        gmb: GMB parameters
        
    Returns:
        Complex shear modulus
    """
    
    # high-frequency (unrelaxed) modulus
    mu_u      = rho * beta**2
        
    # Modulus defect (using constant Q approximation)
    Q_inv     = 1.0 / Qs
    # The modulus defect is approximated as: delta_mu = (2/pi*Q) * mu_u
    delta_mu  = mu_u * (2.0 / np.pi) * Q_inv
    
    # Relaxed modulus
    mu_R      = mu_u - delta_mu
    
    # GMB summation: μ(ω) = μ_R + Σ(iω*δμ*a_l / (ω_l + iω))
    mu_complex = mu_R + 0j
    
    for a_l, f_l in zip(gmb.weights, gmb.relax_freqs):
        omega_l      = 2 * np.pi * f_l
        mu_complex += (1j * omega * delta_mu * a_l) / (omega_l + 1j * omega)
            
    return mu_complex

# ================================================================
# 3. Vertical wavenumber (handles complex phase velocity)
# ================================================================
def compute_vertical_wavenumber(omega, c, layer, gmb):
    """
    Compute vertical wavenumber ν_j and complex modulus for a layer
    
    ν_j = k * sqrt((c/β_j)² - 1)
    where k = ω/c is the horizontal wavenumber (can be complex)
    
    Args:
        omega: Angular frequency (rad/s)
        c: Phase velocity (m/s) - can be complex
        layer: Layer parameters
        gmb: GMB parameters
        
    Returns:
        (nu, mu): Vertical wavenumber and complex modulus
    """
            
    k = omega / c  # Complex if c is complex
    
    # Complex shear modulus
    mu = complex_modulus_gmb(layer.beta, layer.rho, layer.Qs, omega, gmb)
    
    # Complex shear velocity: β(ω) = sqrt(μ(ω)/ρ)
    beta_complex = np.sqrt(mu / layer.rho)
    
    # ν = k * sqrt((c/β)² - 1)
    ratio = c / beta_complex
    nu    = k * np.sqrt(ratio**2 - 1)
    
    # Enforce Im(nu) >= 0 (decay direction into the layer)
    # This is a critical step for stability and ensuring energy decays away.
    if np.isfinite(nu) and (np.imag(nu) < 0):
        nu = -nu
    
    return nu, mu

# ================================================================
# 4. Thomson-Haskell layer propagator (handles complex c)
# ================================================================
def calculate_global_propagator(omega, c, layers, gmb):
    """
    Thomson–Haskell global propagator with support for complex phase velocity.
    Multiplies layer propagators T_j for all layers (top to bottom).
    """
    G = np.eye(2, dtype=complex)

    for layer in layers:
        nu, mu = compute_vertical_wavenumber(omega, c, layer, gmb)

        # Skip invalid values
        if not np.isfinite(nu) or not np.isfinite(mu):
            continue

        # Compute trig terms safely
        try:
            arg = nu * layer.h
            # For large imaginary parts (high frequency, deep evanescent wave),
            # use exponential form to avoid cmath.cos/sin overflow.
            if abs(arg.imag) > 50:
                # e^(i*nu*h) = e^(-Im(nu)*h) * e^(i*Re(nu)*h)
                exp_pos = np.exp(1j * arg)
                exp_neg = np.exp(-1j * arg)
                C = 0.5 * (exp_pos + exp_neg)
                S = 0.5j * (exp_neg - exp_pos)
            else:
                C = np.cos(arg)
                S = np.sin(arg)
        except (OverflowError, ValueError):
            C = 1.0 + 0j
            S = 0.0 + 0j

        # Compute mu*nu with protection against zero division
        mu_nu = mu * nu
        if np.isnan(mu_nu) or abs(mu_nu) < 1e-20:
            mu_nu = 1e-20 + 0j
            
        # force all to scalar complex
        C     = complex(C)
        S     = complex(S)
        mu_nu = complex(mu_nu)

        # Layer propagator matrix (Love wave SH)
        T = np.array([[         C, S / mu_nu],
                      [-mu_nu * S,        C]], dtype=complex)

        # Multiply into the global propagator
        G = G @ T

    return G

# ================================================================
# 5. Dispersion function (handles complex c)
# ================================================================
def dispersion_function(c, omega, layers, halfspace, gmb):
    """
    Evaluate dispersion function D(ω,c) for complex phase velocity (Love wave)
    
    D(ω,c) = G_21 - iμ_6*ν_6*G_22 = 0
    
    Args:
        c: Phase velocity (m/s) - can be complex
        omega: Angular frequency (rad/s)
        layers: List of layers
        halfspace: Half-space parameters
        gmb: GMB parameters
        
    Returns:
        (D_real, D_imag, D_abs): Components of dispersion function
    """
    
    # Compute global propagator
    G = calculate_global_propagator(omega, c, layers, gmb)
    
    # Half-space parameters
    nu6, mu6 = compute_vertical_wavenumber(omega, c, halfspace, gmb)
    
    # Ensure Im(ν_6) > 0 for proper decay into half-space (radiation condition)
    if nu6.imag < 0:
        nu6 = -nu6
    
    # Dispersion function: D(ω,c) = G_21 - iμ_6*ν_6*G_22
    D = G[1, 0] - 1j * mu6 * nu6 * G[1, 1]
    
    return np.real(D), np.imag(D), abs(D)

# ================================================================
# 6. Find real phase velocity roots
#    Robust root-finding (for each frequency) - real roots search
# ================================================================
def find_roots_dispersion(omega, layers, halfspace, gmb, n_modes=3):
    """
    Scan the interval [β1, β6] for sign changes in Re(D)
    and refine roots with brentq.
    Returns up to n_modes real roots (m/s).
    """
    beta_min = layers[0].beta
    beta_max = halfspace.beta
    # Scan grid slightly outside min/max shear velocities for safety
    c_grid   = np.linspace(beta_min * 1.001, beta_max * 0.999, 2000)

#     D_real = []
#     for c in c_grid:
#         val_r, _, _ = dispersion_function(c, omega, layers, halfspace, gmb)
#         D_real.append(val_r)
#     D_real = np.array(D_real)

    # Calculate Real part of the dispersion function
    D_real = np.array([dispersion_function(c, omega, layers, halfspace, gmb)[0] for c in c_grid])
    
    roots = []
    for i in range(len(c_grid) - 1):
        # Look for sign changes between adjacent grid points
        if (np.sign(D_real[i]) != np.sign(D_real[i + 1]) and 
            np.isfinite(D_real[i]) and np.isfinite(D_real[i + 1])):
            try:
                # Use Brent's method to find the precise root
                root = brentq(
                    lambda cc: dispersion_function(cc, omega, layers, halfspace, gmb)[0],
                    c_grid[i], c_grid[i + 1], xtol=1e-6, rtol=1e-6, maxiter=1000 )
                roots.append(root)
            except Exception:
                continue

    roots = np.unique(np.round(roots, 3))
    roots = np.sort(roots)
    return roots[:n_modes]

# ================================================================
# 7. ATTENUATION: Complex phase velocity method
# ================================================================
def compute_attenuation_complex(omega, c_real, layers, halfspace, gmb, prev_c_imag=None):
    """
    Compute attenuation coefficient α using complex phase velocity method
    by minimizing |D(c_real + i*c_imag)|.
    with mode tracking and stable bounds.
    
    Args:
        omega: Angular frequency (rad/s)
        c_real: Real phase velocity from dispersion (m/s)
        layers: List of layers
        halfspace: Half-space parameters
        gmb: GMB parameters
        c_imag_prev: Previous c_imag for continuity (optional)
    
    Returns:
        alpha: attenuation [1/m]
        c_imag: imaginary part of complex velocity [m/s]
    """
    
    # Conservative bounds
    c_imag_max = c_real / 10.0 # Upper bound for Im(c), typically much smaller

    # Initial guess for c_imag
    if prev_c_imag is not None and prev_c_imag > 0:
        c_imag_init = prev_c_imag
    else:
        # Initial guess based on average Q
        total_h = sum(layer.h for layer in layers if np.isfinite(layer.h))
        if total_h > 0:
            Q_avg = sum(layer.Qs * layer.h for layer in layers if np.isfinite(layer.h)) / total_h
            c_imag_init = max(c_real / (2 * Q_avg), 1e-6)
        else:
            c_imag_init = 1e-6

    # Objective function: minimize the absolute value of the dispersion function
    # Minimize |D|
    def objective(c_imag):
        c_complex = c_real + 1j * c_imag
        _, _, D_abs = dispersion_function(c_complex, omega, layers, halfspace, gmb)
        return D_abs

    try:
        # Use bounded minimization
        result = minimize_scalar(
            objective,
            bounds=(0, c_imag_max),
            method='bounded',
            options={'xatol':1e-10}
        )
        c_imag = result.x
    except:
        c_imag = 0.0

    # Exact formula for attenuation α = -Im[k] = ω * c_i / (|c|²)
    alpha = omega * c_imag / (c_real**2 + c_imag**2)
    return alpha, c_imag

# ================================================================
# 8. Build models (from Yuan et al. 2024)
# ================================================================
def build_model(model_id=1):
    if model_id == 1:
        # Table 4 of the paper
        layers = [
            Layer(5, 2000, 180, 18),
            Layer(5, 2000, 300, 30),
            Layer(5, 2000, 420, 42) ]
        halfspace = Layer(np.inf, 2000, 500, 50)
    elif model_id == 2:
        # Table S2 of the paper
        layers = [
            Layer(10, 1800, 200,  20),
            Layer(20, 2200, 600,  40),
            Layer(25, 2400, 800,  50),
            Layer(30, 2600, 1000, 60) ]
        halfspace = Layer(np.inf, 3000, 1500, 100)
    elif model_id == 3:
        # Table S1 of the paper
        layers = [
            Layer(10e3, 2500, 2500, 250),
            Layer(10e3, 2600, 3300, 330),
            Layer(10e3, 2700, 3600, 360),
            Layer(10e3, 2800, 4000, 400),
            Layer(10e3, 2900, 4300, 430) ]
        halfspace = Layer(np.inf, 3000, 4500, 450)
    else:
        layers = [
            Layer(3, 1820, 190, 19),
            Layer(4, 1860, 320, 32),
            Layer(5, 1910, 280, 28),
            Layer(5, 1950, 460, 46),
            Layer(6, 2000, 630, 63) ]
        halfspace = Layer(np.inf, 2100, 750, 75)

    return layers, halfspace

# ================================================================
# 9. Main computation with rigorous attenuation
# ================================================================
def run_dispersion(model_id=1, n_modes=2, freq_min=0.1, freq_max=80.0, n_freq=200):
    """
    Run dispersion analysis with rigorous complex-velocity attenuation calculation.
    Mode tracking ensures smooth α(f) curves.
    """
    layers, halfspace = build_model(model_id)
    
    # MODIFIED GMB PARAMETERS (L=5, logarithmic spread, equal weights)
    # This change better approximates Constant-Q behavior over the wide frequency band [0.1, 80] Hz,
    # thereby stabilizing the low-frequency attenuation for M0 and M1.
#     gmb = GMBParameters(
#         L=5, 
#         weights=[0.2, 0.2, 0.2, 0.2, 0.2], 
#         relax_freqs=[0.05, 0.5, 5.0, 50.0, 500.0]
#     )
    
    L     = 200
    f_min = 0.05
    f_max = 500
    
#     L     = 10
#     f_min = 0.05
#     f_max = 200

    gmb = GMBParameters(
        L=L,
        weights=[1 / L] * L,
        relax_freqs=np.logspace(np.log10(f_min), np.log10(f_max), L).tolist()
    )

    freq_range  = np.linspace(freq_min, freq_max, n_freq)
    phase_vel   = np.full((n_modes, n_freq), np.nan)
    attenuation = np.full_like(phase_vel, np.nan)

    # Track previous c_i for each mode for better mode continuity (important!)
    prev_c_imag_modes = [None]*n_modes

    print("="*70)
    print("Computing Love-wave dispersion with rigorous attenuation...")
    print("="*70)

    for i, f in enumerate(freq_range):
        omega = 2*np.pi*f
        roots = find_roots_dispersion(omega, layers, halfspace, gmb, n_modes)

        # Iterate through the found roots
        for m, c in enumerate(roots):
            phase_vel[m, i] = c

            # Mode tracking: use previous c_i as initial guess
            prev_c_imag = prev_c_imag_modes[m]  if i > 0 else None
            alpha, c_imag = compute_attenuation_complex(omega, c, layers, halfspace, gmb, prev_c_imag)

            attenuation[m, i] = alpha
            prev_c_imag_modes[m] = c_imag # store for next frequency

        if i % 5 == 0:
            print(f"f = {f:6.3f} Hz → {len(roots)} mode(s): c = {np.round(roots,1)}")
            if len(roots) > 0:
                vals = attenuation[:len(roots), i] * 1e6
                formatted = ", ".join(f"{v:.3f}" for v in vals)
                print(f"              → α = [{formatted}] × 10⁻⁶ [1/m]")

    print("="*70)
    return freq_range, phase_vel, attenuation, layers, halfspace, gmb

# ================================================================
# 10. Plotting section (Re-enabled)
# ================================================================
def plot_results(freq_range, phase_vel, attenuation, layers, halfspace, gmb, model_id):
    n_modes = phase_vel.shape[0]
    colors  = plt.cm.plasma(np.linspace(0, 1, n_modes))

    # ===========================
    # (a) Phase velocity 
    # ===========================
    plt.figure(figsize=(7, 5), dpi=600) 

    for m in range(n_modes):
        plt.plot(freq_range, phase_vel[m, :], color=colors[m], lw=1.3, label=f"M{m}")

    plt.axhline(layers[0].beta, color='black', ls=':', alpha=0.7, lw=1., 
                 label=f"$\\beta_1$ = {layers[0].beta:.0f} m/s")
    plt.axhline(halfspace.beta, color='black', ls='--', alpha=0.7, lw=1., 
                 label=f"$\\beta_6$ = {halfspace.beta:.0f} m/s")

    plt.xlabel("Frequency (Hz)", fontsize=11)
    plt.ylabel("Phase velocity (m/s)", fontsize=11)
    plt.xlim(freq_range.min(), freq_range.max())
    plt.ylim(150,525)
    plt.title(f"Model {model_id}: Love Wave in Isotropic ViscoElastic Media", fontsize=10, fontweight="bold")
    plt.grid(True, ls=':', alpha=0.3)
    plt.legend(ncol=len(plt.gca().lines), fontsize='x-small', loc='upper left', 
               frameon=True, framealpha=0.9, labelspacing=0.2, columnspacing=0.2,
               handlelength=1.0, handletextpad=0.4, edgecolor='black')
    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_PhaseVelocity.png", dpi=1200, bbox_inches='tight')
    plt.show()

    # ===========================
    # (b) Attenuation
    # ===========================
    plt.figure(figsize=(7, 5), dpi=600)

    for m in range(n_modes):
        # Scale for visualization if needed, but plotting in 1/m as requested
        plt.plot(freq_range, attenuation[m, :], color=colors[m], lw=1.5, label=f"M{m}")

    plt.xlabel("Frequency (Hz)", fontsize=11)
    plt.ylabel("Attenuation coefficient, $\\alpha$ (1/m)", fontsize=11)
    plt.xlim(freq_range.min(), freq_range.max())
#     plt.ylim(0, 0.1)
    plt.title(f"Model {model_id}: Love Wave Attenuation (Complex Velocity Method)", fontsize=11, fontweight="bold")
    plt.grid(True, ls=':', alpha=0.3)
    plt.legend(ncol=len(plt.gca().lines), fontsize='x-small', loc='upper left', 
               frameon=True, framealpha=0.9, labelspacing=0.2, columnspacing=0.2,
               handlelength=1.0, handletextpad=0.4, edgecolor='black')
    
    # Use scientific notation for the y-axis
#     plt.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_Attenuation.png", dpi=1200, bbox_inches='tight')
    plt.show()

    # ===========================
    # (c) 3D plot: α–c–f 
    # ===========================
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = plt.figure(figsize=(9, 6), dpi=600)
    ax  = fig.add_subplot(111, projection='3d')

    for m in range(n_modes):
        # Filter valid data
        mask = np.isfinite(phase_vel[m, :]) & np.isfinite(attenuation[m, :])
        ax.plot(phase_vel[m, mask], freq_range[mask], attenuation[m, mask], color=colors[m], lw=1.5, label=f"M{m}")
        ax.scatter(phase_vel[m, mask][::10], freq_range[mask][::10], 
                   attenuation[m, mask][::10], color=colors[m], s=10, alpha=0.5)
        
    ax.set_xlabel("Phase Velocity (m/s)", fontsize=10)
    ax.set_ylabel("Frequency (Hz)", fontsize=10)
    ax.set_zlabel("Attenuation $\\alpha$ (1/m)", fontsize=10)
    ax.set_title(f"Model {model_id}: 3D α–c–f Relationship", fontsize=11, fontweight="bold")
    ax.set_xlim(550, 100)
    ax.set_ylim(0, 80)
#     ax.set_zlim(0, 0.01)
    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_3D_α–c–f.png", dpi=1200, bbox_inches='tight')
    plt.show()    
    
# ================================================================
# Run all
# ================================================================
if __name__ == "__main__":
    
    # Configuration
    freq_min = 0.1   
    freq_max = 80
    n_freqs  = 200    # High resolution for smooth curves
    model_id = 1
    n_modes  = 7      # Enough modes to reproduce Fig 9c

    # Run dispersion solver with MODIFIED GMB parameters and RIGOROUS attenuation
    freq_range, phase_vel, attenuation, layers, halfspace, gmb = run_dispersion(model_id=model_id, 
        n_modes=n_modes,
        freq_min=freq_min, 
        freq_max=freq_max, 
        n_freq=n_freqs
    )
    
    # Plot results
    plot_results(freq_range, phase_vel, attenuation, layers, halfspace, gmb, model_id)
    
    print("\n✓ Analysis complete with rigorous attenuation calculation!")
    


# In[ ]:





# In[ ]:





# ### Model 3 - Deep Crustal Model

# In[3]:


import numpy as np
import matplotlib.pyplot as plt
import cmath
from scipy.optimize import brentq, minimize_scalar, fsolve
from dataclasses import dataclass
from typing import List, Tuple

# ================================================================
# 1. Layer and viscoelastic model definitions
# ================================================================

@dataclass
class Layer:
    """Represents a single layer in the model"""
    h: float      # thickness (m)
    rho: float    # density (kg/m³)
    beta: float   # shear-wave velocity (m/s)
    Qs: float     # quality factor
    name: str = ""

@dataclass
class GMBParameters:
    """Generalized Maxwell Body parameters"""
    L: int                    # Number of Maxwell bodies
    weights: List[float]      # Weighting coefficients a_l (must sum to 1)
    relax_freqs: List[float]  # Relaxation frequencies (Hz)        
        
# ================================================================
# 2. Complex shear modulus using Generalized Maxwell Body (GMB)
# ================================================================
def complex_modulus_gmb(beta, rho, Qs, omega, gmb):
    """
    Compute complex shear modulus using Generalized Maxwell Body model
    
    μ(ω) = μ_R + δμ * Σ(iω*a_l / (ω_l + iω))
    
    Args:
        beta: Shear wave velocity (m/s)
        rho: Density (kg/m³)
        Qs: Quality factor
        omega: Angular frequency (rad/s)
        gmb: GMB parameters
        
    Returns:
        Complex shear modulus
    """
    
    # high-frequency (unrelaxed) modulus
    mu_u     = rho * beta**2
       
    # Modulus defect (using constant Q approximation)
    Q_inv    = 1.0 / Qs
    delta_mu = mu_u * (2.0 / np.pi) * Q_inv
    
    # Relaxed modulus
    mu_R     = mu_u - delta_mu
    
    # GMB summation: μ(ω) = μ_R + Σ(iω*δμ*a_l / (ω_l + iω))
    mu_complex = mu_R + 0j
    
    for a_l, f_l in zip(gmb.weights, gmb.relax_freqs):
        omega_l     = 2 * np.pi * f_l
        mu_complex += (1j * omega * delta_mu * a_l) / (omega_l + 1j * omega)
        
    return mu_complex

# ================================================================
# 3. Vertical wavenumber (handles complex phase velocity)
# ================================================================
def compute_vertical_wavenumber(omega, c, layer, gmb):
    """
    Compute vertical wavenumber ν_j and complex modulus for a layer
    
    ν_j = k * sqrt((c/β_j)² - 1)
    where k = ω/c is the horizontal wavenumber (can be complex)
    
    Args:
        omega: Angular frequency (rad/s)
        c: Phase velocity (m/s) - can be complex
        layer: Layer parameters
        gmb: GMB parameters
        
    Returns:
        (nu, mu): Vertical wavenumber and complex modulus
    """
    
    k = omega / c  # Complex if c is complex
    
    # Complex shear modulus
    mu = complex_modulus_gmb(layer.beta, layer.rho, layer.Qs, omega, gmb)
    
    # Complex shear velocity: β(ω) = sqrt(μ(ω)/ρ)
    beta_complex = np.sqrt(mu / layer.rho)
    
    # ν = k * sqrt((c/β)² - 1)
    ratio = c / beta_complex
    nu    = k * np.sqrt(ratio**2 - 1)
    
    # Enforce Im(nu) >= 0 (decay direction)
    if np.isfinite(nu) and (np.imag(nu) < 0):
        nu = -nu
    
    return nu, mu

# ================================================================
# 4. Thomson-Haskell layer propagator (handles complex c)
# ================================================================
def calculate_global_propagator(omega, c, layers, gmb):
    """
    Thomson–Haskell global propagator with support for complex phase velocity.
    Multiplies layer propagators T_j for all layers (top to bottom).
    """
    G = np.eye(2, dtype=complex)

    for layer in layers:
        nu, mu = compute_vertical_wavenumber(omega, c, layer, gmb)

        # Ensure correct decay direction for evanescent/propagating
        if np.isfinite(nu) and nu.imag < 0:
            nu = -nu

        # Skip invalid values
        if not np.isfinite(nu) or not np.isfinite(mu):
            continue

        # Compute trig terms safely
        try:
            arg = nu * layer.h
            # For large imaginary parts, use exponential form to avoid overflow
            if abs(arg.imag) > 50:
                # e^(i*nu*h) = e^(-Im(nu)*h) * e^(i*Re(nu)*h)
                exp_pos = np.exp(1j * arg)
                exp_neg = np.exp(-1j * arg)
                C = 0.5 * (exp_pos + exp_neg)
                S = 0.5j * (exp_neg - exp_pos)
            else:
                C = cmath.cos(arg)
                S = cmath.sin(arg)
        except (OverflowError, ValueError):
            C = 1.0 + 0j
            S = 0.0 + 0j

        # Compute mu*nu with protection against zero division
        mu_nu = mu * nu
        if np.isnan(mu_nu) or abs(mu_nu) < 1e-20:
            mu_nu = 1e-20 + 0j
        
        # force all to scalar complex ----
        C     = complex(C)
        S     = complex(S)
        mu_nu = complex(mu_nu)

        # Layer propagator matrix
        T = np.array([[C, S / mu_nu],
                      [-mu_nu * S, C]], dtype=complex)

        # Multiply into the global propagator
        G = G @ T

    return G

# ================================================================
# 5. Dispersion function (handles complex c)
# ================================================================
def dispersion_function(c, omega, layers, halfspace, gmb):
    """
    Evaluate dispersion function D(ω,c) for complex phase velocity
    
    D(ω,c) = G_21 - iμ_6*ν_6*G_22
    
    Args:
        c: Phase velocity (m/s) - can be complex
        omega: Angular frequency (rad/s)
        layers: List of layers
        halfspace: Half-space parameters
        gmb: GMB parameters
        
    Returns:
        (D_real, D_imag, D_abs): Components of dispersion function
    """
    
    # Compute global propagator
    G = calculate_global_propagator(omega, c, layers, gmb)
    
    # Half-space parameters
    nu6, mu6 = compute_vertical_wavenumber(omega, c, halfspace, gmb)
    
    # Ensure Im(ν_6) > 0 for proper decay into half-space
    if nu6.imag < 0:
        nu6 = -nu6
        
    # Dispersion function: D(ω,c) = G_21 - iμ_6*ν_6*G_22
    D = G[1, 0] - 1j * mu6 * nu6 * G[1, 1]
    
    return np.real(D), np.imag(D), abs(D)

# ================================================================
# 6. Find real phase velocity roots
#    Robust root-finding (for each frequency) - real roots search
# ================================================================
def find_roots_dispersion(omega, layers, halfspace, gmb, n_modes=3):
    """
    Scan the interval [β1, β6] for sign changes in Re(D)
    and refine roots with brentq.
    Returns up to n_modes real roots (m/s).
    """
    beta_min = layers[0].beta
    beta_max = halfspace.beta
    c_grid   = np.linspace(beta_min * 1.001, beta_max * 0.999, 5000)

    D_real = []
    for c in c_grid:
        val_r, _, _ = dispersion_function(c, omega, layers, halfspace, gmb)
        D_real.append(val_r)
    D_real = np.array(D_real)
#     D_real = np.array([dispersion_function(c, omega, layers, halfspace, gmb)[0] for c in c_grid])

    roots = []
    for i in range(len(c_grid) - 1):
        if (np.sign(D_real[i]) != np.sign(D_real[i + 1]) and 
            np.isfinite(D_real[i]) and np.isfinite(D_real[i + 1])):
            try:
                root = brentq(
                    lambda cc: dispersion_function(cc, omega, layers, halfspace, gmb)[0],
                    c_grid[i], c_grid[i + 1], 
                    xtol=1e-6, rtol=1e-6, maxiter=500
                )
                roots.append(root)
            except Exception:
                continue

    roots = np.unique(np.round(roots, 3))
    roots = np.sort(roots)
    return roots[:n_modes]

# ================================================================
# 7. ATTENUATION: Complex phase velocity method (with Mode Tracking)
# ================================================================
# def compute_attenuation_complex(omega, c_real, layers, halfspace, gmb):
def compute_attenuation_complex(omega, c_real, layers, halfspace, gmb, prev_c_imag=None):
    """
    Compute attenuation coefficient using complex phase velocity approach by minimizing |D(c_real + i*c_imag)|.
    
    Uses prev_c_imag for mode tracking/initial guess stability.
    
    The complex phase velocity is: c = c_r + i*c_i
    where c_r is the real phase velocity from dispersion analysis.
    
    We find c_i by requiring that the complex dispersion function
    D(ω, c_r + i*c_i) = 0 + 0i
    
    The attenuation coefficient is:
        α = ω * c_i / c_r²  [1/m]
    
    Or equivalently, using complex wavenumber k = ω/c:
        k = k_r + i*α
        α = Im(k) = Im(ω/c) ≈ ω*c_i/c_r²  (for small c_i)
    
    Args:
        omega: Angular frequency (rad/s)
        c_real: Real phase velocity from dispersion (m/s)
        layers: List of layers
        halfspace: Half-space parameters
        gmb: GMB parameters
        
    Returns:
        alpha: Attenuation coefficient [1/m]
    
    Find c_i such that D(ω, c_r + i*c_i) = 0
    Then: α = ω * c_i / c_r²
        
    """
    
    # Initial guess for imaginary part (small perturbation)
    if prev_c_imag is not None and prev_c_imag > 0:
        c_imag_init = prev_c_imag
    else:
        # Estimate based on average Q
        total_h = sum(layer.h for layer in layers if np.isfinite(layer.h))
        if total_h > 0:
            Q_avg = sum(layer.Qs * layer.h for layer in layers if np.isfinite(layer.h)) / total_h
            c_imag_init = max(c_real / (2 * Q_avg), 1e-8) # Ensure small positive start # Initial guess
        else:
            c_imag_init = c_real * 1e-4
    
    # Conservative bounds
    c_imag_max = c_real / 10.0
    
    # Method 1: Minimize |D(ω, c_r + i*c_i)|
    def objective(c_imag):
        c_complex   = c_real + 1j * c_imag
        _, _, D_abs = dispersion_function(c_complex, omega, layers, halfspace, gmb)
        return D_abs
    
    # Use scipy minimize_scalar for robustness
    try:
        # Search in reasonable range
#         c_imag_max = c_real / 10.0  # Upper bound
        result = minimize_scalar(
            objective, 
            bounds=(0, c_imag_max),
            method='bounded',
            options={'xatol': 1e-10}
        )
        c_imag = result.x
        
        # Verify solution quality
        if result.fun > 1e-3:  # If residual is large, try fsolve
            # Method 2: Solve both Re(D) = 0 and Im(D) = 0
            def equations(c_imag):
                c_complex = c_real + 1j * c_imag
                D_re, D_im, _ = dispersion_function(c_complex, omega, layers, halfspace, gmb)
                return D_im  # Since Re(D) ≈ 0 at c_real, minimize Im(D)
            
            c_imag_sol = fsolve(equations, c_imag_init, full_output=False)
            if np.isfinite(c_imag_sol[0]) and c_imag_sol[0] > 0:
                c_imag = c_imag_sol[0]
    
    except Exception as e:
        print(f"Warning: Optimization failed at ω={omega:.4e}, using initial guess. Error: {e}")
        c_imag = c_imag_init
    
    # Ensure c_imag is positive and reasonable
    c_imag = abs(c_imag)
    if c_imag > c_real * 0.5:  # Sanity check: c_i shouldn't exceed 50% of c_r
        c_imag = c_real * 0.01
    
    # Calculate attenuation coefficient
    # α = Im(k) where k = ω/c = ω/(c_r + i*c_i)
    # For small c_i: α ≈ ω * c_i / c_r²
#     alpha = omega * c_imag / (c_real**2)
    
    # The exact formula
    alpha = omega * c_imag / (c_real**2 + c_imag**2)
    
    return alpha, c_imag

# ================================================================
# 8. Build models (from Yuan et al. 2024)
# ================================================================
def build_model(model_id=1):
    if model_id == 1:
        # Table 4 of the paper
        layers = [
            Layer(5, 2000, 180, 18),
            Layer(5, 2000, 300, 30),
            Layer(5, 2000, 420, 42) ]
        halfspace = Layer(np.inf, 2000, 500, 50)
    elif model_id == 2:
        # Table S2 of the paper
        layers = [
            Layer(10, 1800, 200,  20),
            Layer(20, 2200, 600,  40),
            Layer(25, 2400, 800,  50),
            Layer(30, 2600, 1000, 60) ]
        halfspace = Layer(np.inf, 3000, 1500, 100)
    elif model_id == 3:
        # Table S1 of the paper
        layers = [
            Layer(10e3, 2500, 2500, 250),
            Layer(10e3, 2600, 3300, 330),
            Layer(10e3, 2700, 3600, 360),
            Layer(10e3, 2800, 4000, 400),
            Layer(10e3, 2900, 4300, 430) ]
        halfspace = Layer(np.inf, 3000, 4500, 450)
    else:
        layers = [
            Layer(3, 1820, 190, 19),
            Layer(4, 1860, 320, 32),
            Layer(5, 1910, 280, 28),
            Layer(5, 1950, 460, 46),
            Layer(6, 2000, 630, 63) ]
        halfspace = Layer(np.inf, 2100, 750, 75)

    return layers, halfspace

# ================================================================
# 9. Main computation with rigorous attenuation
# ================================================================
def run_dispersion(model_id=1, n_modes=2, freq_min=0.01, freq_max=50.0, n_freq=100):
    """
    Run dispersion analysis with rigorous complex-velocity attenuation calculation.
    """
    layers, halfspace = build_model(model_id)
    
    # MODIFIED GMB PARAMETERS: Increased L and spread relaxation frequencies 
    # logarithmically to better approximate Constant-Q in the low-frequency range.
    L     = 20
    f_min = 0.001
    f_max = 10

    gmb = GMBParameters(
        L=L,
        weights=[1 / L] * L,
        relax_freqs=np.logspace(np.log10(f_min), np.log10(f_max), L).tolist()
    )
    
    # Frequency range (avoid exact zero)
    freq_range  = np.linspace(freq_min, freq_max, n_freq)
    phase_vel   = np.full((n_modes, len(freq_range)), np.nan)
    attenuation = np.full_like(phase_vel, np.nan)
    
    # CRITICAL: Track previous c_imag for stability (Mode Tracking)
    prev_c_imag_modes = [None]*n_modes
    
    print("=" * 80)
    print("Computing Love-wave dispersion with rigorous attenuation and Mode Tracking...")
    print("=" * 80)
    
    for i, f in enumerate(freq_range):
        omega = 2 * np.pi * f
        roots = find_roots_dispersion(omega, layers, halfspace, gmb, n_modes)
        
        for m, c in enumerate(roots):
            phase_vel[m, i] = c
            
            # Use previous imaginary part for guess
            prev_c_imag = prev_c_imag_modes[m]
            alpha, c_imag = compute_attenuation_complex(omega, c, layers, halfspace, gmb, prev_c_imag)
            
            # RIGOROUS attenuation using complex phase velocity
            attenuation[m, i] = alpha
            prev_c_imag_modes[m] = c_imag # Store for next frequency step
                
        if i % 10 == 0:
            print(f"f = {f:6.3f} Hz → {len(roots)} mode(s): c = {np.round(roots, 1)}")
            if len(roots) > 0:
                vals = attenuation[:len(roots), i] * 1e6
                formatted = ", ".join(f"{v:.3f}" for v in vals)
                print(f"              → α = [{formatted}] × 10⁻⁶ [1/m]")

    print("=" * 70)
    return freq_range, phase_vel, attenuation, layers, halfspace, gmb

# ================================================================
# 10. Plotting section
# ================================================================
def plot_results(freq_range, phase_vel, attenuation, layers, halfspace, gmb, model_id):
    n_modes = phase_vel.shape[0]
    colors  = plt.cm.viridis(np.linspace(0, 1, n_modes))

    # ===========================
    # (a) Phase velocity
    # ===========================
    plt.figure(figsize=(7, 5), dpi=600) 

    for m in range(n_modes):
        plt.plot(freq_range, phase_vel[m, :], color=colors[m], lw=1.3, label=f"M{m}")

    plt.axhline(layers[0].beta, color='black', ls=':', alpha=0.7, lw=1., 
                label=f"$\\beta_1$ = {layers[0].beta:.0f} m/s")
    plt.axhline(halfspace.beta, color='black', ls='--', alpha=0.7, lw=1., 
                label=f"$\\beta_6$ = {halfspace.beta:.0f} m/s")

    plt.xlabel("Frequency (Hz)", fontsize=11)
    plt.ylabel("Phase velocity (m/s)", fontsize=11)
    plt.xlim(freq_range.min(), freq_range.max())
    plt.ylim(2450, 4650)
    plt.title(f"Model {model_id}: Love Wave in Isotropic ViscoElastic Media", fontsize=10, fontweight="bold")
    plt.grid(True, ls=':', alpha=0.3)

    # Create legend
    leg = plt.legend(
        ncol=len(plt.gca().lines),  # one row with all items
        fontsize='x-small',
        loc='upper left',
        frameon=True,           # Ensure frame is visible
        fancybox=False,         # Use rectangular box instead of rounded
        framealpha=0.8,         # Transparency (0-1)
        handlelength=1.0,       # Shorter line handles
        handletextpad=0.4,      # Less padding between handle and text
        labelspacing=0.2,
        columnspacing=0.3 )

    # Add border color
    frame = leg.get_frame()
    frame.set_edgecolor('black')

    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_PhaseVelocity.png", dpi=1200, bbox_inches='tight')
    plt.show()

    # ===========================
    # (b) Attenuation 
    # ===========================
    plt.figure(figsize=(7, 5), dpi=600)

    for m in range(n_modes):
        plt.plot(freq_range, attenuation[m, :], color=colors[m], lw=1.5, label=f"M{m}")

    plt.xlabel("Frequency (Hz)", fontsize=11)
    plt.ylabel("Attenuation coefficient $\\alpha$ (1/m)", fontsize=11)
    plt.xlim(freq_range.min(), freq_range.max())
    plt.title(f"Model {model_id}: Love Wave Attenuation (Complex Velocity Method)", fontsize=11, fontweight="bold")
    plt.grid(True, ls=':', alpha=0.3)
    plt.legend(
        ncol=len(plt.gca().lines),  # one row with all items
        fontsize='x-small',
        loc='upper left',
        frameon=True,           # Ensure frame is visible
        fancybox=False,         # Use rectangular box instead of rounded
        framealpha=0.8,         # Transparency (0-1)
        handlelength=1.0,       # Shorter line handles
        handletextpad=0.4,      # Less padding between handle and text
        labelspacing=0.2,
        columnspacing=0.2, edgecolor='black')  
    
    plt.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_Attenuation.png", dpi=1200, bbox_inches='tight')
    plt.show() 
    
    # ===========================
    # (c) 3D plot: α–c–f 
    # ===========================
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = plt.figure(figsize=(9, 6), dpi=600)
    ax  = fig.add_subplot(111, projection='3d')

    for m in range(n_modes):
        # Filter valid data
        mask = np.isfinite(phase_vel[m, :]) & np.isfinite(attenuation[m, :])
        ax.plot(freq_range[mask], phase_vel[m, mask], attenuation[m, mask], 
                color=colors[m], lw=1.5, label=f"M{m}")
        ax.scatter(freq_range[mask][::5], phase_vel[m, mask][::5], 
                   attenuation[m, mask][::5], color=colors[m], s=3, alpha=0.8)
        
    ax.set_xlabel("Frequency (Hz)", fontsize=10)
    ax.set_ylabel("Phase Velocity (m/s)", fontsize=10)
    ax.set_zlabel("Attenuation $\\alpha$ (1/m)", fontsize=10)
    ax.set_title(f"Model {model_id}: 3D α–c–f Relationship", fontsize=11, fontweight="bold")
    ax.set_xlim(0.5, 0.)
    ax.set_ylim(4500, 2500)
    ax.set_zlim(0, 3e-7)
    
    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_3D_α–c–f.png", dpi=1200, bbox_inches='tight')
    plt.show()


# ================================================================
# Run all
# ================================================================
if __name__ == "__main__":
    # Configuration
    freq_min = 0.01   # Start from small positive value
    freq_max = 0.5
    n_freqs  = 200     
    model_id = 3
    n_modes  = 8 

    # Run dispersion solver with RIGOROUS attenuation
    freq_range, phase_vel, attenuation, layers, halfspace, gmb = run_dispersion(
        model_id=model_id, 
        n_modes=n_modes,
        freq_min=freq_min, 
        freq_max=freq_max, 
        n_freq=n_freqs
    )
    
    # Plot results
    plot_results(freq_range, phase_vel, attenuation, layers, halfspace, gmb, model_id)
    
    print("\n✓ Analysis complete with rigorous attenuation calculation and mode tracking.!")
    


# In[ ]:





# In[ ]:





# ### Model 4 - Near Surface Model

# In[4]:


import numpy as np
import matplotlib.pyplot as plt
import cmath
from scipy.optimize import brentq, minimize_scalar, fsolve
from scipy.signal import medfilt
from dataclasses import dataclass
from typing import List, Tuple, Optional

# ================================================================
# 1. Layer and viscoelastic model definitions
# ================================================================

@dataclass
class Layer:
    """Represents a single layer in the model"""
    h: float      # thickness (m)
    rho: float    # density (kg/m³)
    beta: float   # shear-wave velocity (m/s)
    Qs: float     # quality factor
    name: str = ""

@dataclass
class GMBParameters:
    """Generalized Maxwell Body parameters"""
    L: int                    # Number of Maxwell bodies
    weights: List[float]      # Weighting coefficients a_l (must sum to 1)
    relax_freqs: List[float]  # Relaxation frequencies (Hz)        
        
# ================================================================
# 2. Complex shear modulus using Generalized Maxwell Body (GMB)
# ================================================================
def complex_modulus_gmb(beta, rho, Qs, omega, gmb):
    """
    Compute complex shear modulus using Generalized Maxwell Body model
    
    μ(ω) = μ_R + δμ * Σ(iω*a_l / (ω_l + iω))
    """
    
    # high-frequency (unrelaxed) modulus
    mu_u = rho * beta**2
       
    # Modulus defect (using constant Q approximation)
    Q_inv = 1.0 / Qs
    # delta_mu = mu_u * (2/pi) * (1/Q)
    delta_mu = mu_u * (2.0 / np.pi) * Q_inv
    
    # Relaxed modulus
    mu_R = mu_u - delta_mu
    
    # GMB summation
    mu_complex = mu_R + 0j
    
    for a_l, f_l in zip(gmb.weights, gmb.relax_freqs):
        omega_l = 2 * np.pi * f_l
        mu_complex += (1j * omega * delta_mu * a_l) / (omega_l + 1j * omega)
        
    return mu_complex

# ================================================================
# 3. Vertical wavenumber (handles complex phase velocity)
# ================================================================
def compute_vertical_wavenumber(omega, c, layer, gmb):
    """
    Compute vertical wavenumber ν_j and complex modulus for a layer
    """
    
    k = omega / c  # Complex if c is complex
    
    # Complex shear modulus
    mu = complex_modulus_gmb(layer.beta, layer.rho, layer.Qs, omega, gmb)
    
    # Complex shear velocity: β(ω) = sqrt(μ(ω)/ρ)
    beta_complex = np.sqrt(mu / layer.rho)
    
    # ν = k * sqrt((c/β)² - 1)
    ratio = c / beta_complex
    nu = k * np.sqrt(ratio**2 - 1)
    
    # Enforce Im(nu) >= 0 (decay direction)
    if np.isfinite(nu) and (np.imag(nu) < 0):
        nu = -nu
    
    return nu, mu


# ================================================================
# 4. Thomson-Haskell layer propagator (handles complex c)
# ================================================================
def calculate_global_propagator(omega, c, layers, gmb):
    """
    Thomson–Haskell global propagator with support for complex phase velocity.
    """
    G = np.eye(2, dtype=complex)

    for layer in layers:
        nu, mu = compute_vertical_wavenumber(omega, c, layer, gmb)

        # Ensure correct decay direction
        if np.isfinite(nu) and nu.imag < 0:
            nu = -nu

        # Skip invalid values
        if not np.isfinite(nu) or not np.isfinite(mu):
            continue

        # Compute trig terms safely
        try:
            arg = nu * layer.h
            # For large imaginary parts, use exponential form to avoid overflow
            if abs(arg.imag) > 50:
                exp_pos = np.exp(1j * arg)
                exp_neg = np.exp(-1j * arg)
                C       = 0.5 * (exp_pos + exp_neg)
                S       = 0.5j * (exp_neg - exp_pos)
            else:
                C = cmath.cos(arg)
                S = cmath.sin(arg)
        except (OverflowError, ValueError):
            C = 1.0 + 0j
            S = 0.0 + 0j

        # Compute mu*nu with protection against zero division
        mu_nu = mu * nu
        if np.isnan(mu_nu) or abs(mu_nu) < 1e-20:
            mu_nu = 1e-20 + 0j
        
        # Force all to scalar complex
        C     = complex(C)
        S     = complex(S)
        mu_nu = complex(mu_nu)

        # Layer propagator matrix (Love wave SH)
        T = np.array([[         C, S / mu_nu],
                      [-mu_nu * S,         C]], dtype=complex)

        # Multiply into the global propagator
        G = G @ T

    return G


# ================================================================
# 5. Dispersion function (handles complex c)
# ================================================================
def dispersion_function(c, omega, layers, halfspace, gmb):
    """
    Evaluate dispersion function D(ω,c) for complex phase velocity (Love wave)
    
    D(ω,c) = G_21 - iμ_6*ν_6*G_22
    """
    
    # Compute global propagator
    G = calculate_global_propagator(omega, c, layers, gmb)
    
    # Half-space parameters
    nu6, mu6 = compute_vertical_wavenumber(omega, c, halfspace, gmb)
    
    # Ensure Im(ν_6) > 0 for proper decay into half-space
    if nu6.imag < 0:
        nu6 = -nu6
        
    # Dispersion function: D(ω,c) = G_21 - iμ_6*ν_6*G_22
    D = G[1, 0] - 1j * mu6 * nu6 * G[1, 1]
    
    return np.real(D), np.imag(D), abs(D)

# ================================================================
# 6. Find real phase velocity roots
# ================================================================
def find_roots_dispersion(omega, layers, halfspace, gmb, n_modes=3):
    """
    Scan the interval [β1, β6] for sign changes in Re(D)
    and refine roots with brentq.
    """
    beta_min = layers[0].beta
    beta_max = halfspace.beta
    # Use higher grid resolution for more robust root finding
    c_grid = np.linspace(beta_min * 1.001, beta_max * 0.999, 5000)

    D_real = []
    for c in c_grid:
        val_r, _, _ = dispersion_function(c, omega, layers, halfspace, gmb)
        D_real.append(val_r)
    D_real = np.array(D_real)

    roots = []
    for i in range(len(c_grid) - 1):
        if (np.sign(D_real[i]) != np.sign(D_real[i + 1]) and 
            np.isfinite(D_real[i]) and np.isfinite(D_real[i + 1])):
            try:
                root = brentq(
                    lambda cc: dispersion_function(cc, omega, layers, halfspace, gmb)[0],
                    c_grid[i], c_grid[i + 1], 
                    xtol=1e-6, rtol=1e-6, maxiter=500
                )
                roots.append(root)
            except Exception:
                continue

    roots = np.unique(np.round(roots, 3))
    roots = np.sort(roots)
    return roots[:n_modes]

# ================================================================
# 7. ATTENUATION: Complex phase velocity method (CLEANED)
# ================================================================
def compute_attenuation_complex(omega, c_real, layers, halfspace, gmb, 
                                 c_imag_prev: Optional[float] = None):
    """
    Compute attenuation coefficient using complex phase velocity
    by minimizing |D(c_real + i*c_imag)|.
    
    Simplified and robust implementation.
    """
    
    if omega < 1e-10:
        return 0.0, 0.0
    
    # Determine safe maximum bound for c_imag
    Q_min = min(layer.Qs for layer in layers) if layers else 10.0
    c_imag_max = c_real / (2.0 * Q_min) * 2.0  # Conservative bound (e.g., Q=10 gives c_i_max = c_r/10)
    
    # Initial guess - use previous value if available for continuity
    c_imag_init = c_imag_prev if c_imag_prev is not None and np.isfinite(c_imag_prev) and c_imag_prev > 0 else 1e-6
    
    # Ensure initial guess is within bounds
    c_imag_init = np.clip(c_imag_init, 1e-12, c_imag_max * 0.5)

    # Objective function: minimize the absolute value of the dispersion function
    def objective(c_imag):
        c_complex = c_real + 1j * c_imag
        try:
            _, _, D_abs = dispersion_function(c_complex, omega, layers, halfspace, gmb)
            return D_abs if np.isfinite(D_abs) else 1e10
        except:
            return 1e10
    
    # Optimization with bounded search (0 up to c_imag_max)
    try:
        result = minimize_scalar(
            objective, 
            bounds=(0, c_imag_max),
            method='bounded',
            options={'xatol': 1e-12, 'maxiter': 500}
        )
        c_imag = result.x
    except Exception:
        c_imag = 0.0 # Fallback
    
    # Final bounds check
    c_imag = max(0.0, c_imag)
    
    # Exact formula for attenuation: α = ω * c_i / (c_r² + c_i²)
    alpha = omega * c_imag / (c_real**2 + c_imag**2)
    
    return alpha, c_imag


# ================================================================
# 8. Build models (from Yuan et al. 2024)
# ================================================================
def build_model(model_id=1):
    # (Leaving original model definitions, but focusing on model_id=4)
    if model_id == 1:
        # Table 4 of the paper
        layers = [
            Layer(5, 2000, 180, 18),
            Layer(5, 2000, 300, 30),
            Layer(5, 2000, 420, 42) ]
        halfspace = Layer(np.inf, 2000, 500, 50)
    # ... other models
    else: # model_id = 4 (The one used in __main__)
        layers = [
            Layer(3, 1820, 190, 19),
            Layer(4, 1860, 320, 32),
            Layer(5, 1910, 280, 28),
            Layer(5, 1950, 460, 46),
            Layer(6, 2000, 630, 63) ]
        halfspace = Layer(np.inf, 2100, 750, 75)

    return layers, halfspace


# ================================================================
# 9. Main computation with fixed GMB parameters
# ================================================================
def run_dispersion(model_id=1, n_modes=2, freq_min=0.01, freq_max=50.0, 
                   n_freq=100, apply_smoothing=True):
    """
    Run dispersion analysis with MODIFIED GMB parameters.
    """
    layers, halfspace = build_model(model_id)
    
    # <--- FIX APPLIED HERE: MODIFIED GMB PARAMETERS --->
    # The original L=3 parameters caused low-frequency instability.
    # Increasing L and widening the log-spaced frequency spectrum 
    # provides a much smoother approximation of constant-Q behavior.
    L_new = 9
#     f_l = np.logspace(np.log10(1e-3), np.log10(1e3), L_new) # 10^-3 Hz to 10^3 Hz
#     a_l = [1.0 / L_new] * L_new
#     gmb = GMBParameters(L=L_new, weights=a_l, relax_freqs=f_l.tolist())
    # <--- END FIX --->
    
    
    L     = L_new
    f_min = 0.001
    f_max = 100

    gmb = GMBParameters(
        L=L,
        weights=[1 / L] * L,
        relax_freqs=np.logspace(np.log10(f_min), np.log10(f_max), L).tolist()
    )
    
    # Ensure freq_min > 0
    if freq_min < 1e-6:
        freq_min = 1e-6
    
    # Frequency range
    freq_range  = np.linspace(freq_min, freq_max, n_freq)
    phase_vel   = np.full((n_modes, len(freq_range)), np.nan)
    attenuation = np.full_like(phase_vel, np.nan)
    
    # Track previous c_imag for continuity
    c_imag_prev = np.full(n_modes, np.nan)

    print("=" * 75)
    print(f"Computing Love-wave dispersion for Model {model_id} with fixed GMB (L={L_new})...")
    print("=" * 75)
    
    for i, f in enumerate(freq_range):
        omega = 2 * np.pi * f
        roots = find_roots_dispersion(omega, layers, halfspace, gmb, n_modes)
        
        # Sort roots to ensure mode order is preserved (Mode 0, Mode 1, ...)
        roots = np.sort(roots)
        
        for m, c in enumerate(roots):
            phase_vel[m, i] = c
            
            # Mode tracking: Use previous c_imag for continuity.
            prev_val = c_imag_prev[m] if np.isfinite(c_imag_prev[m]) else None
            
            attenuation[m, i], c_imag_prev[m] = compute_attenuation_complex(
                omega, c, layers, halfspace, gmb, c_imag_prev=prev_val
            )
                
        if i % 10 == 0:
            print(f"f = {f:6.2f} Hz → {len(roots)} mode(s)")
            if len(roots) > 0:
                for m in range(len(roots)):
                    if np.isfinite(attenuation[m, i]):
                        print(f"  M{m}: c={roots[m]:6.1f} m/s, "
                              f"α={attenuation[m,i]*1e6:6.3f}×10⁻⁶ 1/m")

    print("=" * 75)
       
    return freq_range, phase_vel, attenuation, layers, halfspace, gmb


# ================================================================
# 10. Plotting section
# ================================================================
def plot_results(freq_range, phase_vel, attenuation, layers, halfspace, gmb, model_id):
    n_modes = phase_vel.shape[0]
    colors  = plt.cm.viridis(np.linspace(0, 1, n_modes))

    # ===========================
    # (a) Phase velocity
    # ===========================
    plt.figure(figsize=(7, 5), dpi=600) 

    for m in range(n_modes):
        plt.plot(freq_range, phase_vel[m, :], color=colors[m], lw=1.3, label=f"M{m}")

    plt.axhline(layers[0].beta, color='black', ls=':', alpha=0.7, lw=1., 
                label=f"$\\beta_1$ = {layers[0].beta:.0f} m/s")
    plt.axhline(halfspace.beta, color='black', ls='--', alpha=0.7, lw=1., 
                label=f"$\\beta_6$ = {halfspace.beta:.0f} m/s")

    plt.xlabel("Frequency (Hz)", fontsize=11)
    plt.ylabel("Phase velocity (m/s)", fontsize=11)
    plt.xlim(freq_range.min(), freq_range.max())
    plt.ylim(100,800)
    plt.title(f"Model {model_id}: Love Wave Dispersion (Fixed GMB)", fontsize=11, fontweight="bold")
    plt.grid(True, ls=':', alpha=0.3)
    plt.legend(fontsize='x-small', loc='upper left', columnspacing=0.2, handletextpad=0.4, ncol=min(n_modes, 4), framealpha=0.9, edgecolor='black')
    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_PhaseVelocity.png", dpi=1200, bbox_inches='tight')
    plt.show()

    # ===========================
    # (b) Attenuation 
    # ===========================
    plt.figure(figsize=(8, 5), dpi=600)

    for m in range(n_modes):
        plt.plot(freq_range, attenuation[m, :], color=colors[m], lw=1.3, label=f"M{m}")

    plt.xlabel("Frequency (Hz)", fontsize=11)
    plt.ylabel("Attenuation coefficient $\\alpha$ (1/m)", fontsize=11)
    plt.xlim(freq_range.min(), freq_range.max())
    plt.title(f"Model {model_id}: Love Wave Attenuation (Fixed GMB)", fontsize=11, fontweight="bold")
    plt.grid(True, ls=':', alpha=0.3)
    plt.legend(fontsize='x-small', loc='best', framealpha=0.9, edgecolor='black', columnspacing=0.2)
    plt.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_Attenuation.png", dpi=1200, bbox_inches='tight')
    plt.show() 
    
    # ===========================
    # (c) 3D plot: α–c–f 
    # ===========================
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = plt.figure(figsize=(10, 6), dpi=600)
    ax  = fig.add_subplot(111, projection='3d')

    for m in range(n_modes):
        # Filter valid data
        mask = np.isfinite(phase_vel[m, :]) & np.isfinite(attenuation[m, :])
        ax.plot(freq_range[mask],phase_vel[m, mask], attenuation[m, mask], 
                color=colors[m], lw=1.3, label=f"M{m}")
        
    ax.set_ylabel("Phase Velocity (m/s)", fontsize=10)
    ax.set_xlabel("Frequency (Hz)", fontsize=10)
    ax.set_zlabel("Attenuation $\\alpha$ (1/m)", fontsize=10)
    ax.set_title(f"Model {model_id}: 3D α–c–f Relationship", fontsize=11, fontweight="bold")
    ax.set_xlim(80, 0)
    ax.set_ylim(800, 100)
    plt.tight_layout()
    plt.savefig(f"Model_{model_id}_3D_α–c–f.png", dpi=1200, bbox_inches='tight')
    plt.show()


# ================================================================
# Run all
# ================================================================
if __name__ == "__main__":
    print("\n" + "="*90)
    print("LOVE WAVE ANALYSIS - GMB MODEL FIX APPLIED")
    print("="*90 + "\n")
    
    # Configuration
    freq_min = 0.01    # Start from small positive value
    freq_max = 80.0
    n_freqs  = 200     # High resolution for smooth curves
    model_id = 4
    n_modes  = 9

    # Run dispersion solver with FIXED GMB parameters
    freq_range, phase_vel, attenuation, layers, halfspace, gmb = run_dispersion(
        model_id=model_id, 
        n_modes=n_modes,
        freq_min=freq_min, 
        freq_max=freq_max, 
        n_freq=n_freqs,
        apply_smoothing=True
    )
    
    # Plot results
    plot_results(freq_range, phase_vel, attenuation, layers, halfspace, gmb, model_id)
    
    print("\n" + "="*75)
    print("✓ Analysis complete with fixed GMB parameters for smooth low-frequency attenuation.")
    print("="*75)
    


