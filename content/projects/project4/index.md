---
title: "Dispersion and Attenuation of Love Waves in an Isotropic Viscoelastic Layers over a Half-Space"

weight: 1

date: 2025-12-12

tags: ["Love Waves", "Phase Velocity Dispersion Curve", "Attenuation Coefficient", "Half-Space"]

author: ["Oyekan Hammed"]

description: "Dispersion and Attenuation of Love Waves"

summary: "Dispersion and Attenuation of Love Waves in a Stack of N Isotropic Viscoelastic Layers over a Half-Space: A Thomson-Haskell Propagator Matrix Approach."

cover:
    image: "Dispersion.png"
    alt: "Dispersion and Attenuation"
    relative: true

---

---

##### Download

<!-- + [Code 1](/projects/project2/Acoustic_Wave_Modelling.py)
+ [Code 2](/projects/project2/2D%20Acoustic%20Wave%20Equation.py) <br>
  The velocity models (`model1.npy` and `model2.npy`) used in the modelling can be downloaded [here](https://github.com/hoyekan/hoyekan.github.io/tree/main/content/projects/project2) -->

---

##### Abstract

This project presents a comprehensive modelling approach for analyzing Love wave propagation in a stratified medium consisting of $N$ isotropic viscoelastic shear layers overlying a semi-infinite half-space. The theoretical formulation employs the Thomson-Haskell propagator matrix method ([Haskell (1953)](#haskell-1953) and [Thomson (1950)](#thomson-1950)) as implemented in [Chen k. et al. (2025)](#Chenetal2025) to relate state vectors—comprising displacement and shear stress—across layer interfaces, effectively linking the half-space to the free surface.To account for realistic material damping, the study incorporates the Generalized Maxwell Body (GMB) model (Maxwell-Wiechert model), which utilizes a superposition of multiple Maxwell elements to simulate a frequency-independent quality factor ($Q$) over the seismic frequency bandwidth.The project derives the complex dispersion equation by applying traction-free boundary conditions at the surface and radiation conditions in the half-space. Finally, a numerical algorithm is implemented to solve the dispersion function, generating phase velocity and attenuation coefficient curves that characterize the dispersive and dissipative properties of the viscoelastic multilayered system.

The theoretical framework follows the classic work of [Haskell (1953)](#haskell-1953) and [Thomson (1950)](#thomson-1950) as implemented in [Chen k. et al. (2025)](#Chenetal2025). Model 1, as well as model 3 and 4, used in my project are from [Yuan S. et al (2024)](#Yuanetal). The attanuation coefficient is computed via complex velocity approach 

---

##### Figure 1: Phase Velocity Dispersion Curve
![](Model_1_PhaseVelocity.png)

##### Figure 2: Attenuation Coefficient Curve (Using Complex Velocity Method)
![](Model_1_Attenuation.png)

##### Figure 3: Attenuation Coefficient Curve (Simplified Q averaging)
This assumes that attenuation is very small (Q>>1) and its computed using the Kolsky-Futterman model (see attached code for derivation).

For small attenuation ![Q >> 1](https://latex.codecogs.com/svg.latex?Q%20%3E%3E%201), the standard formulation is:

<img src="https://latex.codecogs.com/svg.image?\frac{1}{c^*}\approx\frac{1}{c}\left(1+\frac{i}{2Q}\right)" />

where ![c*](https://latex.codecogs.com/svg.latex?c)* is the real-valued phase velocity at low loss.

Starting from the wavenumber:

<img src="https://latex.codecogs.com/svg.image?k=\frac{\omega}{c^*}\approx\frac{\omega}{c}\left(1+\frac{i}{2Q}\right)=\frac{\omega}{c}+i\frac{\omega}{2cQ}" />

Comparing this to ![k = k_r - i\alpha](https://latex.codecogs.com/svg.latex?k%20%3D%20k_r%20-%20i%5Calpha), we immediately identify:

<img src="https://latex.codecogs.com/svg.image?\boxed{\alpha\approx\frac{\omega}{2cQ}}" />

![](Model_1_Attenuation_Simplified.png)

##### Figure 4: The Love-wave modal solutions in the frequency-phase velocity-attenuation coefficient domain.
![](Model_1_3D.png)

---

##### Key Reference

`Yuan, S., Pan, L., Shi, C., Song, X., & Chen, X.` ***Computation and analysis of surface wave dispersion and attenuation in layered viscoelastic–vertical transversely isotropic media by the generalized R/T coefficient method.*** *Geophysical Journal International* 238, no. 3 (2024): 1505–1529.

`Chen, K., Li, Z., Wang, M., & Sacchi, M. D.` ***Theoretical calculation of dispersion and attenuation curves of deep-guided wave in viscoelastic media.*** *Geophysical Journal International* 243, no. 3 (2025): ggaf393.

`Haskell, N. A.` ***The dispersion of surface waves on multilayered media.*** *Bulletin of the Seismological Society of America* 43 (1953): 17–34.

`Thomson, W. T.` ***Transmission of elastic waves through a stratified solid medium.*** *Journal of Applied Physics* 21, no. 2 (1950): 89–93.



