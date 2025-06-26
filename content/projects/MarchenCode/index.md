---
title: "MarchenCode"

weight: 2

date: 2025-05-05

tags: ["Marchenko Equation", "Green's Functions", "Upgoing Green Function", "Downgoing Green Function", "Focusing Functions"]

author: ["Oyekan Hammed"]

description: "MarchenCode"

summary: " *MarchenCode* contains Python code to compute the Green's functions from single-sided surface reflectivity data. There's also a code to perform Marchenko-based imaging using decomposed Green’s functions (upgoing). I embarked on this project in an effort to understand the Marchenko method. For learners proficient in MATLAB, you can download the original version of the code in MATLAB [here](https://wiki.seg.org/wiki/Software:Marchenko_for_imaging) by Angus Lomas and Andrew Curtis "


cover:
    image: "marchenko.png"
    alt: "Marchenko Imaging"
    relative: true

---

---

##### Download

+ [Green Function Estimation](/projects/project2/Acoustic_Wave_Modelling.py)
+ [Marchenko Imaging](/projects/project2/2D%20Acoustic%20Wave%20Equation.py) <br>
  The velocity models (`model1.npy` and `model2.npy`) used in the modelling can be downloaded [here](https://github.com/hoyekan/hoyekan.github.io/tree/main/content/projects/project2)

---

##### Abstract

<!-- This project involves modeling groundwater flow around pumping wells in confined and unconfined aquifers using discharge potential theory. Building on the analytical framework presented by [Korkmaz (2017)](https://www.ewra.net/ew/pdf/EW_2017_57_52.pdf), the study focuses on steady-state flow in confined, unconfined, and combined aquifer systems by leveraging the linearity of Laplace’s equation and the principle of superposition. The analytical solutions are extended to heterogeneous and anisotropic aquifers through coordinate transformations. Additionally, the classical Theis equation [(Theis, C.V. 1935)](https://doi.org/10.1029/TR016i002p00519) is employed to model transient (time-dependent) flow around two pumping wells in an unconfined aquifer system.

The implementation begins with a 1D solution for groundwater flow using discharge potential in both confined and unconfined zones, followed by a 2D radial flow model around a single well. Discharge potential and stream function solutions for a single well are generalized to two-well systems using superposition for both homogeneous isotropic and heterogeneous anisotropic aquifers. The final segment of the code implements transient flow modeling to compute the evolution of hydraulic head and discharge potential over time. -->

---

##### Figure 1: Input Shot Gather

![](marchenko1.png)

##### Figure 2: Focusing Functions (Initial and After 5 Iterations)

![](marchenko2.png)

##### Figure 3: Estimated Green's Function

![](marchenko3.png)

##### Figure 4: Marchenko Imaging

![](marchenko4.png)


---

