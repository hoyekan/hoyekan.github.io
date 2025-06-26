---
title: "Modeling Groundwater Flow Around Wells with Discharge Potential"

venue: "For GEOP 629 in Term 242"

weight: 3

date: 2025-05-05

tags: ["Discharge Potential", "Groundwater Modelling", "Aquifer", "Groundwater Flow", "Confined Aquifer", "Unconfined Aquifer", "Pumping Wells", "Steady Flow", "Transient (Unsteady) Flow", "Heterogeneous", "Anisotropic Flow", "Theis Equation", "Wolfram Mathematica"]

author: ["Oyekan Hammed"]

description: "Analytical and numerical exploration of groundwater flow in confined, unconfined, and combined aquifer systems using discharge potential theory."

summary: "A Mathematica-based study of groundwater flow around wells using analytical discharge potential formulations for steady-state and transient conditions in both isotropic and anisotropic aquifers."

cover:
    image: "project.png"
    alt: "Equipotential and streamlines around wells"
    relative: true

---

---

##### Download

+ [Project](final_project.pdf)
+ [Code](/projects/project1/final_project.nb)

---

##### Abstract

This project involves modeling groundwater flow around pumping wells in confined and unconfined aquifers using discharge potential theory. Building on the analytical framework presented by [Korkmaz (2017)](https://www.ewra.net/ew/pdf/EW_2017_57_52.pdf), the study focuses on steady-state flow in confined, unconfined, and combined aquifer systems by leveraging the linearity of Laplace’s equation and the principle of superposition. The analytical solutions are extended to heterogeneous and anisotropic aquifers through coordinate transformations. Additionally, the classical Theis equation [(Theis, C.V. 1935)](https://doi.org/10.1029/TR016i002p00519) is employed to model transient (time-dependent) flow around two pumping wells in an unconfined aquifer system.

The implementation begins with a 1D solution for groundwater flow using discharge potential in both confined and unconfined zones, followed by a 2D radial flow model around a single well. Discharge potential and stream function solutions for a single well are generalized to two-well systems using superposition for both homogeneous isotropic and heterogeneous anisotropic aquifers. The final segment of the code implements transient flow modeling to compute the evolution of hydraulic head and discharge potential over time.

---

##### Figure 1: Flownet around a single well (Discharge potential, Stream function, and Hydraulic head)

![](pic1.png)

##### Figure 2: Contour plot of equipotentials and streamlines for 2 wells in an unconfined aquifer (Homogeneous and Isotropic)

![](pic2.png)

##### Figure 3: Contour plot of equipotentials and streamlines for 2 wells in an unconfined aquifer with a uniform flow (Homogeneous and Isotropic)

![](pic3.png)

##### Figure 4: Flow Net: Discharge and Stream Functions (Heterogeneous and Anisotropic unconfined aquifer)

![](pic4.png)

##### Figure 5: Hydraulic Head distribution (Heterogeneous and Anisotropic unconfined aquifer)

![](pic5.png)

---

<!-- ##### Related material

+ [Presentation slides](presentation2.pdf)
+ [Wikipedia entry](https://en.wikipedia.org/wiki/The_Finer_Points_of_Sausage_Dogs) -->
