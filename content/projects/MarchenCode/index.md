---
title: "MarchenCode"

weight: 2

date: 2025-04-10

tags: ["Marchenko Equation", "Green's Functions", "Upgoing Green Function", "Downgoing Green Function", "Focusing Functions"]

author: ["Oyekan Hammed"]

description: "MarchenCode"

summary: "MarchenCode contains Python code to compute the Green's functions from single-sided surface reflectivity data. It also contains code to perform Marchenko-based imaging using decomposed Green’s functions (upgoing)." 

cover:
    image: "marchenko.png"
    alt: "Marchenko Imaging"
    relative: true

# showToc: true
---

---

##### Download

+ [Green Function Estimation 1](/projects/project2/Acoustic_Wave_Modelling.py)
+ [Green Function Estimation 2](/projects/project2/Acoustic_Wave_Modelling.py)
+ [Marchenko Imaging](projects/MarchenCode/Marchenko_Imaging.py)

---

##### Abstract

*MarchenCode* contains Python code to compute the Green's functions from single-sided surface reflectivity data. Learners will also find code to perform Marchenko-based imaging using decomposed Green’s functions (upgoing). I embarked on this project in an effort to understand the Marchenko method. For learners proficient in MATLAB, you can download the original version of the code in MATLAB [here](https://wiki.seg.org/wiki/Software:Marchenko_for_imaging) by [`Angus Lomas and Andrew Curtis`](https://www.geos.ed.ac.uk/~acurtis/assets/Lomas_Curtis_Geop_March_2019.pdf)

---

##### Figure 1: Input Shot Gather

![](marchenko1.png)

##### Figure 2: Focusing Functions (Initial and After 5 Iterations)

![](marchenko2.png)

##### Figure 3: Estimated Green's Function

![](marchenko3.png)

##### Figure 4: Marchenko Imaging

![](marchenko4.png)

A look at the output from the tqdm (**taqaddum**, which mean progress in arabic) progress bar `132/132 [36:28:28<00:00, 994.76s/it]` indicates:

1. There are 132 iterations..

2. **Total Runtime**: 36 hours, 28 minutes, and 28 seconds (extremely slow).

3. **Per-Iteration Time**: ~994.76 seconds (16.6 min) per point, indicating inefficiency.

The progress bar reveals the computation is **`super slow`**. It should be noted that the imaging grid does not include all the full coordinates. The imaging grid ranges from 400m - 2500m (x position) and 200m - 1600m (z position) at 16m spacing. So, interested learners should consider improving the performance of the code significantly.


##### Key Reference

`Lomas, Angus, and Andrew Curtis`. "An introduction to Marchenko methods for imaging." *Geophysics* 84, no. 2 (2019): F35-F45.

---

