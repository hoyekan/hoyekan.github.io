---
title: "MarchenCode"

weight: 2

date: 2025-04-10

tags: ["Marchenko Equation", "Green's Functions", "Upgoing Green Function", "Downgoing Green Function", "Focusing Functions"]

author: ["Oyekan Hammed"]

description: "MarchenCode"

summary: "MarchenCode contains GPU-optimized Python code to compute the Green's functions from single-sided surface reflectivity data. It also contains code to perform Marchenko-based imaging using decomposed Green’s functions (upgoing)." 

cover:
    image: "marchenko.png"
    alt: "Marchenko Imaging"
    relative: true

# showToc: true
---

---

##### Download <br>

The data used in the code below can be downloaded [here](https://drive.google.com/drive/folders/1SOQUvALC9upk8hrCzLSN_JkUvkq5ezh_?usp=sharing). The directory has two sub-folders:  
- `DAT` folder containing files with `.dat` extensions  
- `MAT` folder containing files with `.mat` extensions  

The code below uses data from the `MAT` folder.

---

<!--+ [Green Function Estimation 1](/projects/project3/Marchenko.ipynb)-->
<!--+ [Green Function Estimation 1](https://raw.githubusercontent.com/hoyekan/hoyekan.github.io/main/content/projects/project3/Marchenko.ipynb)-->
+ [Green Function Estimation 1](https://github.com/hoyekan/hoyekan.github.io/blob/main/content/projects/project3/marchenko_gpu_optimized.ipynb)
+ [Green Function Estimation 1](https://github.com/hoyekan/hoyekan.github.io/blob/main/content/projects/project3/marchenko_imaging_gpu_optimized.ipynb)
<!--+ [Green Function Estimation 2](/projects/project3/marchenko_modular.py)-->
<!--+ [Marchenko Imaging](/projects/project3/marchenko_imaging.py)-->

---

##### Abstract

*MarchenCode* contains GPU-optimized Python code to compute the Green's functions from single-sided surface reflectivity data. Learners will also find code to perform Marchenko-based imaging using decomposed Green’s functions (upgoing). I embarked on this project in an effort to understand the Marchenko methods. For learners proficient in MATLAB, you can download the original version of the code in MATLAB [here](https://wiki.seg.org/wiki/Software:Marchenko_for_imaging) by [`Angus Lomas and Andrew Curtis`](https://www.geos.ed.ac.uk/~acurtis/assets/Lomas_Curtis_Geop_March_2019.pdf)

---

##### Figure 1: Input Shot Gather

![](marchenko1.png)

##### Figure 2: Focusing Functions (Initial and After 5 Iterations)

![](marchenko2.png)

##### Figure 3: Estimated Green's Function

![](marchenko3.png)

##### Figure 4: Marchenko Imaging

![](marchenko4.png)

A look at the output from the tqdm (**taqaddum**, which mean progress in arabic) progress bar `64/64 [23:06<00:00, 22.04s/it]` indicates:

1. There are 64 iterations..

2. **Total Runtime**: 0 hours, 23 minutes, and 06 seconds.

3. **Per-Iteration Time**: ~22.04 seconds per point.

<!--The progress bar reveals the computation is **`super slow`**. It should be noted that the imaging grid does not include all the full coordinates. The imaging grid ranges from 400m - 2500m (x position) and 200m - 1600m (z position) at 16m spacing. So, interested learners should consider improving the performance of the code significantly.-->


##### Key Reference

`Lomas, Angus, and Andrew Curtis`. "***An introduction to Marchenko methods for imaging.***" *Geophysics* 84, no. 2 (2019): F35-F45.

---
