---
title: "Goal-Space Planning with Subgoal Models" 
venue: "JMLR 2024"
date: 2024-10-22
lastmod: 2024-10-28
tags: ["Model-Based Reinforcement Learning","Temporal Abstraction","Planning"]
author: ["Chunlok Lo*","Kevin Roice*", "Parham Mohammad Panahi*", "Scott M. Jordan", "Adam White", "Gabor Michuz", "Farzane Aminmansour", "Martha White"]
description: "" 
summary: "We propose a long-horizon background-planning algorithm for online RL. This used subgoal
models (abstract in state & time) for faster long-term decision making & smarter value propagation." 
cover:
    image: "paper1.png"
    alt: "Abstract Models for Planning"
    relative: false
editPost:
    URL: "http://jmlr.org/papers/v25/24-0040.html"
    Text: "JMLR 2024"

---

---

##### Download

+ [Paper](paper1.pdf)
<!-- + [Code and data](https://github.com/pmichaillat/feru) -->

---

##### Abstract

Geophysical DC resistivity imaging is crucial in subsurface exploration, environmental studies, and resource assessment. However, traditional inversion techniques face challenges in accurately resolving complex subsurface features because of the inherent nonlinearities in geophysical data. To overcome these challenges, we propose a hybrid optimization approach that combines incomprehensible but intelligible-in-time (IbI) logic with the interior point method (IPM). The IbI logic framework leverages complexity and temporal intelligibility, allowing for a dynamic interpretation of subsurface phenomena. By integrating IbI with IPM, our approach benefits from global exploration and local refinement, leading to improved convergence speed and solution quality. The objective function formulated for the inversion process includes data misfit and model regularization terms, which promote accurate and smooth solutions. Our methodology involves the use of the IbI logic algorithm (ILA) for initial global search, which identifies promising regions in the search space. This is followed by the application of IPM for local optimization. This synergy between the two algorithms ensures robustness and efficiency in handling large datasets and complex geological models. We conducted tests using synthetic and real DC resistivity data to validate our approach. The synthetic test demonstrated the accurate reconstruction of subsurface anomalies, whereas the real data test successfully identified fault zones, which is consistent with previous studies. The hybrid optimization algorithm significantly improves the resolution of subsurface structures and enhances geophysical data inversion practices. It balances the exploration and refinement phases effectively, optimizing the computation time and ensuring precise model delineation.

---

##### Figure 5: Display and comparison of the (a) synthetic positive anomaly DC resistivity model, (b) positive anomaly inverted resistivity model, (c) synthetic negative anomaly DC resistivity model, and (d) negative anomaly inverted resistivity model.

![](paper1.png)

---

##### Citation

Edigbue, Paul, Hammed Oyekan, Abdullatif Al-Shuhail, and Sherif Hanafy. "Hybrid optimization for DC resistivity imaging via intelligible-in-time logic and the interior point method." *Scientific Reports 14*, no. 1 (2024): 27558.

```BibTeX
@article{edigbue2024hybrid,
  title={Hybrid optimization for DC resistivity imaging via intelligible-in-time logic and the interior point method},
  author={Edigbue, Paul and Oyekan, Hammed and Al-Shuhail, Abdullatif and Hanafy, Sherif},
  journal={Scientific Reports},
  volume={14},
  number={1},
  pages={27558},
  year={2024},
  publisher={Nature Publishing Group UK London}
}
```

---

<!-- ##### Related material

+ [Presentation slides](presentation1.pdf)
+ [Summary of the paper](https://www.penguinrandomhouse.com/books/110403/unusual-uses-for-olive-oil-by-alexander-mccall-smith/) -->
