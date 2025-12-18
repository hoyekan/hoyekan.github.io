---
title: "Projects"
date: 2025-06-01
summary: "My software and past projects."

_build:
  list: false
---

<style>
/* Container for the projects */
.projects-grid {
  display: grid;
  grid-template-columns: repeat(2, 1fr); /* 2 columns */
  gap: 2rem;
  margin-top: 2rem;
  width: 100%;
  grid-auto-rows: auto; /* adjust row height automatically */
}

/* Responsive: 1 column on smaller screens */
@media (max-width: 900px) {
  .projects-grid {
    grid-template-columns: 1fr;
  }
}

/* Individual project card */
.project-card {
  border: 1px solid #ddd;
  border-radius: 12px;
  padding: 1.5rem;
  background: #ffffff;
  width: 100%;
  box-sizing: border-box;
  display: flex;
  flex-direction: column;
  justify-content: space-between;
}

/* Project images */
.project-card img {
  width: 100%;
  border-radius: 8px;
  margin-bottom: 1rem;
  object-fit: cover;
}

/* Project title */
.project-card h3 {
  margin: 0.5rem 0;
}

/* Project links */
.project-card a {
  margin-top: 1rem;
  color: #0070f3;
  text-decoration: none;
  font-weight: bold;
}

.project-card a:hover {
  text-decoration: underline;
}
</style>

## Projects

Here are some of the projects I have embarked on.

<div class="projects-grid">

  <div class="project-card">
    <img src="/images/Dispersion.png" alt="Love Waves Dispersion">
    <h3>Dispersion and Attenuation of Love Waves in an Isotropic Viscoelastic Layers over Half-Space</h3>
    <p>
      Dispersion and Attenuation of Love Waves in a Stack of N Isotropic
      Viscoelastic Layers over a Half-Space using the
      Thomson–Haskell propagator matrix approach.
    </p>
    <p><strong>December, 2025 · Oyekan Hammed</strong></p>
    <a href="/projects/project4/">View project →</a>
  </div>

  <div class="project-card">
    <img src="/images/marchenko.png" alt="Marchenko Imaging">
    <h3>MarchenCode</h3>
    <p>
      MarchenCode contains Python code to compute the Green's functions from single-sided surface reflectivity data. 
      It also contains code to perform Marchenko-based imaging using decomposed Green’s functions (upgoing).
    </p>
    <p><strong>April, 2025 · Oyekan Hammed</strong></p>
    <a href="/projects/project3/">View project →</a>
  </div>

  <div class="project-card">
    <img src="/images/acoustic.png" alt="Acoustic Wave Modelling">
    <h3>2D Acoustic Wave Modelling</h3>
    <p>
      Modelling of 2D acoustic wave equation in a homogeneous model and two-layers velocity model.
    </p>
    <p><strong>February, 2025 · Oyekan Hammed</strong></p>
    <a href="/projects/project2/">View project →</a>
  </div>

  <div class="project-card">
    <img src="/images/groundwater.png" alt="Equipotential and Streamlines around Wells">
    <h3>Modelling groundwater flow around wells in a confined, unconfined, and combined aquifer systems using discharge potential theory.</h3>
    <p>
      A Mathematica-based study of groundwater flow around wells using analytical discharge potential formulations 
      for steady-state and transient conditions in both isotropic and anisotropic aquifers.
    </p>
    <p><strong>May, 2025 · Oyekan Hammed</strong></p>
    <a href="/projects/project1/">View project →</a>
  </div>

</div>
