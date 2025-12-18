---
title: "Projects"
date: 2025-06-01
summary: "My software and past projects."
---

{{< rawhtml >}}

<style>
.projects-grid {
  display: grid;
  grid-template-columns: repeat(2, 1fr);
  gap: 2rem;
  margin-top: 2rem;
  width: 100%;
}

@media (max-width: 900px) {
  .projects-grid {
    grid-template-columns: 1fr;
  }
}

.project-card {
  border: 1px solid #ddd;
  border-radius: 12px;
  padding: 1.5rem;
  background: #ffffff;
  box-sizing: border-box;
}

.project-card img {
  width: 100%;
  border-radius: 8px;
  margin-bottom: 1rem;
}
</style>

<h2>Projects</h2>
<p>Here are some of the projects I have embarked on.</p>

<div class="projects-grid">

  <div class="project-card">
    <img src="/images/Dispersion.png" alt="Love Waves Dispersion">
    <h3>Dispersion and Attenuation of Love Waves</h3>
    <p>Thomson–Haskell propagator matrix approach.</p>
    <p><strong>December, 2025 · Oyekan Hammed</strong></p>
    <a href="/projects/project4/">View project →</a>
  </div>

  <div class="project-card">
    <img src="/images/marchenko.png" alt="Marchenko Imaging">
    <h3>MarchenCode</h3>
    <p>Green’s functions and Marchenko-based imaging.</p>
    <p><strong>April, 2025 · Oyekan Hammed</strong></p>
    <a href="/projects/project3/">View project →</a>
  </div>

  <div class="project-card">
    <img src="/images/acoustic.png" alt="Acoustic Wave Modelling">
    <h3>2D Acoustic Wave Modelling</h3>
    <p>Homogeneous and two-layer velocity models.</p>
    <p><strong>February, 2025 · Oyekan Hammed</strong></p>
    <a href="/projects/project2/">View project →</a>
  </div>

  <div class="project-card">
    <img src="/images/groundwater.png" alt="Groundwater Flow">
    <h3>Groundwater Flow Around Wells</h3>
    <p>Analytical discharge potential theory.</p>
    <p><strong>May, 2025 · Oyekan Hammed</strong></p>
    <a href="/projects/project1/">View project →</a>
  </div>

</div>

{{< /rawhtml >}}
