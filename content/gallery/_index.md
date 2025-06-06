---
title: "Gallery"
author: ["My Gallery"]
# date: 2025-06-01
# hidemeta: true  # Hides date/author if enabled by your theme
type: "page"
layout: "gallery"
description: "My photo album"

# showToc: true
# disableAnchoredHeadings: false
---

---

# Gallery

A collection of photos and video highlights...


<style>
.gallery-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
  gap: 15px;
  margin-top: 20px;
}

.gallery-grid img {
  width: 100%; /* Full container width */
  height: auto; /* Maintain aspect ratio */
  border-radius: 8px; /* Rounded corners */
  box-shadow: 0 2px 6px rgba(0,0,0,0.15); /* Subtle shadow */
  transition: transform 0.2s ease-in-out; /* Smooth hover animation */
}

.gallery-grid img:hover {
  transform: scale(1.03); /* Slight zoom on hover */
}

.gallery-grid video {
  width: 100%;
  border-radius: 8px; /* Matches image styling */
  box-shadow: 0 2px 6px rgba(0,0,0,0.15); 
}

</style>



<div class="gallery-grid">

  <figure>
    <img src="/images/gallery/pic1.jpeg" alt="Fieldwork 1" />
    <figcaption>Geophone layout during seismic fieldwork.</figcaption>
  </figure>

  <figure>
    <img src="/images/gallery/pic2.jpg" alt="Lab Work 1" />
    <figcaption>Rock sample analysis under microscope.</figcaption>
  </figure>

  <figure>
    <img src="/images/gallery/pic3.JPG" alt="Vibroseis" />
    <figcaption>Vibroseis truck in operation.</figcaption>
  </figure>

  
  <!--<img src="/images/gallery/pic4.jpg" alt="Description 4" />
  <img src="/images/gallery/pic5.jpg" alt="Description 5" />
  <img src="/images/gallery/pic6.jpg" alt="Description 6" />
  <img src="/images/gallery/pic7.jpg" alt="Description 7" />
  <img src="/images/gallery/pic8.jpg" alt="Description 8" />
  <img src="/images/gallery/pic9.jpg" alt="Description 9" />
  <img src="/images/gallery/pic10.jpg" alt="Description 10" />-->

  <!--<video controls>
    <source src="/images/gallery/video1.mp4" type="video/mp4">
    Your browser does not support the video tag.
  </video>-->

</div>

