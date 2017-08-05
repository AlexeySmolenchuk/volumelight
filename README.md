# volumelight
## Overview
Fog shader for Houdini Mantra based on the first part of [this](https://www.solidangle.com/research/egsr2012_volume.pdf) paper.

It works similar to Lit Fog sahder. Create **Atmosphere** object an assign **volumelight** shader to it.
Shader works only for primary rays, it simple, adjustable, animatable and has a simple noise to modulate albedo (not density).
![beauty.jpg](https://github.com/somesanctus/volumelight/blob/master/img/beauty.jpg)
Unlike **Lit Fog** it physicly correct, takes into account [Beer–Lambert](https://en.wikipedia.org/wiki/Beer%E2%80%93Lambert_law) law. It use modern algorithms such as Exponential sampling (works well in most cases, but noisy near the lights)

**Equiangular sampling** (well in lights area, but noisy in shadows)

**Multiple Importance Sampling** (select the preferable sampling stategy [Eric Veach, Leonidas J. Guibas](https://graphics.stanford.edu/courses/cs348b-03/papers/veach-chapter9.pdf). I use the power heuristic). Heat map of distributions green - exponential, red - equiangular.

The beauty channel unclamped and alpha channel not depend on beauty intensity.

You can simply export distributions on per-light basis.



## Limitations
Work well only with decaying lights (excluding Environment, Sun, Distant)