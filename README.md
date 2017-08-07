# volumelight
## Overview
Fog shader for Houdini Mantra based on the first part of [this](https://www.solidangle.com/research/egsr2012_volume.pdf) paper.

![beauty](https://github.com/somesanctus/volumelight/blob/master/img/beauty.jpg)
![tunnel](https://github.com/somesanctus/volumelight/blob/master/img/tunnel.jpg)

It works similar to **Lit Fog** shader. Create **Atmosphere** object and assign **volumelight** shader to it.
Shader works only for primary rays, it simple, adjustable, animatable and has a built-in noise to modulate albedo (not density).

![noise](https://github.com/somesanctus/volumelight/blob/master/img/noise.jpg)

Unlike **Lit Fog** it physically correct and takes into account [Beer–Lambert](https://en.wikipedia.org/wiki/Beer%E2%80%93Lambert_law) law. It uses modern algorithms such as **Exponential sampling** (works well in most cases, but noisy near the lights) **Equiangular sampling** (well in areas around lights, but noisy in shadows)

![compare](https://github.com/somesanctus/volumelight/blob/master/img/compare.jpg)

**Multiple Importance Sampling** (select the preferable sampling stategy [Eric Veach, Leonidas J. Guibas](https://graphics.stanford.edu/courses/cs348b-03/papers/veach-chapter9.pdf). I use the power heuristic). Heat map of distributions (green - exponential, red - equiangular).

![mis_weights](https://github.com/somesanctus/volumelight/blob/master/img/mis_weights.jpg)

The beauty channel unclamped and alpha channel doesn’t depend on beauty intensity, only density.
You can simply export distributions on per-light basis.

![auto_passes](https://github.com/somesanctus/volumelight/blob/master/img/auto_passes.jpg)
![tunnel_passes](https://github.com/somesanctus/volumelight/blob/master/img/tunnel_passes.jpg)

For better sampling set **Max Distance** carefully. I recommend to render **volumelight** separately by turning objects to matte shading. You can set **Color** parameter greater then 1. Although it’s not physically correct sometimes it help to achive needed visual effect. Use IES profiles for better look.

## Limitations
- Work well only with decaying lights (excluding **Environment**, **Sun**, **Distant**)
- Rendering time grows linearly with number of lights.
