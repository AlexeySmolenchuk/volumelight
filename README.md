# volumelight
## Overview
Fog shader for Houdini Mantra based on the first part of [this](https://www.solidangle.com/research/egsr2012_volume.pdf) paper.
![beauty](https://github.com/somesanctus/volumelight/blob/master/img/beauty.jpg)

It works similar to **Lit Fog** sahder. Create **Atmosphere** object an assign **volumelight** shader to it.
Shader works only for primary rays, it simple, adjustable, animatable and has a simple noise to modulate albedo (not density).
![noise](https://github.com/somesanctus/volumelight/blob/master/img/noise.jpg)

Unlike **Lit Fog** it physicly correct, takes into account [Beer–Lambert](https://en.wikipedia.org/wiki/Beer%E2%80%93Lambert_law) law. It use modern algorithms such as **Exponential sampling** (works well in most cases, but noisy near the lights) *(Time in min:sec for equal numer of samples.)* 
![exponential](https://github.com/somesanctus/volumelight/blob/master/img/exponential.jpg) *1:32*

**Equiangular sampling** (well in lights area, but noisy in shadows)
![equiangular](https://github.com/somesanctus/volumelight/blob/master/img/equiangular.jpg) *2:22*

**Multiple Importance Sampling** (select the preferable sampling stategy [Eric Veach, Leonidas J. Guibas](https://graphics.stanford.edu/courses/cs348b-03/papers/veach-chapter9.pdf). I use the power heuristic)
![mis](https://github.com/somesanctus/volumelight/blob/master/img/mis.jpg) *1:56*

Heat map of distributions green - exponential, red - equiangular.
![mis_weights](https://github.com/somesanctus/volumelight/blob/master/img/mis_weights.jpg)

The beauty channel unclamped and alpha channel not depend on beauty intensity.
![vol_opacity](https://github.com/somesanctus/volumelight/blob/master/img/vol_opacity.jpg)

You can simply export distributions on per-light basis.

![spotlight](https://github.com/somesanctus/volumelight/blob/master/img/spotlight.jpg)
![lead](https://github.com/somesanctus/volumelight/blob/master/img/lead.jpg)
![lightin](https://github.com/somesanctus/volumelight/blob/master/img/lightin.jpg)

For better samplig set Max Distance carefuly. I recommend to render **volumelight** separately by turning objects to matte shading. You can set **Color** parameter grater then 1, although it's not physicaly correct, sometimes help to achive needed visual effect.
## Limitations
Work well only with decaying lights (excluding **Environment**, **Sun**, **Distant**)
