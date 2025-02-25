# Generator of a periodic random self-affine field

## Information

+ Author: Vladislav A. Yastrebov
+ Affiliation: MINES Paris, Centre des Materiaux, CNRS
+ Date: 2015-2023
+ Licence: BSD 3-Clause 

## Description

It is a C++ implementation of the method proposed in
Y.Z. Hu and K. Tonder, "Simulation of 3-D random rough surface by 2-D digital filter in fourier space", Int. J. Mach. Tools Manufact. 32:83-90 (1992): [DOI](https://doi.org/10.1016/0890-6955%2892%2990064-N)
One could also check our paper where this method was used to generate rough surfaces for contact analysis between rough surfaces: V.A. Yastrebov, G. Anciaux, J.F. Molinari.
"The role of the roughness spectral breadth in elastic contact of rough surfaces", Journal of the Mechanics and Physics of Solids, 107:469-493 (2017): [DOI](https://doi.org/10.1016/j.jmps.2017.07.016), [arXiv](https://arxiv.org/abs/1704.05650)

## Input: it requires 8 arguments: 

+ [1] lower cutoff wavenumber k1 (int)
+ [2] upper cutoff wavenumber k2 (int)
+ [3] Hurst exponent, H (double) // 0 < H < 1
+ [4] number of points per size,  L (int) // powers of 2, i.e. 128, 256, 512, etc
+ [5] seed for the random number generator, s (int)
+ [6] standard deviation of heights, rms (double)
+ [7] number of bins in pdf data, Npdf (int)
+ [8] if the power spectral density has a plateau up to k1? (0 - non, 1 -yes)


## Short help:

+ to compile `src/*`<br>
`$ make`
+ the executable SURFACE_GENERATOR is put in `bin/`
to execute run<br> 
`$./SURFACE_GENERATOR` <br>
if you do not provide arguments, it will print the arguments it needs.

## Python implementation

Python implementation is available under CC0 license, you can find it in `python/RandomField.py` and a test in `python/test.py`.

Usage:
```python
import numpy as np
import RandomField as rf

N0 = 1024
random_field = rf.periodic_gaussian_random_field(dim = 2, N = N0, Hurst = 0.5, k_low = 4 / N0, k_high = 128 / N0)
```


