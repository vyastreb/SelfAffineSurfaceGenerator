# Generator of self-affine roughness
## Vladislav A. Yastrebov
## MINES Paris, Centre des Materiaux, CNRS
## 2015-2023
## Licence: BSD 3-Clause 

## Description

It is a C++ implementation of the method proposed in
Y.Z. Hu and K. Tonder, "Simulation of 3-D random rough surface by 2-D digital filter in fourier space", Int. J. Mach. Tools Manufact. 32:83-90 (1992) DOI: https://doi.org/10.1016/0890-6955(92)90064-N
(see doc/*, I also put our papers where this method was used to generate rough surfaces for contact analysis between rough surfaces)

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

+ to compile src/*
`$ make`
+ the executable SURFACE_GENERATOR is put in bin/
to execute run 
`$./SURFACE_GENERATOR` 
if you do not provide arguments, it will print the arguments it needs.
