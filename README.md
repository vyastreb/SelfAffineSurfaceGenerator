# Generator of a periodic random self-affine field

Generate **periodic** 1D/2D/3D self‑affine Gaussian random fields with a prescribed power‑law spectrum in Fourier space. Includes two variants: a "noisy" spectrum (filtered white noise) and an **exact‑magnitude** spectrum with random phases only. C++ and Python implementations.

---

## Project information

+ Author: Vladislav A. Yastrebov
+ Affiliation: CNRS, Mines Paris, Centre des Materiaux, Versailles, Paris
+ Date: 2015-2025
+ Licence: BSD 3-Clause 

---

## Description & references

This repository implements (in C++ and Python) the Fourier‑space filtering approach introduced by:

- **Hu, Y.Z.; Tonder, K.** Simulation of 3‑D random rough surface by 2‑D digital filter in Fourier space. *Int. J. Mach. Tools Manufact.* **32** (1992) 83–90. DOI: [10.1016/0890-6955(92)90064-N](https://doi.org/10.1016/0890-6955%2892%2990064-N)

It has been used, e.g., in:

- **Yastrebov, V.A.; Anciaux, G.; Molinari, J.F.** The role of the roughness spectral breadth in elastic contact of rough surfaces. *J. Mech. Phys. Solids* **107** (2017) 469–493. DOI: [10.1016/j.jmps.2017.07.016](https://doi.org/10.1016/j.jmps.2017.07.016), [arXiv:1704.05650](https://arxiv.org/abs/1704.05650)

**Python** also includes a generator with an *idealized spectrum*: Fourier magnitudes match the target exactly; only phases are random.

---

## Features

- Periodic Gaussian random fields in **1D/2D/3D**
- Self‑affine spectrum with Hurst exponent 
- **Two methods**:
  - *Filtered noise:* spectrum follows the target on average
  - *Prescribed magnitudes:* exact magnitudes, random phases (real field guaranteed via Hermitian symmetry)
- Configurable spectral band: $k_{\text{low}}, k_{\text{high}}$
- Optional **PSD plateau** for $0 \le k < k_{\text{low}}$
- C++ implementation and pure‑NumPy Python implementation



## C++ usage

### Build

```bash 
make # build src/* into bin/SURFACE_GENERATOR
```

### Run

```bash
./bin/SURFACE_GENERATOR k1 k2 H N seed rms Npdf plateau
```
if you do not provide arguments, it will print the arguments it needs.

**Arguments**

+ [1] `k1 (int)` -- lower cutoff wavenumber  
+ [2] `k2 (int)` -- upper cutoff wavenumber
+ [3] `H (double)` -- Hurst exponent, $0 < H < 1$
+ [4] `L (int)` -- number of points per size (prefer powers of 2, i.e. 128, 256, 512, etc.)
+ [5] `s (int)` -- seed for the random number generator 
+ [6] `rms (double)` -- standard deviation of heights
+ [7] `Npdf (int)` -- number of bins in pdf data
+ [8] `if_plateau (bool)` -- boolean argument determining whether the power spectral density has a plateau up to k1? (0 - non, 1 -yes)

## Python usage

The Python code lives in `python/RandomField.py` with a simple test in `python/test.py`.
Import the module (add the folder to `PYTHONPATH` or place it next to your script):

```python
import numpy as np
import RandomField as rf

np.random.seed(42)
N = 1024

random_field = rf.periodic_gaussian_random_field(dim = 2, N = N, Hurst = 0.5, k_low = 4 / N, k_high = 128 / N)

ideal_random_field = rf.ideal_periodic_gaussian_random_field(dim = 2, N = N, Hurst = 0.5, k_low = 4 / N, k_high = 128 / N)

ideal_random_field_with_plateau = rf.ideal_periodic_gaussian_random_field(dim = 2, N = N, Hurst = 0.5, k_low = 4 / N, k_high = 128 / N, plateau = True)
```

## API(Python)

```
periodic_gaussian_random_field(dim, N, Hurst, k_low, k_high, plateau=False, verbose=False)
    -> ndarray
  - Filtered white noise. The radially averaged PSD follows the target power law.

ideal_periodic_gaussian_random_field(dim, N, Hurst, k_low, k_high, plateau=False, verbose=False)
    -> ndarray
  - Exact Fourier magnitudes per bin (target spectrum), random phases only.
```

**Parameters:**

+ `dim` (int): Dimension of the field (1, 2, or 3).
+ `N` (int): number of points per side in the field (N x N for 2D, N x N x N for 3D).
+ `Hurst` (float): Hurst exponent, in the range [0, 1].
+ `k_low` (float): Lower bound of the wavenumber range (k_low > 0), given in $k1/L$.
+ `k_high` (float): Upper bound of the wavenumber range (k_high < 0.5, which represents Nyquist frequency), given in $k2/L$.
+ `plateau` (bool): If True, the power spectrum is flat up to k_low, otherwise all wavenumbers below k_low are set to zero.
+ `verbose` (bool): If True, print the parameters used for generating the random field.


## License

+ C++ repository: BSD 3‑Clause.
+ Python module: CC0 





