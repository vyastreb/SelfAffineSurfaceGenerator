"""
Generate a 1D/2D/3D Gaussian random field.

Method 1: generation from white noise filtered in Fourier space.
The method is described in the following paper:
    Hu, Y.Z. and Tonder, K., 1992. Simulation of 3-D random rough surface by 2-D digital filter and Fourier analysis. International journal of machine tools and manufacture, 32(1-2), pp.83-90.
    DOI: 10.1016/0890-6955(92)90064-N

Method 2: generation with a perfect self-affine spectrum and a random phase.
    
Author: Vladislav Yastrebov, CNRS,  Mines Paris - PSL, Centre des matériaux, Versailles, France 
Help: GPT4, Copilot, Claude Sonnet 4.0
Date: Dec 2023 - Jul 2025
License: CC0
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft2, ifft2, fftfreq, fftshift, ifftshift, fftn, ifftn

def periodic_gaussian_random_field(dim = 2, N = 256, Hurst = .5, k_low = 0.03, k_high = 0.3, plateau = False, verbose = False):
    """
    Generate a periodic Gaussian random field with a self-affine spectrum.
    Randomness is introduced by a random noise generated in the real space and then filtered in Fourier space.

    Parameters:
    dim : int
        Dimension of the field (1, 2, or 3).
    N : int
        Size of the field (N x N for 2D, N x N x N for 3D).
    Hurst : float
        Hurst exponent, in the range [0, 1].
    k_low : float
        Lower bound of the wavenumber range (k_low > 0).
    k_high : float
        Upper bound of the wavenumber range (k_high < 0.5, which represents Nyquist frequency).
    plateau : bool
        If True, the power spectrum is flat up to k_low, otherwise all wavenumbers below k_low are set to zero.
    verbose : bool
        If True, print the parameters used for generating the random field.

    Returns:
    z : ndarray
        The generated random field.
    """
    if k_low < 0 or k_high < 0 or k_low > k_high:
        print("k_low and k_high must be positive and k_low <= k_high")
        exit(1)
    if k_high > 0.5:
        print("k_high must be less than 0.5, which represents Nyquist frequency")
        exit(1)

    if verbose == True:
        print("Random Field Generation parameters:")
        print("/    dim = ", dim)
        print("/    N = ", N)
        print("/    Hurst = ", Hurst)
        print("/    k_low = ", k_low)
        print("/    k_high = ", k_high)

    # Define the power spectrum
    k = fftfreq(N)
    power = -(0.5*dim+Hurst)
    noise = []
    if dim == 1:
        noise = np.fft.fft(np.random.normal(size=(N)))
    if dim == 2:
        k = np.sqrt(k[:,None]**2 + k[None,:]**2)
        noise = np.fft.fft2(np.random.normal(size=(N, N)))
    elif dim == 3:
        k = np.sqrt(k[:,None,None]**2 + k[None,:,None]**2 + k[None,None,:]**2)
        noise = fftn(np.random.normal(size=(N, N, N)))
    else:
        print("Unknown dimension ", dim)
        exit(1)
    k = fftshift(k)

    # First create mask, then compute power spectrum only for valid k values
    mask = (k >= k_low) & (k <= k_high)
    sqrt_power_spectrum = np.zeros_like(k)
    sqrt_power_spectrum[mask] = (k[mask]/k_low)**power
    if plateau:
        sqrt_power_spectrum[k < k_low] = 1.0

    # Make sure the mean power is finite
    sqrt_power_spectrum[np.isinf(sqrt_power_spectrum)] = 0

    # Take the product of the white noise field with the power spectrum (filter)
    grf_fourier = noise * sqrt_power_spectrum

    # Transform back to real space
    if dim == 2:
            z = np.real(ifft2(ifftshift(grf_fourier)))
    else:
        z = np.real(ifftn(ifftshift(grf_fourier)))
    return z


def ideal_periodic_gaussian_random_field(dim = 2, N = 256, Hurst = .5, k_low = 0.03, k_high = 0.3, plateau = False, verbose = False):
    """
    Generate a periodic Gaussian random field with a self-affine spectrum.
    Every "power" in the spectrum follows the self-affine scaling law:
      Phi(kx, ky) = FFT(z)*FFT(z)' = (kx^2 + ky^2)^(-0.5*dim-Hurst)
    but the repartition between Re(FFT(z)) and Im(FFT(z)) is random.

    Parameters:
    dim : int
        Dimension of the field (1, 2, or 3).
    N : int
        Size of the field (N x N for 2D, N x N x N for 3D).
    Hurst : float
        Hurst exponent, in the range [0, 1].
    k_low : float
        Lower bound of the wavenumber range (k_low > 0).
    k_high : float
        Upper bound of the wavenumber range (k_high < 0.5, which represents Nyquist frequency).
    plateau : bool
        If True, the power spectrum is flat up to k_low, otherwise all wavenumbers below k_low are set to zero.
    verbose : bool
        If True, print the parameters used for generating the random field.

    Returns:
    z : ndarray
        The generated random field.
    """
    if not (0 < k_low <= k_high <= 0.5):
        raise ValueError("Require 0 < k_low <= k_high <= 0.5 (Nyquist)")

    if verbose:
        print("Random Phase Field parameters:",
              f"dim={dim}, N={N}, H={Hurst}, k_low={k_low}, k_high={k_high}")

    # Frequency grid (unshifted, consistent with numpy.fft)
    k1 = fftfreq(N)
    if dim == 1:
        k = np.abs(k1)
        shape = (N,)
    elif dim == 2:
        kx, ky = np.meshgrid(k1, k1, indexing='ij')
        k = np.sqrt(kx**2 + ky**2)
        shape = (N, N)
    elif dim == 3:
        kx, ky, kz = np.meshgrid(k1, k1, k1, indexing='ij')
        k = np.sqrt(kx**2 + ky**2 + kz**2)
        shape = (N, N, N)
    else:
        raise ValueError(f"Unknown dim={dim}")

    # Self-affine exponent for the amplitude (not the PSD):
    power = -(0.5 * dim + Hurst)

    # Prescribed Fourier magnitude A(k) (zero outside [k_low, k_high])
    sqrt_power_spectrum = np.zeros_like(k, dtype=float)
    mask = (k >= k_low) & (k <= k_high)
    # Avoid division by zero; k==0 is excluded by mask anyway
    sqrt_power_spectrum[mask] = (k[mask] / k_low) ** power
    sqrt_power_spectrum[k == 0] = 0.0
    if plateau:
        sqrt_power_spectrum[k < k_low] = 1.0

    # Get random phases that obey Hermitian symmetry automatically
    # (from FFT of a real white-noise field)
    w = np.random.normal(size=shape)
    W = fftn(w)
    phase = W / (np.abs(W) + 1e-30)  # unit-magnitude complex numbers

    # Impose the target magnitudes with those phases
    F = sqrt_power_spectrum * phase

    # Back to real space (Hermitian symmetry is preserved by construction)
    z = np.real(ifftn(F))

    return z
