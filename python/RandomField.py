"""
Generate a 2D Gaussian random field from white noise filtered in Fourier space.
The method is described in the following paper:
    Hu, Y.Z. and Tonder, K., 1992. Simulation of 3-D random rough surface by 2-D digital filter and Fourier analysis. International journal of machine tools and manufacture, 32(1-2), pp.83-90.
    DOI: 10.1016/0890-6955(92)90064-N
Author: Vladislav Yastrebov, MINES Paris - PSL, Centre des matériaux, CNRS UMR 7633, BP 87, 91003 Evry, France (with help of GPT4 and Copilot)
Date: Dec 2023
License: CC0
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft2, ifft2, fftfreq, fftshift, ifftshift, fftn, ifftn


def periodic_gaussian_random_field(dim = 2, N = 256, Hurst = .5, k_low = 0.03, k_high = 0.3):
    print("parameters:")
    print("dim = ", dim)
    print("N = ", N)
    print("Hurst = ", Hurst)
    print("k_low = ", k_low)
    print("k_high = ", k_high)    
    if k_low < 0 or k_high < 0 or k_low > k_high:
        print("k_low and k_high must be positive and k_low <= k_high")
        exit(1)
    if k_high > 0.5:
        print("k_high must be less than 0.5, which represents Nyquist frequency")
        exit(1)

    # Define the power spectrum
    k = fftfreq(N)
    power = -(0.5*dim+Hurst)
    noise = []
    if dim == 1:
        noise = np.fft.fft(np.random.normal(size=(N)))
    if dim == 2:
        k = np.sqrt(k[:,None]**2 + k[None,:]**2)
        print("k.shape = ", k.shape)
        print("kmin = ", np.min(k))
        print("kmax = ", np.max(k))
        noise = np.fft.fft2(np.random.normal(size=(N, N)))
    elif dim == 3:
        k = np.sqrt(k[:,None,None]**2 + k[None,:,None]**2 + k[None,None,:]**2)
        noise = fftn(np.random.normal(size=(N, N, N)))
    else:
        print("Unknown dimension ", dim)
        exit(1)
    k = fftshift(k)

    # Define the power spectrum for k values in the range [k_low, k_high]
    power_spectrum = np.where((k >= k_low) & (k <= k_high), (k/k_low)**power, 0)

    # Make sure the mean power is finite
    power_spectrum[np.isinf(power_spectrum)] = 0

    # Take the product of the white noise field with the power spectrum (filter)
    grf_fourier = noise * power_spectrum

    # Transform back to real space
    if dim == 2:
            grf = np.real(ifft2(ifftshift(grf_fourier)))
    else:
        grf = np.real(ifftn(ifftshift(grf_fourier)))
    return grf


######################
# Test 2D version
######################

N0 = 1024           # Size of the random field
k_low = 4 / N0      # Lower cutoff of the power spectrum
k_high = 256 / N0   # Upper cutoff of the power spectrum
Hurst = 0.6         # Hurst exponent
dim = 2             # Dimension of the random field
seed = 124          # Seed for the random number generator
np.random.seed(seed)

# Generate the random field
random_field = periodic_gaussian_random_field(dim = dim, N = N0, Hurst = Hurst, k_low = k_low, k_high = k_high)

# Normalize
random_field /= np.std(random_field)

# Plot
plt.rcParams['figure.figsize'] = [5,5]
fig,ax = plt.subplots()
ax.axis('off')
im = plt.imshow(random_field, cmap='RdYlBu_r', interpolation='bicubic')
plt.colorbar(im, orientation='horizontal', shrink=0.82, pad=0.02)
plt.tight_layout()
plt.show()
fig.savefig("RandomField.png", dpi=300)
