######################
# Test 2D version
######################

import numpy as np
import matplotlib.pyplot as plt
import RandomField as rf

def plot_surface_and_spectrum(random_field, file_name="RandomField.png"):
    plt.rcParams['figure.figsize'] = [10,5]
    fig,ax = plt.subplots(1,2)

    # Colorbar width control
    colorbar_shrink = 0.7  # Adjust this value to control colorbar width (0.1 to 1.0)

    # Surface
    # ax[0].axis('off')
    im = ax[0].imshow(random_field, cmap='RdYlBu_r', interpolation='bicubic')
    ax[0].set_title(r"Random Field (2D), $z/\sigma$")
    ax[0].set_aspect('equal')
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_xlabel("$x$")
    ax[0].set_ylabel("$y$")
    plt.colorbar(im, ax=ax[0], orientation='horizontal', shrink=colorbar_shrink, pad=0.05, label=r'Field Value, $z/\sigma$')

    # Power spectrum
    k = np.fft.fftfreq(random_field.shape[0])
    k = np.sqrt(k[:,None]**2 + k[None,:]**2)
    power_spectrum = np.abs(np.fft.fft2(random_field))**2
    # Center the power spectrum
    k = np.fft.fftshift(k)
    power_spectrum = np.fft.fftshift(power_spectrum)
    # Apply the cutoff
    power_spectrum[power_spectrum < 1e-20] = np.nan 
    ax[1].set_aspect('equal')
    # ax[1].axis('off')
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[1].set_xlabel("$k_x$")
    ax[1].set_ylabel("$k_y$")

    im = ax[1].imshow(np.log10(power_spectrum), cmap='viridis', interpolation='none')
    ax[1].set_title(r"Power Spectrum (log scale), $\Phi(k_x, k_y)$")
    plt.colorbar(im, ax=ax[1], orientation='horizontal', shrink=colorbar_shrink, pad=0.05, label=r'Power Spectrum, $log_{10}(\Phi)$')

    plt.tight_layout()
    plt.show()
    fig.savefig(file_name, dpi=300)


def main():
    N = 1024            # Size of the random field
    k_low = 64 / N      # Lower cutoff of the power spectrum
    k_high = 256 / N    # Upper cutoff of the power spectrum
    Hurst = 0.6         # Hurst exponent
    dim = 2             # Dimension of the random field
    seed = 124          # Seed for the random number generator
    plateau = True      # Use plateau in the power spectrum
    np.random.seed(seed)

    # Generate the random field
    random_field = rf.periodic_gaussian_random_field(dim = dim, N = N, Hurst = Hurst, k_low = k_low, k_high = k_high, plateau = plateau)

    # Normalize
    random_field /= np.std(random_field)

    # Plot the surface and its power spectrum
    plot_surface_and_spectrum(random_field, file_name="RandomField_2D_filtering.png")

    # Generate the random field with an ideal spectrum
    ideal_random_field = rf.ideal_periodic_gaussian_random_field(dim = dim, N = N, Hurst = Hurst, k_low = k_low, k_high = k_high, plateau = plateau)

    # Normalize
    ideal_random_field /= np.std(ideal_random_field)

    # Plot the surface and its power spectrum
    plot_surface_and_spectrum(ideal_random_field, file_name="RandomField_2D_ideal.png")

if __name__ == "__main__":
    main()
