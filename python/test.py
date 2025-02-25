######################
# Test 2D version
######################

import numpy as np
import matplotlib.pyplot as plt
import RandomField as rf

N0 = 1024           # Size of the random field
k_low = 4 / N0      # Lower cutoff of the power spectrum
k_high = 256 / N0   # Upper cutoff of the power spectrum
Hurst = 0.6         # Hurst exponent
dim = 2             # Dimension of the random field
seed = 124          # Seed for the random number generator
np.random.seed(seed)

# Generate the random field
random_field = rf.periodic_gaussian_random_field(dim = dim, N = N0, Hurst = Hurst, k_low = k_low, k_high = k_high)

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
