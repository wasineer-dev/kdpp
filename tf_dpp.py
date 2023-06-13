import numpy as np
import itertools
import tensorflow as tf
import tensorflow_probability as tfp
import matplotlib.pyplot as plt

def elementary_sym_poly(N, k, e_vals):
    polyn = np.zeros((N, k))
    for i in np.arange(N):
       polyn[i][0] = 1
    for j in np.arange(k):
       polyn[0][j] = 0

    for l in np.arange(k):
       for n in np.arange(N):
          polyn[n][l] = polyn[n-1][l] + e_vals[n]*polyn[n-1][l-1]

    return polyn
 
tfd = tfp.distributions
tfpk = tfp.math.psd_kernels

grid_size = 16
# Generate grid_size**2 pts on the unit square.
grid = np.arange(0, 1, 1./grid_size)
points = np.array(list(itertools.product(grid, grid)))

# Create the kernel L that parameterizes the DPP.
kernel_amplitude = 2.
kernel_lengthscale = 2. / grid_size
kernel = tfpk.ExponentiatedQuadratic(kernel_amplitude, kernel_lengthscale)
kernel_matrix = kernel.matrix(points, points)

eigenvalues, eigenvectors = tf.linalg.eigh(kernel_matrix)
dpp = tfd.DeterminantalPointProcess(eigenvalues, eigenvectors)

# The inner-most dimension of the result of `dpp.sample` is a multi-hot
# encoding of a subset of {1, ..., ground_set_size}.

plt.figure(figsize=(6, 6))
for i, samp in enumerate(dpp.sample(4, seed=(1, 2))):  # 4 x grid_size**2
  plt.subplot(221 + i)
  plt.scatter(*points[np.where(samp)].T)
  plt.xticks([])
  plt.yticks([])
plt.tight_layout()
plt.show()

# Like any TFP distribution, the DPP supports batching and shaped samples.

kernel_amplitude = [2., 3, 4]  # Build a batch of 3 PSD kernels.
kernel_lengthscale = 2. / grid_size
kernel = tfpk.ExponentiatedQuadratic(kernel_amplitude, kernel_lengthscale)
kernel_matrix = kernel.matrix(points, points)  # 3 x 256 x 256

eigenvalues, eigenvectors = tf.linalg.eigh(kernel_matrix)
dpp = tfd.DeterminantalPointProcess(eigenvalues, eigenvectors)
print(dpp)  # batch shape: [3], event shape: [256]
samps = dpp.sample(2, seed=(10, 20))
print(samps.shape)  # shape: [2, 3, 256]
print(dpp.log_prob(samps))  # tensor with shape [2, 3]