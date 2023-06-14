# from numpy import sqrt
import argparse
import numpy as np
import pandas as pd

import itertools
import tensorflow as tf
import tensorflow_probability as tfp

from scipy.stats import gamma

tfd = tfp.distributions
tfpk = tfp.math.psd_kernels

from numpy.random import rand, randn
from numpy.linalg import eigh
from Finite_kDPP import Finite_kDPP

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Testing DPPy ....')
parser.add_argument('filename')
args = parser.parse_args()
distance_file = args.filename

N_SAMPLES = 16
df = pd.read_table(distance_file, header=None)

ids = ["S01",
       "S02",
       "S03",
       "S04",
       "S05",
       "S06",
       "S07",
       "S08",
       "S09",
       "S10",
       "S11",
       "S12",
       "S13",
       "S14",
       "S15",
       "S16"]

def convert_to_tf(sampl):
    assert(len(sampl) <= N_SAMPLES)
    q = np.zeros(N_SAMPLES)
    for k in sampl:
        q[k] = 1
    return tf.constant(q)

def sampleIndex(sampleName):
    return ids.index(sampleName)

def getSampleId(basename):
    id = sampleIndex(basename.split("_")[0])
    rndId = basename.split("_")[-1]
    return (id, rndId)

dmatrix = np.zeros((N_SAMPLES, N_SAMPLES))
for source, target, dt in zip(df[1], df[2], df[3]):
    #print(sampleIndex(source.split("_")[0]), sampleIndex(target.split("_")[0]))
    i, rndId = getSampleId(source.split(".")[0])
    j, rndId = getSampleId(target.split(".")[0].split("/")[1])
    if (i != j):
        dmatrix[i,j] = float(dt)

beta = 0.2
similarity = 0.01*np.exp(-np.square(dmatrix/beta))
e_vals_K, e_vecs = eigh(similarity)
print(e_vals_K)

#dpp_K = FiniteDPP('correlation', **{'K_eig_dec': (e_vals_K/np.max(e_vals_K), e_vecs)})
# or
#K = (e_vecs * e_vals_K).dot(e_vecs.T)
#K = similarity
#dpp_K = FiniteDPP('correlation', **{'K': K})
#dpp_K.plot_kernel()
#plt.show()

dpp_K = Finite_kDPP('correlation', **{'K_eig_dec': (e_vals_K, e_vecs)})
dpp_K.plot_kernel()
plt.show()

rng = np.random.RandomState(457)

e_vals_K = e_vals_K
e_vals_L = e_vals_K / (1.0 - e_vals_K)
print(e_vals_L)
dpp_L = Finite_kDPP('likelihood', **{'L_eig_dec': (e_vals_L, e_vecs)})
dpp_L.plot_kernel()
plt.show()

k = 10

def sampler(number_of_samples):
    dpp_L.flush_samples()
    for _ in range(number_of_samples):
        dpp_L.sample_exact_k_dpp(size=k)


kernel_matrix = 0.1*tf.exp(-np.square(dmatrix/beta))
eigenvalues, eigenvectors = tf.linalg.eigh(kernel_matrix)
tf_dpp = tfd.DeterminantalPointProcess(eigenvalues, eigenvectors)

def sample_prob(n_iter):
    a = 7.7
    loc = 53
    scale = 0.423
    sampler(n_iter)
    result_log_prob = np.abs(list(map(tf_dpp.log_prob, map( convert_to_tf, dpp_L.list_of_samples ))))
    max_ind = np.argmax(result_log_prob)
    x = np.max(result_log_prob)
    ts = convert_to_tf(dpp_L.list_of_samples[max_ind])
    if gamma.sf(x, a, loc, scale) < 1e-05:
        print(np.sort(dpp_L.list_of_samples[max_ind]))
    #print(gamma.fit(result_log_prob))

for _ in range(10):
    sample_prob(1000)

#plt.hist(result_log_prob)
#plt.show()

