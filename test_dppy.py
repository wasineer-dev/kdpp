# from numpy import sqrt
import argparse
import numpy as np
import pandas as pd

from numpy.random import rand, randn
from numpy.linalg import eigh
from scipy.linalg import qr
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

similarity = np.exp(-np.square(dmatrix/0.2))
e_vals_K, e_vecs = eigh(similarity)
print(e_vals_K)

#dpp_K = FiniteDPP('correlation', **{'K_eig_dec': (e_vals_K/np.max(e_vals_K), e_vecs)})
# or
#K = (e_vecs * e_vals_K).dot(e_vecs.T)
#K = similarity
#dpp_K = FiniteDPP('correlation', **{'K': K})
#dpp_K.plot_kernel()
#plt.show()

rng = np.random.RandomState(457)

e_vals_K = e_vals_K/10.0
e_vals_L = e_vals_K / (1.0 - e_vals_K)
dpp_L = Finite_kDPP('likelihood', **{'L_eig_dec': (e_vals_L, e_vecs)})
k = 10

dpp_L.flush_samples()

for _ in range(1000):
    dpp_L.sample_exact_k_dpp(size=k)

nb_samples = 100
for _ in range(nb_samples):
    dpp_L.flush_samples()
    sampl = dpp_L.sample_mcmc_k_dpp(size=8, nb_iter=10000)
    print(sampl)

sizes = list(map(len, dpp_L.list_of_samples))
print('E[|X|]:\n emp={:.3f}, theo={:.3f}'
      .format(np.mean(sizes), np.sum(e_vals_K)))
print('Var[|X|]:\n emp={:.3f}, theo={:.3f}'
      .format(np.var(sizes), np.sum(e_vals_K*(1-e_vals_K))))
