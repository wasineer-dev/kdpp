import numpy as np
from dppy.finite_dpps import FiniteDPP
from dppy.exact_sampling import elementary_symmetric_polynomials, proj_dpp_sampler_eig

from dppy.utils import (check_random_state,
                        is_symmetric,
                        is_projection,
                        is_orthonormal_columns,
                        is_full_row_rank,
                        is_in_01,
                        is_geq_0,
                        is_equal_to_O_or_1)

class Finite_kDPP(FiniteDPP):

    ###############
    # Constructor #
    ###############
    
    def k_dpp_eig_vecs_selector(self, eig_vals, eig_vecs, size,
                            E_poly=None, random_state=None):
    
        rng = check_random_state(random_state)

        # Size of: ground set / sample
        N, k = eig_vecs.shape[0], size

        # as in np.linalg.matrix_rank
        tol = np.max(eig_vals) * N * np.finfo(np.float64).eps
        rank = np.count_nonzero(eig_vals > tol)
        if k > rank:
            raise ValueError('size k={} > rank={}'.format(k, rank))

        if E_poly is None:
            E_poly = elementary_symmetric_polynomials(eig_vals, k)

        ind_selected = np.zeros(k, dtype=int)
        for n in range(eig_vals.size, 0, -1):

            if rng.rand() < eig_vals[n - 1] * E_poly[k - 1, n - 1] / E_poly[k, n]:
                k -= 1
                ind_selected[k] = n - 1
                if k == 0:
                    break

        return eig_vecs[:, ind_selected]

    def sample_exact_k_dpp(self, size, mode='GS', **params):
         
        rng = check_random_state(params.get('random_state', None))

        self.sampling_mode = mode
        self.size_k_dpp = size
    
        assert(self.L_eig_vals is not None)

        # Phase 1
        # Precompute elementary symmetric polynomials
        if self.E_poly is None or self.size_k_dpp < size:
            self.E_poly = elementary_symmetric_polynomials(self.L_eig_vals,
                                                            size)
        # Select eigenvectors
        V = self.k_dpp_eig_vecs_selector(self.L_eig_vals, self.eig_vecs,
                                    self.size_k_dpp,
                                    self.E_poly,
                                    rng)
        # Phase 2
        self.size_k_dpp = size
        sampl = proj_dpp_sampler_eig(V, self.sampling_mode,
                                        random_state=rng)
        
        self.list_of_samples.append(sampl)
        return sampl
    