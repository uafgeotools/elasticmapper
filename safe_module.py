import concurrent.futures
import os
from functools import partial

import numpy as np
from scipy.optimize import differential_evolution
from scipy.optimize import minimize

from change_of_basis import tmat_of_cmat
from find_symmetry_groups import betaT
from find_symmetry_groups import DistToVSigmaofU
from find_symmetry_groups import objective
from find_symmetry_groups import proj_to_vsig_of_u_new
from find_symmetry_groups import UsHat
from utilities import v2sm

def distance(c_vec, sigma, tracker=None, method="differential_evolution", popsize=15, number_of_runs=10, sample_size=250,
             number_of_minima=10, use_parallel_processing=True):

    # wraps around GetTempAndT0S0P0
    # function designed to call GetTempAndT0S0P0 within a user defined parallelizatin routine

    # To keep track of how much of the parallelized job is completed
    if tracker:
        print(f'{tracker} \n')

    c_mat = v2sm(c_vec)
    t_mat = tmat_of_cmat(c_mat)

    temp = GetTempAndT0S0P0(t_mat, sigma, method=method, popsize=popsize, number_of_runs=number_of_runs,
                            sample_size=sample_size, number_of_minima=number_of_minima,
                            use_parallel_processing=use_parallel_processing)

    return [betaT(t_mat, temp[0]), temp[1]['theta'], temp[1]['sigma'], temp[1]['phi']]

def closest(Tmat, Sigma):
    # returns closest Tmat (6x6 matrix)
    temp = GetTempAndT0S0P0(Tmat, Sigma)
    return proj_to_vsig_of_u_new(Tmat, UsHat([temp[1]['theta'], temp[1]['sigma'], temp[1]['phi']]), Sigma)

def GetTempAndT0S0P0(Tmat, Sigma, method="differential_evolution", popsize=15, number_of_runs=10, sample_size=250,
                     number_of_minima=10, use_parallel_processing=True):

    if use_parallel_processing:
        workers = os.cpu_count()
    else:
        workers = 1

    objective_new = partial(objective, Tmat=Tmat, Sigma=Sigma)

    if Sigma == "ISO":

        id   = np.eye(3)
        temp = (DistToVSigmaofU(Tmat, id, "ISO"), {'theta': 0, 'sigma': 0, 'phi': 0})

    else:

        theta_bounds = (0, 2 * np.pi)
        sigma_bounds = (-np.pi, np.pi)
        phi_bounds   = (0, np.pi)
        bounds       = [theta_bounds, sigma_bounds, phi_bounds]

        if method=="differential_evolution":

            results = []
            for i in range(number_of_runs):
                result = differential_evolution(objective_new, bounds, maxiter=1000, popsize=popsize, tol=0.01,
                                                workers=workers)
                results.append(result)

        elif method in ["gradient", "random_search_and_gradient"]:

            thetas = np.random.uniform(theta_bounds[0], theta_bounds[1], sample_size)
            sigmas = np.random.uniform(sigma_bounds[0], sigma_bounds[1], sample_size)
            phis   = np.arccos(np.random.uniform(np.cos(phi_bounds[0]), np.cos(phi_bounds[1]), sample_size))

            minimize_new = partial(minimize, bounds=bounds, method='L-BFGS-B')

            if method == "gradient":

                if use_parallel_processing:
                    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                        results = list(executor.map(minimize_new, [objective_new]*sample_size, zip(thetas, sigmas, phis)))

                else:
                    results = []
                    for i in range(sample_size):
                        result = minimize_new(objective_new, [thetas[i], sigmas[i], phis[i]])
                        results.append(result)

            elif method == "random_search_and_gradient":

                if use_parallel_processing:
                    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                        betas = list(executor.map(objective_new, zip(thetas, sigmas, phis)))

                    indices = np.argsort(betas)[:number_of_minima]

                    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                        results = list(executor.map(minimize_new, [objective_new] * number_of_minima,
                                                    zip(thetas[indices], sigmas[indices], phis[indices])))

                else:
                    betas   = []
                    results = []

                    for i in range(sample_size):
                        beta = objective_new([thetas[i], sigmas[i], phis[i]])
                        betas.append(beta)

                    indices = np.argsort(betas)[:number_of_minima]

                    for index in indices:
                        result = minimize_new(objective_new, [thetas[index], sigmas[index], phis[index]])
                        results.append(result)

        else:

            print("ERROR: selected optimization method does not exist!")

        minima = [result.fun for result in results]
        result = results[np.argmin(minima)]

        temp = (result.fun, {'theta': result.x[0], 'sigma': result.x[1], 'phi': result.x[2]})

    return temp