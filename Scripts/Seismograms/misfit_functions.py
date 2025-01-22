# misfit function definitions
# after https://github.com/adjtomo/seisflows (commit 6afdd56)

import numpy as np
from scipy.signal import hilbert as analytic

def waveform(st_syn, st_obs):
    misfit = 0
    nt = st_obs[0].stats.npts
    dt = st_obs[0].stats.delta
    misfit_components = ['Z', 'Y', 'X']
    for component in misfit_components:
        syn = st_syn.select(component=component)[0].data
        obs = st_obs.select(component=component)[0].data
        wrsd = syn - obs
        misfit += np.sqrt(np.sum(wrsd * wrsd * dt))
    return misfit

def envelope(st_syn, st_obs):
    misfit = 0
    nt = st_obs[0].stats.npts
    dt = st_obs[0].stats.delta
    misfit_components = ['Z', 'Y', 'X']
    for component in misfit_components:
        syn = st_syn.select(component=component)[0].data
        obs = st_obs.select(component=component)[0].data
        env_syn = abs(analytic(syn))
        env_obs = abs(analytic(obs))
        env_rsd = env_syn - env_obs
        misfit += np.sqrt(np.sum(env_rsd * env_rsd * dt))
    return misfit

def instantaneous_phase(st_syn, st_obs):
    misfit = 0
    nt = st_obs[0].stats.npts
    dt = st_obs[0].stats.delta
    misfit_components = ['Z', 'Y', 'X']
    for component in misfit_components:
        syn = st_syn.select(component=component)[0].data
        r = np.real(analytic(syn))
        i = np.imag(analytic(syn))
        phi_syn = np.arctan2(i, r)
        obs = st_obs.select(component=component)[0].data
        r = np.real(analytic(obs))
        i = np.imag(analytic(obs))
        phi_obs = np.arctan2(i, r)
        phi_rsd = phi_syn - phi_obs
        misfit += np.sqrt(np.sum(phi_rsd * phi_rsd * dt))
    return misfit

def traveltime(st_syn, st_obs):
    misfit = 0
    nt = st_obs[0].stats.npts
    dt = st_obs[0].stats.delta
    misfit_components = ['Z', 'Y', 'X']
    for component in misfit_components:
        syn = st_syn.select(component=component)[0].data
        obs = st_obs.select(component=component)[0].data
        cc = abs(np.convolve(obs, np.flipud(syn)))
        misfit += (np.argmax(cc) - nt + 1) * dt
    return misfit