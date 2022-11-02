"""
@author:Ivan Munoz-Rodriguez, Antonis Georgakakis and Angel Ruiz Camunas

This script contain the functions to implement selection effects of observations onto simulations
"""
# from astropy import units as u
# from astropy import constants as const
import logging
from colossus.halo import mass_defs

# from colossus.cosmology import cosmology
import numpy as np
import pandas as pd

# from mockUniv.config import *
# import mag_mass_convert
# import calc_k_correct

logger = logging.getLogger(__name__)


def get_mask_distFace(dist_spec, cube):
    min_dist_face = dist_spec
    max_dist_face = cube.size - dist_spec
    mask_x = np.logical_and(
        np.array(cube["x"] > min_dist_face), np.array(cube["x"] < max_dist_face)
    )
    mask_y = np.logical_and(
        np.array(cube["y"] > min_dist_face), np.array(cube["y"] < max_dist_face)
    )
    mask_xy = np.logical_and(mask_x, mask_y)
    mask_z = np.logical_and(
        np.array(cube["z"] > min_dist_face), np.array(cube["z"] < max_dist_face)
    )

    return np.logical_and(mask_xy, mask_z)


def mask_magnitude(mag_lim, filter, cube):
    return np.array(cube[f"{filter}"] < mag_lim)


def mask_mass(lgM_smCut, cube):
    return np.array(cube["lgobs_sm"] > lgM_smCut)


def mask_lumi(lglumi_cut, model, cube):
    return np.array(cube["{}".format(model)] > lglumi_cut)


def _get_masks_mag_sm_lumi_lc_list(mag_lim, sm_lim, lc_list, mag_filter, lx_lim):
    """ """
    if mag_lim != False:
        mask_mag = [mask_magnitude(mag_lim, mag_filter, lc) for lc in lc_list]
    else:
        mask_mag = [False]

    if sm_lim != False:
        mask_sm = mask_sm = [mask_mass(sm_lim, lc) for lc in lc_list]
    else:
        mask_sm = [False]

    mask_xrlumi_dict = []
    for lc in lc_list:
        mask_xrlumi_dict.append(
            {
                "aird": np.array(lc["lgLxAMedian"] > lx_lim),
                "age": np.array(lc["lgLxGMedian"] > lx_lim),
            }
        )
    return mask_mag, mask_sm, mask_xrlumi_dict


def get_spatial_mask(
    id_cluster, ang, lc, delta_z, field=False, ang_max=6 / 60, r200=False
):
    """
    Function that returns the two masks, one that limits the objects around the cluster in redshift
    for a given delta(redshift). The other mask objects at a certain distance to the center of the cluster. Note that
    everything is done respect to the cluster and this avoid to project actively using light cones. This means that redshift
    and angle are translated in cosmological distances, comoving and angular diameter respectively.

    Input:
    ------
        - id_cluster (int): cluster id of the cluster that is the center of the observation
        - ang (float): angle that defines the mask around the clusted, in deg!
        - cube (cube object): box cube object, id_cluster should be of this box.
        - delta_z (float): delta redshift around the cluster limit
        - field (bool): set if the observation is cluster (False) or the Field (True)
        - ang_max (bool or float): maximum angle for field observation, if field==False this should be also False
    Output:
    -------
        - mask_spec, mask_ang
    """
    mask_clust_parent = lc["id"] == id_cluster

    RA_clust = float(lc[mask_clust_parent]["RA"])
    dec_clust = float(lc[mask_clust_parent]["dec"])

    mask_spec = get_mask_spec(lc)
    if field == False:
        r200 = float(lc[mask_clust_parent]["r"])
        mask_ang = get_mask_ang(lc, ang, RA_clust, dec_clust, r200)

    if field == True:
        ang = ang
        ang_max = ang_max
        mask_ang = get_mask_ang_field(lc, ang, ang_max)
    return mask_spec, mask_ang


def _get_spatial_mask_join(id_cluster, angle, delta_z, cube, field):
    """
    angle in deg
    """
    mask_spec, mask_ang = get_spatial_mask(id_cluster, angle, cube, delta_z)  # , field)
    return np.logical_and(mask_spec, mask_ang)


def get_mask_spec(lc):
    """
    Function that calculate the spatial mask i.e. spectroscopic criteria for a given lc and delta redshift.
    """
    mask_spec = np.array(
        np.abs(lc["final_redshift"] - lc.observation["desire_redshift"])
        < delta_z_spec_crit(lc)
    )
    return mask_spec


def delta_z_spec_crit(lc, cube=False):
    if cube == False:
        if lc.z > 1:
            sigma = 2000
        else:
            c = 10**0.7
            m_clust = float(10 ** (lc[lc["id"] == lc.id_obs]["lgm"]))
            lgm_clust_200 = np.log10(
                mass_defs.changeMassDefinition(m_clust, c, lc.z, "vir", "200c")[0]
            )
            sigma = 3 * get_sigma_cluster(lc.cosmo, lgm_clust_200, lc.z)
    else:
        sigma = 2000
    return sigma * (1 + lc.z) / ((const.c).to(u.km / u.s)).value


def comovingDist_2_proper(com_dist, z):
    return com_dist / (1 + z)


def get_mask_ang(lc, ang, r200):
    """
    angle in arcmin
    """
    if r200 == False:
        dist_ang = lc.cosmo.angularDiameterDistance(lc.z) * (np.pi / 180) * ang
    else:
        dist_ang = _calc_angle_from_r(lc)
    mask_ang = lc["angle_los"] < dist_ang
    return mask_ang


def get_mask_ang_field(cube, ang_min, ang_max, x_clust, y_clust):
    """
    angle in arcmin
    """
    dist_ang_min = cube.cosmo.angularDiameterDistance(cube.z) * (np.pi / 180) * ang_min
    dist_ang_max = cube.cosmo.angularDiameterDistance(cube.z) * (np.pi / 180) * ang_max

    mask_ang_min = np.array(
        np.sqrt((cube["x"] - x_clust) ** 2 + (cube["y"] - y_clust) ** 2) > dist_ang_min
    )
    mask_ang_max = np.array(
        np.sqrt((cube["x"] - x_clust) ** 2 + (cube["y"] - y_clust) ** 2) < dist_ang_max
    )
    return mask_ang_min & mask_ang_max


def get_dyn_mass(cosmo, sigma, redshift):
    """
    Dynamical mass as defined in eq. (1) Munari et al. +13 (https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2638M/abstract)

    Return:
    ------
    M_200,dyn in solar masses
    """
    A = 1177
    alpha = 0.364
    hz = cosmo.Hz(redshift) / 100
    return 15 + (np.log10(sigma / A) / alpha) - np.log(hz)


def get_sigma_cluster(cosmo, lgm_200, redshift):
    """
    Dynamical mass as defined in eq. (1) Munari et al. +13 (https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2638M/abstract)

    Input:
    ------
    cosmo
    lgm_200 (Msun)
    redshift

    Return:
    ------
    sigma in km/s
    """
    A = 1177
    alpha = 0.364
    hz = cosmo.Hz(redshift) / 100
    return 10 ** (np.log10(A) + alpha * (np.log10(hz) + lgm_200 - 15))


def _add_m200(cube):
    """
    Concentration from Ludlow+14: https://ui.adsabs.harvard.edu/abs/2014MNRAS.441..378L/abstract
    """
    c = 10**0.7  # Fig.1 Ludlow+14
    for mass in ["lgmp", "lgm"]:
        mass200m = mass_defs.changeMassDefinition(
            10 ** cube[mass], c, cube.z, "vir", "200m"
        )
        cube[mass + "200m"] = np.log10(mass200m[0] / cube.cosmo.h)
        cube[f"R200m_{mass}"] = ((np.array(mass200m[1]) * u.kpc).to(u.Mpc)).value  #
    cube[["lgmp", "lgm", "lgmp200m", "lgm200m"]]
    return


def _calc_angle_from_r(lc):
    """
    Function that calculates the angle corresponding to a cluster in **PHYSICAL** units in degs
    """
    mask_id = lc["id"].isin([lc.id_obs])
    d_ang_cosmo = lc.cosmo.angularDiameterDistance(
        lc.observation["desire_redshift"]
    ) * (u.Mpc / u.rad)
    r_vir_physUnits = comovingDist_2_proper((float(lc[mask_id]["r"]) * u.kpc), lc.z)
    FOV = float(((r_vir_physUnits / d_ang_cosmo).to(u.rad)).value)
    return FOV


def _calc_angle_from_r_qties(cosmo, redshift, r):
    """
    Function that calculates the angle corresponding to a cluster in **PHYSICAL** units in degs
    """
    d_ang_cosmo = cosmo.angularDiameterDistance(redshift) * (u.Mpc / u.rad)
    r_vir_physUnits = comovingDist_2_proper((r * u.kpc), redshift)
    FOV = float(((r_vir_physUnits / d_ang_cosmo).to(u.deg)).value)
    return FOV
