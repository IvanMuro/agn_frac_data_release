import copy
import logging

#from astropy import units as u
from colossus.halo import mass_defs
import numpy as np
import pandas as pd

from tqdm import tqdm
#import mockUniv.seed.seed_luminosity as seed_lumi
import calc_frac_paper.selection_effects as select_eff

logger = logging.getLogger(__name__)


def _get_spatial_mask(id_cluster, angle, delta_z, cube, field):
    """
    angle in deg
    """
    mask_spec, mask_ang = select_eff.get_spatial_mask(
        id_cluster, angle, cube, delta_z)  # , field)
    return np.logical_and(mask_spec, mask_ang)


def _calc_num_gal(cube, mask_selection, mask_spatial):
    mask_gal = np.logical_and(mask_spatial, mask_selection)
    ngal = len(cube[mask_gal])
    return ngal


def _calc_num_agn(cube, mask_xrlumi, mask_spatial, mask_selection_gal, mask_gal_agn=True):
    mask_agn = np.logical_and(mask_spatial, mask_xrlumi)
    if mask_gal_agn == True:
        mask_agn = np.logical_and(mask_agn, mask_selection_gal)
    nagn = len(cube[mask_agn])
    return nagn

def get_resample_agns(df_num_gal_agn, mask_xrlumi_dict, num_realizations='default', mag_lim=False,
                      sm_lim=False, id_label='id'):
    """
    Function that resample a given data frame with ids, num gal and num agns. Used to calculate error assotiated
    to the mean. From the data frame w number of galaxies/agns for diff. selection effect we generate a sample of
    random number for the range of index in the data frame. We loop through each of these randm numbers and 
    we sabe the ngal_[...], n_agn_[...] in lists that will be the outputs. The random numbers have dimensions
    of (num_realizations, len(df_num_gal_agn)), i.e. create a random sample of clusters w the same size as the input
    sample. 

    Input:
    ------

    Return:
    -------
    n_gal_realizations, n_agn_realizations (list of arrays): each of the list has dimensions (num_realiz, X), where
    X depend on the number of repetitions. Each vector of the list are the number of gal or agns in a particular random
    resampling
    """

    #logger.info(f'Resampling clusters with {num_realizations} realizations')

    rand_realiz = _get_random_realizations(
        df_num_gal_agn, num_realizations, id_label)

    n_gal_realizations_mag = []
    n_gal_realizations_sm = []

    n_agn_realizations_mag = {}
    n_agn_realizations_sm = {}
    for loop_realization, realization in tqdm(enumerate(rand_realiz)):
        # take number galax./agn for given random idx
        df_resample = df_num_gal_agn.loc[realization]
        # MAG
        if mag_lim != False:
            n_gal_realizations_mag.append(np.array(df_resample['n_gal_mag']))
        else:
            n_gal_realizations_mag.append(np.full((len(df_resample)), np.nan))
        # SM
        if sm_lim != False:
            n_gal_realizations_sm.append(np.array(df_resample['n_gal_sm']))
        else:
            n_gal_realizations_sm.append(np.full((len(df_resample)), np.nan))
        # LUMI
        for mask_xrlumi in mask_xrlumi_dict[0]:
            if loop_realization == 0:
                n_agn_realizations_mag.update({f'n_agn_mag_{mask_xrlumi}':
                                               [np.array(df_resample[f'n_agn_mag_{mask_xrlumi}'])]})
                n_agn_realizations_sm.update({f'n_agn_sm_{mask_xrlumi}':
                                              [np.array(df_resample[f'n_agn_sm_{mask_xrlumi}'])]})
            else:
                n_agn_aux = n_agn_realizations_mag[f'n_agn_mag_{mask_xrlumi}']
                n_agn_aux.append(
                    np.array(df_resample[f'n_agn_mag_{mask_xrlumi}']))
                n_agn_realizations_mag.update(
                    {f'n_agn_mag_{mask_xrlumi}': n_agn_aux})

                n_agn_aux = n_agn_realizations_sm[f'n_agn_sm_{mask_xrlumi}']
                n_agn_aux.append(
                    np.array(df_resample[f'n_agn_sm_{mask_xrlumi}']))
                n_agn_realizations_sm.update(
                    {f'n_agn_sm_{mask_xrlumi}': n_agn_aux})
        # Below comment was supposed to reproduce the above loop, it doesnt (:
        # n_agn_realizations = [ _update_agn_num(loop_realization, n_agn_realizations, mask_xrlumi, ids_repetitions, df_resample)
        #                       for mask_xrlumi in mask_xrlumi_dict ]
    return n_gal_realizations_mag, n_gal_realizations_sm, n_agn_realizations_mag, n_agn_realizations_sm


def _get_random_realizations(df_num_gal_agn, num_realizations='default', id_label='id'):
    print(num_realizations)
    if num_realizations == 'default':
        num_realizations = 100 #len(df_num_gal_agn[id_label])

    logger.info(f'Generating {num_realizations} random realizations')
    size_resample = len(df_num_gal_agn)
    # np.random.seed(121517) #DELETE JUST 4 TESTING
    # Use to random generation below to select by index, the commented select by ids
    # seleting by ids sneaky thing to get the number of galaxies
    # rand_realiz = np.random.choice(df_num_gal_agn['id'],
    #                               size=(num_realizations, size_resample), replace=True)
    rand_realiz = np.random.choice(df_num_gal_agn.index,
                                   size=(num_realizations, size_resample), replace=True)
    return rand_realiz

def _create_df_resample(mask_xrlumi_dict, mag_lim=False, sm_lim=False):
    """
    Funciton that creates a pandas DataFrame with the columns correspondant to the 
    different input selection effects. It returns an empty DataFrame with the correspondent
    columns.
    Input:
    ------

    Output:
    ------
    df_resample_sum (pandas DataFrame): empty with the correspondent columns matching different
        selection effects
    """

    df_resample_sum = pd.DataFrame()
    for mask_xrlumi in mask_xrlumi_dict:
        df_resample_sum[f'num_agn_mag_{mask_xrlumi}'] = []
        df_resample_sum[f'num_agn_mag_{mask_xrlumi}'] = df_resample_sum[f'num_agn_mag_{mask_xrlumi}'].astype(
            'i8')
    for mask_xrlumi in mask_xrlumi_dict:  # the loop is repeated to match how the data is written after
        # at the end of get_resample_summed( ... )
        df_resample_sum[f'num_agn_sm_{mask_xrlumi}'] = []
        df_resample_sum[f'num_agn_sm_{mask_xrlumi}'] = df_resample_sum[f'num_agn_sm_{mask_xrlumi}'].astype(
            'i8')

    # if mag_lim!=False:
    df_resample_sum['num_gal_mag'] = []
    df_resample_sum['num_gal_mag'] = df_resample_sum['num_gal_mag'].astype(
        'i8')
    # if sm_lim!=False:
    df_resample_sum['num_gal_sm'] = []
    df_resample_sum['num_gal_sm'] = df_resample_sum['num_gal_sm'].astype('i8')

    for mask_xrlumi in mask_xrlumi_dict:
        df_resample_sum[f'f_mag_{mask_xrlumi}'] = []
        df_resample_sum[f'f_mag_{mask_xrlumi}'] = df_resample_sum[f'f_mag_{mask_xrlumi}'].astype(
            'f4')
    for mask_xrlumi in mask_xrlumi_dict:
        df_resample_sum[f'f_sm_{mask_xrlumi}'] = []
        df_resample_sum[f'f_sm_{mask_xrlumi}'] = df_resample_sum[f'f_sm_{mask_xrlumi}'].astype(
            'f4')
    return df_resample_sum

#############################
###### RE-SEED VERSION ######
#############################


def join_lc_lst(lc_list, BCG=False, only_sat=False):
    """
    Function that concatenate a list of light-cones in a single data-frame. Care with the attributes they are meaningless
    """

    for loop_lc, lc in enumerate(lc_list):
        # print(lc)
        mask_id = lc['id'].isin([lc.id_obs])
        FOV = select_eff._calc_angle_from_r(lc)
        lc['FOV'] = FOV
        lc['r_parent'] = float(lc[mask_id]['r']/(1+lc.z))
        lc['id_obs'] = lc.id_obs
        lc['delta_redshift_lim'] = select_eff.delta_z_spec_crit(lc)
        RA_lc = lc.RA
        dec_lc = lc.dec
        if loop_lc == 0:
            lc_catalogue_concat = copy.deepcopy(lc)
        else:
            lc_catalogue_concat = lc_catalogue_concat .append(lc,
                                                              ignore_index=True)
    lc_catalogue_concat['delta_redshift'] = lc_catalogue_concat['final_redshift'] - lc.z
    lc_catalogue_concat['delta_RA'] = lc_catalogue_concat['RA'] - RA_lc
    lc_catalogue_concat['delta_dec'] = lc_catalogue_concat['dec'] - dec_lc

    if BCG == True:
        mask_bcg = lc_catalogue_concat .id == lc_catalogue_concat .id_obs
        lc_catalogue_concat = lc_catalogue_concat[mask_bcg]
    if only_sat == True:
        mask_sat = lc_catalogue_concat .upid != -1
        lc_catalogue_concat = lc_catalogue_concat[mask_sat]
    else:
        lc_catalogue_concat = lc_catalogue_concat

    return lc_catalogue_concat


def _create_df_reseed_num_gal_agn_df(lc_catalogue, mask_xrlumi_dict):
    """
    Function that creates the dataframe where the number of gal/agn are saved after
    """
    df_num_gal_agn = pd.DataFrame({'id_cluster_resed': lc_catalogue['id_clust_reseed'].unique()})
    df_num_gal_agn[['id', 'resed']] = df_num_gal_agn['id_cluster_resed'].str.split('_', 1, expand=True)

    df_num_gal_agn['n_gal_mag'] = np.nan
    df_num_gal_agn['n_gal_sm'] = np.nan
    for mask_xrlumi in mask_xrlumi_dict:
        df_num_gal_agn[f'n_agn_mag_{mask_xrlumi}'] = np.nan
        df_num_gal_agn[f'n_agn_sm_{mask_xrlumi}'] = np.nan
    return df_num_gal_agn


def _calc_num_gal_agn_lc_catalogue(df_num_gal_agn, lc, id_cluster_resed, m_cut, band, lx_lim, 
                                   lgM_smCut, mask_gal_agn=True, mask_xrlumi_dict=False, extra_frac=False):
    """
    """
    mask_mag = select_eff.mask_magnitude(m_cut, band, lc)
    mask_sm = select_eff.mask_mass(lgM_smCut, lc)

    mask_aux = {}
    for key in mask_xrlumi_dict:
        mask_aux.update({f'{key}': (np.array(lc[f'lgLx{key}'] > lx_lim))})
    mask_xrlumi_dict = mask_aux

    if m_cut != False:
        df_num_gal_agn.loc[df_num_gal_agn['id_cluster_resed'] == id_cluster_resed,'n_gal_mag'] = len(lc[mask_mag])

    if lgM_smCut != False:
        df_num_gal_agn.loc[df_num_gal_agn['id_cluster_resed'] == id_cluster_resed,'n_gal_sm'] = len(lc[mask_sm])

    for mask_xrlumi in mask_xrlumi_dict:
        mask_lx = mask_xrlumi_dict[mask_xrlumi]

        mask_sm_lx = np.logical_and(mask_lx, mask_sm)
        df_num_gal_agn.loc[df_num_gal_agn['id_cluster_resed'] == id_cluster_resed,
                           f'n_agn_sm_{mask_xrlumi}'] = \
            len(lc[mask_sm_lx])
        if mask_gal_agn:
            mask = np.logical_and(mask_mag, mask_lx)
        else:
            mask = np.array(mask_lx)

        if extra_frac == False:
            df_num_gal_agn.loc[df_num_gal_agn['id_cluster_resed'] == id_cluster_resed,
                               f'n_agn_mag_{mask_xrlumi}'] = \
                len(lc[mask])
        else:
            df_num_gal_agn.loc[df_num_gal_agn['id_cluster_resed'] == id_cluster_resed,
                               f'n_agn_mag_{mask_xrlumi}'] = \
                len(lc[mask])/extra_frac
    return df_num_gal_agn


def get_num_galAgn(lc_catalogue, mask_xrlumi_dict, m_cut, band, lx_lim,
                   lgM_smCut, mask_gal_agn=True, extra_frac=False):
    """

    """
    df_num_gal_agn = _create_df_reseed_num_gal_agn_df(lc_catalogue, mask_xrlumi_dict)
    for id_cluster_resed in lc_catalogue['id_clust_reseed'].unique():
        lc = copy.deepcopy(lc_catalogue)
        lc = lc[lc['id_clust_reseed'] == id_cluster_resed]
        df_num_gal_agn = _calc_num_gal_agn_lc_catalogue(df_num_gal_agn, lc, id_cluster_resed,
                                                        m_cut, band, lx_lim, lgM_smCut,
                                                        mask_gal_agn=mask_gal_agn,
                                                        mask_xrlumi_dict=mask_xrlumi_dict,
                                                        extra_frac=extra_frac)
    return df_num_gal_agn


def _get_resample_reseed_sum(df_resed_resample, mask_xrlumi_dict, lgM_smCut, mag_cut):
    """
    """
    n_gal_realizations_mag, n_gal_realizations_sm, n_agn_realizations_mag, n_agn_realizations_sm = \
        get_resample_agns(df_resed_resample,
                          [mask_xrlumi_dict],
                          num_realizations='default',
                          mag_lim=mag_cut,
                          sm_lim=lgM_smCut,
                          id_label='id_cluster_resed')

    df_resample_sum = _create_df_resample(mask_xrlumi_dict,
                                          mag_lim=mag_cut,
                                          sm_lim=lgM_smCut)

    for loop_realiz, ngal in enumerate(n_gal_realizations_mag):
        df_resample_sum = _add_realization_reseed_sum(df_resample_sum, ngal, n_gal_realizations_sm,
                                                      n_agn_realizations_mag, n_agn_realizations_sm,
                                                      mask_xrlumi_dict, loop_realiz)
    return df_resample_sum


def _add_realization_reseed_sum(df_resample_sum, ngal, n_gal_realizations_sm,
                                n_agn_realizations_mag, n_agn_realizations_sm,
                                mask_xrlumi_dict, loop_realiz):
    """
    """
    num_agn_mag = [sum(n_agn_realizations_mag[f'n_agn_mag_{mask_xrlumi}'][loop_realiz])
                   for mask_xrlumi in mask_xrlumi_dict]

    data = copy.deepcopy(num_agn_mag)
    num_agn_sm = [sum(n_agn_realizations_sm[f'n_agn_sm_{mask_xrlumi}'][loop_realiz])
                  for mask_xrlumi in mask_xrlumi_dict]
    data.extend(num_agn_sm)

    n_gal_total_mag = sum(ngal)
    data.extend([n_gal_total_mag])
    n_gal_total_sm = sum(n_gal_realizations_sm[loop_realiz])
    data.extend([n_gal_total_sm])

    data.extend(num_agn_mag/n_gal_total_mag)
    data.extend(num_agn_sm/n_gal_total_sm)
    df_resample_sum.loc[len(df_resample_sum.index)] = data
    return df_resample_sum


def get_resample_reseed_summed(lc_catalogue, mask_xrlumi_dict, mag_cut, band,
                               lx_lim, lgM_smCut, mask_gal_agn=True, extra_frac=False):

    df_num_gal_agn = get_num_galAgn(lc_catalogue, mask_xrlumi_dict, mag_cut, band, lx_lim, lgM_smCut, 
                                    mask_gal_agn=mask_gal_agn, extra_frac=extra_frac)
    df_resample_sum = _get_resample_reseed_sum(df_num_gal_agn, mask_xrlumi_dict, lgM_smCut, mag_cut)
    return df_resample_sum
