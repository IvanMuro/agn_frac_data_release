import copy
import logging

#from astropy import units as u
from colossus.halo import mass_defs
import numpy as np
import pandas as pd

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


#def calc_num_gal_and_agn(id_cluster, angle, delta_z, mask_xrlumi_dict, cube, field,
#                         lc=False, loop_lc=False, z_object=False, mag_lim=False, mask_mag=False,
#                         sm_lim=False, mask_sm=False, mask_gal_agn=True):
#    """
#    Function that calculate the number of galaxies/AGNs for a given selection effects:
#    magnitude cut, stellar mass cut (galaxies + AGNS), luminosity cut (AGNs). 
#    """
#    results_dict = {}
#
#    # GALAXIES
#    if lc == False:
#        mask_spatial = _get_spatial_mask(
#            id_cluster, angle, delta_z, cube, field)
#        n_gal = _calc_num_gal(cube, mask_mag, mask_spatial)
#
#    else:
#        #mask_spatial = np.array(np.abs(lc ['final_redshift'] - z_object) < delta_z)
#        mask_spatial = _get_spatial_mask(id_cluster, angle, delta_z, lc, field)
#        #   (len(lc [np.logical_and(mask_spatial, mask_mag)]['final_redshift']))
#        if mag_lim != False:
#            n_gal = _calc_num_gal(lc, mask_mag, mask_spatial)
#            results_dict.update({'n_gal_mag': n_gal})
#        if sm_lim != False:
#            n_gal = _calc_num_gal(lc, mask_sm, mask_spatial)
#            results_dict.update({'n_gal_sm': n_gal})
#    # AGNS
#    if lc == False:
#        for mask_xrlumi in mask_xrlumi_dict:
#            model = f'{mask_xrlumi}'
#            n_agn = _calc_num_agn(cube, mask_xrlumi_dict[model], mask_spatial,
#                                  mask_gal_agn=mask_gal_agn)
#            results_dict.update({f'n_agn_{model}': n_agn})
#    else:
#        for mask_xrlumi in mask_xrlumi_dict:  # loop through dict. keys
#            model = f'{mask_xrlumi}'
#            if mag_lim != False:
#                n_agn_mag = _calc_num_agn(lc, mask_xrlumi_dict[model], mask_spatial,
#                                          mask_mag, mask_gal_agn=mask_gal_agn)
#                results_dict.update({f'n_agn_mag_{model}': n_agn_mag})
#            if sm_lim != False:
#                n_agn_sm = _calc_num_agn(lc, mask_xrlumi_dict[model], mask_spatial,
#                                         mask_sm, mask_gal_agn=mask_gal_agn)
#                results_dict.update({f'n_agn_sm_{model}': n_agn_sm})
#
#    return results_dict
#
#
#def get_num_gal_agn_df(id_clust_list, mask_xrlumi_dict, angle, delta_z, cube, field=False,
#                       sample_size='default', lc_list=False, z_object=False,
#                       mag_lim=False, mask_mag=False, sm_lim=False, mask_sm=False,
#                       lim_lum=False, z=False, write_ngal_nagn=True,
#                       mask_gal_agn=True):  # delete this last line only to write in the file the num of agns
#    """
#    Function that return a pandas.DataFrame with the number of galaxies and AGNs (for given
#    models) for an input cluster list.
#
#    Input:
#    ------
#    id_clust_list (list of int): list of cluster to be analyzed
#    mask_xrlumi_dict (list of mask): X-ray luminosity masks
#    angle (float, deg):  used to constrain the selection effects. See selection_effects.py
#    delta_z (float): used to constrain the selection effects. See selection_effects.py
#    cube (box object): see cube.py
#    sample_size(int): number of cluster that want to be analyzed, if no value is given by default it will 
#    analyze all the cluster on the id list
#
#    Output:
#    -------
#    df_num_gal_agn (pandas data frame): data frame that constains ids, num galaxies and num ags for different
#                                        luminosity seeding models
#    """
#    logger.info(
#        f'Getting number of galaxies and AGNs for {sample_size} clusters')
#
#    if sample_size == 'default':
#        sample_size = len(id_clust_list)
#
#    # if field == False:
#    df_num_gal_agn = _create_df_num_gal_agn_df(id_clust_list, mask_xrlumi_dict[0], sample_size,
#                                               mag_lim=mag_lim, sm_lim=sm_lim, field=field)
#    # else:
#    #    df_num_gal_agn = _create_df_num_gal_agn_df(np.arange(0,len(lc_list)), mask_xrlumi_dict[0], sample_size,
#    #                                               mag_lim=mag_lim, sm_lim=sm_lim, field=field)
#
#    if lc_list == False:
#        for loop_clust, id_clust in enumerate(df_num_gal_agn['id']):
#            results_dict = calc_num_gal_and_agn(id_clust, angle, delta_z, mask_xrlumi_dict,
#                                                mask_mag, cube, field, mask_gal_agn=mask_gal_agn)
#            _add_ngal_nagn_df(id_clust, df_num_gal_agn,
#                              results_dict, field=field)
#    else:
#        for loop_lc, lc in enumerate(lc_list):
#            mask_mag_aux, mask_sm_aux = _get_mask_mag_sm_aux(
#                loop_lc, mag_lim, sm_lim, mask_mag, mask_sm)
#
#            results_dict = calc_num_gal_and_agn(lc.id_obs, lc.angle, delta_z,
#                                                mask_xrlumi_dict[loop_lc],
#                                                cube, field, lc=lc, loop_lc=loop_lc, z_object=z_object,
#                                                mag_lim=mag_lim, mask_mag=mask_mag_aux, sm_lim=sm_lim,
#                                                mask_sm=mask_sm_aux, mask_gal_agn=mask_gal_agn)
#
#            results_dict_true_members = _calc_num_gal_and_agn_true_members(loop_lc, lc,
#                                                                           mask_xrlumi_dict,
#                                                                           delta_z, cube, field,
#                                                                           z_object, mag_lim,
#                                                                           mask_mag_aux, sm_lim,
#                                                                           mask_sm_aux,
#                                                                           mask_gal_agn=mask_gal_agn)
#            if field == False:
#                _add_ngal_nagn_df(lc.id_obs, df_num_gal_agn,
#                                  results_dict, field)
#                _add_ngal_nagn_df(
#                    lc.id_obs, df_num_gal_agn, results_dict_true_members, field, true_members=True)
#            else:
#                _add_ngal_nagn_df(loop_lc, df_num_gal_agn, results_dict, field)
#    if write_ngal_nagn == True:
#        mass_dm, sigma_clust = _get_mass_sigma_parent_sample(
#            cube, id_clust_list)
#        _write_ngal_nagn_file(df_num_gal_agn, lim_lum,
#                              field, z, mass_dm, sigma_clust)
#    return df_num_gal_agn
#
#
#def _calc_num_gal_and_agn_true_members(loop_lc, lc, mask_xrlumi_dict, delta_z, cube, field,
#                                       z_object, mag_lim, mask_mag_aux, sm_lim, mask_sm_aux,
#                                       mask_gal_agn=True):
#    """
#    Function that calculate the number of galaxies/AGNs for a given selection effects
#    using the main function calc_num_gal_and_agn(...). We is done is to mask the general
#    lc and its mask to contain only cluster members using their upids
#    """
#    mask_true_mem = np.logical_or(np.isin(lc['upid'], lc.id_obs),
#                                  np.isin(lc['id'], lc.id_obs))
#    lc_true_mem = copy.deepcopy(lc)
#    lc_true_mem = lc_true_mem[mask_true_mem]
#    mask_xrlumi_dict_true_mem = copy.deepcopy(mask_xrlumi_dict[loop_lc])
#    for mask_xrlumi in mask_xrlumi_dict_true_mem:
#        model = f'{mask_xrlumi}'
#        mask_xrlumi_dict_true_mem[model] = mask_xrlumi_dict_true_mem[model][mask_true_mem]
#    results_dict_true_members = calc_num_gal_and_agn(lc.id_obs, lc.angle, delta_z,
#                                                     mask_xrlumi_dict_true_mem,
#                                                     cube, field, lc=lc_true_mem,
#                                                     loop_lc=loop_lc, z_object=z_object,
#                                                     mag_lim=mag_lim, mask_mag=mask_mag_aux[mask_true_mem],
#                                                     sm_lim=sm_lim,
#                                                     mask_sm=mask_sm_aux[mask_true_mem],
#                                                     mask_gal_agn=mask_gal_agn)
#    return results_dict_true_members
#
#
#def _get_mask_mag_sm_aux(loop_lc, mag_lim, sm_lim, mask_mag, mask_sm):
#    """
#    Function that returns the corresponding mask of magnitude or stellar mass cut, in the case
#    there is one. In other case it returns False.
#    """
#    if mag_lim != False:
#        mask_mag_aux = mask_mag[loop_lc]
#    else:
#        mask_mag_aux = False
#
#    if sm_lim != False:
#        mask_sm_aux = mask_sm[loop_lc]
#    else:
#        mask_sm_aux = False
#    return mask_mag_aux, mask_sm_aux
#
#
#def _get_mass_sigma_parent_sample(cube, id_clust_list):
#    mask_clust = np.isin(cube['id'], id_clust_list)
#    mass_dm = np.array(cube[mask_clust]['lgm'])
#    mass_dm_200 = np.log10(mass_defs.changeMassDefinition(
#        10**mass_dm, 10**.7, cube.z, 'vir', '200m')[0]/cube.cosmo.h)
#    sigma_clust = select_eff.get_sigma_cluster(cube.cosmo,
#                                               mass_dm_200,
#                                               cube.z)
#    return mass_dm, sigma_clust
#
#
#def _write_ngal_nagn_file(df_num_gal_agn, lim_lum, field, z, mass_dm, sigma_clust):
#    """
#    Function that writes the number of galaxies/agns w different selction effects in the input
#    parent cluster sample.
#    """
#    logger.info(f'Writting ngal, nagn: z={z}, lx>{lim_lum}')
#    df_num_gal_agn['lgm'] = mass_dm
#    df_num_gal_agn['sigma'] = sigma_clust
#    if field == False:
#        type = 'clust'
#    else:
#        type = 'field'
#    df_num_gal_agn.to_csv('/storage/ivan/mock_univ_output/MDP2/martini/num_galagn_lc/' +
#                          f'lc_numgal_numagn_lxgt{lim_lum}_{type}_z_{z}.csv')
#    return
#
#
#def _create_df_num_gal_agn_df(id_clust_list, mask_xrlumi_dict, sample_size,
#                              mag_lim=False, sm_lim=False, field=False):
#    """
#    Creates (and fill w NaNs) a dataFrame w the columns = id, n_gal_mag, n_gal_sm,
#    n_agn_mag_[model], n_agn_sm_[model]. This is used after to save the number of 
#    galaxies and AGNs w different selection effects
#    Input:
#    ------
#    id_clust_list (list): list of cluster ids
#    mask_xrlumi_dict (dict): dictation w keys names of the model and value 
#                             the mask for the particular cube or lc
#    sample_size (int): number of clusters that want to be used, maximum len(id_clust_list)
#    mag_lim (float): magnitude limit
#    sm_lim (float): stellar mass limit 
#    field (bool): if the lc observes a cluster set this to False, True otherwhise 
#
#    Output:
#    -------
#    df_num_gal_agn (pandas DataFrame)
#    """
#
#    if field == False:
#        df_num_gal_agn = pd.DataFrame({'id': id_clust_list[:sample_size]})
#    else:
#        id_clust_list = np.arange(0, len(sample_size))
#        df_num_gal_agn = pd.DataFrame({'id': id_clust_list[:sample_size]})
#
#    df_num_gal_agn['id'] = df_num_gal_agn['id'].astype('int64')
#
#    if mag_lim != False:
#        df_num_gal_agn['n_gal_mag'] = np.nan
#        df_num_gal_agn['n_gal_mag_true_mem'] = np.nan
#    if sm_lim != False:
#        df_num_gal_agn['n_gal_sm'] = np.nan
#        df_num_gal_agn['n_gal_sm_true_mem'] = np.nan
#    for mask_xrlumi in mask_xrlumi_dict:
#        df_num_gal_agn[f'n_agn_mag_{mask_xrlumi}'] = np.nan
#        df_num_gal_agn[f'n_agn_mag_{mask_xrlumi}_true_mem'] = np.nan
#        df_num_gal_agn[f'n_agn_sm_{mask_xrlumi}'] = np.nan
#        df_num_gal_agn[f'n_agn_sm_{mask_xrlumi}_true_mem'] = np.nan
#    return df_num_gal_agn
#
#
#def _add_ngal_nagn_df(id_clust, df, results_dict, field=False, true_members=False):
#    """
#    Function that adds the number of galaxies  for a given pandas.DataFrame
#    """
#    for result in results_dict:
#        if field == False:
#            if true_members == False:
#                df.loc[df['id'] == id_clust,
#                       f'{result}'] = results_dict[result]
#            else:
#                df.loc[df['id'] == id_clust,
#                       f'{result}_true_mem'] = results_dict[result]
#        else:
#            df.loc[df['id'] == id_clust, result] = results_dict[result]
#    return
#
#
#def get_resample_agns_sneaky(df_num_gal_agn, mask_xrlumi_dict, num_realizations='default'):
#    """
#    Function that resample a given data frame with ids, num gal and num agns. Used to calculate error assotiated
#    to the mean. 
#    Input:lim_lum=False, type=False, z=z
#    """
#    logger.info(f'Resampling clusters with {num_realizations} realizations')
#    rand_realiz = _get_random_realizations(df_num_gal_agn, num_realizations)
#    n_gal_realizations = []
#    n_agn_realizations = {}
#
#    for loop_realization, realization in enumerate(rand_realiz):
#        # Sneacky way of selecting by ids
#        unique_realiz = np.unique(sorted(realization), return_counts=True)
#        unique_ids = unique_realiz[0]
#        ids_repetitions = unique_realiz[1]
#
#        df_resample = df_num_gal_agn[df_num_gal_agn['id'].isin(realization)]
#        n_gal_realizations.append(
#            np.array(df_resample['n_gal']*ids_repetitions))
#
#        for mask_xrlumi in mask_xrlumi_dict[0]:
#            if loop_realization == 0:
#                n_agn_realizations.update({f'n_agn_{mask_xrlumi}':
#                                           [np.array(df_resample[f'n_agn_{mask_xrlumi}']*ids_repetitions)]})
#            else:
#                n_agn_aux = n_agn_realizations[f'n_agn_{mask_xrlumi}']
#                n_agn_aux.append(
#                    np.array(df_resample[f'n_agn_{mask_xrlumi}']*ids_repetitions))
#                n_agn_realizations.update({f'n_agn_{mask_xrlumi}': n_agn_aux})
#
#        # Below comment was supposed to reproduce the above loop, it doesnt (:
#        # n_agn_realizations = [ _update_agn_num(loop_realization, n_agn_realizations, mask_xrlumi, ids_repetitions, df_resample)
#        # for mask_xrlumi in mask_xrlumi_dict ]
#    return n_gal_realizations, n_agn_realizations
#
#
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

    logger.info(f'Resampling clusters with {num_realizations} realizations')

    rand_realiz = _get_random_realizations(
        df_num_gal_agn, num_realizations, id_label)

    n_gal_realizations_mag = []
    n_gal_realizations_sm = []

    n_agn_realizations_mag = {}
    n_agn_realizations_sm = {}
    for loop_realization, realization in enumerate(rand_realiz):
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
#
#
#def _update_agn_num(loop_realization, n_agn_realizations, mask_xrlumi, ids_repetitions, df_resample):
#    if loop_realization == 0:
#        n_agn_realizations.update({f'n_agn_{mask_xrlumi}':
#                                  [np.array(df_resample[f'n_agn_{mask_xrlumi}']*ids_repetitions)]})
#    else:
#        n_agn_aux = n_agn_realizations[f'n_agn_{mask_xrlumi}']
#        n_agn_aux.append(
#            np.array(df_resample[f'n_agn_{mask_xrlumi}']*ids_repetitions))
#        n_agn_realizations.update({f'n_agn_{mask_xrlumi}': n_agn_aux})
#
#    return n_agn_realizations
#
#

def _get_random_realizations(df_num_gal_agn, num_realizations='default', id_label='id'):
    if num_realizations == 'default':
        num_realizations = len(df_num_gal_agn[id_label])

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
#
#
#def get_resample_summed(id_clust_list, cube, delta_z, angle, field=False, sample_size='default',
#                        num_realizations='default', lc_list=False, lx_lim=False, sm_lim=False,
#                        z_object=False, mag_lim=False, mag_filter='IRAC1_mine', mask_xrlumi_dict=False,
#                        mask_mag=False, mask_gal_agn=True):  # moved to utils
#    """
#    If there is a list of lc first create a list of masks (magnitude, sm and luminosity). 
#    Then this masks are filled for each lc.
#    Input:
#    ------
#    id_clust_list (int or list of int)
#    cube (cube object): see cube.py
#    delta_z (float or list)
#    angle (float): angle in deg
#    field (bool): False is the observation is a cluster, True otherwhise
#    sample_size (int): Number of elements in the resample
#    num_realizations (int): Number of times re-sample is repeated
#    lc_list (list of objects type LightCone): see light_cone.py
#    lx_lim (False or float or list of float): **LOG10** luminosity limit
#    sm_lim (False or float or list of float): **Log10** stellar mass limit
#    z_object (float): redshift of the cluster:
#    mag_lim (False or float or list of float):
#    mag_filter (str or list):
#    mask_xrlumi_dict (False or dictation):
#    mask_mag (False or list of bools):
#
#    Output:
#    -------
#    df_resample_sum (pandas DataFrame): table with len=num_realization summed over the sample_size elements
#    """
#    if lc_list != False:
#        mask_mag, mask_sm_list, mask_xrlumi_dict = _get_masks_mag_sm_lumi(
#            mag_lim, sm_lim, lc_list, mag_filter, lx_lim)
#    df_num_gal_agn = get_num_gal_agn_df(id_clust_list, mask_xrlumi_dict, angle,
#                                        delta_z, cube, sample_size=sample_size,
#                                        lc_list=lc_list, z_object=z_object, field=field,
#                                        mag_lim=mag_lim, mask_mag=mask_mag,
#                                        sm_lim=sm_lim, mask_sm=mask_sm_list,
#                                        lim_lum=lx_lim, z=z_object,
#                                        mask_gal_agn=mask_gal_agn)
#
#    n_gal_realizations_mag, n_gal_realizations_sm, n_agn_realizations_mag, n_agn_realizations_sm = \
#        get_resample_agns(df_num_gal_agn, mask_xrlumi_dict,
#                          num_realizations=num_realizations,
#                          mag_lim=mag_lim, sm_lim=sm_lim)
#
#    df_resample_sum = _get_df_resample_sum(mask_xrlumi_dict, mag_lim, sm_lim, n_gal_realizations_mag, n_agn_realizations_mag,
#                                           n_gal_realizations_sm, n_agn_realizations_sm)
#
#    return df_resample_sum
#
#
#def _get_masks_mag_sm_lumi(mag_lim, sm_lim, lc_list, mag_filter, lx_lim):
#    mask_xrlumi_dict = []
#
#    if mag_lim != False:
#        mask_mag = []
#    else:
#        mask_mag = [False]
#
#    if sm_lim != False:
#        mask_sm = []
#    else:
#        mask_sm = [False]
#
#    for loop_lc, lc in enumerate(lc_list):  # fill masks
#        if mag_lim != False:
#            # print(loop_lc)
#            mask_mag.append(select_eff.mask_magnitude(mag_lim, mag_filter, lc))
#        if sm_lim != False:
#            mask_sm.append(select_eff.mask_mass(sm_lim, lc))
#
#        # mask_xrlumi_dict.append({'aird': np.array( lc ['lgLxAMedian'] > lx_lim ),
#        #                         'age' : np.array( lc ['lgLxGMedian'] > lx_lim )})
#        mask_aux = {}
#        for key in mask_xrlumi_dict:
#            mask_aux.update({f'{key}': (np.array(lc[f'lgLx{key}'] > lx_lim))})
#        mask_xrlumi_dict = mask_aux
#    return mask_mag, mask_sm, mask_xrlumi_dict
#
#
#def _get_df_resample_sum(mask_xrlumi_dict, mag_lim, sm_lim, n_gal_realizations_mag, n_agn_realizations_mag,
#                         n_gal_realizations_sm, n_agn_realizations_sm):
#    """
#    *** NEED TO IMPROVE FUTURE VERSION ***
#    These needs to be improved including data types, the problem is that I use extend so I cannot give types 
#    as it is now because at some point I divide a list of ints between error, but list cannot be divided. Probably
#    I need to use numpyarray, but cannot use extend anymore since np.arrays are fixed structures. The idea would be 
#    to create an array w the proper dimension since the beggining and fill it by positions.
#    """
#
#    df_resample_sum = _create_df_resample(
#        mask_xrlumi_dict[0], mag_lim=mag_lim, sm_lim=sm_lim)
#    logger.info('Summing resample')
#    for loop_realiz, ngal in enumerate(n_gal_realizations_mag):
#        num_agn_mag = [sum(n_agn_realizations_mag[f'n_agn_mag_{mask_xrlumi}'][loop_realiz])
#                       for mask_xrlumi in mask_xrlumi_dict[0]]
#        data = copy.deepcopy(num_agn_mag)
#        num_agn_sm = [sum(n_agn_realizations_sm[f'n_agn_sm_{mask_xrlumi}'][loop_realiz])
#                      for mask_xrlumi in mask_xrlumi_dict[0]]
#        data.extend(num_agn_sm)
#
#        n_gal_total_mag = sum(ngal)
#        data.extend([n_gal_total_mag])
#        n_gal_total_sm = sum(n_gal_realizations_sm[loop_realiz])
#        data.extend([n_gal_total_sm])
#
#        data.extend(num_agn_mag/n_gal_total_mag)
#        data.extend(num_agn_sm/n_gal_total_sm)
#        df_resample_sum.loc[len(df_resample_sum.index)] = data
#
#    return df_resample_sum
#
#
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


#def _get_lc_concat_filter_spt(lc_catalogue):
#    """
#    Function that filter spatially a lc concat catalogue
#    """
#    d_ang_cosmo = lc_catalogue.cosmo.angularDiameterDistance(
#        lc_catalogue.z) * (u.Mpc/u.rad)
#    dist_ang = (
#        ((np.array(lc_catalogue['r_parent'])*u.kpc)/d_ang_cosmo).to(u.deg)).value
#    #dist_ang = select_eff._calc_angle_from_r_qties(lc_catalogue.cosmo, lc_catalogue.z, lc_catalogue ['r'])
#    #mask_ang = np.sqrt(np.array(lc_catalogue['delta_RA'])**2 +
#    #                   np.array(lc_catalogue['delta_dec'])**2) < dist_ang
#    mask_ang = lc_catalogue['angle_los']< dist_ang
#    mask_spec = np.array(np.abs(lc_catalogue['delta_redshift']) <
#                         lc_catalogue['delta_redshift_lim'])
#    mask_spat = np.logical_and(mask_ang, mask_spec)
#    lc_catalogue_concat_mask_spt = copy.deepcopy(lc_catalogue)
#    lc_catalogue_concat_mask_spt = lc_catalogue_concat_mask_spt[mask_spat]
#
#    return lc_catalogue_concat_mask_spt


#def duplicate_lc_lst(lc_catalogue, n_duplication):
#    """
#    Function that duplicates a data frame an input number times
#    """
#
#    lc_catalogue_concat_mask_spt = _get_lc_concat_filter_spt(lc_catalogue)
#
#    for reseed_loop in range(n_duplication):
#        if reseed_loop == 0:
#            lc_catalogue_concat_mask_spt_duplicate = copy.deepcopy(
#                lc_catalogue_concat_mask_spt)
#            lc_catalogue_concat_mask_spt_duplicate['reseed'] = reseed_loop
#            id_clust_reseed = [f'{id_obs}_{reseed_loop}'
#                               for id_obs in lc_catalogue_concat_mask_spt_duplicate['id_obs']]
#            lc_catalogue_concat_mask_spt_duplicate['id_clust_reseed'] = id_clust_reseed
#
#        else:
#            lc_catalogue_aux = copy.deepcopy(lc_catalogue_concat_mask_spt)
#            lc_catalogue_aux['reseed'] = reseed_loop
#            id_clust_reseed = [f'{id_obs}_{reseed_loop}'
#                               for id_obs in lc_catalogue_aux['id_obs']]
#            lc_catalogue_aux['id_clust_reseed'] = id_clust_reseed
#
#            lc_catalogue_concat_mask_spt_duplicate = \
#                lc_catalogue_concat_mask_spt_duplicate.append(
#                    lc_catalogue_aux,
#                    ignore_index=True)
#    return lc_catalogue_concat_mask_spt_duplicate


def _create_df_reseed_num_gal_agn_df(lc_catalogue, mask_xrlumi_dict):
    """
    Function that creates the dataframe where the number of gal/agn are saved after
    """
    df_num_gal_agn = pd.DataFrame(
        {'id_cluster_resed': lc_catalogue['id_clust_reseed'].unique()})
    df_num_gal_agn[['id', 'resed']] = df_num_gal_agn['id_cluster_resed'].str.split(
        '_', 1, expand=True)

    df_num_gal_agn['n_gal_mag'] = np.nan
    df_num_gal_agn['n_gal_sm'] = np.nan
    for mask_xrlumi in mask_xrlumi_dict:
        df_num_gal_agn[f'n_agn_mag_{mask_xrlumi}'] = np.nan
        df_num_gal_agn[f'n_agn_sm_{mask_xrlumi}'] = np.nan
    return df_num_gal_agn


def _calc_num_gal_agn_lc_catalogue(df_num_gal_agn, lc, id_cluster_resed, m_cut,
                                   band, lx_lim, lgM_smCut, mask_gal_agn=True,
                                   mask_xrlumi_dict=False,
                                   extra_frac=False):
    """
    """
    mask_mag = select_eff.mask_magnitude(m_cut, band, lc)
    mask_sm = select_eff.mask_mass(lgM_smCut, lc)

    mask_aux = {}
    for key in mask_xrlumi_dict:
        mask_aux.update({f'{key}': (np.array(lc[f'lgLx{key}'] > lx_lim))})
    mask_xrlumi_dict = mask_aux

    if m_cut != False:
        df_num_gal_agn.loc[df_num_gal_agn['id_cluster_resed'] == id_cluster_resed,
                           'n_gal_mag'] = len(lc[mask_mag])

    if lgM_smCut != False:
        df_num_gal_agn.loc[df_num_gal_agn['id_cluster_resed'] == id_cluster_resed,
                           'n_gal_sm'] = len(lc[mask_sm])

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
    df_num_gal_agn = _create_df_reseed_num_gal_agn_df(
        lc_catalogue, mask_xrlumi_dict)
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

    df_num_gal_agn = get_num_galAgn(lc_catalogue, mask_xrlumi_dict, mag_cut, band,
                                    lx_lim, lgM_smCut, mask_gal_agn=mask_gal_agn,
                                    extra_frac=extra_frac)
    df_resample_sum = _get_resample_reseed_sum(
        df_num_gal_agn, mask_xrlumi_dict, lgM_smCut, mag_cut)
    return df_resample_sum
