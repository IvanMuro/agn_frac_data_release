# Cosmic evolution of the incidence of Active Galactic Nuclei in massive clusters: Simulations versus observations

This repository provides the code that reproduces the results presented in the paper "Cosmic evolution of the incidence of Active Galactic Nuclei in massive clusters: Simulations versus observations" by Munoz-Rodriguez et al. (2022, MNRAS accepted, https://doi.org/10.1093/mnras/stac3114, https://doi.org/10.48550/arXiv.2211.00032). The notebook paper_results.ipynb demonstrates how to reproduce the plots on the fraction of AGN in clusters of galaxies presented in the above paper. Running this code requires the simulated light-cones of massive clusters stored in the hdf5 files


lc_resample_reseed_0.4505.h5

lc_resample_reseed_0.5747.h5

lc_resample_reseed_0.8376.h5


These are accessible via zenodo https://zenodo.org/record/7266395. These files should be copied to the "./input" directory of the repository. Each of these files contains a collection of light cones of massive halos (typically few times 10^14 solar) selected from the MDPL2 (MultiDark PLanck-2; https://ui.adsabs.harvard.edu/abs/2016MNRAS.457.4340K/abstract) N-body simulation boxes with scale parameters 0.4505, 0.5747 and 0.8376 that approximately correspond to redshifts 0.2, 0.75 and 1.25. These halos and their associated galaxies (via UniverseMachine; https://ui.adsabs.harvard.edu/abs/2020MNRAS.499.5702B/abstract) are populated with AGN following the prescriptions described in Munoz-Rodriguez et al. (2022). They are then projected on the sky to produce mock observations that mimic the selection of massive clusters, their member galaxies and AGN described in Martini et al. 2009, 2013. For each MDPL2 massive halo multiple (total 10) realisations of the AGN seeding are produced. The realisations are then used to estimate the AGN fraction averaged over all massive halos and their corresponding realisations in each light-cone.  

The README.pdf document provides information on the columns included in the light cone hdf5 tables.

Contribution to the code: Ivan Munoz Rodriguez, Antonis Georgakakis and Angel Ruiz Camunas.
