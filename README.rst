Code for HST M33 UV extinction project
======================================

Routines for HST M33 UV extinction curves program.
PI: G. Clayton

In Development!
---------------

Active development.
Everything changing.
Use at your own risk.

Contributors
------------
Karl Gordon

License
-------

This code is licensed under a 3-clause BSD style license (see the
``LICENSE`` file).

Extinction Curves
-----------------

Extinction curves created by running the `fitstars` bash script.  This fits the
STIS spectra and the available photometry (mostly HST) using `utils/fit_model.py`.
The foreground extinction is inlcuded in the fitting using the MW velocity integrated
HI columns converted to A(V) using the high-latitude N(HI)/A(V) ratio and a R(V) = 3.1
is assumed.

Most of the extinction curves are created relative to the F475W band as V band photometry
does not exist for them.   For these curves, `utils/convert_to_av_norm.py` is used to 
convert the curves to E(lambda - V) from E(lambda - F475W) using the fit R(V) value and 
the dust_extinction G23 R(V) dependent extinction model.

Figures
------- 

1. Locations of the stars on MIPS 24um image: `plotting/plot_positions_mips24.py --names`

2. Spectra stack: `plotting/plot_spec_stack.py`

3. Example fit: `meplot_model m33_e5_j013344.59+304436 --obspath ~/Python/extstar_data/M33/ --picmodname _modinfo.p``

4. Extinction stack: `plotting/plot_uvoptir_ext.py --rebin_fac 5`

5. A(V) vs R(V) and FM90 parmaeters versus each other: uses the the extinction_ensemble package
   with `plot_many_param_vs_param.py --fm90_noc1 --datasets GCC09 G03_lmc G24_smc C25_m31 G25_m33`

6. A(V) and main FM90 parameters versus N(HI)/A(V): uses the the extinction_ensemble package
   with `plot_many_param_vs_param.py --gdprops --datasets GCC09 G03_lmc G24_smc C25_m31 G25_m33`

7. Average M33 extinction curve: `plotting/plot_m33ave.py`

Tables
------

1. Names, positions, spectral types, and ids: By hand

2. Photometry: By Hand

3. Fitting parameters: By hand

4. Stellar parameters: utils/create_param_table.py

5. Column parameters: utils/create_param_table.py

6. FM90 parameters: utils/create_param_table.py. 
   For M33 average values, `utils/fit_ave_fm90`.