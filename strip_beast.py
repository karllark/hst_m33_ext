# strip the BEAST SED file down to just want is needed
#  photometry SEDs and associated physical parameters

import numpy as np
import h5py
from astropy.table import Table, Column

if __name__ == '__main__':

    # define the filters
    basefilters = ['F275W', 'F336W', 'F475W',
                   'F814W', 'F110W', 'F160W']

    # read in the full BEAST sed file
    fname = 'beast_example_phat_seds.grid.hd5'
    hdf_cache = h5py.File(fname, 'r')

    # for sname in hdf_cache.keys():
    #     print(sname)
    #     if hdf_cache[sname].dtype.fields is None:
    #         print(hdf_cache[sname].value)
    #     else:
    #         print(hdf_cache[sname].dtype.fields.keys())

    # build the new table
    a = Table()

    # only save SEDs for only hot, "lightly" reddened stars
    indxs, = np.where(np.logical_and(
        hdf_cache['grid'].value['logT'] >= 4.0,
        hdf_cache['grid'].value['Av'] <= 3.0))
    print(len(indxs))

    # model parameters
    model_parameters = ['logT', 'logg', 'Z', 'logA', 'M_ini', 'logL',
                        'Av', 'Rv', 'f_A']
    for cparam in model_parameters:
        a[cparam] = Column(hdf_cache['grid'].value[cparam][indxs])

    # SEDs
    for k, cname in enumerate(basefilters):
        a[cname] = Column(hdf_cache['seds'].value[:, k][indxs])

    # write
    a.write('m33_beast_sed_grid_glogT4_lAv3.fits.gz', overwrite=True)
