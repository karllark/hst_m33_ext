import numpy as np

from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u


def get_postage_catalog(radecstr,
                        size,
                        catalog,
                        realcatalog):
    """
    Extract a postage stamp catalog from a larger catalog

    Parameters
    ----------
    radecstr : str
        SkyCoord pareseable string for ra, dec coordinates

    size : 2 element tuple
        size in ra and dec in arcsec

    catalog : SkyCoord catalog
        all the ra, decs from the catalog in SkyCoord format

    realcatalog : astropy Table
        full catalog (all columns)
    """
    position = SkyCoord(radecstr)

    d2d = position.separation(catalog)
    catalogmsk = d2d < size
    idxcatalog = np.where(catalogmsk)[0]

    if len(idxcatalog):
        out_table = realcatalog[catalogmsk]
    else:
        out_table = None

    return out_table


if __name__ == '__main__':

    # get the large fits file and wcs info
    filenames = ['/home/kgordon/HST/M33-Ext/14610_M33-B01_1.phot.fits.gz',
                 '/home/kgordon/HST/M33-Ext/14610_M33-B01_2.phot.fits.gz',
                 '/home/kgordon/HST/M33-Ext/merged-6filt-bricks23.fits.gz']
    ftags = ['B01_2', 'B01_1', 'B23']

    # h, m, s have to be lower case for radecstr.  That is specific
    # maybe check in astropy why this is the case
    starnames = ['J013415.71+303341_V01',
                 'J013358.07+303308_V02',
                 'J013352.13+304319_V03',
                 'J013344.59+304436_V04',
                 'J013406.63+304147_V05',
                 'J013339.52+304540_V06',
                 'J013410.59+304616_V07',
                 'J013334.26+303327_V08']
    radecstrs = ['01h34m15.6800s +30d33m40.90s',
                 '01h33m58.0700s +30d33m08.20s',
                 '01h33m52.1300s +30d43m19.30s',
                 '01h33m44.5900s +30d44m36.90s',
                 '01h34m06.6300s +30d41m47.80s',
                 '01h33m39.5200s +30d45m40.50s',
                 '01h34m10.5900s +30d46m16.10s',
                 '01h33m34.2600s +30d33m27.60s']

    # new set
    starnames = ['J013341.93+304728_VX',
                 'J013306.53+303010_VX',
                 'J013340.10+304138_VX',
                 'J013313.76+305240_VX']

    radecstrs = ['01h33m41.93s +30d47m28.3s',
                 '01h33m06.53s +30d20m10.8s',
                 '01h33m40.10s +30d41m38.4s',
                 '01h33m13.76s +30d52m40.1s']

    # size to extract in (ra, dec) arcsec
    size = (10./3600.)*u.degree

    for k, cfile in enumerate(filenames):

        bigcat = Table.read(cfile)

        bigcat_skycoord = SkyCoord(ra=bigcat['ra']*u.degree,
                                   dec=bigcat['dec']*u.degree)

        for i, cradec in enumerate(radecstrs):
            pcat = get_postage_catalog(cradec, size, bigcat_skycoord, bigcat)
            print(starnames[i])
            print(pcat)

            if pcat is not None:
                pcat.write('%s_%s_pcatalog.fits' % (starnames[i], ftags[k]),
                           overwrite=True)
