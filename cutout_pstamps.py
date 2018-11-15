from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import wcs
from astropy.wcs.utils import proj_plane_pixel_scales


def get_postage_stamp(radecstr,
                      size,
                      image,
                      image_wcs):
    """
    Extract a postage stamp image from a larger image

    Parameters
    ----------
    radecstr : str
        SkyCoord pareseable string for ra, dec coordinates

    size : 2 element tuple
        size in ra and dec in arcsec

    image : 2D numpy.ndarray
        image to cut postage stamps from

    image_wcs : astropy.wcs
        WCS for image
    """
    position = SkyCoord(radecstr)
    pix_scales = 3600.*proj_plane_pixel_scales(image_wcs)
    pix_size = (size[0]/pix_scales[0], size[1]/pix_scales[1])

    return Cutout2D(image, position, pix_size, wcs=image_wcs)


if __name__ == '__main__':

    # get the large fits file and wcs info
    filename = 'F475W_mosaic.fits.gz'
    hdulist = fits.open(filename)
    image_wcs = wcs.WCS(hdulist[0].header)
    image = hdulist[0].data

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

    # size to extract in (ra, dec) arcsec
    size = (10., 10.)

    for k, cradec in enumerate(radecstrs):
        pimage = get_postage_stamp(cradec, size, image, image_wcs)
        print(pimage.wcs)

        ohdu = fits.PrimaryHDU(pimage.data, header=pimage.wcs.to_header())
        ohdul = fits.HDUList([ohdu])
        ohdul.writeto('%s_pstamp.fits' % starnames[k], overwrite=True)
