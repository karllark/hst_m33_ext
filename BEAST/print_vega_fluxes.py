from beast.observationmodel.vega import Vega


if __name__ == '__main__':
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']

    with Vega() as v:
        _, vega_flux, _ = v.getFlux(filters)

    for cfilt, cvflux in zip(filters, vega_flux):
        print(cfilt, cvflux)
