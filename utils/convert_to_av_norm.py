import astropy.units as u
from dust_extinction.parameter_averages import G23
from measure_extinction.extdata import ExtData


if __name__ == "__main__":

    filename = "data/m33_exts.dat"

    f = open(filename, "r")
    file_lines = list(f)
    starnames = []
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            starnames.append(name)

    emod = G23()
    emod.Rv_range = [1.5, 7.0]
    for cname in starnames:
        cfile = f"exts/{cname}_mefit_ext.fits"
        edata = ExtData(filename=cfile)

        emod.Rv = edata.columns["RV"][0]

        if edata.type_rel_band != "V":
            print(cname, edata.columns["RV"][0], edata.type_rel_band)
            if edata.type_rel_band == "ACS_F814W":
                refwave_av = emod(0.802932 * u.micron)
            else:
                refwave_av = emod(0.477217 * u.micron)
            av = edata.columns["AV"][0]

            for ckey in edata.exts.keys():
                edata.exts[ckey] += (refwave_av - 1.0) * av

            edata.type_rel_band = "V"

        edata.save(cfile.replace("_ext.fits", "_ext_elv.fits"))