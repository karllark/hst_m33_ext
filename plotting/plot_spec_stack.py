import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import astropy.units as u

from measure_extinction.stardata import StarData


def plot_set(
    ax,
    starnames,
    extra_off_val=0.0,
    norm_wave_range=[0.26, 0.30] * u.micron,
    col_vals=["b", "g", "r", "m", "c", "y"],
    ann_wave_range=[0.3, 0.5] * u.micron,
    ann_rot=20.0,
    ann_offset=0.0,
    fontsize=18,
    path="/home/kgordon/Python/extstar_data/M33/",
    subpath="",
    rebin_res=None,
):
    """
    Plot a set of spectra
    """
    only_bands = None

    n_col = len(col_vals)
    for i in range(len(starnames)):
        stardata = StarData(
            subpath + starnames[i] + ".dat", path=path, only_bands=only_bands
        )

        # rebin spectra if desired
        if rebin_res is not None:
            gkeys = list(stardata.data.keys())
            gkeys.remove("BAND")
            for src in gkeys:
                stardata.data[src].rebin_constres(stardata.data[src].wave_range, rebin_res)

        # remove "bad" data
        gkeys = list(stardata.data.keys())
        gkeys.remove("BAND")
        for src in gkeys:
            medval = np.nanmedian(stardata.data[src].fluxes)
            bvals = stardata.data[src].fluxes < medval * 1e-1
            stardata.data[src].fluxes[bvals] = np.nan

        stardata.plot(
            ax,
            norm_wave_range=norm_wave_range,
            yoffset=extra_off_val + 0.5 * i,
            yoffset_type="multiply",
            pcolor=col_vals[i % n_col],
            annotate_key="BAND",
            annotate_wave_range=ann_wave_range,
            annotate_text=starnames[i].split("_")[1],
            fontsize=fontsize * 0.8,
            annotate_rotation=ann_rot,
            annotate_yoffset=ann_offset,
        )


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    fontsize = 18

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    # fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=(10, 13))
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 11))

    all = "m33_e2_j013334.26+303327  m33_e3_j013339.52+304540 m33_e4_j013341.93+304728 m33_e5_j013344.59+304436 m33_e6_j013406.63+304147 m33_e8_j013416.10+303344"

    starnames = np.array(all.split())

    # avs = [4.76, 5.53, 4.49, 5.29, 4.99, 4.45, 5.50, 4.97, 6.24, 6.79, 6.18, 6.17]
    # sindxs = np.argsort(avs)
    # starnames = starnames[sindxs]

    # read in the data and determine the optical spectral slope
    # then sort by that
    path = "/home/kgordon/Python/extstar_data/M33/"
    sslope = []
    for cstar in starnames:
        stardata = StarData(f"{cstar}.dat", path=path)
        gvals = (
            np.absolute(stardata.data["STIS"].waves - 0.15 * u.micron) < 0.02 * u.micron
        )
        bval = np.nanmedian(stardata.data["STIS"].fluxes[gvals])
        ctype = "STIS"
        gvals = (
            np.absolute(stardata.data[ctype].waves - 0.3 * u.micron) < 0.02 * u.micron
        )
        rval = np.nanmedian(stardata.data[ctype].fluxes[gvals])
        sslope.append(bval / rval)

    sindxs = np.argsort(sslope)
    starnames = starnames[sindxs]

    col_vals = ["b", "g"]
    plot_set(
        ax,
        starnames,
        ann_rot=-20.0,
        ann_offset=0.0,
        col_vals=col_vals,
        rebin_res=200,
    )

    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1e6)
    ax.set_xscale("log")
    # ax.set_xlim(kxrange)
    ax.set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)
    ax.set_ylabel(r"$F(\lambda)/F(0.28 \mu m)$ + offset", fontsize=1.3 * fontsize)

    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    # ax.spines["right"].set_visible(False)
    # ax.spines["top"].set_visible(False)

    ax.xaxis.set_minor_formatter(ScalarFormatter())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xticks([0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2], minor=True)

    ax.set_ylim(1e-2,2e3)

    fig.tight_layout()

    save_file = "figs/m33_ext_spec"
    if args.png:
        fig.savefig(save_file + ".png")
    elif args.pdf:
        fig.savefig(save_file + ".pdf")
    else:
        plt.show()
