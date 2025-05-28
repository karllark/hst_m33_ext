import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from astropy.wcs import WCS
from astropy.io import fits
from astropy.table import QTable
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.visualization import SqrtStretch, ImageNormalize


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--names", help="annotate with star names", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filename = (
        "~/Spitzer/M33/final_data/m33_mips24_map1234_kdg_19nov08.fits"
    )

    hdu = fits.open(filename)[0]
    data = hdu.data
    wcs = WCS(hdu.header)

    # make a smaller image
    coord = SkyCoord("01:33:49", "+30:38:35", unit=(u.hourangle, u.deg))
    cutout = Cutout2D(data, coord, (1.0 * u.deg, 0.75 * u.deg), wcs=wcs)
    data = cutout.data
    wcs = cutout.wcs

    fontsize = 14

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    ax = plt.subplot(projection=wcs)
    fig = plt.gcf()
    fig.set_size_inches(7.0, 8.5)

    norm = ImageNormalize(vmin=0.0, vmax=3.0, stretch=SqrtStretch())
    norm2 = ImageNormalize(vmin=0.0, vmax=20.0, stretch=SqrtStretch())

    # ax.set_figsize((10, 6))
    # ax.imshow(data, vmin=-0.1, vmax=2, origin='lower', cmap="binary")
    ax.imshow(data, origin="lower", norm=norm, cmap="binary")
    ax.grid(color="black", ls="dotted", alpha=0.5)
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")

    fnames = ["all"]
    fsyms = ["o"]
    fcols = ["b"]

    names = []
    for cname, csym, ccol in zip(fnames, fsyms, fcols):
        ptab = QTable.read(
            f"data/m33_positions_{cname}.dat", format="ascii.commented_header"
        )
        for k in range(len(ptab)):
            coord = SkyCoord(
                ptab["ra"][k],
                ptab["dec"][k],
                unit=(u.hourangle, u.deg),
            )
            ax.scatter(
                coord.ra.degree,
                coord.dec.degree,
                transform=ax.get_transform("fk5"),
                s=30,
                edgecolor=ccol,
                facecolor="none",
                linewidth=2,
                alpha=0.75,
                marker=csym,
            )

            # print(ptab["name"][k], wcs.world_to_pixel(coord))
            if args.names:
                ax.annotate(f" {ptab["name"][k]}", wcs.world_to_pixel(coord), color=ccol, alpha=0.75)

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            markerfacecolor="none",
            color="none",
            label="M33",
            markeredgecolor="b",
            markersize=10,
        ),
    ]
    # ax.legend(handles=legend_elements, loc="upper left", fontsize=0.8 * fontsize)

    ax.annotate(r"MIPS 24 $\mu$m", (100, 100))

    #ax.annotate(r"Wing", (3000, 500), fontsize=1.5*fontsize, alpha=0.5)
    #ax.annotate(r"Bar", (5000, 3800), fontsize=1.5*fontsize, alpha=0.5)

    plt.tight_layout()

    save_str = "figs/m33_mips24_positions"
    if args.png:
        plt.savefig(f"{save_str}.png")
    elif args.pdf:
        plt.savefig(f"{save_str}.pdf")
    else:
        plt.show()
