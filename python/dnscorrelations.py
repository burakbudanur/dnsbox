#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.lib.npyio import save
from tqdm.auto import tqdm

import dns

cmap = "coolwarm"


def main():

    parser = argparse.ArgumentParser(
        description="Produce 2D visualizations of the velocity field.",
    )
    parser.add_argument(
        "statesPath",
        type=str,
        help="path to the folder containing states of interest.",
    )
    parser.add_argument(
        "si", type=int, help="initial state to average.",
    )
    parser.add_argument(
        "sf", type=int, help="final state to average.",
    )
    parser.add_argument(
        "--tex", action="store_true", dest="tex", help="use LaTeX to render text."
    )
    parser.add_argument(
        "--noshow", action="store_true", dest="noshow", help="do not display the plots."
    )
    parser.add_argument(
        "--mirror_y",
        action="store_true",
        dest="mirror_y",
        help="display the fundamental domain of mirror_y.",
    )
    parser.add_argument(
        "--mirror_z",
        action="store_true",
        dest="mirror_z",
        help="display the fundamental domain of mirror_z.",
    )
    parser.add_argument(
        "--savetxt",
        action="store_true",
        dest="savetxt",
        help="save results as text files",
    )
    parser.add_argument(
        "--loadtxt",
        action="store_true",
        dest="loadtxt",
        help="load results.",
    )
    args = vars(parser.parse_args())

    statesPath = Path(args["statesPath"])
    si = args["si"]
    sf = args["sf"]
    tex = args["tex"]
    noshow = args["noshow"]
    mirror_y = args["mirror_y"]
    mirror_z = args["mirror_z"]
    savetxt = args["savetxt"]
    loadtxt = args["loadtxt"]

    if loadtxt:
        corrs_avg = np.loadtxt("corrs_avg.gp")

    else:

        nstates = sf - si + 1

        print("Finding the average state.")
        state_avg = None
        for i in tqdm(range(si, sf + 1)):
            state_i = statesPath / f"state.{str(i).zfill(6)}"
            state, header = dns.readState(state_i)
            if state_avg is None:
                state_avg = state
            else:
                state_avg += state
        state_avg /= nstates

        print("Averaging over fluctuations.")
        corrs_avg = None

        for i in tqdm(range(si, sf + 1)):
            state_i = statesPath / f"state.{str(i).zfill(6)}"
            state, header = dns.readState(state_i)
            forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = header

            corrs, midy, nz_display = dnscorrelations(
                state - state_avg, header, mirror_y=mirror_y, mirror_z=mirror_z
            )

            if corrs_avg is None:
                corrs_avg = corrs
            else:
                corrs_avg += corrs

        corrs_avg = corrs_avg / nstates
        if savetxt:
            np.savetxt(statesPath / "corrs_avg.gp", corrs_avg)

    R_uu = corrs_avg / corrs_avg[0, 0]

    state_0 = statesPath / f"state.{str(si).zfill(6)}"
    state, header = dns.readState(state_0)
    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = header

    if not mirror_y:
        midy = ny // 2
    else:
        midy = ny // 4
    if not mirror_z:
        nz_display = nz
    else:
        nz_display = nz // 2 + 1

    title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $N=({nx},{ny},{nz})$"

    xs = np.array([i * (Lx / nx) for i in range(0, nx)])
    zs = np.array([k * (Lz / nz) for k in range(0, nz_display)])

    dns.setPlotDefaults(tex=tex)
    figuresDir = dns.createFiguresDir(statesPath.parent)

    fig, ax = plt.subplots()
    cbar = ax.imshow(
        R_uu.T,
        cmap=cmap,
        aspect="equal",
        origin="lower",
        vmin=-1,
        vmax=1,
        interpolation="spline16",
        extent=[xs[0], xs[-1], zs[0], zs[-1]],
    )
    ax.set_xlabel("$r_x$")
    ax.set_ylabel("$r_z$")
    ax.set_xlim(left=xs[0], right=xs[-1])
    ax.set_ylim(bottom=zs[0], top=zs[-1])
    ax.set_title(title)
    colorbar(ax, cbar, label=f"$R_{{uu}}(r_x,0,r_z)$")
    fig.savefig(figuresDir / f"R_uu.png", bbox_inches="tight")

    if not noshow:
        plt.show()


def dnscorrelations(
    stateIn, header, mirror_y=False, mirror_z=False,
):

    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = header
    ny_half = ny // 2

    vel_phys = dns.fftSpecToPhysAll(stateIn)

    velx = vel_phys[:, :, :, 0]

    if not mirror_y:
        midy = ny // 2
    else:
        midy = ny // 4
    if not mirror_z:
        nz_display = nz
    else:
        nz_display = nz // 2 + 1

    velx_midy = velx[:, midy, :nz_display]

    corrs = np.zeros((nx, nz_display))

    for ix in range(nx):
        for iz in range(nz_display):
            u_roll = np.roll(velx_midy, shift=(-ix, -iz), axis=(0, 1))
            corrs[ix, iz] = np.average(velx_midy * u_roll)

    return corrs, midy, nz_display


def colorbar(ax, im, label=None):

    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, label=label)

    return cbar


if __name__ == "__main__":
    main()
