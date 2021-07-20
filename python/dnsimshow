#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import dns

cmap = "Spectral"


def main():

    parser = argparse.ArgumentParser(
        description="Produce 2D visualizations of velocity and vorticity fields.",
    )
    parser.add_argument("state", type=str, help="path to the state.")
    parser.add_argument(
        "--tex", action="store_true", dest="tex", help="use LaTeX to render text."
    )
    parser.add_argument(
        "--noshow", action="store_true", dest="noshow", help="do not display the plots."
    )
    parser.add_argument(
        "--undotilt",
        action="store_true",
        dest="undotilt",
        help="rotate u and w to be parallel/orthogonal to the forcing.",
    )
    parser.add_argument(
        "--sublam",
        action="store_true",
        dest="sublam",
        help="remove laminar part from states.",
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
    args = vars(parser.parse_args())
    dnsimshow(**args)


def dnsimshow(
    state,
    tex=False,
    noshow=False,
    si=None,
    sf=None,
    undotilt=False,
    sublam=False,
    mirror_y=False,
    mirror_z=False,
):
    state = Path(state).resolve()
    figuresDir = dns.createFiguresDir(state.parent)

    stateIn, headers = dns.readState(state)

    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers
    ny_half = ny // 2

    if sublam:
        stateIn = stateIn - dns.laminar(forcing, nx, ny_half, nz, tilt_angle=tilt_angle)

    uw_untilted = False
    if abs(tilt_angle) > 0 and undotilt:
        uw_untilted = True
        stateIn = dns.tilt_state(stateIn, tilt_angle)

    vel_phys = dns.fftSpecToPhysAll(stateIn)
    vor_spec = dns.vorticity(stateIn, Lx, Lz)
    vor_phys = dns.fftSpecToPhysAll(vor_spec)

    velx = vel_phys[:, :, :, 0]
    vorx = vor_phys[:, :, :, 0]

    xs = np.array([i * (Lx / nx) for i in range(0, nx)])
    ys = np.array([j * (dns.Ly / ny) for j in range(0, ny)])
    zs = np.array([k * (Lz / nz) for k in range(0, nz)])

    xLabel = "$x$"
    zLabel = "$z$"
    yLabel = "$y$"

    if not mirror_y:
        yMid = "2"
        ny_display = ny
        midy = ny // 2
    else:
        yMid = "1"
        ny_display = ny // 2 + 1
        midy = ny // 4
    if not mirror_z:
        zMid = "L_z/2"
        nz_display = nz
        midz = nz // 2
    else:
        zMid = "L_z/4"
        nz_display = nz // 2 + 1
        midz = nz // 4

    # mid-y
    velx_midy = velx[:, midy, :]
    vorx_midy = vorx[:, midy, :]

    # mid-z
    velx_midz = velx[:, :, midz]
    vorx_midz = vorx[:, :, midz]

    if abs(tilt_angle) > 0:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $\\theta={tilt_angle:.1f}$, $N=({nx},{ny},{nz})$"
    else:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $N=({nx},{ny},{nz})$"

    if sublam:
        title = "sublam, " + title
    if uw_untilted:
        title = "uw-untilted, " + title

    dns.setPlotDefaults(tex=tex)

    # streamwise velocity, shearwise midplane
    figVelMid, axVelMid = plt.subplots()
    cVelMid = axVelMid.imshow(
        velx_midy[:, :nz_display].T,
        cmap=cmap,
        aspect="equal",
        origin="lower",
        interpolation="spline16",
        extent=[xs[0], xs[-1], zs[0], zs[nz_display - 1]],
    )
    axVelMid.set_xlabel(xLabel)
    axVelMid.set_ylabel(zLabel)
    axVelMid.set_xlim(left=xs[0], right=xs[-1])
    axVelMid.set_ylim(bottom=zs[0], top=zs[nz_display - 1])
    axVelMid.set_title(title)
    colorbar(axVelMid, cVelMid, label=f"$u(y={yMid})$")
    figVelMid.savefig(
        figuresDir / f"{state.name[-6:]}_velocity_midy.png", bbox_inches="tight"
    )

    # streamwise vorticity, shearwise midplane
    figVorMid, axVorMid = plt.subplots()
    cVorMid = axVorMid.imshow(
        vorx_midy[:, :nz_display].T,
        cmap=cmap,
        aspect="equal",
        origin="lower",
        interpolation="spline16",
        extent=[xs[0], xs[-1], zs[0], zs[nz_display - 1]],
    )
    axVorMid.set_xlabel(xLabel)
    axVorMid.set_ylabel(zLabel)
    axVorMid.set_xlim(left=xs[0], right=xs[-1])
    axVorMid.set_ylim(bottom=zs[0], top=zs[nz_display - 1])
    axVorMid.set_title(title)
    colorbar(axVorMid, cVorMid, label=f"$\\omega(y={yMid})$")
    figVorMid.savefig(
        figuresDir / f"{state.name[-6:]}_vorticity_midy.png", bbox_inches="tight"
    )

    # streamwise velocity, spanwise midplane
    figVelMidZ, axVelMidZ = plt.subplots()
    cVelMidZ = axVelMidZ.imshow(
        velx_midz[:, :ny_display].T,
        cmap=cmap,
        aspect="equal",
        origin="lower",
        interpolation="spline16",
        extent=[xs[0], xs[-1], ys[0], ys[ny_display - 1]],
    )
    axVelMidZ.set_xlabel(xLabel)
    axVelMidZ.set_ylabel(yLabel)
    axVelMidZ.set_xlim(left=xs[0], right=xs[-1])
    axVelMidZ.set_ylim(bottom=ys[0], top=ys[ny_display - 1])
    axVelMidZ.set_title(title)
    colorbar(axVelMidZ, cVelMidZ, label=f"$u(z={zMid})$")
    figVelMidZ.savefig(
        figuresDir / f"{state.name[-6:]}_velocity_midz.png", bbox_inches="tight"
    )

    # streamwise vorticity, spanwise midplane
    figVorMidZ, axVorMidZ = plt.subplots()
    cVorMidZ = axVorMidZ.imshow(
        vorx_midz[:, :ny_display].T,
        cmap=cmap,
        aspect="equal",
        origin="lower",
        interpolation="spline16",
        extent=[xs[0], xs[-1], ys[0], ys[ny_display - 1]],
    )
    axVorMidZ.set_xlabel(xLabel)
    axVorMidZ.set_ylabel(yLabel)
    axVorMidZ.set_xlim(left=xs[0], right=xs[-1])
    axVorMidZ.set_ylim(bottom=ys[0], top=ys[ny_display - 1])
    axVorMidZ.set_title(title)
    colorbar(axVorMidZ, cVorMidZ, label=f"$\\omega(z={zMid})$")
    figVorMidZ.savefig(
        figuresDir / f"{state.name[-6:]}_vorticity_midz.png", bbox_inches="tight"
    )

    if not noshow:
        plt.show()


def colorbar(ax, im, label=None):

    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, label=label)

    return cbar


if __name__ == "__main__":
    main()
