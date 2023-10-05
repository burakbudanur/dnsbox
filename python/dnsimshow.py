#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import dnsbox as dns

cmap = "Greys"


def main():

    parser = argparse.ArgumentParser(
        description="Produce 2D visualizations of the velocity field.",
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

    state = Path(args["state"])
    tex = args["tex"]
    noshow = args["noshow"]
    undotilt = args["undotilt"]
    sublam = args["sublam"]
    mirror_y = args["mirror_y"]
    mirror_z = args["mirror_z"]

    results, headers, data = dnsimshow(
        state, undotilt=undotilt, sublam=sublam, mirror_y=mirror_y, mirror_z=mirror_z
    )
    # velx_midy, vorx_midy, velx_midz, vorx_midz = results
    velx_midy, velx_midz = results
    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers
    uw_untilted, ny_display, nz_display = data

    if abs(tilt_angle) > 0:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $\\theta={tilt_angle:.1f}$, $N=({nx},{ny},{nz})$"
    else:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $N=({nx},{ny},{nz})$"

    if sublam:
        title = "sublam, " + title
    if uw_untilted:
        title = "uw-untilted, " + title

    xs = np.array([i * (Lx / nx) for i in range(0, nx)])
    ys = np.array([j * (dns.Ly / ny) for j in range(0, ny_display)])
    zs = np.array([k * (Lz / nz) for k in range(0, nz_display)])

    xLabel = "$x$"
    zLabel = "$z$"
    yLabel = "$y$"

    if not mirror_y:
        yMid = "2"
    else:
        yMid = "1"
    if not mirror_z:
        zMid = "L_z/2"
    else:
        zMid = "L_z/4"

    dns.setPlotDefaults(tex=tex)
    figuresDir = dns.createFiguresDir(state.parent)

    # get the color scales

    # min_velx_midy, min_vorx_midy, min_velx_midz, min_vorx_midz = (
    #     np.amin(velx_midy),
    #     np.amin(vorx_midy),
    #     np.amin(velx_midz),
    #     np.amin(vorx_midz),
    # )

    # max_velx_midy, max_vorx_midy, max_velx_midz, max_vorx_midz = (
    #     np.amax(velx_midy),
    #     np.amax(vorx_midy),
    #     np.amax(velx_midz),
    #     np.amax(vorx_midz),
    # )

    min_velx_midy, min_velx_midz = (
        np.amin(velx_midy),
        np.amin(velx_midz),
    )

    max_velx_midy, max_velx_midz = (
        np.amax(velx_midy),
        np.amax(velx_midz),
    )

    scale_velx_midy = max(abs(min_velx_midy), abs(max_velx_midy))
    scale_velx_midz = max(abs(min_velx_midz), abs(max_velx_midz))

    # streamwise velocity, shearwise midplane
    figVelMid, axVelMid = plt.subplots()
    cVelMid = axVelMid.imshow(
        velx_midy.T,
        cmap=cmap,
        aspect="equal",
        origin="lower",
        vmin=-scale_velx_midy,
        vmax=scale_velx_midy,
        interpolation="spline16",
        extent=[xs[0], xs[-1], zs[0], zs[-1]],
    )
    axVelMid.set_xlabel(xLabel)
    axVelMid.set_ylabel(zLabel)
    axVelMid.set_xlim(left=xs[0], right=xs[-1])
    axVelMid.set_ylim(bottom=zs[0], top=zs[-1])
    axVelMid.set_title(title)
    colorbar(axVelMid, cVelMid, label=f"$u(y={yMid})$")
    figVelMid.savefig(
        figuresDir / f"{state.name}_velocity_midy.png", bbox_inches="tight"
    )

    # # streamwise vorticity, shearwise midplane
    # figVorMid, axVorMid = plt.subplots()
    # cVorMid = axVorMid.imshow(
    #     vorx_midy.T,
    #     cmap=cmap,
    #     aspect="equal",
    #     origin="lower",
    #     vmin=min_vorx_midy,
    #     vmax=max_vorx_midy,
    #     interpolation="spline16",
    #     extent=[xs[0], xs[-1], zs[0], zs[-1]],
    # )
    # axVorMid.set_xlabel(xLabel)
    # axVorMid.set_ylabel(zLabel)
    # axVorMid.set_xlim(left=xs[0], right=xs[-1])
    # axVorMid.set_ylim(bottom=zs[0], top=zs[-1])
    # axVorMid.set_title(title)
    # colorbar(axVorMid, cVorMid, label=f"$\\omega(y={yMid})$")
    # figVorMid.savefig(
    #     figuresDir / f"{state.name}_vorticity_midy.png", bbox_inches="tight"
    # )

    # streamwise velocity, spanwise midplane
    figVelMidZ, axVelMidZ = plt.subplots()
    cVelMidZ = axVelMidZ.imshow(
        velx_midz.T,
        cmap=cmap,
        aspect="equal",
        origin="lower",
        vmin=-scale_velx_midz,
        vmax=scale_velx_midz,
        interpolation="spline16",
        extent=[xs[0], xs[-1], ys[0], ys[-1]],
    )
    axVelMidZ.set_xlabel(xLabel)
    axVelMidZ.set_ylabel(yLabel)
    axVelMidZ.set_xlim(left=xs[0], right=xs[-1])
    axVelMidZ.set_ylim(bottom=ys[0], top=ys[-1])
    axVelMidZ.set_title(title)
    colorbar(axVelMidZ, cVelMidZ, label=f"$u(z={zMid})$")
    figVelMidZ.savefig(
        figuresDir / f"{state.name}_velocity_midz.png", bbox_inches="tight"
    )

    # # streamwise vorticity, spanwise midplane
    # figVorMidZ, axVorMidZ = plt.subplots()
    # cVorMidZ = axVorMidZ.imshow(
    #     vorx_midz.T,
    #     cmap=cmap,
    #     aspect="equal",
    #     origin="lower",
    #     vmin=min_vorx_midz,
    #     vmax=max_vorx_midz,
    #     interpolation="spline16",
    #     extent=[xs[0], xs[-1], ys[0], ys[-1]],
    # )
    # axVorMidZ.set_xlabel(xLabel)
    # axVorMidZ.set_ylabel(yLabel)
    # axVorMidZ.set_xlim(left=xs[0], right=xs[-1])
    # axVorMidZ.set_ylim(bottom=ys[0], top=ys[-1])
    # axVorMidZ.set_title(title)
    # colorbar(axVorMidZ, cVorMidZ, label=f"$\\omega(z={zMid})$")
    # figVorMidZ.savefig(
    #     figuresDir / f"{state.name}_vorticity_midz.png", bbox_inches="tight"
    # )

    if not noshow:
        plt.show()


def dnsimshow(
    state, undotilt=False, sublam=False, mirror_y=False, mirror_z=False,
):
    state = Path(state)

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
    # vor_spec = dns.vorticity(stateIn, Lx, Lz)
    # vor_phys = dns.fftSpecToPhysAll(vor_spec)

    velx = vel_phys[:, :, :, 0]
    # vorx = vor_phys[:, :, :, 0]

    if not mirror_y:
        ny_display = ny
        midy = ny // 2
    else:
        ny_display = ny // 2 + 1
        midy = ny // 4
    if not mirror_z:
        nz_display = nz
        midz = nz // 2
    else:
        nz_display = nz // 2 + 1
        midz = nz // 4

    # mid-y
    velx_midy = velx[:, midy, :nz_display]
    # vorx_midy = vorx[:, midy, :nz_display]

    # mid-z
    velx_midz = velx[:, :ny_display, midz]
    # vorx_midz = vorx[:, :ny_display, midz]

    # results = velx_midy, vorx_midy, velx_midz, vorx_midz
    results = velx_midy, velx_midz
    data = uw_untilted, ny_display, nz_display

    return results, headers, data


def colorbar(ax, im, label=None):

    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, label=label)

    return cbar


if __name__ == "__main__":
    main()
