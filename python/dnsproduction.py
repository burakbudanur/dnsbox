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
        "--mirror_y",
        action="store_true",
        dest="mirror_y",
        help="display the fundamental domain of mirror_y.",
    )
    args = vars(parser.parse_args())

    state = Path(args["state"])
    tex = args["tex"]
    noshow = args["noshow"]
    mirror_y = args["mirror_y"]

    results, headers, data = dnsproduction(state, mirror_y=mirror_y)
    input_field, input_y = results
    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers
    ny_display = data

    if abs(tilt_angle) > 0:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $\\theta={tilt_angle:.1f}$, $N=({nx},{ny},{nz})$"
    else:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $N=({nx},{ny},{nz})$"

    xs = np.array([i * (Lx / nx) for i in range(0, nx)])
    ys = np.array([j * (dns.Ly / ny) for j in range(0, ny_display)])

    xLabel = "$x$"
    yLabel = "$y$"

    if not mirror_y:
        yMid = "2"
    else:
        yMid = "1"

    dns.setPlotDefaults(tex=tex)
    figuresDir = dns.createFiguresDir(state.parent)

    # get the color scales
    min_input_field = np.inf
    max_input_field = -np.inf
    min_input_field = min(min_input_field, np.amin(input_field))

    max_input_field = max(max_input_field, np.amax(input_field))

    # streamwise velocity, shearwise midplane
    figInput, axInput = plt.subplots()
    cInput = axInput.imshow(
        input_field.T,
        cmap=cmap,
        aspect="equal",
        origin="lower",
        vmin=min_input_field,
        vmax=max_input_field,
        interpolation="spline16",
        extent=[xs[0], xs[-1], ys[0], ys[-1]],
    )
    axInput.set_xlabel(xLabel)
    axInput.set_ylabel(yLabel)
    axInput.set_xlim(left=xs[0], right=xs[-1])
    axInput.set_ylim(bottom=ys[0], top=ys[-1])
    axInput.set_title(title)
    colorbar(axInput, cInput, label=f"$I(x,y)$")
    figInput.savefig(figuresDir / f"{state.name}_input.png", bbox_inches="tight")

    fig_input_y, ax_input_y = plt.subplots()
    ax_input_y.plot(ys, input_y)
    ax_input_y.set_xlabel(yLabel)
    ax_input_y.set_ylabel("I(y)")
    fig_input_y.savefig(figuresDir / f"{state.name}_input_y.png", bbox_inches="tight")

    if not noshow:
        plt.show()


def dnsproduction(
    state,
    mirror_y=False,
):
    state = Path(state)

    stateIn, headers = dns.readState(state)

    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers
    ny_half = ny // 2

    vel_phys = dns.fftSpecToPhysAll(stateIn)
    forcing_phys = (np.pi ** 2 / (4 * Re)) * dns.fftSpecToPhysAll(
        dns.laminar(forcing, nx, ny_half, nz, tilt_angle=tilt_angle)
    )

    if not mirror_y:
        ny_display = ny
        midy = ny // 2
    else:
        ny_display = ny // 2 + 1
        midy = ny // 4

    # dot product + integrate over z
    input_field = np.sum(vel_phys * forcing_phys, axis=(2, 3)) / (2 * nx * ny * nz)

    # integrate over x
    input_y = np.sum(input_field, axis=0)

    input_field = input_field[:, :ny_display]
    input_y = input_y[:ny_display]

    results = input_field, input_y
    data = ny_display

    return results, headers, data


def colorbar(ax, im, label=None):

    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, label=label)

    return cbar


if __name__ == "__main__":
    main()
