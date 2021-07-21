#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import dns

cmap = "Spectral"

# Example ffmpeg command to merge png files to an mp4:
# ffmpeg -framerate 5 -i %06d.png -c:v libx264 -r 60 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" out.mp4

# This script removes the laminar part


def main():

    parser = argparse.ArgumentParser(
        description="Batch produce 2D visualizations of streamwise velocity on the mid-y plane of many states.",
    )
    parser.add_argument(
        "statesDir", type=str, help="path to the directory containing states."
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
    dnsimshowmany(**args)


def dnsimshowmany(statesDir, mirror_y=False, mirror_z=False):

    statesDir = Path(statesDir)
    states = sorted(list(statesDir.glob("state.*")))

    velmax = 0
    # get the scales first
    print("Finding the maximum speed in the snapshots.")
    for i, state in enumerate(states):
        print(f"m {i} / {len(states)}")
        state = Path(state).resolve()

        stateIn, headers = dns.readState(state)

        forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers
        ny_half = ny // 2

        stateIn = stateIn - dns.laminar(forcing, nx, ny_half, nz, tilt_angle=tilt_angle)

        velx = dns.fftSpecToPhys(stateIn[:, :, :, 0])

        xs = np.array([i * (Lx / nx) for i in range(0, nx)])
        zs = np.array([k * (Lz / nz) for k in range(0, nz)])

        if not mirror_y:
            midy = ny // 2
        else:
            midy = ny // 4
        if not mirror_z:
            nz_display = nz
        else:
            nz_display = nz // 2 + 1

        # mid-y
        velx_midy = velx[:, midy, :]

        velmax = max(velmax, np.amax(velx_midy))

    dns.setPlotDefaults()
    figuresDir = dns.createFiguresDir(statesDir)

    # Now plot
    print("Plotting the snapshots.")
    for i, state in enumerate(states):
        print(f"p {i} / {len(states)}")
        state = Path(state).resolve()

        stateIn, headers = dns.readState(state)

        forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers
        ny_half = ny // 2
        stateIn = stateIn - dns.laminar(forcing, nx, ny_half, nz, tilt_angle=tilt_angle)
        velx = dns.fftSpecToPhys(stateIn[:, :, :, 0])

        xs = np.array([i * (Lx / nx) for i in range(0, nx)])
        zs = np.array([k * (Lz / nz) for k in range(0, nz)])

        if not mirror_y:
            midy = ny // 2
        else:
            midy = ny // 4
        if not mirror_z:
            nz_display = nz
        else:
            nz_display = nz // 2 + 1

        # mid-y
        velx_midy = velx[:, midy, :]

        # raw, without the colorbar
        figVelMid, axVelMid = plt.subplots()
        axVelMid.imshow(
            velx_midy[:, :nz_display].T,
            cmap=cmap,
            aspect="equal",
            origin="lower",
            vmin=-velmax,
            vmax=velmax,
            interpolation="spline16",
            extent=[xs[0], xs[-1], zs[0], zs[nz_display - 1]],
        )
        axVelMid.set_xlim(left=xs[0], right=xs[-1])
        axVelMid.set_ylim(bottom=zs[0], top=zs[nz_display - 1])
        axVelMid.axis("off")
        figVelMid.savefig(
            figuresDir / f"{state.name[-6:]}.png", bbox_inches="tight", pad_inches=0,
        )
        plt.close(figVelMid)


if __name__ == "__main__":
    main()
