#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import psutil
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from tqdm import tqdm

import dns
from dnsimshow import dnsimshow

cmap = "Spectral"
def_n_jobs = len(psutil.Process().cpu_affinity())
if def_n_jobs > 1:
    def_n_jobs = def_n_jobs - 1

def_joblib_verbosity = 0
def_joblib_backend = "loky"

# Example ffmpeg command to merge png files to an mp4:
# ffmpeg -framerate 24 -i %06d.png -c:v libx264 -r 60 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" out.mp4

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
    dnsbatchimshow(**args)


def dnsbatchimshow(
    statesDir, mirror_y=False, mirror_z=False, n_jobs=def_n_jobs, print_messages=True
):

    statesDir = Path(statesDir)
    states = sorted(list(statesDir.glob("state.*")))

    min_velx_midy, min_vorx_midy, min_velx_midz, min_vorx_midz = (
        np.inf,
        np.inf,
        np.inf,
        np.inf,
    )
    max_velx_midy, max_vorx_midy, max_velx_midz, max_vorx_midz = (
        -np.inf,
        -np.inf,
        -np.inf,
        -np.inf,
    )
    # get the scales first
    print("Finding the maxima in the snapshots.")
    with tqdm(total=len(states), disable=not print_messages) as pbar:
        for i, state in enumerate(states):
            state = Path(state)

            results, headers, data = dnsimshow(
                state, undotilt=False, sublam=True, mirror_y=mirror_y, mirror_z=mirror_z
            )
            velx_midy, vorx_midy, velx_midz, vorx_midz = results
            forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers
            uw_untilted, ny_display, nz_display = data

            min_velx_midy, min_vorx_midy, min_velx_midz, min_vorx_midz = (
                min(min_velx_midy, np.amin(velx_midy)),
                min(min_vorx_midy, np.amin(vorx_midy)),
                min(min_velx_midz, np.amin(velx_midz)),
                min(min_vorx_midz, np.amin(vorx_midz)),
            )

            max_velx_midy, max_vorx_midy, max_velx_midz, max_vorx_midz = (
                max(max_velx_midy, np.amax(velx_midy)),
                max(max_vorx_midy, np.amax(vorx_midy)),
                max(max_velx_midz, np.amax(velx_midz)),
                max(max_vorx_midz, np.amax(vorx_midz)),
            )
            pbar.update()

    print("Minima:")
    print(min_velx_midy, min_vorx_midy, min_velx_midz, min_vorx_midz)
    print("Maxima:")
    print(max_velx_midy, max_vorx_midy, max_velx_midz, max_vorx_midz)

    xs = np.array([i * (Lx / nx) for i in range(0, nx)])
    ys = np.array([j * (dns.Ly / ny) for j in range(0, ny_display)])
    zs = np.array([k * (Lz / nz) for k in range(0, nz_display)])

    dns.setPlotDefaults()
    figuresDir = dns.createFiguresDir(statesDir)

    # Now plot
    print("Plotting the snapshots.")

    def plot_state_i(i):
        state = Path(states[i])

        results, headers, data = dnsimshow(
            state, undotilt=False, sublam=True, mirror_y=mirror_y, mirror_z=mirror_z
        )
        velx_midy, vorx_midy, velx_midz, vorx_midz = results

        # streamwise velocity, shearwise midplane
        figVelMid, axVelMid = plt.subplots()
        axVelMid.imshow(
            velx_midy.T,
            cmap=cmap,
            aspect="equal",
            origin="lower",
            vmin=min_velx_midy,
            vmax=max_velx_midy,
            interpolation="spline16",
            extent=[xs[0], xs[-1], zs[0], zs[-1]],
        )
        axVelMid.set_xlim(left=xs[0], right=xs[-1])
        axVelMid.set_ylim(bottom=zs[0], top=zs[-1])
        axVelMid.axis("off")
        figVelMid.savefig(
            figuresDir / f"velocity_midy_{state.name}.png",
            bbox_inches="tight",
            pad_inches=0,
        )
        plt.close(figVelMid)

        # streamwise vorticity, shearwise midplane
        figVorMid, axVorMid = plt.subplots()
        axVorMid.imshow(
            vorx_midy.T,
            cmap=cmap,
            aspect="equal",
            origin="lower",
            vmin=min_vorx_midy,
            vmax=max_vorx_midy,
            interpolation="spline16",
            extent=[xs[0], xs[-1], zs[0], zs[-1]],
        )
        axVorMid.set_xlim(left=xs[0], right=xs[-1])
        axVorMid.set_ylim(bottom=zs[0], top=zs[-1])
        axVorMid.axis("off")
        figVorMid.savefig(
            figuresDir / f"vorticity_midy_{state.name}.png",
            bbox_inches="tight",
            pad_inches=0,
        )
        plt.close(figVorMid)

        # streamwise velocity, spanwise midplane
        figVelMidZ, axVelMidZ = plt.subplots()
        axVelMidZ.imshow(
            velx_midz.T,
            cmap=cmap,
            aspect="equal",
            origin="lower",
            vmin=min_velx_midz,
            vmax=max_velx_midz,
            interpolation="spline16",
            extent=[xs[0], xs[-1], ys[0], ys[-1]],
        )
        axVelMidZ.set_xlim(left=xs[0], right=xs[-1])
        axVelMidZ.set_ylim(bottom=ys[0], top=ys[-1])
        axVelMidZ.axis("off")
        figVelMidZ.savefig(
            figuresDir / f"velocity_midz_{state.name}.png",
            bbox_inches="tight",
            pad_inches=0,
        )
        plt.close(figVelMidZ)

        # streamwise vorticity, spanwise midplane
        figVorMidZ, axVorMidZ = plt.subplots()
        axVorMidZ.imshow(
            vorx_midz.T,
            cmap=cmap,
            aspect="equal",
            origin="lower",
            vmin=min_vorx_midz,
            vmax=max_vorx_midz,
            interpolation="spline16",
            extent=[xs[0], xs[-1], ys[0], ys[-1]],
        )
        axVorMidZ.set_xlim(left=xs[0], right=xs[-1])
        axVorMidZ.set_ylim(bottom=ys[0], top=ys[-1])
        axVorMidZ.axis("off")
        figVorMidZ.savefig(
            figuresDir / f"vorticity_midz_{state.name}.png",
            bbox_inches="tight",
            pad_inches=0,
        )
        plt.close(figVorMidZ)

    Parallel(n_jobs=n_jobs, backend=def_joblib_backend, verbose=def_joblib_verbosity)(
        delayed(plot_state_i)(i) for i in tqdm(range(len(states)))
    )


if __name__ == "__main__":
    main()
