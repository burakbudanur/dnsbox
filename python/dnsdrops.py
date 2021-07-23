#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

import dns


def main():
    parser = argparse.ArgumentParser("Computes direction-dependent dropoffs.")
    parser.add_argument(
        "statePath",
        type=str,
        help="path to the state of interest.",
    )
    parser.add_argument(
        "--tex", action="store_true", dest="tex", help="use LaTeX to render text."
    )
    parser.add_argument(
        "--noshow", action="store_true", dest="noshow", help="do not display the plots."
    )

    args = vars(parser.parse_args())
    statePath = Path(args["statePath"])
    tex = args["tex"]
    noshow = args["noshow"]

    state, header = dns.readState(statePath)
    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = header
    (
        drop_x,
        drop_y,
        drop_z,
    ) = dnsdrop(state, header)

    wavenums_x = np.arange(drop_x.shape[0])
    wavenums_y = np.arange(drop_y.shape[0])
    wavenums_z = np.arange(drop_z.shape[0])

    dns.setPlotDefaults(tex=tex)
    figuresDir = dns.createFiguresDir(statePath.parent)

    if abs(tilt_angle) > 0:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $\\theta={tilt_angle:.1f}$, $N=({nx},{ny},{nz})$"
    else:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $N=({nx},{ny},{nz})$"

    fig, ax = plt.subplots()
    ax.plot(
        wavenums_x[1:],
        drop_x[1:],
        label="$\\max |{{\\bf u}}(n_x, n_y \\neq 0, n_z \\neq 0)|$",
    )
    ax.plot(
        wavenums_y[1:],
        drop_y[1:],
        label="$\\max |{{\\bf u}}(n_x \\neq 0, n_y, n_z \\neq 0)|$",
    )
    ax.plot(
        wavenums_z[1:],
        drop_z[1:],
        label="$\\max |{{\\bf u}}(n_x \\neq 0, n_y \\neq 0, n_z)|$",
    )
    ax.grid(True, which="both")
    ax.set_xlabel("$n$")
    ax.xaxis.get_major_locator().set_params(integer=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_title(title)
    fig.savefig(figuresDir / f"{statePath.name}_drops.png")

    if not noshow:
        plt.show()


def dnsdrop(state, header):

    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = header
    nxp, nyp, nzp = nx // 2 - 1, ny // 2 - 1, nz // 2 - 1

    drop_x = np.zeros((nxp + 1))
    drop_y = np.zeros((nyp + 1))
    drop_z = np.zeros((nzp + 1))

    norm2 = np.sum((np.conj(state) * state).real, axis=3)
    norm2[:, 0, :] = 0.5 * norm2[:, 0, :]
    norm = np.sqrt(norm2)

    for i in range(1, nxp + 1):
        drop_x[i] = max(np.amax(norm[i, 1:, 1:]), np.amax(norm[-i, 1:, 1:]))

    for j in range(nyp + 1):
        drop_y[j] = np.amax(norm[1:, j, 1:])

    for k in range(nzp + 1):
        drop_z[k] = max(np.amax(norm[1:, 1:, k]), np.amax(norm[1:, 1:, -k]))

    return (
        drop_x,
        drop_y,
        drop_z,
    )


if __name__ == "__main__":
    main()
