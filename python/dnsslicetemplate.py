#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import dnsbox as dns


def main():
    parser = argparse.ArgumentParser(description="Produce the optimal slice template.")
    parser.add_argument("rundir", type=str, help="path to the run folder.")
    parser.add_argument("savedir", type=str, help="folder to save results.")
    parser.add_argument(
        "-cutpercent",
        type=int,
        default=20,
        help="how much percent of data data to cut from the start and the end of the trajectory ",
    )
    parser.add_argument(
        "--noshow", action="store_true", dest="noshow", help="do not display the plots."
    )
    parser.add_argument(
        "--tex", action="store_true", dest="tex", help="use LaTeX to render text."
    )
    args = vars(parser.parse_args())

    dnsslicetemplate(**args)


def dnsslicetemplate(rundir, savedir, cutpercent=20, tex=False, noshow=False):

    dns.setPlotDefaults(tex=tex)

    rundir = Path(rundir)
    savedir = Path(savedir)

    # non-dimensionalize based on laminar velocity and its characteristic
    # length scale
    # assumed here qF = 1
    nml = dns.readParameters(rundir / "parameters.in")
    Re = nml["physics"]["Re"]
    try:
        tilt_angle = nml["physics"]["tilt_angle"]
    except:
        tilt_angle = 0

    Lx = nml["grid"]["Lx"]
    Lz = nml["grid"]["Lz"]
    nx = nml["grid"]["nx"]
    ny = nml["grid"]["ny"]
    nz = nml["grid"]["nz"]
    Ry = nml["symmetries"]["Ry"]

    ny_half = ny // 2

    if Ry and ny_half % 2 != 0:
        exit("Ry but ny_half is not even.")

    if not Ry:
        exit("non-Ry not implemented yet.")

    if abs(tilt_angle) > 0:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $\\theta={tilt_angle:.1f}$, $N=({nx},{ny},{nz})$"
    else:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $N=({nx},{ny},{nz})$"

    data_x = np.loadtxt(rundir / "slice_projections_x.gp")
    data_z = np.loadtxt(rundir / "slice_projections_z.gp")

    if cutpercent > 0:
        data_x = data_x[
            cutpercent * len(data_x) // 100 : -cutpercent * len(data_x) // 100
        ]
        data_z = data_z[
            cutpercent * len(data_z) // 100 : -cutpercent * len(data_z) // 100
        ]

    projections_x_rview = data_x[:, 2:]
    projections_z_rview = data_z[:, 2:]

    projections_x = (
        projections_x_rview[:, ::2] + 1j * projections_x_rview[:, 1:][:, ::2]
    )
    projections_z = (
        projections_z_rview[:, ::2] + 1j * projections_z_rview[:, 1:][:, ::2]
    )

    times_x = data_x[:, 1]
    times_z = data_z[:, 1]

    opt_maximize_projection_amplitudes(
        savedir,
        np.conj(projections_x),
        "x",
        times_x,
        nx,
        ny_half,
        nz,
        Lx / 2,
        title=title,
    )
    opt_maximize_projection_amplitudes(
        savedir,
        np.conj(projections_z),
        "z",
        times_z,
        nx,
        ny_half,
        nz,
        Lz / 2,
        title=title,
    )

    if not noshow:
        plt.show()


def opt_maximize_projection_amplitudes(
    savedir,
    projections,
    str_projections,
    times,
    nx,
    ny_half,
    nz,
    L,
    title=None,
):

    u, _, _ = np.linalg.svd(
        np.einsum("ki,kj", np.conj(projections), projections).real, hermitian=True
    )
    coeffs = u[:, 0] / np.sqrt(np.sum(u[:, 0] * u[:, 0]))

    # save coefficients
    np.savetxt(savedir / f"u_{str_projections}p.gp", coeffs)

    # save template
    template = np.zeros((nx, ny_half, nz, 3), dtype=np.complex128)
    if str_projections == "x":
        template[1, 1:, 0, 0] = coeffs[: ny_half - 1] / 4
        template[-1, 1:, 0, 0] = coeffs[: ny_half - 1] / 4

        template[1, 1:, 0, 1] = -1j * coeffs[ny_half - 1 : 2 * ny_half - 2] / 4
        template[-1, 1:, 0, 1] = -1j * coeffs[ny_half - 1 : 2 * ny_half - 2] / 4

        template[1, 0, 0, 2] = coeffs[2 * ny_half - 2] / 2 / np.sqrt(2)
        template[-1, 0, 0, 2] = coeffs[2 * ny_half - 2] / 2 / np.sqrt(2)

        template[1, 1:, 0, 2] = coeffs[2 * ny_half - 1 :] / 4
        template[-1, 1:, 0, 2] = coeffs[2 * ny_half - 1 :] / 4

    elif str_projections == "z":
        template[0, 0, 1, 0] = coeffs[0] / 2 / np.sqrt(2)
        template[0, 0, -1, 0] = coeffs[0] / 2 / np.sqrt(2)

        template[0, 1:, 1, 0] = coeffs[1:ny_half] / 4
        template[0, 1:, -1, 0] = coeffs[1:ny_half] / 4

        template[0, 1:, 1, 1] = -1j * coeffs[ny_half : 2 * ny_half - 1] / 4
        template[0, 1:, -1, 1] = -1j * coeffs[ny_half : 2 * ny_half - 1] / 4

        template[0, 1:, 1, 2] = coeffs[2 * ny_half - 1 :] / 4
        template[0, 1:, -1, 2] = coeffs[2 * ny_half - 1 :] / 4

    dns.writeState_nocompact(template, outFile=savedir / f"u_{str_projections}p.000000")

    projections_opt = np.einsum("ij,j", projections, coeffs, dtype=np.complex128)
    dphases = find_dphases(projections_opt, times)

    figuresdir = dns.createFiguresDir(savedir)

    # slice templates have norm 1 / sqrt(8)
    fig_proj, ax_proj = plt.subplots()
    ax_proj.plot(times, np.abs(projections_opt) * np.sqrt(8))
    ax_proj.set_xlabel("$t$")
    ax_proj.set_ylabel(f"$|p_{str_projections}|$")
    ax_proj.set_yscale("log")
    ax_proj.set_xlim(left=times[0], right=times[-1])
    fig_proj.savefig(figuresdir / f"p{str_projections}.png", bbox_inches="tight")

    fig_phasevel, ax_phasevel = plt.subplots()
    ax_phasevel.plot(times[:-1], dphases * L / (2 * np.pi))
    ax_phasevel.set_xlabel("$t$")
    ax_phasevel.set_ylabel(
        f"$\\dot{{\\phi_{str_projections}}} L_{str_projections} / (2\\pi)$"
    )
    ax_phasevel.set_xlim(left=times[0], right=times[-2])
    fig_phasevel.savefig(
        figuresdir / f"dphase_{str_projections}.png", bbox_inches="tight"
    )

    return fig_proj, ax_proj, fig_phasevel, ax_phasevel


def find_shift(projection):
    shift = (np.arctan2(projection.imag, projection.real) / (2 * np.pi)) % 1
    return shift


def find_dphases(projections, times):
    shifts = find_shift(projections)
    phases = np.unwrap(shifts * 2 * np.pi)
    dts = times[1:] - times[:-1]
    dphases = (phases[1:] - phases[:-1]) / dts
    return dphases


if __name__ == "__main__":
    main()
