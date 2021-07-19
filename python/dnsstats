#!/usr/bin/env python3
"""Produces plots of statistics of a run.
Figures are saved in a seperate directory within the run directory, and by
default they are shown.

"""

import argparse
from pathlib import Path
from sys import argv

import numpy as np
from matplotlib import pyplot as plt

import dns


def main():
    parser = argparse.ArgumentParser(description="Produce plots of time series.")
    parser.add_argument("runDir", type=str, help="path to the run folder.")
    parser.add_argument(
        "Ni",
        type=int if "--t" not in argv else float,
        help="initial line number in stat.gp, or if --tfilter, initial time.",
    )
    parser.add_argument(
        "Nf",
        type=int if "--t" not in argv else float,
        help="final line number if stat.gp, or if --tfilter, final time.",
    )
    parser.add_argument(
        "--noshow", action="store_true", dest="noshow", help="do not display the plots."
    )
    parser.add_argument(
        "--tfilter",
        action="store_true",
        dest="tfilter",
        help="use initial final time instead of line numbers.",
    )
    parser.add_argument(
        "--tex", action="store_true", dest="tex", help="use LaTeX to render text."
    )
    parser.add_argument(
        "--diet",
        action="store_true",
        dest="diet",
        help="plot only energy and input-dissipation.",
    )
    args = vars(parser.parse_args())

    dnsstats(**args)


def dnsstats(
    runDir,
    Ni,
    Nf,
    tfilter=False,
    noshow=False,
    tex=False,
    diet=False,
):

    dns.setPlotDefaults(tex=tex)

    runDir = Path(runDir).resolve()
    figuresDir = dns.createFiguresDir(runDir)

    statsfile = "stat.gp"
    stepsfile = "steps.gp"
    stats = np.loadtxt(runDir / statsfile, ndmin=2)

    # non-dimensionalize based on laminar velocity and its characteristic
    # length scale
    # assumed here qF = 1
    nml = dns.readParameters(runDir / "parameters.in")
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

    if abs(tilt_angle) > 0:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $\\theta={tilt_angle:.1f}$, $N=({nx},{ny},{nz})$"
    else:
        title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $N=({nx},{ny},{nz})$"

    # no non-dimensionalization whatsoever, use stat.gp as is
    if tfilter:
        Ni = np.transpose(np.nonzero(stats[:, 1] > Ni))[0][0]
        Nf = np.transpose(np.nonzero(stats[:, 1] < Nf))[-1][0]

    Elam = 1 / 4
    Edotlam = np.pi ** 2 / (8 * Re)
    stats = stats[Ni:Nf]
    time = stats[:, 1]
    KineticEnergy = stats[:, 2] / Elam
    Production = stats[:, 3] / Edotlam
    Dissipation = stats[:, 4] / Edotlam
    normRHS = stats[:, 5]
    # Courant = stats[:, 6]
    # ncorr = stats[:, 7]
    # err_corr = stats[:, 8]
    # dt = stats[:, 9]

    timeLabel = "$t$"
    ekinLabel = "$E / E_L$"
    prodLabel = "$I / I_L$"
    dissLabel = "$\\epsilon / \\epsilon_L$"
    normLabel = "$|{{\\bf F}}({{\\bf u}})|$"

    # Turbulent kinetic energy
    figKin, axKin = plt.subplots()
    axKin.set_xlabel(timeLabel)
    axKin.set_ylabel(ekinLabel)
    axKin.plot(time, KineticEnergy)
    axKin.set_title(title)
    figKin.savefig(figuresDir / "ekin.png")

    # Production-Dissipation plot
    factor = 0.1
    figProdDis, axProdDis = plt.subplots()
    axProdDis.set_xlabel(prodLabel)
    axProdDis.set_ylabel(dissLabel)
    axProdDis.plot(Production, Dissipation, color="gray", alpha=0.5)

    minDP = min(0, min(min(Production), min(Dissipation)))
    maxDP = max(max(Production), max(Dissipation))
    unitline = np.linspace(minDP, (1.0 + factor) * maxDP, 100)
    axProdDis.plot(unitline, unitline, "g")
    axProdDis.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    axProdDis.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    axProdDis.set_xlim(min(0, min(Production)), (1.0 + factor) * max(Production))
    axProdDis.set_ylim(min(0, min(Dissipation)), (1.0 + factor) * max(Dissipation))
    axProdDis.set_title(title)
    figProdDis.savefig(figuresDir / "peps.png")

    if not diet:

        # |RHS|
        figRHS, axRHS = plt.subplots()
        axRHS.set_xlabel(timeLabel)
        axRHS.set_ylabel(normLabel)
        axRHS.plot(time, normRHS)
        axRHS.set_title(title)
        figRHS.savefig(figuresDir / "normrhs.png")

        if Path.is_file(runDir / stepsfile):
            steps = np.loadtxt(runDir / stepsfile, ndmin=2)
            time = steps[:, 1]
            dt = steps[:, 2]
            courant = steps[:, 3]
            err_corr = steps[:, 4]
            ncorr = steps[:, 5]

            # time stepping stats
            figts, (axts1, axts2, axts4, axts3) = plt.subplots(
                nrows=4, ncols=1, sharex=True
            )
            axts3.plot(time, courant)
            axts3.set_xlabel(timeLabel)
            axts3.set_ylabel("Courant number")

            axts4.plot(time, dt)
            axts4.set_xlabel(timeLabel)
            axts4.set_ylabel("$dt$")

            axts2.plot(time[:-1], err_corr[1:])
            axts2.set_ylabel("Corrector rel. err.")

            axts1.plot(time[:-1], ncorr[1:])
            axts1.set_ylabel("Corrector iterations")
            axts1.yaxis.get_major_locator().set_params(integer=True)
            axts1.set_title(title)

            figts.savefig(figuresDir / "steps.png")

    if not noshow:
        plt.show()


if __name__ == "__main__":
    main()
