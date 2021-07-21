#!/usr/bin/env python3
"""Produces plots of statistics of a run.
Figures are saved in a seperate directory within the run directory, and by
default they are shown.

"""

import argparse
from os import makedirs
from pathlib import Path
from sys import argv

import numpy as np
from matplotlib import pyplot as plt

import dnsCommon as dns

figuresDirName = "figures"


def main():
    parser = argparse.ArgumentParser(
        description="Produce plots of time series.", prog="dnsbox analysis"
    )
    parser.add_argument("runDir", type=str, help="path to the run folder.")
    parser.add_argument(
        "Ni",
        type=int if "--t" not in argv else float,
        help="initial line number in stat.gp, or if --t, initial time.",
    )
    parser.add_argument(
        "Nf",
        type=int if "--t" not in argv else float,
        help="final line number if stat.gp, or if --t, final time.",
    )
    parser.add_argument(
        "--noShow", action="store_true", dest="noShow", help="do not display the plots."
    )
    parser.add_argument(
        "--t",
        action="store_true",
        dest="t",
        help="use initial final time instead of line numbers.",
    )
    args = vars(parser.parse_args())

    noShow = args["noShow"]
    timeFilter = args["t"]

    runDir = args["runDir"]
    Ni = args["Ni"]
    Nf = args["Nf"]

    dnsPlot(
        runDir,
        Ni,
        Nf,
        timeFilter=timeFilter,
        noShow=noShow,
        setDefaults=True,
    )


def dnsPlot(
    runDir,
    Ni,
    Nf,
    timeFilter=False,
    noShow=False,
    setDefaults=False,
    plots=["ekin", "peps", "courant", "normrhs"],
):

    if setDefaults:
        # Set plotting defaults
        dns.setPlotDefaults()

    runDir = dns.checkDir(runDir)
    parameters = dns.readParameters(runDir / "parameters.in")

    # Save figures in the figures dir in run dir
    figuresDir = runDir / figuresDirName
    # Create it if it doesn't exist
    if not Path.is_dir(figuresDir):
        makedirs(figuresDir)

    dt = parameters["dt"]

    statsfile = "stat.gp"
    stats = np.loadtxt(runDir / statsfile, ndmin=2)

    if Path.is_dir(runDir / "continue"):
        stats2 = np.loadtxt(runDir / "continue/" / statsfile, ndmin=2)
        stats = np.append(stats, stats2, axis=0)

    if Path.is_dir(runDir / "continue" / "continue"):
        stats2 = np.loadtxt(runDir / "continue" / "continue/" / statsfile, ndmin=2)
        stats = np.append(stats, stats2, axis=0)

    if timeFilter:
        Ni = np.transpose(np.nonzero(stats[:, 1] > Ni - dt / 2))[0][0]
        Nf = np.transpose(np.nonzero(stats[:, 1] < Nf + dt / 2))[-1][0]

    stats = stats[Ni:Nf]
    time = np.copy(stats[:, 1])

    KineticEnergy = stats[:, 2]
    Production = stats[:, 3]
    Dissipation = stats[:, 4]
    normRHS = stats[:, 6]
    Courant = stats[:, 7]

    timeLabel = "$t$"

    # Turbulent kinetic energy
    if "ekin" in plots:
        figKin, axKin = plt.subplots()
        axKin.set_xlabel(timeLabel)
        axKin.set_ylabel("$k$")
        axKin.plot(time, KineticEnergy)
        figKin.savefig(figuresDir / "ekin.png")

    # Production-Dissipation plot
    if "peps" in plots:
        factor = 0.1
        figProdDis, axProdDis = plt.subplots()
        axProdDis.set_xlabel("$\\mathcal{P}$")
        axProdDis.set_ylabel("$\\epsilon$")
        axProdDis.plot(Production, Dissipation, color="gray", alpha=0.5)
        # unitline = np.linspace((1.0 - factor) * min(Dissipation),
        # 		       (1.0 + factor) * max(Dissipation), 100)

        minDP = min(0, min(min(Production), min(Dissipation)))
        maxDP = max(max(Production), max(Dissipation))
        # unitline = np.linspace(0.0, (1.0 + factor) * max(Dissipation), 100)
        unitline = np.linspace(minDP, (1.0 + factor) * maxDP, 100)
        axProdDis.plot(unitline, unitline, "g")
        axProdDis.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
        axProdDis.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
        # axProdDis.set_xlim((1.0 - factor) * min(Dissipation),
        #                    (1.0 + factor) * max(Dissipation))
        # axProdDis.set_xlim(0.0, (1.0 + factor) * max(Dissipation))
        axProdDis.set_xlim(min(0, min(Production)), (1.0 + factor) * max(Production))
        # axProdDis.set_ylim((1.0 - factor) * min(Dissipation),
        #                    (1.0 + factor) * max(Dissipation))
        # axProdDis.set_ylim(0.0, (1.0 + factor) * max(Dissipation))
        axProdDis.set_ylim(min(0, min(Dissipation)), (1.0 + factor) * max(Dissipation))
        figProdDis.savefig(figuresDir / "peps.png")

    # Courant number
    if "courant" in plots:
        figCourant, axCourant = plt.subplots()
        # axCourant.set_xlabel("$t\,(\\tau)$")
        axCourant.set_xlabel(timeLabel)
        axCourant.set_ylabel("$C$")
        axCourant.plot(time, Courant)
        figCourant.savefig(figuresDir / "courant.png")

    if "normrhs" in plots:
        # |RHS| / |u|
        figRHS, axRHS = plt.subplots()
        axRHS.set_xlabel(timeLabel)
        axRHS.set_ylabel("$|\\bf{F}(\\bf{u})| / |\\bf{u}|$")
        axRHS.plot(time, normRHS / np.sqrt(KineticEnergy))
        figRHS.savefig(figuresDir / "normrhs.png")

    if not noShow:
        plt.show()


if __name__ == "__main__":
    main()
