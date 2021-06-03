#!/usr/bin/env python3
"""
This is a module that keeps functions common to all utilities.

"""

from os import makedirs
from pathlib import Path
from struct import pack, unpack
from sys import exit

import numpy as np
from matplotlib import rcParams

figuresDirName = "figures"


def createFiguresDir(mainDir):
    mainDir = Path(mainDir).resolve()
    figuresDir = mainDir / figuresDirName
    # Create it if it doesn't exist
    if not Path.is_dir(figuresDir):
        makedirs(figuresDir)

    return figuresDir


def checkFile(file):
    # Check that file at the path exists, return a Path instance if so, exit
    # otherwise

    fileIn = Path(file).resolve()
    if not Path.is_file(fileIn):
        exit(f"File not found: {str(fileIn.resolve())}")

    return fileIn


def checkDir(directory):
    # Check that the folder at the path exists, return a Path instance if so,
    # exit otherwise

    dirIn = Path(directory).resolve()
    if not Path.is_dir(dirIn):
        exit(f"Directory not found: {str(dirIn.resolve())}")

    return dirIn


def checkGrid(Nz, Ny, NxAll):
    # Check to see if the grid sizes are all even and Ny=NxAll

    if (Nz % 2 != 0) or (Ny % 2 != 0) or (NxAll % 2 != 0):
        exit("Some grid lengths are not even.")

    if Ny != NxAll:
        exit("Ny != NxAll.")


def setPlotDefaults(usetex=False):
    # Used to apply some sane plotting defaults.

    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    rcParams["axes.labelsize"] = 36
    rcParams["figure.autolayout"] = True
    rcParams["figure.dpi"] = 100
    rcParams["figure.figsize"] = [8, 8]

    if usetex:
        # This should only be called if LaTeX is in path
        # And the default font choices look problematic, that needs work
        rcParams["text.usetex"] = True
        rcParams["text.latex.preamble"] = [
            r"\usepackage[T1]{fontenc}",
            r"\usepackage[utf8]{inputenc}",
            r"\usepackage{lmodern}",
            r"\usepackage{mathtools}",
            r"\usepackage{amssymb}",
            r"\usepackage{physics}",
            r"\usepackage{bm}",
        ]


def readParameters(parametersFile):
    """Reads parameters.in, returns parameters in a dictionary."""

    # Convert to Path
    parametersFile = checkFile(parametersFile)

    parameters = dict()

    with open(parametersFile) as readingParameters:
        for lineNo, line in enumerate(readingParameters):
            if lineNo == 3:
                parameters["nz"], parameters["ny"], parameters["nx_all"] = [
                    int(i) for i in line.split()[:3]
                ]
            elif lineNo == 4:
                parameters["Lz"] = float(line.split()[0])
            elif lineNo == 6:
                parameters["ITMIN"] = int(line.split()[0])
            elif lineNo == 7:
                parameters["ITMAX"] = int(line.split()[0])
            elif lineNo == 8:
                parameters["IPRINT1"] = int(line.split()[0])
            elif lineNo == 9:
                parameters["IPRINT2"] = int(line.split()[0])

            elif lineNo == 11:
                parameters["nu"] = float(line.split()[0])

            elif lineNo == 13:
                parameters["dt"] = float(line.split()[0])

            # Skipping min/max time steps

            elif lineNo == 15:
                parameters["IC"] = int(line.split()[0])

            # Skipping nonlinear term form

            elif lineNo == 17:
                parameters["Rz"] = int(line.split()[0])
            elif lineNo == 18:
                parameters["Sx"] = int(line.split()[0])
            # Skipping time stepping parameters and parsowrite

            elif lineNo == 24:
                parameters["lamSwitch"] = int(line.split()[0])
            elif lineNo == 25:
                parameters["lamThreshold"] = float(line.split()[0])

            elif lineNo == 27:
                parameters["seedNo"] = int(line.split()[0])

    return parameters


def readState(stateFilePath, symReduced=False):
    # Read dnsbox state file and return the array.

    # Convert to Path
    stateFilePath = checkFile(stateFilePath)
    # Open the file
    stateFile = open(stateFilePath, "rb")

    # Read grid size
    # They are kept as MPI_INTEGER4, and there are 3 many to read
    nz, ny, nxAll = unpack("=3i", stateFile.read(3 * 4))
    nx = nxAll
    time, dt, gamma, nu = unpack("=4d", stateFile.read(4 * 8))

    header = [nz, ny, nxAll, time, dt, gamma, nu]

    # Create the wrk array
    if not symReduced:
        state = np.zeros((nz + 2, ny, nx, 3), dtype=np.float64, order="F")
    else:
        state = np.zeros((2 * nz, ny // 2 + 1, nx, 3), dtype=np.float64, order="F")

    # Loop over the file

    for i in range(0, 3):
        if not symReduced:
            nCells = (nz + 2) * ny * nx
        else:
            nCells = (2 * nz) * (ny // 2 + 1) * nx
        toReshape = np.asarray(
            unpack("={}d".format(nCells), stateFile.read(nCells * 8))
        )

        # These were written in Fortran order, so read back in that order
        if not symReduced:
            state[:, :, :, i] = np.reshape(toReshape, (nz + 2, ny, nx), order="F")
        else:
            state[:, :, :, i] = np.reshape(
                toReshape, (2 * nz, ny // 2 + 1, nx), order="F"
            )

    stateFile.close()

    return state, header


def writeState(
    state,
    time=0.0,
    dt=0.005,
    gamma=1.0,
    nu=0.05,
    outFile="state.000000",
    symReduced=False,
):
    # Write array to a dnsbox state file.

    # Open the file
    outFile = Path(outFile).resolve()
    outDir = outFile.parent
    checkDir(outDir)
    stateFile = open(outFile, "wb")

    nz_, ny_, nxAll, _ = state.shape
    if not symReduced:
        nz = nz_ - 2
        ny = ny_
    else:
        nz = nz_ // 2
        ny = (ny_ - 1) * 2
    nx = nxAll

    # Write grid size
    # They are kept as MPI_INTEGER4, and there are 3 many to write
    gridDims = pack("=3i", nz, ny, nxAll)
    stateFile.write(gridDims)

    # Write TIME, DT, gamma and nu
    # They are MPI_REAL8
    header = pack("=4d", time, dt, gamma, nu)
    stateFile.write(header)

    # Loop over the file

    # Convert the array in
    # We will write each dimension in order
    for i in range(0, 3):
        if not symReduced:
            nCells = (nz + 2) * ny * nx
        else:
            nCells = (2 * nz) * (ny // 2 + 1) * nx
        # Careful about having this in Fortran order
        dataPack = pack("={}d".format(nCells), *state[:, :, :, i].flatten(order="F"))
        stateFile.write(dataPack)

    stateFile.close()


def inprodCoefficients(Nz, Ny, NxAll):
    # Return the inner product metric for dnsbox.

    checkGrid(Nz, Ny, NxAll)
    Nx = NxAll

    # Instead of summing the squares of all the coefficients, this goes over
    # the first half of the coefficients, paying care to the two non-hermiticity
    # related modes

    saveName = f"inprodCoefficients-{Nz}-{Ny}-{NxAll}-1-0.npz"
    savePath = Path(__file__).resolve().parent / saveName
    if Path.is_file(savePath):
        saveContent = np.load(savePath)
        coefficients = saveContent["coefficients"]
    else:
        fac = 1 / (Nz * Ny * NxAll) ** 2
        coefficients = fac * np.ones((Nz + 2, Ny, Nx, 3), dtype=np.float64, order="F")
        for i in range(0, Nz + 2):
            # Instead of summing the squares of all the coefficients, this goes over
            # the first half of the coefficients, paying care to the two non-hermiticity
            # related modes
            if i == 0 or i == 1 or i == Nz or i == Nz + 1:
                coefficients[i, :, :, :] = coefficients[i, :, :, :] / 2
        np.savez(savePath, coefficients=coefficients)

    return coefficients


def inprod(state1, state2, coefficients=None, normalized=False):
    # Inner product metric for dnsbox.
    # Set normalized to True if input arrays have unit metric.
    # (Such as symmetry reduced states.)
    # Input coefficients if you want to give the metric yourself.

    Nz_, Ny, NxAll, _ = state1.shape
    Nz = Nz_ - 2

    if state1.shape != state2.shape:
        exit("u1 and u2 have different shapes.")

    if coefficients is not None and normalized:
        exit("normalized mode doesn't need coefficients.")

    # Instead of summing the squares of all the coefficients, this goes over
    # the first half of the coefficients, paying care to the two non-hermiticity
    # related modes
    # This part should also ideally be vectorized with the appropriate
    # coefficient array
    if not normalized:
        if coefficients is None:
            coefficients = inprodCoefficients(Nz, Ny, NxAll)
    else:
        # We'll have absorbed coefficients for symmetry reduction
        coefficients = np.ones(state1.shape, dtype=np.float64, order="F")
    kinAll = np.sum(coefficients * state1 * state2)

    return kinAll


def wavenumbers(Nz, Ny, NxAll, Lz=1):
    # Return dnsbox wavenumber arrays.

    checkGrid(Nz, Ny, NxAll)

    # The wave number corresponding to u[i,j,k,:] is akx[j], akz[i], aky[k]
    # This sets up the akx, aky and akz arrays

    Nx = NxAll

    saveName = f"wavenumbers-{Lz}-{Nz}-{Ny}-{NxAll}-1-0.npz"
    savePath = Path(__file__).resolve().parent / saveName
    if Path.is_file(savePath):
        saveContent = np.load(savePath)
        akx = saveContent["akx"]
        aky = saveContent["aky"]
        akz = saveContent["akz"]
    else:
        akx = np.zeros((Ny), dtype=np.float64)
        # For parallelization reasons, there's a flip between the y and z indices
        # in the Fourier space
        aky = np.zeros((Nx), dtype=np.float64)
        akz = np.zeros((Nz + 2), dtype=np.float64)

        # akz is filled with [0, 0, 1, 1, ..., Nx/2, Nx/2] * (NxAll/Nz)
        for i in range(0, Nz + 1, 2):
            akz[i] = (i // 2) * (2 / Lz)
            akz[i + 1] = akz[i]

        # aky is filled with [0, 1, ..., Nx / 2 - 1, Nx / 2, -(Nx / 2) + 1, ..., -1]
        # for a single process run.
        # Parallelized runs distribute this between processes.

        for j in range(0, Nx):
            roll = 0 + j
            if roll <= NxAll // 2:
                aky[j] = roll
            else:
                aky[j] = roll - NxAll

        # akx is filled with [0, 1, ..., Ny / 2 - 1, Ny / 2, -(Ny / 2) + 1, ..., -1]
        for k in range(0, Ny):
            if k <= Ny // 2:
                akx[k] = k
            else:
                akx[k] = k - Ny

        np.savez(savePath, akx=akx, aky=aky, akz=akz)

    return akx, aky, akz
