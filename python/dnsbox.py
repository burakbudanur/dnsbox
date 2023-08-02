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
from numpy.lib import NumpyVersion

figuresDirName = "figures"

Ly = 4.0
qF = 1
nStretch = 3


def isEven(N):
    # Return True if N is even, False otherwise

    if N % 2 == 0:
        return True
    else:
        return False


def createFiguresDir(mainDir):
    mainDir = Path(mainDir)
    figuresDir = mainDir / figuresDirName
    # Create it if it doesn't exist
    if not Path.is_dir(figuresDir):
        makedirs(figuresDir)

    return figuresDir


def setPlotDefaults(tex=False):
    # Used to apply some sane plotting defaults.

    rcParams.update(
        {
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "axes.labelsize": 13,
            "axes.titlesize": 14,
            "legend.fontsize": 10,
            "figure.autolayout": True,
            "figure.dpi": 100,
            "figure.figsize": [8, 8],
        }
    )

    if tex:

        rcParams.update(
            {
                "text.usetex": True,
                "text.latex.preamble": r"\usepackage[p,osf]{scholax}"
                r"\usepackage[T1]{fontenc}"
                r"\usepackage[utf8]{inputenc}"
                r"\usepackage[varqu,varl]{inconsolata}"
                r"\usepackage{mathtools}"
                r"\usepackage[scaled=1.075,ncf,vvarbb]{newtxmath}"
                r"\usepackage{bm}",
            }
        )


def readParameters(parametersFile):
    """Reads parameters.in, returns parameters in a dictionary."""
    import f90nml

    # Convert to Path
    parametersFile = Path(parametersFile)

    params = f90nml.read(parametersFile)

    return params


def writeParameters(params, parametersFile):
    parametersFile = Path(parametersFile)
    with open(parametersFile, "w") as f:
        params.write(f)


def readState_xcompact(stateFilePath):

    stateFilePath = Path(stateFilePath)
    stateFile = open(stateFilePath, "rb")

    (forcing,) = unpack("=1i", stateFile.read(1 * 4))
    nx, ny, nz = unpack("=3i", stateFile.read(3 * 4))
    Lx, Lz = unpack("=2d", stateFile.read(2 * 8))
    (Re,) = unpack("=1d", stateFile.read(1 * 8))
    (tilt_angle,) = unpack("=1d", stateFile.read(1 * 8))
    (dt,) = unpack("=1d", stateFile.read(1 * 8))
    (itime,) = unpack("=1i", stateFile.read(1 * 4))
    (time,) = unpack("=1d", stateFile.read(1 * 8))

    ny_half = ny // 2
    nxp, nzp = nx // 2 - 1, nz // 2 - 1
    header = (forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time)

    # Read into y, z, x
    state = np.zeros((2 * ny_half, nz, nx, 3), dtype=np.float64)

    # read only the zero mode for u
    nCells = (2 * ny_half) * (nz - 1)
    toReshape = np.asarray(unpack("={}d".format(nCells), stateFile.read(nCells * 8)))

    buff = np.reshape(toReshape, (2 * ny_half, nz - 1, 1), order="F")
    state[:, : nzp + 1, 0, 0] = buff[:, : nzp + 1, 0]
    state[:, nzp + 2 :, 0, 0] = buff[:, nzp + 1 :, 0]

    # read v and w in full
    for n in range(1, 3):
        nCells = (nx - 1) * (2 * ny_half) * (nz - 1)
        toReshape = np.asarray(
            unpack("={}d".format(nCells), stateFile.read(nCells * 8))
        )

        # These were written in Fortran order, so read back in that order
        buff = np.reshape(toReshape, (2 * ny_half, nz - 1, nx - 1), order="F")
        state[:, : nzp + 1, : nxp + 1, n] = buff[:, : nzp + 1, : nxp + 1]
        state[:, : nzp + 1, nxp + 2 :, n] = buff[:, : nzp + 1, nxp + 1 :]
        state[:, nzp + 2 :, : nxp + 1 :, n] = buff[:, nzp + 1 :, : nxp + 1]
        state[:, nzp + 2 :, nxp + 2 :, n] = buff[:, nzp + 1 :, nxp + 1 :]

    stateFile.close()

    # compute the rest from
    # kx u + ky v + kz w = 0
    kx, ky, kz = wavenumbers(Lx, Lz, nx, ny_half, nz)
    for j in range(ny_half):
        for k in range(nz):
            for i in range(1, nx):
                state[2 * j, k, i, 0] = (
                    -(kz[k] * state[2 * j, k, i, 2] + ky[j] * state[2 * j, k, i, 1])
                    / kx[i]
                )
                state[2 * j + 1, k, i, 0] = (
                    -(
                        kz[k] * state[2 * j + 1, k, i, 2]
                        + ky[j] * state[2 * j + 1, k, i, 1]
                    )
                    / kx[i]
                )

    # Pack into a complex array
    state_ = np.zeros((ny_half, nz, nx, 3), dtype=np.complex128)
    for j in range(ny_half):
        state_[j, :, :, :] = state[2 * j, :, :, :] + 1j * state[2 * j + 1, :, :, :]

    # Rotate to x, y, z
    state__ = np.moveaxis(state_, [0, 1, 2, 3], [1, 2, 0, 3])

    return state__, header


def readState_ycompact(stateFilePath):

    stateFilePath = Path(stateFilePath)
    stateFile = open(stateFilePath, "rb")

    (forcing,) = unpack("=1i", stateFile.read(1 * 4))
    nx, ny, nz = unpack("=3i", stateFile.read(3 * 4))
    Lx, Lz = unpack("=2d", stateFile.read(2 * 8))
    (Re,) = unpack("=1d", stateFile.read(1 * 8))
    (tilt_angle,) = unpack("=1d", stateFile.read(1 * 8))
    (dt,) = unpack("=1d", stateFile.read(1 * 8))
    (itime,) = unpack("=1i", stateFile.read(1 * 4))
    (time,) = unpack("=1d", stateFile.read(1 * 8))

    ny_half = ny // 2
    nxp, nzp = nx // 2 - 1, nz // 2 - 1
    header = (forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time)

    # Read into y, z, x
    state = np.zeros((2 * ny_half, nz, nx, 3), dtype=np.float64)

    # read u in full
    nCells = (nx - 1) * (2 * ny_half) * (nz - 1)
    toReshape = np.asarray(unpack("={}d".format(nCells), stateFile.read(nCells * 8)))

    buff = np.reshape(toReshape, (2 * ny_half, nz - 1, nx - 1), order="F")
    state[:, : nzp + 1, : nxp + 1, 0] = buff[:, : nzp + 1, : nxp + 1]
    state[:, : nzp + 1, nxp + 2 :, 0] = buff[:, : nzp + 1, nxp + 1 :]
    state[:, nzp + 2 :, : nxp + 1 :, 0] = buff[:, nzp + 1 :, : nxp + 1]
    state[:, nzp + 2 :, nxp + 2 :, 0] = buff[:, nzp + 1 :, nxp + 1 :]

    # read only the zero mode for v
    nCells = 2 * (nx - 1) * (nz - 1)
    toReshape = np.asarray(unpack("={}d".format(nCells), stateFile.read(nCells * 8)))

    buff = np.reshape(toReshape, (2, nz - 1, nx - 1), order="F")
    state[:2, : nzp + 1, : nxp + 1, 1] = buff[:2, : nzp + 1, : nxp + 1]
    state[:2, : nzp + 1, nxp + 2 :, 1] = buff[:2, : nzp + 1, nxp + 1 :]
    state[:2, nzp + 2 :, : nxp + 1, 1] = buff[:2, nzp + 1 :, : nxp + 1]
    state[:2, nzp + 2 :, nxp + 2 :, 1] = buff[:2, nzp + 1 :, nxp + 1 :]

    # read w in full
    nCells = (nx - 1) * (2 * ny_half) * (nz - 1)
    toReshape = np.asarray(unpack("={}d".format(nCells), stateFile.read(nCells * 8)))

    buff = np.reshape(toReshape, (2 * ny_half, nz - 1, nx - 1), order="F")
    state[:, : nzp + 1, : nxp + 1, 2] = buff[:, : nzp + 1, : nxp + 1]
    state[:, : nzp + 1, nxp + 2 :, 2] = buff[:, : nzp + 1, nxp + 1 :]
    state[:, nzp + 2 :, : nxp + 1 :, 2] = buff[:, nzp + 1 :, : nxp + 1]
    state[:, nzp + 2 :, nxp + 2 :, 2] = buff[:, nzp + 1 :, nxp + 1 :]

    stateFile.close()

    # compute the rest from
    # kx u + ky v + kz w = 0
    kx, ky, kz = wavenumbers(Lx, Lz, nx, ny_half, nz)
    for j in range(1, ny_half):
        for k in range(nz):
            for i in range(nx):
                state[2 * j, k, i, 1] = (
                    -(kx[i] * state[2 * j, k, i, 0] + kz[k] * state[2 * j, k, i, 2])
                    / ky[j]
                )
                state[2 * j + 1, k, i, 1] = (
                    -(
                        kx[i] * state[2 * j + 1, k, i, 0]
                        + kz[k] * state[2 * j + 1, k, i, 2]
                    )
                    / ky[j]
                )

    # Pack into a complex array
    state_ = np.zeros((ny_half, nz, nx, 3), dtype=np.complex128)
    for j in range(ny_half):
        state_[j, :, :, :] = state[2 * j, :, :, :] + 1j * state[2 * j + 1, :, :, :]

    # Rotate to x, y, z
    state__ = np.moveaxis(state_, [0, 1, 2, 3], [1, 2, 0, 3])

    return state__, header


def readState_zcompact(stateFilePath):

    stateFilePath = Path(stateFilePath)
    stateFile = open(stateFilePath, "rb")

    (forcing,) = unpack("=1i", stateFile.read(1 * 4))
    nx, ny, nz = unpack("=3i", stateFile.read(3 * 4))
    Lx, Lz = unpack("=2d", stateFile.read(2 * 8))
    (Re,) = unpack("=1d", stateFile.read(1 * 8))
    (tilt_angle,) = unpack("=1d", stateFile.read(1 * 8))
    (dt,) = unpack("=1d", stateFile.read(1 * 8))
    (itime,) = unpack("=1i", stateFile.read(1 * 4))
    (time,) = unpack("=1d", stateFile.read(1 * 8))

    ny_half = ny // 2
    nxp, nzp = nx // 2 - 1, nz // 2 - 1
    header = (forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time)

    # Read into y, z, x
    state = np.zeros((2 * ny_half, nz, nx, 3), dtype=np.float64)

    # read u and v in full
    for n in range(0, 2):
        nCells = (nx - 1) * (2 * ny_half) * (nz - 1)
        toReshape = np.asarray(
            unpack("={}d".format(nCells), stateFile.read(nCells * 8))
        )

        # These were written in Fortran order, so read back in that order
        buff = np.reshape(toReshape, (2 * ny_half, nz - 1, nx - 1), order="F")
        state[:, : nzp + 1, : nxp + 1, n] = buff[:, : nzp + 1, : nxp + 1]
        state[:, : nzp + 1, nxp + 2 :, n] = buff[:, : nzp + 1, nxp + 1 :]
        state[:, nzp + 2 :, : nxp + 1 :, n] = buff[:, nzp + 1 :, : nxp + 1]
        state[:, nzp + 2 :, nxp + 2 :, n] = buff[:, nzp + 1 :, nxp + 1 :]

    # read only the zero mode for w
    nCells = (nx - 1) * (2 * ny_half)
    toReshape = np.asarray(unpack("={}d".format(nCells), stateFile.read(nCells * 8)))

    buff = np.reshape(toReshape, (2 * ny_half, 1, nx - 1), order="F")
    state[:, 0, : nxp + 1, 2] = buff[:, 0, : nxp + 1]
    state[:, 0, nxp + 2 :, 2] = buff[:, 0, nxp + 1 :]

    stateFile.close()

    # compute the rest from
    # kx u + ky v + kz w = 0
    kx, ky, kz = wavenumbers(Lx, Lz, nx, ny_half, nz)
    for j in range(ny_half):
        for k in range(1, nz):
            for i in range(nx):
                state[2 * j, k, i, 2] = (
                    -(kx[i] * state[2 * j, k, i, 0] + ky[j] * state[2 * j, k, i, 1])
                    / kz[k]
                )
                state[2 * j + 1, k, i, 2] = (
                    -(
                        kx[i] * state[2 * j + 1, k, i, 0]
                        + ky[j] * state[2 * j + 1, k, i, 1]
                    )
                    / kz[k]
                )

    # Pack into a complex array
    state_ = np.zeros((ny_half, nz, nx, 3), dtype=np.complex128)
    for j in range(ny_half):
        state_[j, :, :, :] = state[2 * j, :, :, :] + 1j * state[2 * j + 1, :, :, :]

    # Rotate to x, y, z
    state__ = np.moveaxis(state_, [0, 1, 2, 3], [1, 2, 0, 3])

    return state__, header


def readState_nocompact(stateFilePath):

    stateFilePath = Path(stateFilePath)
    stateFile = open(stateFilePath, "rb")

    (forcing,) = unpack("=1i", stateFile.read(1 * 4))
    nx, ny, nz = unpack("=3i", stateFile.read(3 * 4))
    Lx, Lz = unpack("=2d", stateFile.read(2 * 8))
    (Re,) = unpack("=1d", stateFile.read(1 * 8))
    (tilt_angle,) = unpack("=1d", stateFile.read(1 * 8))
    (dt,) = unpack("=1d", stateFile.read(1 * 8))
    (itime,) = unpack("=1i", stateFile.read(1 * 4))
    (time,) = unpack("=1d", stateFile.read(1 * 8))

    ny_half = ny // 2
    nxp, nzp = nx // 2 - 1, nz // 2 - 1
    header = (forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time)

    # Read into y, z, x
    state = np.zeros((2 * ny_half, nz, nx, 3), dtype=np.float64)

    for n in range(3):
        nCells = (nx - 1) * (2 * ny_half) * (nz - 1)
        toReshape = np.asarray(
            unpack("={}d".format(nCells), stateFile.read(nCells * 8))
        )

        # These were written in Fortran order, so read back in that order
        buff = np.reshape(toReshape, (2 * ny_half, nz - 1, nx - 1), order="F")
        state[:, : nzp + 1, : nxp + 1, n] = buff[:, : nzp + 1, : nxp + 1]
        state[:, : nzp + 1, nxp + 2 :, n] = buff[:, : nzp + 1, nxp + 1 :]
        state[:, nzp + 2 :, : nxp + 1 :, n] = buff[:, nzp + 1 :, : nxp + 1]
        state[:, nzp + 2 :, nxp + 2 :, n] = buff[:, nzp + 1 :, nxp + 1 :]

    stateFile.close()

    # Pack into a complex array
    state_ = np.zeros((ny_half, nz, nx, 3), dtype=np.complex128)
    for j in range(ny_half):
        state_[j, :, :, :] = state[2 * j, :, :, :] + 1j * state[2 * j + 1, :, :, :]

    # Rotate to x, y, z
    state__ = np.moveaxis(state_, [0, 1, 2, 3], [1, 2, 0, 3])

    return state__, header


def readState(stateFilePath, nocompact=False):
    stateFilePath = Path(stateFilePath)
    stateFile = open(stateFilePath, "rb")

    (forcing,) = unpack("=1i", stateFile.read(1 * 4))
    nx, ny, nz = unpack("=3i", stateFile.read(3 * 4))
    stateFile.close()

    ny_half = ny // 2

    if not nocompact:
        if nx - 1 >= ny_half and nx - 1 >= nz - 1:
            return readState_xcompact(stateFilePath)
        elif ny_half >= nx - 1 and ny_half >= nz - 1:
            return readState_ycompact(stateFilePath)
        else:
            return readState_zcompact(stateFilePath)
    else:
        return readState_nocompact(stateFilePath)


def writeState_xcompact(
    state,
    forcing=1,
    Lx=4,
    Lz=2,
    Re=628.3185307179584,
    tilt_angle=0.0,
    dt=0.03183098861837907,
    itime=0,
    time=0.0,
    outFile="state.000000",
):
    # Write array to a dnsbox state file.
    outFile = Path(outFile)
    stateFile = open(outFile, "wb")

    nx, ny_half, nz, _ = state.shape
    ny = ny_half * 2
    nxp, nzp = nx // 2 - 1, nz // 2 - 1

    stateFile.write(pack("=1i", forcing))
    stateFile.write(pack("=3i", nx, ny, nz))
    stateFile.write(pack("=2d", Lx, Lz))
    stateFile.write(pack("=1d", Re))
    stateFile.write(pack("=1d", tilt_angle))
    stateFile.write(pack("=1d", dt))
    stateFile.write(pack("=1i", itime))
    stateFile.write(pack("=1d", time))

    # convert to y, z, x
    state_ = np.moveaxis(state, [0, 1, 2], [2, 0, 1])

    # Convert to a real valued array
    stateOut = np.zeros((2 * ny_half, nz, nx, 3), dtype=np.float64, order="F")
    for j in range(ny_half):
        stateOut[2 * j, :, :, :] = state_[j, :, :, :].real
        stateOut[2 * j + 1, :, :, :] = state_[j, :, :, :].imag

    # write only the zero mode for u
    nCells = 2 * ny_half * (nz - 1)
    buff = np.zeros((2 * ny_half, nz - 1, 1), dtype=np.float64)
    buff[:, : nzp + 1, 0] = stateOut[:, : nzp + 1, 0, 0]
    buff[:, nzp + 1 :, 0] = stateOut[:, nzp + 2 :, 0, 0]
    dataPack = pack("={}d".format(nCells), *buff.flatten(order="F"))
    stateFile.write(dataPack)

    # write everything for v and w
    for n in range(1, 3):
        nCells = (nx - 1) * 2 * ny_half * (nz - 1)
        buff = np.zeros((2 * ny_half, nz - 1, nx - 1), dtype=np.float64)
        buff[:, : nzp + 1, : nxp + 1] = stateOut[:, : nzp + 1, : nxp + 1, n]
        buff[:, : nzp + 1, nxp + 1 :] = stateOut[:, : nzp + 1, nxp + 2 :, n]
        buff[:, nzp + 1 :, : nxp + 1] = stateOut[:, nzp + 2 :, : nxp + 1, n]
        buff[:, nzp + 1 :, nxp + 1 :] = stateOut[:, nzp + 2 :, nxp + 2 :, n]
        dataPack = pack("={}d".format(nCells), *buff.flatten(order="F"))
        stateFile.write(dataPack)

    stateFile.close()


def writeState_ycompact(
    state,
    forcing=1,
    Lx=4,
    Lz=2,
    Re=628.3185307179584,
    tilt_angle=0.0,
    dt=0.03183098861837907,
    itime=0,
    time=0.0,
    outFile="state.000000",
):
    # Write array to a dnsbox state file.
    outFile = Path(outFile)
    stateFile = open(outFile, "wb")

    nx, ny_half, nz, _ = state.shape
    ny = ny_half * 2
    nxp, nzp = nx // 2 - 1, nz // 2 - 1

    stateFile.write(pack("=1i", forcing))
    stateFile.write(pack("=3i", nx, ny, nz))
    stateFile.write(pack("=2d", Lx, Lz))
    stateFile.write(pack("=1d", Re))
    stateFile.write(pack("=1d", tilt_angle))
    stateFile.write(pack("=1d", dt))
    stateFile.write(pack("=1i", itime))
    stateFile.write(pack("=1d", time))

    # convert to y, z, x
    state_ = np.moveaxis(state, [0, 1, 2], [2, 0, 1])

    # Convert to a real valued array
    stateOut = np.zeros((2 * ny_half, nz, nx, 3), dtype=np.float64, order="F")
    for j in range(ny_half):
        stateOut[2 * j, :, :, :] = state_[j, :, :, :].real
        stateOut[2 * j + 1, :, :, :] = state_[j, :, :, :].imag

    # write everything for u
    nCells = (nx - 1) * 2 * ny_half * (nz - 1)
    buff = np.zeros((2 * ny_half, nz - 1, nx - 1), dtype=np.float64)
    buff[:, : nzp + 1, : nxp + 1] = stateOut[:, : nzp + 1, : nxp + 1, 0]
    buff[:, : nzp + 1, nxp + 1 :] = stateOut[:, : nzp + 1, nxp + 2 :, 0]
    buff[:, nzp + 1 :, : nxp + 1] = stateOut[:, nzp + 2 :, : nxp + 1, 0]
    buff[:, nzp + 1 :, nxp + 1 :] = stateOut[:, nzp + 2 :, nxp + 2 :, 0]
    dataPack = pack("={}d".format(nCells), *buff.flatten(order="F"))
    stateFile.write(dataPack)

    # write only the zero mode for v
    nCells = 2 * (nx - 1) * (nz - 1)
    buff = np.zeros((2, nz - 1, nx - 1), dtype=np.float64)
    buff[:2, : nzp + 1, : nxp + 1] = stateOut[:2, : nzp + 1, : nxp + 1, 1]
    buff[:2, : nzp + 1, nxp + 1 :] = stateOut[:2, : nzp + 1, nxp + 2 :, 1]
    buff[:2, nzp + 1 :, : nxp + 1] = stateOut[:2, nzp + 2 :, : nxp + 1, 1]
    buff[:2, nzp + 1 :, nxp + 1 :] = stateOut[:2, nzp + 2 :, nxp + 2 :, 1]
    dataPack = pack("={}d".format(nCells), *buff.flatten(order="F"))
    stateFile.write(dataPack)

    # write everything for w
    nCells = (nx - 1) * 2 * ny_half * (nz - 1)
    buff = np.zeros((2 * ny_half, nz - 1, nx - 1), dtype=np.float64)
    buff[:, : nzp + 1, : nxp + 1] = stateOut[:, : nzp + 1, : nxp + 1, 2]
    buff[:, : nzp + 1, nxp + 1 :] = stateOut[:, : nzp + 1, nxp + 2 :, 2]
    buff[:, nzp + 1 :, : nxp + 1] = stateOut[:, nzp + 2 :, : nxp + 1, 2]
    buff[:, nzp + 1 :, nxp + 1 :] = stateOut[:, nzp + 2 :, nxp + 2 :, 2]
    dataPack = pack("={}d".format(nCells), *buff.flatten(order="F"))
    stateFile.write(dataPack)

    stateFile.close()


def writeState_zcompact(
    state,
    forcing=1,
    Lx=4,
    Lz=2,
    Re=628.3185307179584,
    tilt_angle=0.0,
    dt=0.03183098861837907,
    itime=0,
    time=0.0,
    outFile="state.000000",
):
    # Write array to a dnsbox state file.
    outFile = Path(outFile)
    stateFile = open(outFile, "wb")

    nx, ny_half, nz, _ = state.shape
    ny = ny_half * 2
    nxp, nzp = nx // 2 - 1, nz // 2 - 1

    stateFile.write(pack("=1i", forcing))
    stateFile.write(pack("=3i", nx, ny, nz))
    stateFile.write(pack("=2d", Lx, Lz))
    stateFile.write(pack("=1d", Re))
    stateFile.write(pack("=1d", tilt_angle))
    stateFile.write(pack("=1d", dt))
    stateFile.write(pack("=1i", itime))
    stateFile.write(pack("=1d", time))

    # convert to y, z, x
    state_ = np.moveaxis(state, [0, 1, 2], [2, 0, 1])

    # Convert to a real valued array
    stateOut = np.zeros((2 * ny_half, nz, nx, 3), dtype=np.float64, order="F")
    for j in range(ny_half):
        stateOut[2 * j, :, :, :] = state_[j, :, :, :].real
        stateOut[2 * j + 1, :, :, :] = state_[j, :, :, :].imag

    # write everything for u and v
    for n in range(2):
        nCells = (nx - 1) * 2 * ny_half * (nz - 1)
        buff = np.zeros((2 * ny_half, nz - 1, nx - 1), dtype=np.float64)
        buff[:, : nzp + 1, : nxp + 1] = stateOut[:, : nzp + 1, : nxp + 1, n]
        buff[:, : nzp + 1, nxp + 1 :] = stateOut[:, : nzp + 1, nxp + 2 :, n]
        buff[:, nzp + 1 :, : nxp + 1] = stateOut[:, nzp + 2 :, : nxp + 1, n]
        buff[:, nzp + 1 :, nxp + 1 :] = stateOut[:, nzp + 2 :, nxp + 2 :, n]
        dataPack = pack("={}d".format(nCells), *buff.flatten(order="F"))
        stateFile.write(dataPack)

    # write only the zero mode for w
    nCells = (nx - 1) * 2 * ny_half
    buff = np.zeros((2 * ny_half, 1, nx - 1), dtype=np.float64)
    buff[:, 0, : nxp + 1] = stateOut[:, 0, : nxp + 1, 2]
    buff[:, 0, nxp + 1 :] = stateOut[:, 0, nxp + 2 :, 2]
    dataPack = pack("={}d".format(nCells), *buff.flatten(order="F"))
    stateFile.write(dataPack)

    stateFile.close()


def writeState_nocompact(
    state,
    forcing=1,
    Lx=4,
    Lz=2,
    Re=628.3185307179584,
    tilt_angle=0.0,
    dt=0.03183098861837907,
    itime=0,
    time=0.0,
    outFile="state.000000",
):
    # Write array to a dnsbox state file.
    outFile = Path(outFile)
    stateFile = open(outFile, "wb")

    nx, ny_half, nz, _ = state.shape
    ny = ny_half * 2
    nxp, nzp = nx // 2 - 1, nz // 2 - 1

    stateFile.write(pack("=1i", forcing))
    stateFile.write(pack("=3i", nx, ny, nz))
    stateFile.write(pack("=2d", Lx, Lz))
    stateFile.write(pack("=1d", Re))
    stateFile.write(pack("=1d", tilt_angle))
    stateFile.write(pack("=1d", dt))
    stateFile.write(pack("=1i", itime))
    stateFile.write(pack("=1d", time))

    # convert to y, z, x
    state_ = np.moveaxis(state, [0, 1, 2], [2, 0, 1])

    # Convert to a real valued array
    stateOut = np.zeros((2 * ny_half, nz, nx, 3), dtype=np.float64, order="F")
    for j in range(ny_half):
        stateOut[2 * j, :, :, :] = state_[j, :, :, :].real
        stateOut[2 * j + 1, :, :, :] = state_[j, :, :, :].imag

    for n in range(3):
        nCells = (nx - 1) * 2 * ny_half * (nz - 1)
        buff = np.zeros((2 * ny_half, nz - 1, nx - 1), dtype=np.float64)
        buff[:, : nzp + 1, : nxp + 1] = stateOut[:, : nzp + 1, : nxp + 1, n]
        buff[:, : nzp + 1, nxp + 1 :] = stateOut[:, : nzp + 1, nxp + 2 :, n]
        buff[:, nzp + 1 :, : nxp + 1] = stateOut[:, nzp + 2 :, : nxp + 1, n]
        buff[:, nzp + 1 :, nxp + 1 :] = stateOut[:, nzp + 2 :, nxp + 2 :, n]
        dataPack = pack("={}d".format(nCells), *buff.flatten(order="F"))
        stateFile.write(dataPack)

    stateFile.close()


def writeState(
    state,
    forcing=1,
    Lx=4,
    Lz=2,
    Re=628.3185307179584,
    tilt_angle=0.0,
    dt=0.03183098861837907,
    itime=0,
    time=0.0,
    outFile="state.000000",
):

    nx, ny_half, nz, _ = state.shape

    if nx - 1 >= ny_half and nx - 1 >= nz - 1:
        writeState_xcompact(
            state,
            forcing=forcing,
            Lx=Lx,
            Lz=Lz,
            Re=Re,
            tilt_angle=tilt_angle,
            dt=dt,
            itime=itime,
            time=time,
            outFile=outFile,
        )
    elif ny_half >= nx - 1 and ny_half >= nz - 1:
        writeState_ycompact(
            state,
            forcing=forcing,
            Lx=Lx,
            Lz=Lz,
            Re=Re,
            tilt_angle=tilt_angle,
            dt=dt,
            itime=itime,
            time=time,
            outFile=outFile,
        )
    else:
        writeState_zcompact(
            state,
            forcing=forcing,
            Lx=Lx,
            Lz=Lz,
            Re=Re,
            tilt_angle=tilt_angle,
            dt=dt,
            itime=itime,
            time=time,
            outFile=outFile,
        )


def inprod(state1, state2):
    # outputs only the real part of the inner product

    if state1.shape != state2.shape:
        exit("u1 and u2 have different shapes.")

    # DO: account for zeroed modes

    # account for double-counting due to hermiticity
    kin = np.sum(np.conj(state1[:, 0, :, :]) * state2[:, 0, :, :]) / 2
    kin += np.sum(np.conj(state1[:, 1:, :, :]) * state2[:, 1:, :, :])
    # the nyquist mode is not there after shrinking the grid

    kin = kin.real

    return kin


def wavenumbers(Lx, Lz, nx, ny_half, nz):

    kx = np.zeros((nx), dtype=np.float64)
    ky = np.zeros((ny_half), dtype=np.float64)
    kz = np.zeros((nz), dtype=np.float64)

    for i in range(nx):
        if i <= nx // 2:
            kx[i] = i * (2 * np.pi / Lx)
        else:
            kx[i] = i * (2 * np.pi / Lx) - nx * (2 * np.pi / Lx)

    for j in range(ny_half):
        ky[j] = j * (2 * np.pi / Ly)

    for k in range(nz):
        if k <= nz // 2:
            kz[k] = k * (2 * np.pi / Lz)
        else:
            kz[k] = k * (2 * np.pi / Lz) - nz * (2 * np.pi / Lz)

    return kx, ky, kz

def positions(Lx,Lz,nx,ny,nz):
    xs = np.array([i * (Lx / nx) for i in range(0, nx)])
    ys = np.array([j * (Ly / ny) for j in range(0, ny)])
    zs = np.array([k * (Lz / nz) for k in range(0, nz)])

    return xs, ys, zs


def derivative(subState, axis, Lx, Lz):
    nx, ny_half, nz = subState.shape
    nxp = nx // 2 - 1
    nzp = nz // 2 - 1
    kx, ky, kz = wavenumbers(Lx, Lz, nx, ny_half, nz)

    out = np.zeros((nx, ny_half, nz), dtype=np.complex128)
    if axis == 0:
        for i in range(nx):
            if i == nxp + 1:
                continue
            out[i, :, :] = 1j * kx[i] * subState[i, :, :]
    elif axis == 1:
        for j in range(ny_half):
            out[:, j, :] = 1j * ky[j] * subState[:, j, :]
    elif axis == 2:
        for k in range(nz):
            if k == nzp + 1:
                continue
            out[:, :, k] = 1j * kz[k] * subState[:, :, k]

    return out


def vorticity(state, Lx, Lz):
    wy = derivative(state[:, :, :, 2], 1, Lx, Lz)
    vz = derivative(state[:, :, :, 1], 2, Lx, Lz)
    uz = derivative(state[:, :, :, 0], 2, Lx, Lz)
    wx = derivative(state[:, :, :, 2], 0, Lx, Lz)
    vx = derivative(state[:, :, :, 1], 0, Lx, Lz)
    uy = derivative(state[:, :, :, 0], 1, Lx, Lz)

    out = np.zeros(state.shape, dtype=np.complex128)
    out[:, :, :, 0] = wy - vz
    out[:, :, :, 1] = uz - wx
    out[:, :, :, 2] = vx - uy

    return out


def dissipation(state, Lx, Lz, Re):
    vor = vorticity(state, Lx, Lz)
    diss = 2 * inprod(vor, vor) / Re

    return diss

def powerin(state, forcing, Re):
    # assumes no tilting
    if forcing == 1:
        return -(np.pi**2 / (4*Re)) * state[0,1,0,0].imag
    elif forcing == 2:
        return (np.pi**2 / (4*Re)) * state[0,1,0,0].real


def fftSpecToPhys(subState, supersample=False):
    nx, ny_half, nz = subState.shape
    nxp, nzp = nx // 2 - 1, nz // 2 - 1

    if supersample:
        nxx, nyy, nzz = nStretch * (nx // 2), nStretch * ny_half, nStretch * (nz // 2)
        nyy_half_pad1 = nyy // 2 + 1

        expandZ = np.zeros((nx, ny_half, nzz), dtype=np.complex128)
        for k in range(nzp + 1):
            expandZ[:, :, k] = subState[:, :, k]
            if 0 < k <= nzp:
                expandZ[:, :, -k] = subState[:, :, -k]
        expandZ = np.fft.ifft(expandZ, axis=2, norm="forward", n=nzz)

        expandX = np.zeros((nxx, ny_half, nzz), dtype=np.complex128)
        for i in range(nxp + 1):
            expandX[i, :, :] = expandZ[i, :, :]
            if 0 < i <= nxp:
                expandX[-i, :, :] = expandZ[-i, :, :]
        expandX = np.fft.ifft(expandX, axis=0, norm="forward", n=nxx)

        expandY = np.zeros((nxx, nyy_half_pad1, nzz), dtype=np.complex128)
        expandY[:, :ny_half, :] = expandX
        stateOut = np.fft.irfft(expandY, axis=1, norm="forward", n=nyy)
    else:
        nxx, nyy, nzz = nx, 2 * ny_half, nz
        stateOut = np.swapaxes(
            np.fft.irfftn(
                np.swapaxes(subState, 1, 2), s=(nxx, nzz, nyy), norm="forward"
            ),
            1,
            2,
        )

    return stateOut


def fftSpecToPhysAll(state, supersample=False):

    nx, ny_half, nz, _ = state.shape

    if supersample:
        nxx, nyy, nzz = nStretch * (nx // 2), nStretch * (ny_half), nStretch * (nz // 2)
    else:
        nxx, nyy, nzz = nx, 2 * ny_half, nz

    stateOutAll = np.zeros((nxx, nyy, nzz, 3), dtype=np.float64)

    for comp in range(3):
        stateOutAll[:, :, :, comp] = fftSpecToPhys(
            state[:, :, :, comp], supersample=supersample
        )

    return stateOutAll


def fftPhysToSpec(stateOut, supersample=False):

    nxx, nyy, nzz = stateOut.shape

    if supersample:
        nx, ny, nz = 2 * (nxx // nStretch), 2 * (nyy // nStretch), 2 * (nzz // nStretch)
        ny_half = ny // 2
        nxp, nzp = nx // 2 - 1, nz // 2 - 1

        # FFT in y
        expandY = np.fft.rfft(stateOut, axis=1, norm="forward", n=nyy)
        # Shrink in y
        expandX = expandY[:, :ny_half, :]

        # FFT in x
        expandX = np.fft.fft(expandX, axis=0, norm="forward", n=nxx)
        # Shrink in x
        expandZ = np.zeros((nx, ny_half, nzz), dtype=np.complex128)
        for i in range(nxp + 1):
            expandZ[i, :, :] = expandX[i, :, :]
            if 0 < i <= nxp:
                expandZ[-i, :, :] = expandX[-i, :, :]

        # FFT in z
        expandZ = np.fft.fft(expandZ, axis=2, norm="forward", n=nzz)
        # Shrink in z
        subState = np.zeros((nx, ny_half, nz), dtype=np.complex128)
        for k in range(nzp + 1):
            subState[:, :, k] = expandZ[:, :, k]
            if 0 < k <= nzp:
                subState[:, :, -k] = expandZ[:, :, -k]
    else:
        nx, ny, nz = nxx, nyy, nzz
        ny_half = ny // 2
        nxp, nzp = nx // 2 - 1, nz // 2 - 1

        subState = np.swapaxes(
            np.fft.rfftn(
                np.swapaxes(stateOut, 1, 2), s=(nxx, nzz, nyy), norm="forward"
            ),
            1,
            2,
        )
        # Zero the Nyquist modes
        subState = subState[:, :ny_half, :]
        subState[nxp + 1, :, :] = 0
        subState[:, :, nzp + 1] = 0

    return subState


def fftPhysToSpecAll(stateOutAll, supersample=False):

    nxx, nyy, nzz, _ = stateOutAll.shape

    if supersample:
        nx, ny, nz = 2 * (nxx // nStretch), 2 * (nyy // nStretch), 2 * (nzz // nStretch)
    else:
        nx, ny, nz = nxx, nyy, nzz

    ny_half = ny // 2

    state = np.zeros((nx, ny_half, nz, 3), dtype=np.complex128)

    for comp in range(3):
        state[:, :, :, comp] = fftPhysToSpec(
            stateOutAll[:, :, :, comp], supersample=supersample
        )

    return state


def tilt_state(state, tilt_angle):

    nx, ny_half, nz, _ = state.shape
    tilt_radians = tilt_angle * np.pi / 180.0

    stateout = np.zeros((nx, ny_half, nz, 3), dtype=np.complex128)
    stateout[:, :, :, 1] = state[:, :, :, 1]
    stateout[:, :, :, 0] = (
        np.cos(tilt_radians) * state[:, :, :, 0]
        - np.sin(tilt_radians) * state[:, :, :, 2]
    )
    stateout[:, :, :, 2] = (
        np.sin(tilt_radians) * state[:, :, :, 0]
        + np.cos(tilt_radians) * state[:, :, :, 2]
    )

    return stateout


def laminar(forcing, nx, ny_half, nz, tilt_angle=0.0):

    state = np.zeros((nx, ny_half, nz, 3), dtype=np.complex128)
    if forcing == 1:
        state[0, qF, 0, 0] = -1j * 0.5
    elif forcing == 2:
        state[0, qF, 0, 0] = 0.5

    if abs(tilt_angle) > 0:
        state = tilt_state(state, tilt_angle)

    return state
