#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

import dnsbox as dns


def main():

    parser = argparse.ArgumentParser(
        description="Change the resolution of a state file by padding with zeros or removing modes.",
        prog="dnsbox change resolution",
    )
    parser.add_argument("stateFileIn", type=str, help="path to the input state.")
    parser.add_argument("stateFileOut", type=str, help="path to the output state.")
    parser.add_argument("nx", type=int, help="new nx.")
    parser.add_argument("ny", type=int, help="new ny.")
    parser.add_argument("nz", type=int, help="new nz.")

    args = vars(parser.parse_args())
    dnschangeres(**args)


def dnschangeres(stateFileIn, stateFileOut, nx, ny, nz):
    """Zero pads the velocity field to a larger grid for upscaling
    or removes modes for downscaling.

    """

    stateFileIn = Path(stateFileIn)
    stateFileOut = Path(stateFileOut)
    u, header = dns.readState(stateFileIn)
    forcing, nx0, ny0, nz0, Lx, Lz, Re, tilt_angle, dt, itime, time = header

    nx_ = min(nx, nx0)
    ny_ = min(ny, ny0)
    nz_ = min(nz, nz0)

    nxp_ = nx_ // 2 - 1
    ny_half = ny // 2
    ny_half_ = ny_ // 2
    nzp_ = nz_ // 2 - 1

    u_pad = np.zeros((nx, ny_half, nz, 3), dtype=np.complex128)

    # Both x and z positive
    u_pad[: nxp_ + 1, :ny_half_, : nzp_ + 1] = u[: nxp_ + 1, :ny_half_, : nzp_ + 1]

    # Both x and z negative
    u_pad[-nxp_:, :ny_half_, -nzp_:] = u[-nxp_:, :ny_half_, -nzp_:]

    # x positive, z negative
    u_pad[: nxp_ + 1, :ny_half_, -nzp_:] = u[: nxp_ + 1, :ny_half_, -nzp_:]

    # x negative, z positive
    u_pad[-nxp_:, :ny_half_, : nzp_ + 1] = u[-nxp_:, :ny_half_, : nzp_ + 1]

    dns.writeState(
        u_pad,
        forcing=forcing,
        Lx=Lx,
        Lz=Lz,
        Re=Re,
        tilt_angle=tilt_angle,
        dt=dt,
        itime=itime,
        time=time,
        outFile=stateFileOut,
    )


if __name__ == "__main__":
    main()
