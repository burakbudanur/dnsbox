#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

import dns
import dns_symmetries as dnss


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("statefile", type=str, help="path to the state to slice.")
    parser.add_argument("slicedir", type=str, help="folder containing slice templates.")
    parser.add_argument("savedir", type=str, help="directory to save the file to.")
    args = vars(parser.parse_args())
    dnssymred(**args)


def dnssymred(statefile, slicedir, savedir):
    statefile = Path(statefile)
    slicedir = Path(savedir)
    savedir = Path(savedir)
    txfile = slicedir / "u_xp.000000"
    tzfile = slicedir / "u_zp.000000"

    state, header = dns.readState(statefile)
    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = header
    tx, header = dns.readState_nocompact(txfile)
    tz, header = dns.readState_nocompact(tzfile)

    tx_shifted = dnss.Tx(Lx / 4, tx, Lx, Lz)
    tz_shifted = dnss.Tz(Lz / 4, tz, Lx, Lz)

    px_r = dns.inprod(state, tx)
    px_i = dns.inprod(state, tx_shifted)
    px = np.arctan2(px_i, px_r)
    shiftx = px / (2 * np.pi) * Lx
    state_ = dnss.Tx(-shiftx, state, Lx, Lz)

    pz_r = dns.inprod(state_, tz)
    pz_i = dns.inprod(state_, tz_shifted)
    pz = np.arctan2(pz_i, pz_r)
    shiftz = pz / (2 * np.pi) * Lz
    state_ = dnss.Tz(-shiftz, state_, Lx, Lz)

    statefile_ = savedir / f"sliced_{statefile.name}"
    dns.writeState(
        state_,
        forcing=forcing,
        Lx=Lx,
        Lz=Lz,
        Re=Re,
        tilt_angle=tilt_angle,
        dt=dt,
        itime=itime,
        time=time,
        outFile=statefile_,
    )


if __name__ == "__main__":
    main()
