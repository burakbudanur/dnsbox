#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
from tqdm import tqdm

import dns


def main():
    parser = argparse.ArgumentParser(
        "Rate of strain and vorticity tensors."
    )
    parser.add_argument(
        "fstate",
        type=str,
        help="path of the state file",
    )

    args = vars(parser.parse_args())
    fstate = args['fstate']
    fstate = Path(fstate)
    state, header = dns.readState(fstate)

    S_ij, Om_ij = dnstensors(fstate)

    fstatenum = str(fstate)[-6:]
    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = header

    for i in range(3):
        dns.writeState(dns.fftPhysToSpec(S_ij[:,:,:,i,:]),
            forcing=forcing,
            Lx=Lx,
            Lz=Lz,
            Re=Re,
            tilt_angle=tilt_angle,
            dt=dt,
            itime=0,
            time=0,
            outfile=fstate.parent / f"S{i}.{fstatenum}"
        )

        dns.writeState(dns.fftPhysToSpec(Om_ij[:,:,:,i,:]),
            forcing=forcing,
            Lx=Lx,
            Lz=Lz,
            Re=Re,
            tilt_angle=tilt_angle,
            dt=dt,
            itime=0,
            time=0,
            outfile=fstate.parent / f"Om{i}.{fstatenum}"
        )


def dnstensors(fstate):

    fstate = Path(fstate)

    state, header = dns.readState(fstate)

    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = header
    ny_half = ny // 2

    d_ui_d_xj = np.zeros((nx, ny_half, nz, 3, 3), dtype=np.complex128)
    for i in range(3):
        for j in range(3):
            d_ui_d_xj[:, :, :, i, j] = dns.derivative(
                state[:, :, :, i], j, Lx, Lz
            )

    S_ij = np.zeros((nx, ny, nz, 3, 3), dtype=np.float64)
    for i in range(3):
        # this can be made faster by only using the symmetric component
        for j in range(3):
            S_ij[:, :, :, i, j] = (
                0.5
                * dns.fftSpecToPhys(
                    d_ui_d_xj[:, :, :, i, j] + d_ui_d_xj[:, :, :, j, i]
                )
            )

    Om_ij = np.zeros((nx, ny, nz, 3, 3), dtype=np.float64)
    for i in range(3):
        # this can be made faster by only using the symmetric component
        for j in range(3):
            
            if i == j:
                continue

            Om_ij[:, :, :, i, j] = (
                0.5
                * dns.fftSpecToPhys(
                    d_ui_d_xj[:, :, :, i, j] - d_ui_d_xj[:, :, :, j, i]
                )
            )
    
    return S_ij, Om_ij

if __name__ == "__main__":
    main()
