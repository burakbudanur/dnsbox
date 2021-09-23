#!/usr/bin/env python3
import argparse
from os import wait
from pathlib import Path

import numpy as np

import dns


def main():

    parser = argparse.ArgumentParser(
        description="Save to disk the physical space version of the input snapshot.",
    )
    parser.add_argument("state", type=str, help="path to the state.")
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

    state = Path(args["state"])
    mirror_y = args["mirror_y"]
    mirror_z = args["mirror_z"]

    stateIn, headers = dns.readState(state)

    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers

    if not mirror_y:
        ny_display = ny
    else:
        ny_display = ny // 2 + 1
    if not mirror_z:
        nz_display = nz
    else:
        nz_display = nz // 2 + 1

    vel_phys = dns.fftSpecToPhysAll(stateIn, supersample=True)

    # Compute pressure
    nxx, nyy, nzz = 3 * (nx // 2), 3 * (ny // 2), 3 * (nz // 2)
    di_uj = np.zeros((nxx, nyy, nzz, 3, 3), dtype=np.float64)
    for i in range(3):
        for j in range(3):
            di_uj[:, :, :, i, j] = dns.fftSpecToPhys(
                dns.derivative(stateIn[:, :, :, j], i, Lx, Lz), supersample=True
            )
    # Nonlinear term
    nj = np.zeros((nx, ny // 2, nz, 3), dtype=np.complex128)
    for j in range(3):
        nj[:, :, :, j] = -dns.fftPhysToSpec(
            np.sum(vel_phys * di_uj[:, :, :, :, j], axis=3), supersample=True
        )
    kx, ky, kz = dns.wavenumbers(Lx, Lz, nx, ny // 2, nz)
    p_spec = np.zeros((nx, ny // 2, nz), dtype=np.complex128)
    for i, kx_i in enumerate(kx):
        for j, ky_j in enumerate(ky):
            for k, kz_k in enumerate(kz):
                if not (i == 0 and j == 0 and k == 0):
                    p_spec[i, j, k] = (
                        -1j
                        * (
                            kx_i
                            + nj[i, j, k, 0]
                            + ky_j
                            + nj[i, j, k, 1]
                            + kz_k
                            + nj[i, j, k, 2]
                        )
                        / (kx_i ** 2 + ky_j ** 2 + kz_k ** 2)
                    )
    p_phys = dns.fftSpecToPhys(p_spec)
    out = np.zeros((nx, ny_display, nz_display, 4), dtype=np.float64)
    out[:, :, :, :3] = dns.fftSpecToPhysAll(stateIn)[:, :ny_display, :nz_display, :]
    out[:, :, :, 3] = p_phys[:, :ny_display, :nz_display]
    np.save(state.parent / f"{state.name}_phys", out)


if __name__ == "__main__":
    main()
