#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
from tqdm import tqdm

import dns


def main():
    parser = argparse.ArgumentParser(
        "Compute all nonzero Reynolds stresses and terms that go into their budget."
    )
    parser.add_argument(
        "statesPath",
        type=str,
        help="path to the folder containing states of interest.",
    )
    parser.add_argument(
        "si", type=int, help="initial state to average.",
    )
    parser.add_argument(
        "sf", type=int, help="final state to average.",
    )
    parser.add_argument(
        "--savetxt",
        action="store_true",
        dest="savetxt",
        help="save results as text files",
    )
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
    dnsstresses(**args)


def dnsstresses(statesPath, si, sf, savetxt=False, mirror_y=False, mirror_z=False):

    statesPath = Path(statesPath)

    nstates = sf - si + 1

    print("Finding the average state.")
    state_avg = None
    for i in tqdm(range(si, sf + 1)):
        state_i = statesPath / f"state.{str(i).zfill(6)}"
        state, header = dns.readState(state_i)
        if state_avg is None:
            state_avg = state
        else:
            state_avg += state
    state_avg /= nstates

    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = header
    ny_half = ny // 2
    if not mirror_y:
        ny_display = ny
    else:
        ny_display = ny // 2 + 1
    if not mirror_z:
        nz_display = nz
    else:
        nz_display = nz // 2 + 1

    dns.writeState(
        state_avg,
        forcing=forcing,
        Lx=Lx,
        Lz=Lz,
        Re=Re,
        tilt_angle=tilt_angle,
        outFile=statesPath / "mean_state.000000",
    )
    # Velocity profile
    mean_profile = np.average(
        dns.fftSpecToPhys(state_avg[:, :, :, 0])[:, :ny_display, :nz_display],
        axis=(0, 2),
    )
    np.savetxt(statesPath / "mean_profile.gp", mean_profile)
    mean_shear = np.average(
        dns.fftSpecToPhys(dns.derivative(state_avg[:, :, :, 0], 1, Lx, Lz))[
            :, :ny_display, :nz_display
        ],
        axis=(0, 2),
    )
    np.savetxt(statesPath / "mean_shear.gp", mean_shear)

    print("Averaging over fluctuations.")
    u2_y_avg = None
    v2_y_avg = None
    w2_y_avg = None
    uv_y_avg = None
    diss_y_avg = None

    for i_state in tqdm(range(si, sf + 1)):
        state_i = statesPath / f"state.{str(i_state).zfill(6)}"
        state, header = dns.readState(state_i)

        velphys = dns.fftSpecToPhysAll(state - state_avg)[
            :, :ny_display, :nz_display, :
        ]
        u2_y = np.average(velphys[:, :, :, 0] ** 2, axis=(0, 2))
        v2_y = np.average(velphys[:, :, :, 1] ** 2, axis=(0, 2))
        w2_y = np.average(velphys[:, :, :, 2] ** 2, axis=(0, 2))
        uv_y = np.average(velphys[:, :, :, 0] * velphys[:, :, :, 1], axis=(0, 2))

        d_ui_d_xj = np.zeros((nx, ny_half, nz, 3, 3), dtype=np.complex128)
        for i in range(3):
            for j in range(3):
                d_ui_d_xj[:, :, :, i, j] = dns.derivative(
                    (state - state_avg)[:, :, :, i], j, Lx, Lz
                )

        S_ij = np.zeros((nx, ny_display, nz_display, 3, 3), dtype=np.float64)
        for i in range(3):
            # this can be made faster by only using the symmetric component
            for j in range(3):
                S_ij[:, :, :, i, j] = (
                    0.5
                    * dns.fftSpecToPhys(
                        d_ui_d_xj[:, :, :, i, j] + d_ui_d_xj[:, :, :, j, i]
                    )[:, :ny_display, :nz_display]
                )
        diss_y = (2 / Re) * np.average(np.sum(S_ij ** 2, axis=(3, 4)), axis=(0, 2))

        if u2_y_avg is None:
            u2_y_avg = u2_y
        else:
            u2_y_avg += u2_y

        if v2_y_avg is None:
            v2_y_avg = v2_y
        else:
            v2_y_avg += v2_y

        if w2_y_avg is None:
            w2_y_avg = w2_y
        else:
            w2_y_avg += w2_y

        if uv_y_avg is None:
            uv_y_avg = uv_y
        else:
            uv_y_avg += uv_y

        if diss_y_avg is None:
            diss_y_avg = diss_y
        else:
            diss_y_avg += diss_y

    u2_y_avg /= nstates
    v2_y_avg /= nstates
    w2_y_avg /= nstates
    uv_y_avg /= nstates
    diss_y_avg /= nstates

    np.savetxt(statesPath / "mean_u2.gp", u2_y_avg)
    np.savetxt(statesPath / "mean_v2.gp", v2_y_avg)
    np.savetxt(statesPath / "mean_w2.gp", w2_y_avg)
    np.savetxt(statesPath / "mean_uv.gp", uv_y_avg)
    np.savetxt(statesPath / "mean_diss.gp", diss_y_avg)

    return mean_profile, mean_shear, u2_y_avg, v2_y_avg, w2_y_avg, uv_y_avg, diss_y_avg


if __name__ == "__main__":
    main()
