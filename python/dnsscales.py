#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
from tqdm import tqdm

import dns


def main():
    parser = argparse.ArgumentParser(
        "Compute various scales, averaging over many states."
    )
    parser.add_argument(
        "statesPath",
        type=str,
        help="path to the folder containing states of interest.",
    )
    parser.add_argument(
        "si",
        type=int,
        help="initial state to average.",
    )
    parser.add_argument(
        "sf",
        type=int,
        help="final state to average.",
    )
    parser.add_argument(
        "--savetxt",
        action="store_true",
        dest="savetxt",
        help="save results as text files",
    )

    args = vars(parser.parse_args())
    statesPath = Path(args["statesPath"])
    si = args["si"]
    sf = args["sf"]
    savetxt = args["savetxt"]

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

    print("Averaging over fluctuations.")
    dissipations = None
    velsquareds = None

    for i in tqdm(range(si, sf + 1)):
        state_i = statesPath / f"state.{str(i).zfill(6)}"
        state, header = dns.readState(state_i)
        forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = header

        dissipation = dns.dissipation(state - state_avg, Lx, Lz, Re)
        velphys = dns.fftSpecToPhysAll(state - state_avg)

        if dissipations is None:
            dissipations = dissipation
        else:
            dissipations += dissipation

        if velsquareds is None:
            velsquareds = velphys ** 2
        else:
            velsquareds += velphys ** 2

    dissipation = dissipations / nstates
    velsquareds = velsquareds / nstates

    # Average over space and the three directions
    velsquared = np.average(velsquareds)
    urms = np.sqrt(velsquared)

    print("Re =", Re)
    print("ϵ =", dissipation)
    print("<u'>_rms = ", urms)

    scales = np.array([dissipation, urms])
    if savetxt:
        np.savetxt(statesPath / "scales.gp", scales)

    # Taylor microscale (\lambda)
    lambda_g = urms * np.sqrt(15 / (Re * dissipation))
    # Re_\lambda
    Re_lambda = urms * lambda_g * Re

    print("λ =", lambda_g)
    print("Re_λ =", Re_lambda)

    # Kolmogorov microscales
    eta = np.power(dissipation * Re ** 3, -1 / 4)
    tau = np.power(Re * dissipation, -1 / 2)
    ueta = np.power(dissipation / Re, 1 / 4)

    print("η =", eta)
    print("τ_η =", tau)
    print("u_η = ", ueta)


if __name__ == "__main__":
    main()