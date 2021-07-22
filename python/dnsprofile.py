#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import dns


def main():

    parser = argparse.ArgumentParser(
        description="Produce spanwise-averaged plots of the streamwise "
        "components of velocity and vorticity.",
    )
    parser.add_argument("state", type=str, help="path to the/a state.")
    parser.add_argument(
        "--tex", action="store_true", dest="tex", help="use LaTeX to render text."
    )
    parser.add_argument(
        "--noshow", action="store_true", dest="noshow", help="do not display the plots."
    )
    parser.add_argument(
        "--average",
        action="store_true",
        dest="average",
        help="average multiple states.",
    )
    parser.add_argument(
        "--si", type=int, help="initial state to average.",
    )
    parser.add_argument(
        "--sf", type=int, help="final state to average.",
    )
    parser.add_argument(
        "--mirror_y",
        action="store_true",
        dest="mirror_y",
        help="display the fundamental domain of mirror_y.",
    )
    parser.add_argument(
        "--savetxt",
        action="store_true",
        dest="savetxt",
        help="save results as a text file.",
    )
    args = vars(parser.parse_args())
    dnsprofile(**args)


def dnsprofile(
    state,
    tex=False,
    noshow=False,
    average=False,
    si=None,
    sf=None,
    mirror_y=False,
    savetxt=False,
):
    state = Path(state)
    figuresDir = dns.createFiguresDir(state.parent)

    if not average:
        # if FFT'ing state files in Python
        stateIn, headers = dns.readState(state)
        forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers

        # measure along the forcing
        if abs(tilt_angle) > 0:
            stateIn = dns.tilt_state(stateIn, tilt_angle)

    else:
        # average many states
        statesDir = state.parent
        states = [f"state.{str(i).zfill(6)}" for i in range(si, sf + 1)]
        stateIn, headers = dns.readState(statesDir / states[0])
        forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers

        for stateName in states:
            stateIn_, headers_ = dns.readState(statesDir / stateName)

            # measure along the forcing
            if abs(tilt_angle) > 0:
                stateIn_ = dns.tilt_state(stateIn_, tilt_angle)
            stateIn += stateIn_

        # profile of the average
        stateIn /= len(states)

    velx_phys = dns.fftSpecToPhys(stateIn[:, :, :, 0])

    # average over z
    velx_phys_zavg = np.average(velx_phys[:, :, :], axis=2)

    # and average over x
    velx_profile = np.average(velx_phys_zavg, axis=0)

    kF = (2 * np.pi / dns.Ly) * dns.qF
    ys = np.array([j * (dns.Ly / ny) for j in range(0, ny)])
    lamvx = np.sin(kF * ys)

    yLabel = "$y$"

    dns.setPlotDefaults(tex=tex)

    ny_display = ny
    if mirror_y:
        ny_display = ny // 2 + 1
    else:
        ny_display = ny

    # profiles
    fig, ax = plt.subplots()
    ax.plot(velx_profile[:ny_display], ys[:ny_display])
    ax.plot(lamvx[:ny_display], ys[:ny_display], linestyle="--")
    ax.set_xlabel(f"$\\langle u \\rangle_{{xz}}$")
    ax.set_ylabel(yLabel)
    ax.set_ylim(bottom=ys[0], top=ys[ny_display - 1])
    fig.savefig(figuresDir / f"{state.name}_velocity_profile.png")
    if savetxt:
        ys_velx = np.zeros((ny_display, 2))
        ys_velx[:, 0] = ys[:ny_display]
        ys_velx[:, 1] = velx_profile[:ny_display]
        np.savetxt(figuresDir / f"{state.name}_velocity_profile.dat", ys_velx)

    if not noshow:
        plt.show()


if __name__ == "__main__":
    main()
