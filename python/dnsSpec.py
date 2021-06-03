#!/usr/bin/env python3
import argparse
from pathlib import Path
from math import floor

import numpy as np
from matplotlib import pyplot as plt

import dnsCommon as dns

figuresDirName = "figures"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("statePath", type=str, help="path to the state of interest.")
    args = vars(parser.parse_args())

    specs_xy, spec_z, speciso, kmax_xy, kmax_z = dnsSpec(**args)
    wavenums_xy = np.arange(specs_xy.shape[-1])
    wavenums_z = np.arange(spec_z.shape[-1])
    wavenumsiso = np.arange(speciso.shape[-1]) + 1

    fig, ax = plt.subplots()
    ax.plot(wavenums_xy, specs_xy[0, :])
    ax.grid(True, which="both")
    ax.set_xlabel("$k_x$")
    ax.set_ylabel("$E_k$")
    ax.set_yscale("log")
    ax.set_xticks(wavenums_xy[::4])
    fig.savefig("spec-x.png")

    fig, ax = plt.subplots()
    ax.plot(wavenums_xy, specs_xy[1, :])
    ax.grid(True, which="both")
    ax.set_xlabel("$k_y$")
    ax.set_ylabel("$E_k$")
    ax.set_yscale("log")
    ax.set_xticks(wavenums_xy[::4])
    fig.savefig("spec-y.png")

    fig, ax = plt.subplots()
    # Skip the odd modes
    ax.plot(wavenums_z[::2], spec_z[::2])
    ax.grid(True, which="both")
    ax.set_xlabel("Even $k_z$")
    ax.set_ylabel("$E_k$")
    ax.set_yscale("log")
    ax.set_xticks(wavenums_z[::4])
    fig.savefig("spec-z.png")

    fig, ax = plt.subplots()
    # Skip the odd modes
    ax.plot(wavenums_z, spec_z)
    ax.grid(True, which="both")
    ax.set_xlabel("$k_z$")
    ax.set_ylabel("$E_k$")
    ax.set_yscale("log")
    ax.set_xticks(wavenums_z[::4])
    fig.savefig("spec-z.png")

    fig, ax = plt.subplots()
    # Skip the kx=ky=kz=0 mode
    ax.loglog(wavenumsiso, speciso[1:])
    ax.grid(True, which="both")
    ax.set_xlabel("$k$")
    ax.set_ylabel("$E_k$")
    ax.set_xticks(wavenumsiso[::4])
    fig.savefig("spec-iso.png")

    plt.show()


def dnsSpec(statePath, Lz=1):

    statePath = Path(statePath).resolve()
    dns.checkFile(statePath)

    state, header = dns.readState(statePath)
    nz, ny, nx = header[:3]
    akx, aky, akz = dns.wavenumbers(nz, ny, nx)
    scales = np.sqrt(dns.inprodCoefficients(nz, ny, nx))

    # Embed the scales into the state for ease
    state *= scales

    rnx = 2 * (ny // 2 - 1) / 3
    rnz = 2 * ((nz // 2) * 2 / Lz - 1) / 3
    kmax_xy = int(floor(rnx))
    kmax_z = int(floor(rnz))
    iso_max = int(round(np.sqrt(2 * kmax_xy ** 2 + kmax_z ** 2)))
    specs_xy = np.zeros((2, kmax_xy + 1))
    spec_z = np.zeros((kmax_z + 1))
    speciso = np.zeros((iso_max + 1))
    for i in range(0, nz + 1, 2):
        kz = int(np.abs(akz[i]))
        if kz <= rnz:
            for k in range(0, ny):
                kx = int(np.abs(akx[k]))
                if kx <= rnx:
                    for j in range(0, nx):
                        ky = int(np.abs(aky[j]))
                        if ky <= rnz:
                            part = np.sum(
                                state[i, k, j, :] ** 2 + state[i + 1, k, j, :] ** 2
                            )
                            specs_xy[0, kx] += part
                            specs_xy[1, ky] += part
                            spec_z[kz] += part

                            iso = int(round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)))
                            speciso[iso] += part

    return specs_xy, spec_z, speciso, kmax_xy, kmax_z


if __name__ == "__main__":
    main()
