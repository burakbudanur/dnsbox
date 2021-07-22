#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import dns

# harcoded dropoff harmonics
dnki = 1
dnkf = 0


def main():
    parser = argparse.ArgumentParser("Computes direction-dependent spectra.")
    parser.add_argument("statePath", type=str, help="path to the state of interest.")
    parser.add_argument(
        "--tex", action="store_true", dest="tex", help="use LaTeX to render text."
    )
    parser.add_argument(
        "--noshow", action="store_true", dest="noshow", help="do not display the plots."
    )
    parser.add_argument(
        "--harmonics",
        action="store_true",
        dest="harmonics",
        help="plot in terms of harmonics.",
    )
    parser.add_argument(
        "--iso",
        action="store_true",
        dest="iso",
        help="compute isotropic spectrum as well.",
    )
    parser.add_argument(
        "--savetxt",
        action="store_true",
        dest="savetxt",
        help="save results as text files"
    )

    args = vars(parser.parse_args())
    statePath = Path(args["statePath"])
    tex = args["tex"]
    noshow = args["noshow"]
    harmonics = args["harmonics"]
    iso = args["iso"]
    savetxt = args["savetxt"]

    figuresDir = dns.createFiguresDir(statePath.parent)

    if iso:
        (spec_x, spec_y, spec_z, spec_iso, Lx, Lz, dissipation,) = dnsspec(
            statePath, iso=iso
        )
    else:
        (spec_x, spec_y, spec_z, Lx, Lz,) = dnsspec(statePath, iso=iso)

    if not harmonics:
        wavenums_x = np.arange(spec_x.shape[0]) * (2 * np.pi / Lx)
        wavenums_y = np.arange(spec_y.shape[0]) * (2 * np.pi / dns.Ly)
        wavenums_z = np.arange(spec_z.shape[0]) * (2 * np.pi / Lz)
    else:
        wavenums_x = np.arange(spec_x.shape[0])
        wavenums_y = np.arange(spec_y.shape[0])
        wavenums_z = np.arange(spec_z.shape[0])

    dns.setPlotDefaults(tex=tex)

    # log log versions
    fig, ax = plt.subplots()
    ax.plot(wavenums_x[1:], spec_x[1:])
    ax.grid(True, which="both")
    if not harmonics:
        ax.set_xlabel(f"$k_x$")
        ax.set_ylabel(f"$E_{{k_x}}$")
    else:
        ax.set_xlabel(f"$n_x$")
        ax.xaxis.get_major_locator().set_params(integer=True)
        ax.set_ylabel(f"$E_{{n_x}}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    fig.savefig(figuresDir / f"{statePath.name[-6:]}_spec_x_log.png")
    # print("dropoff-x", np.sqrt(spec_x[dnki] / spec_x[-1 + dnkf]))
    
    k_spec_x = np.zeros((spec_x.shape[0], 2))
    k_spec_x[:, 0] = wavenums_x
    k_spec_x[:, 1] = spec_x
    if savetxt:
        np.savetxt(figuresDir / f"{statePath.name[-6:]}_spec_x.dat", k_spec_x)

    fig, ax = plt.subplots()
    ax.plot(wavenums_y[1:], spec_y[1:])
    ax.grid(True, which="both")
    if not harmonics:
        ax.set_xlabel(f"$k_y$")
        ax.set_ylabel(f"$E_{{k_y}}$")
    else:
        ax.set_xlabel(f"$n_y$")
        ax.xaxis.get_major_locator().set_params(integer=True)
        ax.set_ylabel(f"$E_{{n_y}}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    fig.savefig(figuresDir / f"{statePath.name[-6:]}_spec_y_log.png")
    
    k_spec_y = np.zeros((spec_y.shape[0], 2))
    k_spec_y[:, 0] = wavenums_y
    k_spec_y[:, 1] = spec_y
    if savetxt:
        np.savetxt(figuresDir / f"{statePath.name[-6:]}_spec_y.dat", k_spec_y)

    # print("dropoff-y", np.sqrt(spec_y[dnki] / spec_y[-1 + dnkf]))

    fig, ax = plt.subplots()
    ax.plot(wavenums_z[1:], spec_z[1:])
    ax.grid(True, which="both")
    if not harmonics:
        ax.set_xlabel(f"$k_z$")
        ax.set_ylabel(f"$E_{{k_z}}$")
    else:
        ax.set_xlabel(f"$n_z$")
        ax.xaxis.get_major_locator().set_params(integer=True)
        ax.set_ylabel(f"$E_{{n_z}}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    fig.savefig(figuresDir / f"{statePath.name[-6:]}_spec_z_log.png")
    # print("dropoff-z", np.sqrt(spec_z[dnki] / spec_z[-1 + dnkf]))

    k_spec_z = np.zeros((spec_z.shape[0], 2))
    k_spec_z[:, 0] = wavenums_z
    k_spec_z[:, 1] = spec_z
    if savetxt:
        np.savetxt(figuresDir / f"{statePath.name[-6:]}_spec_z.dat", k_spec_z)

    if iso:
        fig, ax = plt.subplots()
        ax.plot(spec_iso[:, 0], spec_iso[:, 1])
        ax.plot(
            spec_iso[:, 0],
            1.5 * (dissipation ** (2 / 3)) * np.power(spec_iso[:, 0], -5 / 3),
            color="k",
        )
        ax.grid(True, which="both")
        if not harmonics:
            ax.set_xlabel(f"$k$")
            ax.set_ylabel(f"$E_k$")
        else:
            ax.set_xlabel(f"$n$")
            ax.xaxis.get_major_locator().set_params(integer=True)
            ax.set_ylabel(f"$E_n$")
        ax.set_xscale("log")
        ax.set_yscale("log")
        fig.savefig(figuresDir / f"{statePath.name[-6:]}_spec_iso.png")        
        
        if savetxt:
            np.savetxt(figuresDir / f"{statePath.name[-6:]}_spec_iso.dat", spec_iso)
            np.savetxt(figuresDir / f"{statePath.name[-6:]}_dissipation.dat", np.array([dissipation]))

    if not noshow:
        plt.show()


def dnsspec(statePath, iso=False):

    statePath = Path(statePath)
    state, header = dns.readState(statePath)
    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = header
    nxp, nyp, nzp = nx // 2 - 1, ny // 2 - 1, nz // 2 - 1

    spec_x = np.zeros((nxp + 1))
    spec_y = np.zeros((nyp + 1))
    spec_z = np.zeros((nzp + 1))

    if iso:
        spec_iso = {}
        _kx = 2 * np.pi / Lx
        _ky = 2 * np.pi / dns.Ly
        _kz = 2 * np.pi / Lz
        _kmin = min(_kx, _ky)
        _kmin = min(_kmin, _kz)

    for i in range(0, nx):
        kx = i
        if kx > nx // 2:
            kx = kx - nx
        kx = np.abs(kx)
        # x keeps one zeroed mode (the largest)
        if kx == nxp + 1:
            continue
        for j in range(0, nyp + 1):
            ky = j
            for k in range(0, nz):
                kz = k
                if kz > nz // 2:
                    kz = kz - nz
                kz = np.abs(kz)
                if kz == nzp + 1:
                    continue

                part = (np.sum(np.conj(state[i, j, k, :]) * state[i, j, k, :])).real

                # Correct for double counting (Hermiticity)
                if ky == 0:
                    part = part / 2

                spec_x[kx] += part
                spec_y[ky] += part
                spec_z[kz] += part

                # isotropic spectrum
                if iso and not (kx == 0 and ky == 0 and kz == 0):
                    k_iso = int(
                        round(
                            np.sqrt((_kx * kx) ** 2 + (_ky * ky) ** 2 + (_kz * kz) ** 2)
                            / _kmin
                        )
                    )

                    if k_iso in spec_iso:
                        spec_iso[k_iso] += part
                    else:
                        spec_iso[k_iso] = part

    if iso:
        spec_iso_ = np.zeros((len(spec_iso), 2), dtype=np.float64)
        i = 0
        for key, value in spec_iso.items():
            spec_iso_[i, 0] = key * _kmin
            spec_iso_[i, 1] = value
            i += 1

        sorter = np.argsort(spec_iso_[:, 0])
        spec_iso_ = spec_iso_[sorter]

        vorticity = dns.vorticity(state, Lx, Lz)
        dissipation = 2 * dns.inprod(vorticity, vorticity) / Re

        return (
            spec_x,
            spec_y,
            spec_z,
            spec_iso_,
            Lx,
            Lz,
            dissipation,
        )
    else:
        return (
            spec_x,
            spec_y,
            spec_z,
            Lx,
            Lz,
        )


if __name__ == "__main__":
    main()
