#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import psutil
import pyvista as pv
from joblib import Parallel, delayed
from tqdm import tqdm

import dnsbox as dns


def_n_jobs = len(psutil.Process().cpu_affinity())
if def_n_jobs > 1:
    def_n_jobs = def_n_jobs - 1

def_joblib_verbosity = 0
def_joblib_backend = "loky"

# Example ffmpeg command to merge png files to an mp4:
# ffmpeg -framerate 24 -i %06d.png -c:v libx264 -r 60 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" out.mp4

# This script removes the laminar part


def main():

    parser = argparse.ArgumentParser(
        description="Produce 3D visualizations of a state",
    )
    parser.add_argument("statesDir", type=str, help="path to the state.")
    parser.add_argument(
        "--xvfb", action="store_true", dest="xvfb", help="render to a virtual display."
    )
    parser.add_argument(
        "-cvel",
        default=0.75,
        type=float,
        dest="cvel",
        help="multiplier for velocity isosurfaces",
    )
    parser.add_argument(
        "-cvor",
        default=0.5,
        type=float,
        dest="cvor",
        help="multiplier for vorticity isosurfaces",
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

    print(args)

    dnsbatchvis(**args)


def dnsbatchvis(
    statesDir,
    xvfb=False,
    cvel=0.75,
    cvor=0.5,
    mirror_y=False,
    mirror_z=False,
    print_messages=True,
    n_jobs=def_n_jobs,
):

    statesDir = Path(statesDir)
    states = sorted(list(statesDir.glob("state.*")))
    min_vel, min_vor = np.inf, np.inf
    max_vel, max_vor = -np.inf, -np.inf

    print("Finding the maxima in the snapshots.")
    with tqdm(total=len(states), disable=not print_messages) as pbar:
        for i, state in enumerate(states):
            stateIn, headers = dns.readState(state)

            forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers
            ny_half = ny // 2

            stateIn = stateIn - dns.laminar(
                forcing, nx, ny_half, nz, tilt_angle=tilt_angle
            )

            vel_phys = dns.fftSpecToPhysAll(stateIn)
            vor_spec = dns.vorticity(stateIn, Lx, Lz)
            vor_phys = dns.fftSpecToPhysAll(vor_spec)

            if not mirror_y:
                ny_display = ny
            else:
                ny_display = ny // 2 + 1
            if not mirror_z:
                nz_display = nz
            else:
                nz_display = nz // 2 + 1

            velx = vel_phys[:, :ny_display, :nz_display, 0]
            vorx = vor_phys[:, :ny_display, :nz_display, 0]

            min_vel, min_vor = min(min_vel, np.amin(velx)), min(min_vor, np.amin(vorx))
            max_vel, max_vor = max(max_vel, np.amax(velx)), max(max_vor, np.amax(vorx))
            pbar.update()

    print("Minima:")
    print(min_vel, min_vor)
    print("Maxima:")
    print(max_vel, max_vor)

    figuresDir = dns.createFiguresDir(statesDir.parent)

    if xvfb:
        pv.start_xvfb()

    def render_state_i(i):
        state = Path(states[i])
        stateIn, headers = dns.readState(state)

        forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers
        ny_half = ny // 2

        stateIn = stateIn - dns.laminar(forcing, nx, ny_half, nz, tilt_angle=tilt_angle)

        vel_phys = dns.fftSpecToPhysAll(stateIn)
        vor_spec = dns.vorticity(stateIn, Lx, Lz)
        vor_phys = dns.fftSpecToPhysAll(vor_spec)

        if not mirror_y:
            ny_display = ny
        else:
            ny_display = ny // 2 + 1
        if not mirror_z:
            nz_display = nz
        else:
            nz_display = nz // 2 + 1

        velx = vel_phys[:, :ny_display, :nz_display, 0]
        vorx = vor_phys[:, :ny_display, :nz_display, 0]

        pv.set_plot_theme("document")
        u = pv.wrap(velx)
        om = pv.wrap(vorx)

        vel_levels = cvel * np.array([min_vel, max_vel])
        vor_levels = cvor * np.array([min_vor, max_vor])

        p = pv.Plotter(off_screen=True)

        p.add_mesh(u.outline(), color="k")
        p.add_mesh(
            u.contour(vel_levels),
            smooth_shading=True,
            opacity=0.35,
            cmap=["red", "blue"],
            clim=vel_levels,
            show_scalar_bar=False,
        )
        p.add_mesh(
            om.contour(vor_levels),
            smooth_shading=True,
            opacity=0.35,
            cmap=["green", "purple"],
            clim=vor_levels,
            show_scalar_bar=False,
        )
        p.show_axes()

        #
        p.camera.roll += 90
        p.camera.elevation -= 15
        p.camera.azimuth -= 45
        p.camera.roll += 30
        p.camera.azimuth -= 45
        p.camera.roll -= 10

        p.show(screenshot=figuresDir / f"{state.name}_isosurf.png")

        # hope memory doesn't leak
        p.clear()
        p.deep_clean()

    Parallel(n_jobs=n_jobs, backend=def_joblib_backend, verbose=def_joblib_verbosity)(
        delayed(render_state_i)(i) for i in tqdm(range(len(states)))
    )


if __name__ == "__main__":
    main()
