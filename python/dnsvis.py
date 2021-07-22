#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import dns
import pyvista as pv


def main():

    parser = argparse.ArgumentParser(
        description="Produce 3D visualizations of a state",
    )
    parser.add_argument("state", type=str, help="path to the state.")
    parser.add_argument(
        "--noshow", action="store_true", dest="noshow", help="do not display the plots."
    )
    parser.add_argument(
        "--xvfb", action="store_true", dest="xvfb", help="render to a virtual display."
    )
    parser.add_argument(
        "--undotilt",
        action="store_true",
        dest="undotilt",
        help="rotate u and w to be parallel/orthogonal to the forcing.",
    )
    parser.add_argument(
        "--sublam",
        action="store_true",
        dest="sublam",
        help="remove laminar part from states.",
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

    dnsvis(**args)


def dnsvis(
    state,
    noshow=False,
    xvfb=False,
    undotilt=False,
    sublam=False,
    cvel=0.75,
    cvor=0.5,
    mirror_y=False,
    mirror_z=False,
):

    if xvfb:
        noshow = True
        pv.start_xvfb()

    pv.set_plot_theme("document")
    state = Path(state)
    figuresDir = dns.createFiguresDir(state.parent)
    stateIn, headers = dns.readState(state)

    forcing, nx, ny, nz, Lx, Lz, Re, tilt_angle, dt, itime, time = headers
    ny_half = ny // 2

    if sublam:
        stateIn = stateIn - dns.laminar(forcing, nx, ny_half, nz, tilt_angle=tilt_angle)

    uw_untilted = False
    if abs(tilt_angle) > 0 and undotilt:
        uw_untilted = True
        stateIn = dns.tilt_state(stateIn, tilt_angle)

    vel_phys = dns.fftSpecToPhysAll(stateIn)
    vor_spec = dns.vorticity(stateIn, Lx, Lz)
    vor_phys = dns.fftSpecToPhysAll(vor_spec)

    x_label = "$x$"
    y_label = "$y$"
    z_label = "$z$"

    if not mirror_y:
        yMid = "2"
        ny_display = ny
    else:
        yMid = "1"
        ny_display = ny // 2 + 1
    if not mirror_z:
        zMid = "L_z/2"
        nz_display = nz
    else:
        zMid = "L_z/4"
        nz_display = nz // 2 + 1

    velx = vel_phys[:, :ny_display, :nz_display, 0]
    vorx = vor_phys[:, :ny_display, :nz_display, 0]

    u = pv.wrap(velx)
    om = pv.wrap(vorx)

    vel_levels = cvel * np.array([np.min(velx), np.max(velx)])
    vor_levels = cvor * np.array([np.min(vorx), np.max(vorx)])

    p = pv.Plotter(off_screen=noshow)

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


if __name__ == "__main__":
    main()
