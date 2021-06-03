#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

import dnsCommon as dns
from mayavi import mlab

mlab.options.offscreen = True


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "runDir",
        type=str,
        help="path to the directory of visualization files of interest.",
    )
    parser.add_argument(
        "state", type=int, help="state number of interest. 100 for velocity.000100 etc."
    )
    parser.add_argument(
        "--SxRzUnit",
        action="store_true",
        dest="SxRzUnit",
        help="visualize only a unit region of the SxRz symmetry.",
    )
    parser.add_argument(
        "--showAxes",
        action="store_true",
        dest="showAxes",
        help="show orientation axes.",
    )
    parser.add_argument(
        "--noSubtractLaminar",
        action="store_true",
        dest="noSubtractLaminar",
        help="do not subtract the laminar state.",
    )
    parser.add_argument(
        "--velCC",
        type=float,
        dest="velCC",
        help="coefficient of constant velocity contours.",
        default=0.95,
    )
    parser.add_argument(
        "--vorCC",
        type=float,
        dest="vorCC",
        help="coefficient of constant vorticity contours.",
        default=0.65,
    )
    parser.add_argument(
        "--dimX", type=int, dest="dimX", help="pixels in x.", default=800
    )
    parser.add_argument(
        "--dimY", type=int, dest="dimY", help="pixels in y.", default=800
    )

    args = vars(parser.parse_args())
    dimX = args.pop("dimX")
    dimY = args.pop("dimY")
    dims = (dimX, dimY)

    isosurf(**args, dims=dims)


def isosurf(
    runDir,
    state,
    SxRzUnit=False,
    showAxes=False,
    noSubtractLaminar=False,
    velCC=0.95,
    vorCC=0.65,
    dims=(800, 800),
):

    runDir = Path(runDir).resolve()
    parametersPath = runDir / "parameters.in"
    vorticityPath = runDir / ("vorticity." + str(state).zfill(6))
    velocityPath = runDir / ("velocity." + str(state).zfill(6))
    figurePath = runDir / (str(state).zfill(6) + ".png")

    pars = dns.readParameters(parametersPath)
    gamma = pars["gamma"]
    kF = pars["kF"]
    nu = pars["nu"]
    Nz = pars["nz"]
    Ny = pars["ny"]
    Nx = pars["nx_all"]
    Lz = pars["Lz"]

    if SxRzUnit:
        Ny = Ny // 2
        Nz = Nz // 2

    vorx = dns.readVisualizeState(vorticityPath)
    velx = dns.readVisualizeState(velocityPath)

    if SxRzUnit:
        vorx = np.copy(vorx[0:Nz, Ny // 2 : Ny + Ny // 2, :])
        velx = np.copy(velx[0:Nz, Ny // 2 : Ny + Ny // 2, :])

    # Get things back to x, y, z order in the physical space
    vorx = np.swapaxes(vorx, 0, 2)
    velx = np.swapaxes(velx, 0, 2)

    Lx = 2.0 * np.pi
    Ly = 2.0 * np.pi
    Lz = Lz * np.pi
    if SxRzUnit:
        Ly = Ly / 2
        Lz = Lz / 2

    x = [i * (Lx / Nx) for i in range(0, Nx)]
    y = [i * (Ly / Ny) for i in range(0, Ny)]
    z = [i * (Lz / Nz) for i in range(0, Nz)]

    if SxRzUnit:
        y = [i + np.pi / 2 for i in y]

    xx = np.zeros((Nx, Ny, Nz), dtype="float")
    yy = np.zeros((Nx, Ny, Nz), dtype="float")
    zz = np.zeros((Nx, Ny, Nz), dtype="float")

    if not noSubtractLaminar:
        lamvx = np.zeros(velx.shape)

        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    xx[i, j, k] = x[i]
                    yy[i, j, k] = y[j]
                    zz[i, j, k] = z[k]
                    lamvx[i, j, k] = (gamma / (nu * kF ** 2)) * np.sin(
                        y[j]
                    )  # Laminar sol.

    cmap1 = "PiYG"
    cmap2 = "blue-red"

    surf1 = vorx

    if not noSubtractLaminar:
        surf2 = velx - lamvx  # Subtracting laminar
    else:
        surf2 = velx

    mlab.figure(bgcolor=(1.0, 1.0, 1.0), fgcolor=(0.1, 0.1, 0.1), size=dims)

    src = mlab.pipeline.scalar_field(xx, yy, zz, surf1, colormap=cmap1, name="sfield")

    mlab.pipeline.iso_surface(
        src,
        contours=[surf1.max() * vorCC, surf1.min() * vorCC],
        colormap=cmap1,
        opacity=0.85,
        vmin=surf1.min() * vorCC,
        vmax=surf1.max() * vorCC,
    )

    src = mlab.pipeline.scalar_field(xx, yy, zz, surf2, colormap=cmap2)

    mlab.pipeline.iso_surface(
        src,
        contours=[surf2.max() * velCC, surf2.min() * velCC],
        colormap=cmap2,
        opacity=0.35,
        vmin=surf2.min() * velCC,
        vmax=surf2.max() * velCC,
    )

    mlab.axes(x_axis_visibility=False, y_axis_visibility=False, z_axis_visibility=False)

    scene = mlab.gcf()
    # Disable rendering until everything is ready / g 190430
    scene.scene.disable_render = True

    scene.scene.camera.position = [
        -15.317366231234612,
        10.12151296484298,
        9.414788208003529,
    ]
    scene.scene.camera.focal_point = [
        3.0925052165985107,
        3.0925052165985107,
        3.0925052165985107,
    ]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [
        0.3399799924116968,
        0.9388894854294522,
        -0.05385293780065341,
    ]
    scene.scene.camera.clipping_range = [11.043861252010439, 32.89472646477479]
    scene.scene.camera.compute_view_plane_normal()

    if showAxes:
        mlab.orientation_axes(figure=scene, line_width=5)

    mlab.outline()

    scene.scene.render()

    mlab.savefig(filename=str(figurePath.resolve()))
    # mlab.screenshot()
    # mlab.show()
    mlab.clf()
    mlab.close(all=True)

    # Return the path to the figure
    return figurePath


if __name__ == "__main__":
    main()
