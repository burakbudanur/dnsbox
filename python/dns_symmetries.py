#!/usr/bin/env python3
import numpy as np

from dns import (
    isEven,
    wavenumbers,
)


def TxHalfRepr(nx, ny_half, nz):
    # Representation of T_x(L_x/2) on the states

    # Representative array of the action of T_x(L_x/2) on the Fourier coefficients
    TxHalfRepr = np.ones((nx, ny_half, nz, 3), dtype=np.int)

    # Again, dnsbox flips x and y in the Fourier space
    for i in range(0, nx):
        if not isEven(i):
            TxHalfRepr[i, :, :, :] = -TxHalfRepr[i, :, :, :]

    return TxHalfRepr


def Tx(dx, u, Lx, Lz):
    # Return the T_x(dx) image of a state

    nx, ny_half, nz, _ = u.shape
    kkx, _, _ = wavenumbers(Lx, Lz, nx, ny_half, nz)

    image = np.zeros((nx, ny_half, nz, 3), dtype=np.complex128)
    for ix, kx in enumerate(kkx):
        # Re part
        image[ix, :, :, :] = (
            np.cos(kx * dx) * u[ix, :, :, :].real
            + np.sin(kx * dx) * u[ix, :, :, :].imag
        )
        # Im part
        image[ix, :, :, :] += 1j * (
            -np.sin(kx * dx) * u[ix, :, :, :].real
            + np.cos(kx * dx) * u[ix, :, :, :].imag
        )

    return image


def TyHalfRepr(nx, ny_half, nz):
    # Representation of T_y(L_y/2) on the states

    TyHalfRepr = np.ones((nx, ny_half, nz, 3), dtype=np.int)

    # Again, dnsbox flips z and y in the Fourier space
    for j in range(0, ny_half):
        if not isEven(j):
            TyHalfRepr[:, j, :, :] = -TyHalfRepr[:, j, :, :]

    return TyHalfRepr


def Ty(dy, u, Lx, Lz):
    """
    Return the T_y(dy) image of a state
    Warning: Not a symmetry of the Kolmogorov flow for all dy
    """

    nx, ny_half, nz, _ = u.shape
    _, kky, _ = wavenumbers(Lx, Lz, nx, ny_half, nz)

    image = np.zeros((nx, ny_half, nz, 3), dtype=np.complex128)
    for iy, ky in enumerate(kky):

        # Re part
        image[:, iy, :, :] = (
            np.cos(ky * dy) * u[:, iy, :, :].real
            + np.sin(ky * dy) * u[:, iy, :, :].imag
        )
        # Im part
        image[:, iy, :, :] += 1j * (
            -np.sin(ky * dy) * u[:, iy, :, :].real
            + np.cos(ky * dy) * u[:, iy, :, :].imag
        )

    return image


def TzHalfRepr(nx, ny_half, nz):
    # Representation of T_z(L_z/2) on the states

    TzHalfRepr = np.ones((nx, ny_half, nz, 3), dtype=np.int)
    for k in range(nz):
        if not isEven(k):
            TzHalfRepr[:, :, k, :] = -TzHalfRepr[:, :, k, :]

    return TzHalfRepr


def Tz(dz, u, Lx, Lz):
    # Return the T_z(dz) image of a state

    nx, ny_half, nz, _ = u.shape
    _, _, kkz = wavenumbers(Lx, Lz, nx, ny_half, nz)

    image = np.zeros((nx, ny_half, nz, 3), dtype=np.complex128)
    for iz, kz in enumerate(kkz):

        # Re part
        image[:, :, iz, :] = (
            np.cos(kz * dz) * u[:, :, iz, :].real
            + np.sin(kz * dz) * u[:, :, iz, :].imag
        )
        # Im part
        image[:, :, iz, :] += 1j * (
            -np.sin(kz * dz) * u[:, :, iz, :].real
            + np.cos(kz * dz) * u[:, :, iz, :].imag
        )

    return image


def Rx(u):
    # Return the R_x image of a state

    nx, ny_half, nz, _ = u.shape

    image = np.zeros((nx, ny_half, nz, 3), dtype=u.dtype)

    # Rx maps [u,v,w](kx, ky, kz) to [-u,v,w](-kx, ky, kz)

    for i in range(nx):

        # u component switches sign
        image[i, :, :, 0] = -u[-i, :, :, 0]

        # The rest
        image[i, :, :, 1:3] = u[-i, :, :, 1:3]

    return image


def Ry(u):
    # Return the R_y image of a state
    # Not a symmetry by itself

    nx, ny_half, nz, _ = u.shape

    image = np.zeros((nx, ny_half, nz, 3), dtype=u.dtype)

    # Ry maps [u,v,w](kx, ky, kz) to [u,-v,w](-kx, ky, -kz)*

    for i in range(nx):
        for j in range(ny_half):
            for k in range(nz):

                if j != 0:

                    # v component switches sign
                    # Real part
                    image[i, j, k, 1] = -u[-i, j, -k, 1].real
                    # Imaginary part, complex conjugate
                    image[i, j, k, 1] += 1j * u[-i, j, -k, 1].imag

                    # The rest
                    image[i, j, k, [0, 2]] = u[-i, j, -k, [0, 2]].real
                    # Imaginary part, complex conjugate
                    image[i, j, k, [0, 2]] += -1j * u[-i, j, -k, [0, 2]].imag

                else:

                    # v component switches sign
                    image[i, j, k, 1] = -u[i, j, k, 1].real
                    image[i, j, k, 1] += -1j * u[i, j, k, 1].imag

                    # The rest remain the same
                    image[i, j, k, [0, 2]] = u[i, j, k, [0, 2]].real
                    image[i, j, k, [0, 2]] += 1j * u[i, j, k, [0, 2]].imag

    return image


def Rz(u):
    # Return the R_z image of a state

    nx, ny_half, nz, _ = u.shape

    image = np.zeros((nx, ny_half, nz, 3), dtype=u.dtype)

    # Rz maps [u,v,w](kx, ky, kz) to [u,v,-w](kx, ky, -kz)

    for k in range(nz):

        # w component switches sign
        image[:, :, k, 2] = -u[:, :, -k, 2]

        # The rest
        image[:, :, k, 0:2] = u[:, :, -k, 0:2]

    return image


def Sx(u):
    # Return the S_x image of a state

    nx, ny_half, nz, _ = u.shape

    image = TyHalfRepr(nx, ny_half, nz) * u
    image = Rx(image)

    return image


def RyMid(u):
    # Return the R_y T_y(L_y/2) image of a state

    nx, ny_half, nz, _ = u.shape

    # Reflection about y = L_y/4
    # Equivalent to Ry Ty(Ly/2)

    image = TyHalfRepr(nx, ny_half, nz) * u
    image = Ry(image)

    return image


def Rxy(u):
    # Return the R_xy image of a state.

    nx, ny_half, nz, _ = u.shape

    image = Ry(u)
    image = Rx(image)

    return image
