#!/usr/bin/env python3
import argparse
from collections import deque
from math import floor
from operator import itemgetter
from pathlib import Path

import numpy as np
import psutil
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage import minimum_filter
from tqdm.auto import tqdm

import dns

cmap = "Greys"

def_n_jobs = len(psutil.Process().cpu_affinity())
if def_n_jobs > 1:
    def_n_jobs = def_n_jobs - 1

def_joblib_verbosity = 0
def_joblib_backend = "threading"


def main():

    parser = argparse.ArgumentParser(
        description="Identifies and records local minima in a recurrence data."
    )

    parser.add_argument("rundir", type=str, help="path to the directory of the run.")
    parser.add_argument(
        "savedir", type=str, help="path to the directory to save results."
    )
    parser.add_argument("t_rec", type=float, help="recurrence time.")
    parser.add_argument(
        "dt", type=float, help="dt of states to search recurrences for."
    )
    parser.add_argument(
        "-n_jobs",
        default=def_n_jobs,
        type=int,
        dest="n_jobs",
        help="number of threads to use",
    )
    parser.add_argument(
        "-t_i",
        type=float,
        default=-np.inf,
        dest="t_i",
        help="initial time.",
    )
    parser.add_argument(
        "-t_f",
        type=float,
        default=np.inf,
        dest="t_f",
        help="final time.",
    )
    parser.add_argument(
        "--reprocess",
        action="store_true",
        dest="reprocess",
        help="reprocess existing recurrence data.",
    )
    parser.add_argument(
        "-f_subsamp",
        type=int,
        default=1,
        dest="f_subsamp",
        help="subsampling frequency of the found state files.",
    )
    parser.add_argument(
        "-d_max",
        type=float,
        default=0.5,
        dest="d_max",
        help="maximum recurrence distance.",
    )
    parser.add_argument(
        "-t_radius",
        type=float,
        default=20.0,
        dest="t_radius",
        help="radius of the minimum search area.",
    )

    args = vars(parser.parse_args())
    recurrence(**args)


def recurrence(
    rundir,
    savedir,
    t_rec,
    dt,
    n_jobs=def_n_jobs,
    t_i=-np.inf,
    t_f=np.inf,
    reprocess=False,
    f_subsamp=1,
    d_max=0.15,
    t_radius=20.0,
):
    rundir = Path(rundir).resolve()
    savedir = Path(savedir).resolve()

    statefiles = [
        str(stateFile.resolve()) for stateFile in rundir.glob("sliced_state.*")
    ]
    times = np.array(
        [float(stateFile.name[-6:]) * dt for stateFile in rundir.glob("sliced_state.*")]
    )
    # Sort in time
    sorter = np.argsort(times)
    times = times[sorter]
    statefiles = [statefiles[i] for i in sorter]
    # Filter in time
    t_filter = np.nonzero(np.logical_and(times >= t_i, times < t_f + t_rec))[0]
    times = times[t_filter]
    statefiles = [statefiles[i] for i in t_filter]
    # Subsample
    times = times[::f_subsamp]
    statefiles = statefiles[::f_subsamp]

    dt = dt * f_subsamp
    n_states = len(statefiles)
    n_rec = int(floor(t_rec / dt))
    n_rows = n_rec + 1
    n_cols = n_states - n_rec
    n_radius = int(floor(t_radius / dt))

    print(
        f"Searching from time {times[0]} to time {times[n_cols - 1]} with dt {dt} and T {t_rec}.",
        flush=True,
    )

    if not reprocess:

        phases_data = np.loadtxt(rundir / "phases.gp")
        tphases = phases_data[:, 0]
        xphases = phases_data[:, 1]
        zphases = phases_data[:, 2]

        shifts = np.zeros((n_states, 2))
        states = deque(maxlen=n_rows)
        recs = np.zeros((n_rows, n_cols), dtype=np.float64)

        def fillrecs(j, i):
            deltastate = states[j] - states[0]
            normdelta = np.sqrt(dns.inprod(deltastate, deltastate))
            recs[j, i] = normdelta

        with Parallel(n_jobs=n_jobs, backend=def_joblib_backend) as parallel:

            j_state_new = 0

            for i in tqdm(range(n_cols)):

                if i == 0:
                    for j in range(n_rec + 1):
                        state_new, header = dns.readState(statefiles[j_state_new])
                        time_new = header[-1]
                        it = np.argmin(tphases - time_new)
                        phase_new_x = xphases[it]
                        phase_new_z = zphases[it]

                        shifts[j_state_new, 0] = find_shift_from_phase(phase_new_x)
                        shifts[j_state_new, 1] = find_shift_from_phase(phase_new_z)
                        states.append(state_new)
                        j_state_new += 1
                else:
                    state_new, header = dns.readState(statefiles[j_state_new])
                    time_new = header[-1]
                    it = np.argmin(tphases - time_new)
                    phase_new_x = xphases[it]
                    phase_new_z = zphases[it]

                    shifts[j_state_new, 0] = find_shift_from_phase(phase_new_x)
                    shifts[j_state_new, 1] = find_shift_from_phase(phase_new_z)
                    states.append(state_new)
                    j_state_new += 1

                parallel(delayed(fillrecs)(j, i) for j in range(1, n_rec + 1))

                norm0 = np.sqrt(dns.inprod(states[0], states[0]))
                recs[:, i] = recs[:, i] / norm0

        np.savetxt(savedir / "shifts.gp", shifts)
        np.save(savedir / "recs.npy", recs)

    else:
        shifts = np.loadtxt(savedir / "shifts.gp")
        recs = np.load(savedir / "recs.npy")

    # Apply filters
    mins = minimum_filter(recs, size=n_radius, mode="reflect")  # Run the search
    # Filter the results based on the recurrence distance limits
    # Where they don't fit the filter, return 0
    recs_filtered = np.where(
        np.logical_and(mins == recs, mins < d_max),
        recs,
        0,
    )
    # Get the indices of local minima that are not zero
    # Meaning, we thrash those with zero recurrence distance and those that don't
    # fit the limits
    indices = np.transpose(np.nonzero(recs_filtered))

    signals = list()
    # Loop over minima
    # Transpose is due to how np.nonzero works
    for index in indices:
        # Note that for the recurrence matrix, rows are time delays, columns
        # are times of reference states
        row, column = index
        period = times[column + row] - times[column]
        if period < t_rec:
            d = recs[row, column]
            time = times[column]
            statefile = statefiles[column]
            shiftx = shift_center((shifts[column + row, 0] - shifts[column, 0]) % 1)
            shiftz = shift_center((shifts[column + row, 1] - shifts[column, 1]) % 1)
            shiftd = np.sqrt(shiftx ** 2 + shiftz ** 2) / np.sqrt(2 * 0.5 ** 2)

            signal = [
                d,
                shiftd,
                period,
                shiftx,
                shiftz,
                time,
                statefile,
            ]

            signals.append(signal)

    # Sort the suggestions by relative error
    signals_sorted = sorted(signals, key=itemgetter(0))

    # Save the signals, if any
    if len(signals_sorted) > 0:
        signalsFile = savedir / "signals"
        with open(signalsFile, "w") as f:
            for signal in signals_sorted:
                f.write(" ".join([str(i) for i in signal]) + "\n")

    else:
        print("No recurrences found within set limits.")

    # Plot recurrence data
    fig, ax = plt.subplots()
    cbar = ax.imshow(
        recs,
        cmap=cmap,
        origin="lower",
        vmin=0,
        vmax=np.amax(recs),
        interpolation="spline16",
        extent=[times[0], times[n_cols - 1], 0, t_rec],
    )
    ax.set_xlabel("$t$")
    ax.set_ylabel("$T$")
    ax.set_xlim(left=times[0], right=times[n_cols - 1])
    ax.set_ylim(bottom=0, top=t_rec)
    colorbar(ax, cbar, label=f"$r(t,t+T)$")
    fig.savefig(savedir / f"recs.png", bbox_inches="tight")

    # Plot candidates
    if len(signals_sorted) > 0:
        fig, ax = plt.subplots()
        cbar = ax.imshow(
            recs,
            cmap=cmap,
            origin="lower",
            vmin=0,
            vmax=np.amax(recs),
            interpolation="spline16",
            extent=[times[0], times[n_cols - 1], 0, t_rec],
        )
        for index in indices:
            # Note that for the recurrence matrix, rows are time delays, columns
            # are times of reference states
            row, column = index
            period = times[column + row] - times[column]
            if period < t_rec:
                ax.scatter(times[0] + column * dt, row * dt, c="k", s=75)
        ax.set_xlabel("$t$")
        ax.set_ylabel("$T$")
        ax.set_xlim(left=times[0], right=times[n_cols - 1])
        ax.set_ylim(bottom=0, top=t_rec)
        colorbar(ax, cbar, label=f"$r(t,t+T)$")
        fig.savefig(savedir / f"recs_candidates.png", bbox_inches="tight")


def shift_center(shift):
    if np.abs(shift - 1) < shift:
        shift_ = shift - 1
    else:
        shift_ = shift

    return shift_


def find_shift_from_phase(phase):

    shift = shift_center(((phase / (2 * np.pi)) % 1 / (2 * np.pi)) % 1)

    return shift


def colorbar(ax, im, label=None):

    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, label=label)

    return cbar


if __name__ == "__main__":
    main()
