#!/usr/bin/env python3
import argparse
from math import floor
from operator import itemgetter
from pathlib import Path
from sys import exit

import numpy as np
from scipy.ndimage import minimum_filter

import dnsCommon as dns


def main():

    parser = argparse.ArgumentParser(
        description="Identifies and records local minima in a recurrence data."
    )

    parser.add_argument("runDir", type=str, help="path to the directory of the run.")

    parser.add_argument(
        "--dMin",
        type=float,
        default=0.0,
        dest="dMin",
        help="minimum recurrence distance.",
    )
    parser.add_argument(
        "--dMax",
        type=float,
        default=0.3,
        dest="dMax",
        help="maximum recurrence distance.",
    )
    parser.add_argument(
        "--windowTime",
        type=float,
        default=4.0,
        dest="windowTime",
        help="radius of the minimum search area.",
    )

    args = vars(parser.parse_args())
    recMinima(**args)


def recMinima(runDir, dMin=0.0, dMax=0.3, windowTime=4.0):
    """Identifies and records local minima in recurrence data using
    minimum_filter of SciPy. Boundaries are mirrored, and search regions are
    set based on the maximum recurrence time and windowSize.
    """

    # Convert to Path
    runDir = Path(runDir).resolve()

    # Read parameters
    params = dns.readParameters(runDir / "parameters.in")

    # Read recurrence parameters
    with open(runDir / "recurrence.in") as recParams:
        recInitialTime = float(recParams.readline().strip())
        recFinalTime = float(recParams.readline().strip())

    laminarizedFile = runDir / "LAMINARIZED"
    if not Path.is_file(laminarizedFile):
        exit("Laminarized file is not found.")
    # Read the laminarization time
    with open(laminarizedFile) as f:
        LAMINARIZED = f.readlines()
    # Read the first one, that's what recurrence reads
    time_finalState = float(LAMINARIZED[0].strip())

    # Search from state 0 to the final state
    deltaState = floor((recInitialTime / params["dt"]) / params["IPRINT2"])
    nrecur = floor((recFinalTime / params["dt"]) / params["IPRINT2"])
    # Radius of the minimum searches
    windowSize = floor((windowTime / params["dt"]) / params["IPRINT2"])
    # Set the last state such that there are at least nrecur many states
    # more after it
    initialState = 0
    finalState = floor((time_finalState / params["dt"]) / params["IPRINT2"]) - nrecur

    # This is to convert state file numbers to time
    # Full recurrence is done with respect to state file saving frequency
    timeScale = params["IPRINT2"] * params["dt"]

    recFile = runDir / "recurrence.dat"
    if not Path.is_file(recFile):
        exit("Recurrence data is not found.")

    nColumns = finalState - initialState + 1
    nRows = nrecur - deltaState + 1

    # Let's try and work with stream i/o
    with open(recFile, "rb") as f:
        toReshape = np.fromfile(f, np.float32)

    # Sanity check
    if len(toReshape) != nRows * nColumns:
        print("Length of input data:", len(toReshape))
        print("Length of expected data:", nRows * nColumns)
        exit("Recurrence data does not have the expected size.")

    # Reshape
    # The file was written column-wise
    recs = np.transpose(np.reshape(toReshape, (nColumns, nRows)))

    # Apply filters
    mins = minimum_filter(recs, size=windowSize, mode="reflect")  # Run the search
    # Filter the results based on the recurrence distance limits
    # Where they don't fit the filter, return 0
    recsFiltered = np.where(
        np.logical_and(mins == recs, mins <= dMax, dMin <= mins), recs, 0
    )
    # Get the indices of local minima that are not zero
    # Meaning, we thrash those with zero recurrence distance and those that don't
    # fit the limits
    indices = np.transpose(np.nonzero(recsFiltered))

    signals = list()
    # Loop over minima
    # Transpose is due to how np.nonzero works
    for index in indices:
        # Note that for the recurrence matrix, rows are time delays, columns
        # are times of reference states
        row, column = index
        d = recsFiltered[row, column]
        # Calling this initial time is wrong, but it works for now
        period = (row + deltaState) * timeScale
        time = (initialState + column) * timeScale
        state = floor((time / params["dt"]) / params["IPRINT2"])

        signal = [period, d, time, state]

        signals.append(signal)

    # Sort the suggestions by period
    # I initially had this thinking lower period orbits were more important,
    # but even if that's not generally true, it's still a good idea to get them
    # done first, since their searches should, I suppose, end earlier.
    signalsP = sorted(signals, key=itemgetter(2))

    # Save the signals, if any
    if len(signalsP) > 0:
        signalsFile = runDir / "SIGNALS"
        with open(signalsFile, "w") as f:
            for signal in signalsP:
                f.write(" ".join([str(i) for i in signal]) + "\n")

    else:
        exit("No recurrences found within set limits.")


if __name__ == "__main__":
    main()
