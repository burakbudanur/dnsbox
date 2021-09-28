#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from subprocess import call
import numpy as np

import dns


def main():

    parser = argparse.ArgumentParser(
        description="Clean up a run directory to restart a run from i_start.",
    )
    parser.add_argument("rundir", type=str, help="path to the run")
    parser.add_argument(
        "-script",
        dest="script",
        help="Submission script. If given, submit the job by sbatch script",
    )

    args = vars(parser.parse_args())

    dnscontinue(**args)


def dnscontinue(rundir, script):

    rundir = Path(rundir)
    parameters = dns.readParameters(rundir / "parameters.in")
    i_start =  parameters["initiation"]["i_start"]

    states = sorted(list(rundir.glob("state.*")))
    perturbed_states = sorted(list(rundir.glob("perturb.*")))

    if len(states) > 0:
        for state in states:
            if int(state.name[-6:]) > i_start:
                state.unlink()

    if len(perturbed_states) > 0:
        for state in perturbed_states:
            if int(state.name[-6:]) > i_start:
                state.unlink()

    results = list(rundir.glob("*.gp"))

    for file in results:
        with open(file, "r") as f:
            lines = f.readlines()

        with open(file, "w") as f:
            for line in lines:
                try:
                    i_time = int(re.search(r"\d+", line).group())

                    if i_time > i_start:
                        break
                    else:
                        f.write(line)

                except:
                    f.write(line)

    if not script == None:
        call(["sbatch", script])


if __name__ == "__main__":
    main()
