#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from subprocess import call
import numpy as np

import dnsbox as dns


def main():

    parser = argparse.ArgumentParser(
        description="Continue a run from the last state file saved",
    )
    parser.add_argument("rundir", type=str, help="path to the run")
    parser.add_argument(
        "-i_finish_plus",
        type=int,
        dest="i_finish_plus",
        help="number of time steps to add to i_finish",
    )
    parser.add_argument(
        "-script",
        dest="script",
        help="Submission script. If given, submit the job by sbatch script",
    )

    args = vars(parser.parse_args())

    dnscontinue(**args)


def dnscontinue(rundir, i_finish_plus=None, script=None):

    rundir = Path(rundir)
    states = sorted(list(rundir.glob("state.*")))
    i_final_state = int(states[-1].name[-6:])

    parameters = dns.readParameters(rundir / "parameters.in")
    parameters["initiation"]["ic"] = i_final_state
    itime_final = i_final_state * parameters["output"]["i_save_fields"]
    parameters["initiation"]["i_start"] = itime_final

    if i_finish_plus is not None:
        parameters["termination"]["i_finish"] += i_finish_plus

    stat_file = rundir / "stat.gp"
    if Path.is_file(stat_file):
        stats = np.loadtxt(rundir / "stat.gp")
        times = t_final_state = stats[stats[:, 0] == itime_final]
        if len(times) > 0:
            t_final_state = stats[stats[:, 0] == itime_final][-1][1]
        else:
            t_final_state = itime_final * parameters["time_stepping"]["dt"]
    else:
        t_final_state = itime_final * parameters["time_stepping"]["dt"]

    parameters["initiation"]["t_start"] = t_final_state
    dns.writeParameters(parameters, rundir / "parameters.in")

    files = list(rundir.glob("*.gp"))
    for file in files:
        with open(file, "r") as f:
            lines = f.readlines()

        with open(file, "w") as f:
            for line in lines:
                try:
                    i_time = int(re.search(r"\d+", line).group())

                    if i_time > itime_final:
                        break
                    else:
                        f.write(line)

                except:
                    f.write(line)

    if not script == None:
        call(["sbatch", script])


if __name__ == "__main__":
    main()
