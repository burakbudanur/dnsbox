#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from subprocess import call
import numpy as np

import dnsbox as dns


def main():

    parser = argparse.ArgumentParser(
        description="Delete states in a directory keeping every n. Rename states that are kept so that the remaining are state.000000, state.000001, ..., change parameters.in accordingly.",
    )
    parser.add_argument("rundir", type=str, help="path to the run")
    parser.add_argument("n", type=int, help="states other than 0, n, 2n, ... will be deleted")

    args = vars(parser.parse_args())

    dnsrmstates(**args)

    
def dnsrmstates(rundir, n):

    rundir = Path(rundir)
    parameters = dns.readParameters(rundir / "parameters.in")
    parameters["output"]["i_save_fields"] *= n
    dns.writeParameters(parameters, rundir / "parameters.in")
    
    states = sorted(list(rundir.glob("state.*")))
    perturbed_states = sorted(list(rundir.glob("perturb.*")))

    if len(states) > 0:
        for state in states:
            i_state = int(state.name[-6:]) 
            if i_state%n == 0:
                state.replace(rundir / Path(f'state.{i_state//n:06}'))
            else: 
                state.unlink()

    if len(perturbed_states) > 0:
        for state in perturbed_states:
            i_state = int(state.name[-6:]) 
            if i_state%n == 0:
                state.replace(rundir / Path(f'perturb.{i_state//n:06}'))
            else: 
                state.unlink()

if __name__ == "__main__":
    main()
