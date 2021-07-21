#!/usr/bin/env python3
import argparse
from pathlib import Path
import os
import numpy as np
from pyvista.plotting.renderer import CameraPosition
import dns
import glob
import re
from subprocess import call


def main():

    parser = argparse.ArgumentParser(
        description="Continue a run from the last state file saved",
    )
    parser.add_argument("rundir", type=str, help="path to the run")
    parser.add_argument(
        "-script",
        dest="script",
        help="Submission script. If given, submit the job by sbatch script",
    )

    args = vars(parser.parse_args())

    print(args)

    dnscontinue(**args)


def dnscontinue(rundir, script):

    rundir = Path(rundir).resolve()
    os.chdir(rundir)

    states = list(glob.glob('state.*'))
    i_final_state = 0 
    for state in states:
        i_state = int(state[-6:])
        if i_state  > i_final_state:
            i_final_state = i_state
    
    parameters = dns.readParameters("parameters.in")
    parameters['initiation']['ic'] = i_final_state
    parameters['initiation']['i_start'] = i_final_state \
        * parameters['output']['i_save_fields']
    parameters['initiation']['t_start'] = parameters['initiation']['i_start'] \
                                        * parameters['time_stepping']['dt']
    dns.writeParameters(parameters, "parameters.in")

    files = list(glob.glob('*.gp'))
    print(files)

    for file in files:
        with open(file, "r") as f:
            lines = f.readlines()
            
        with open(file, "w") as f:
            for line in lines:
                try:
                    i_time = int(re.search(r'\d+', line).group())

                    if i_time >= parameters['initiation']['i_start']:
                        break
                    else:
                        print(i_time)
                        f.write(line)

                except:
                    pass 
                    f.write(line)

    if not script == None:
        call(["sbatch", script])


if __name__ == "__main__":
    main()
