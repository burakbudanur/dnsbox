#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from subprocess import call
import numpy as np

import dns


def main():

    parser = argparse.ArgumentParser(
        description="Set ic, modify t_start and i_start accordingly.",
    )
    parser.add_argument("rundir", type=str, help="path to the run")
    parser.add_argument("ic", type=int, help="new ic")

    args = vars(parser.parse_args())

    dnssetic(**args)

    
def dnssetic(rundir, ic):

    rundir = Path(rundir)
    
    parameters = dns.readParameters(rundir / "parameters.in")
    
    i_save_fields = parameters['output']['i_save_fields']
    i_start = ic * i_save_fields
        
    try:
        stats = np.loadtxt(rundir / "stat.gp")
        t_start = stats[np.argwhere(stats[:, 0] == i_start)[0][0], 1]
        print('read off t_start from stat.gp')
    except:
        print('treating as a fixed time-step run')
        t_start = i_start * parameters['time_stepping']['dt']
        
    parameters['initiation']['ic'] = ic
    parameters['initiation']['i_start'] = i_start
    parameters['initiation']['t_start'] = t_start
    
    dns.writeParameters(parameters, rundir / "parameters.in")
    

if __name__ == "__main__":
    main()
