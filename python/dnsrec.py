#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

import dns


def main():
    parser = argparse.ArgumentParser("Computes recurrence among a series of states.")
    parser.add_argument(
        "rundir",
        type=str,
        help="Path to the run directory.",
    )
    parser.add_argument(
        "--si", type=int, help="Initial state.",
    )
    parser.add_argument(
        "--sf", type=int, help="Final state.",
    )
    parser.add_argument(
        "--nrec", type=int, help="Recurrence horizon.",
    )
    parser.add_argument(
        "--sliced", type=bool, default="False", 
        help="Compute for *_sliced states if true.",
    )

    args = parser.parse_args()
    rundir = args.rundir
    si = args.si
    sf = args.sf
    nrec = args.nrec
    sliced = args.sliced

    dnsrec(rundir, si, sf, nrec, sliced)


def dnsrec(rundir, si, sf, nrec, sliced):
    """
    Compute and save the recurrence matrix from si to sf in rundir upto the 
    recurrence horizon nrec.
    Use state*_sliced  
    """

    rundir = Path(rundir)
    rec_mat = np.zeros((sf - si + 1, nrec + 1))

    for i in range(sf + 1 - si):

        i_state = i + si
        fstate_i = f"state:{i_state:06}"

        if sliced:
            fstate_i += "_sliced"
        
        try:
            state_i, _ = dns.readState(rundir / fstate_i)
        except:
            print(f"Cannot read f{str(rundir / fstate_i)}")
            break

        for j in range(1, nrec + 1):

            j_state = i_state + j
            fstate_j = f"state:{j_state:06}"

            if sliced:
                fstate_j += "_sliced"

            try:
                state_j, _ = dns.readState(rundir / fstate_j)
            except:
                print(f"Cannot read f{str(rundir / fstate_j)}")
                break
            
            difference = state_j - state_i
            rec_mat[i, j] = (
                np.sqrt(dns.inprod(difference, difference)) 
                / np.sqrt(dns.inprod(state_i, state_i))
                )

    if sliced:
        outfile = f"rec_sliced_si_{si}_sf_{sf}_nrec_{nrec}"
    else:
        outfile = f"rec_si_{si}_sf_{sf}_nrec_{nrec}"

    outfile = str(rundir / outfile)

    np.savetxt(outfile, rec_mat)

    return 

if __name__ == "__main__":
    main()
