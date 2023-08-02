#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

import dnsbox as dns


def main():

    parser = argparse.ArgumentParser(
        description="Inner product of two states, (state1, state2)."
    )
    parser.add_argument("state1", type=str, help="path to the first state.")
    parser.add_argument("state2", type=str, help="path to the second state.")
    parser.add_argument(
        "--nocompact",
        action="store_true",
        dest="nocompact",
    )
    parser.add_argument(
        "--relerr",
        action="store_true",
        dest="relerr",
    )
    args = vars(parser.parse_args())
    inprod = dnsinprod(**args)
    print(inprod)


def dnsinprod(state1, state2, nocompact=False, relerr=False):
    state1 = Path(state1)
    state2 = Path(state2)

    state1Data, header1 = dns.readState(state1, nocompact=nocompact)
    state2Data, header2 = dns.readState(state2, nocompact=nocompact)

    if relerr:
        delta = state1Data - state2Data
        normdelta = np.sqrt(dns.inprod(delta, delta))
        inprod = normdelta / np.sqrt(dns.inprod(state1Data, state1Data))
    else:
        inprod = dns.inprod(state1Data, state2Data)

    if state1Data.shape != state2Data.shape:
        exit("Grid dimensions of the two states do not match.")

    return inprod


if __name__ == "__main__":
    main()
