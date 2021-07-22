#!/usr/bin/env python3
import argparse
from pathlib import Path

import dns


def main():

    parser = argparse.ArgumentParser(
        description="Inner product of two states, (state1, state2)."
    )
    parser.add_argument("state1", type=str, help="path to the first state.")
    parser.add_argument("state2", type=str, help="path to the second state.")
    args = vars(parser.parse_args())
    inprod = dnsinprod(**args)
    print(f"{inprod:25.16f}")


def dnsinprod(state1, state2):
    state1 = Path(state1)
    state2 = Path(state2)

    state1Data, header1 = dns.readState(state1)
    state2Data, header2 = dns.readState(state2)

    inprod = dns.inprod(state1Data, state2Data)

    if state1Data.shape != state2Data.shape:
        exit("Grid dimensions of the two states do not match.")

    return inprod


if __name__ == "__main__":
    main()
