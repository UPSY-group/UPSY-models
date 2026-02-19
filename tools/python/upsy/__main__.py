"""
Script for diagnosing a UFEMISM or LADDIE run
"""

import argparse

from upsy.run import Run

def main():
    parser = argparse.ArgumentParser(
    description='Diagnose the run'
    )

    parser.add_argument(
        'rundir',
        help='Run directory where output is stored')

    args = parser.parse_args()

    run = Run(args.rundir)

    print('===================')
    print('\033[35m  Vertex variables:\033[0m')
    print(run.vars_vi)
    print('=====================')
    print('\033[35m  Triangle variables:\033[0m')
    print(run.vars_ti)
    print('===========')
    print('\033[35m  Contours:\033[0m')
    print(run.contours)
    print('==============')
    print('\033[35m  Time values:\033[0m')
    for m in run.times:
        print(f"\033[33mMesh {m}\033[0m : {run.times[m]}")
