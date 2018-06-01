#!/usr/bin/env python
'''
Reduce a MUSE cube | G. Tremblay
'''
import os
import sys

import time

import matplotlib.pyplot as plt

import argparse

from aoide import make_sof_files as sof



def main():

    cores = 6

    science_targname = 'Abell 3581'

    parent_dir = '/home/grant/Science/museA3581/Data/MUSE'

    raw_data_dir = parent_dir + '/rawdata'
    reduction_dir = parent_dir + '/reduction'

    static_cal_dir = '/home/grant/Software/ESO/MUSE/calib/muse-2.4.1/cal'
    esorex = '/home/grant/Software/ESO/MUSE/bin/esorex'

    sof.make_sof_files(raw_data_dir, reduction_dir, static_cal_dir, science_targname)

    os.chdir(reduction_dir)

    print("=======  STARTING REDUCTION OF MUSE DATA =======")
    print("Setting number of OMP threads to {} CPU cores".format(cores))


    print("=======  CREATING MASTER BIAS =======")
    os.system("OMP_NUM_THREADS={} esorex --log-file=bias.log muse_bias --nifu=-1 --merge bias.sof".format(cores))


# def parse_args():

#     parser = argparse.ArgumentParser(
#         description="Produce MUSE cube from raw data")

#     parser.add_argument('input_cube', metavar='CUBE_IN', type=str, nargs='?',
#                         help='Input FITS datacube from which sky spectra are selected')

#     args = parser.parse_args()

#     return args



if __name__ == '__main__':

    start_time = time.time()
    main()
    runtime = round((time.time() - start_time), 3)
    print("\n=====================    AOIDE is Finished   ====================\n")
    print("Finished in {} minutes.".format(round(runtime / 60, 3)))
