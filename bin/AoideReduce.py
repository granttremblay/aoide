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


    print("=======  CREATING MASTER DARK =======")
    os.system("OMP_NUM_THREADS={} esorex --log-file=dark.log muse_dark --nifu=-1 --merge dark.sof".format(cores))


    print("=======  CREATING MASTER FLAT =======")
    os.system("OMP_NUM_THREADS={} esorex --log-file=flat.log muse_flat --nifu=-1 --merge flat.sof".format(cores))


    print("=======      WAVELENGTH CALIBRATION     =======")
    os.system("OMP_NUM_THREADS={} esorex --log-file=wavecal.log muse_wavecal --nifu=-1 --resample --residuals --merge arc.sof".format(cores))


    print("=======       LINE SPREAD FUNCTION      =======")
    os.system("OMP_NUM_THREADS={} esorex --log-file=lsf.log muse_lsf --nifu=-1 --merge lsf.sof".format(cores))


    print("=======          TWILIGHT FLATS         =======")
    os.system("OMP_NUM_THREADS={} esorex --log-file=twilight.log muse_twilight twilight.sof".format(cores))


    print("=======    SCIBASIC for SCI FRAMES      =======")
    os.system("OMP_NUM_THREADS={} esorex --log-file=science_scibasic.log muse_scibasic --nifu=-1 --merge science_scibasic.sof".format(cores))



    print("=======    SCIBASIC for STD FRAME      =======")
    os.system("OMP_NUM_THREADS={} esorex --log-file=science_scibasic.log muse_scibasic --nifu=-1 --merge std_scibasic.sof".format(cores))
    


    print("=======        FLUX CALIBRATION         =======")
    os.system("OMP_NUM_THREADS={} esorex --log-file=fluxcal.log muse_standard --filter=white fluxcal.sof".format(cores))
    


    print("=======           SCIPOST              =======")
    os.system("OMP_NUM_THREADS={} esorex --log-file=fluxcal.log muse_standard --filter=white fluxcal.sof".format(cores))    



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
