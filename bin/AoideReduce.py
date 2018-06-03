#!/usr/bin/env python
'''
Reduce a MUSE cube | G. Tremblay
'''
import os
import sys
import time

from distutils import spawn

import argparse

from aoide import make_sof_files as sof


def main():


    print("\n======================    Welcome to AoideReduce   ====================\n")
    print("Prepare Pipeline-reduced MUSE datacubes for Paradise\n")

    args = parse_args()

    skip_existing = args.skip_existing
    cores = args.cores

    ################## TELL USER THEIR SETUP INFO ############

    print("_________________________   Setting Up    _______________________\n")

    raw_data_dir = os.path.abspath(args.rawdata)
    reduction_dir = os.path.abspath(os.path.join(raw_data_dir, '../reduction'))
    esorex_path = spawn.find_executable("esorex")

    print("Raw data directory set to {}".format(raw_data_dir))
    print("Using {} processor cores for reduction.".format(args.cores))

    if os.path.exists(reduction_dir):
        print("Data products will be placed in {}".format(reduction_dir))
    else:
        os.makedirs(reduction_dir)
        print("Creating data products directory {}".format(reduction_dir))

    print("esorex path is {}".format(esorex_path))

    static_cal_dir = os.path.abspath(args.static_cal_dir)
    print("Static Calibration Files directory is {}".format(static_cal_dir))
    print("IF YOU'RE NOT RUNNING THIS ON GRANT TREMBLAY'S MACHINES, this is PROBABLY INCORRECT!")
    print("Specify the correct directory with --static_cal_dir. If you've done so, ignore this. ")

    if args.testsetup is True:
        sys.exit("--testsetup was set, exiting before reduction.")

    #################### CREATE SOF FILES #####################

    print("_________________________   Creating SOF Files    _______________________\n")
    sof.make_sof_files(raw_data_dir, reduction_dir, static_cal_dir, science_targname)

    os.chdir(reduction_dir)

    if skip_existing is True:
        if os.path.isfile('MASTER_BIAS.FITS'):
            print("Found MASTER BIAS")
            skip_bias = True
        else:
            skip_bias = False

        if os.path.isfile('MASTER_DARK.FITS'):
            skip_dark = True
        else:
            skip_dark = False

        if os.path.isfile('MASTER_FLAT.FITS'):
            skip_flat = True
        else:
            skip_flat = False

    print("=======  STARTING REDUCTION OF MUSE DATA =======")
    print("Setting number of OMP threads to {} CPU cores".format(cores))

    if skip_bias is False:
        print("=======  CREATING MASTER BIAS =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=bias.log muse_bias --nifu=-1 --merge bias.sof".format(cores))
    elif skip_bias is True:
        print("SKIPPING MASTER BIAS CREATION")

    if skip_dark is False:
        print("=======  CREATING MASTER DARK =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=dark.log muse_dark --nifu=-1 --merge dark.sof".format(cores))
    elif skip_dark is True:
        print("SKIPPING MASTER DARK CREATION")

    if skip_flat is False:
        print("=======  CREATING MASTER FLAT =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=flat.log muse_flat --nifu=-1 --merge flat.sof".format(cores))
    elif skip_flat is True:
        print("SKIPPING MASTER FLAT CREATION")


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



def parse_args():

    parser = argparse.ArgumentParser(description="Reduce MUSE cube from raw archive files.")

    parser.add_argument('-r', '--rawdata', default='./', metavar='WORKING_DIRECTORY')

    parser.add_argument('--static_cal_dir', default='/home/grant/Software/ESO/MUSE/calib/muse-2.4.1/cal', metavar='STATIC_CALIBRATION_DIRECTORY')

    parser.add_argument('--skip_existing', action="store_true", default=False,
                        help='Flag to skip creation of existing data products.')

    parser.add_argument('-c', '--cores', type=str, nargs='?',
                        default='6', help='Number of cores used for procesing.')

    parser.add_argument('--testsetup', action="store_true", default=False,
                        help='Flag to skip everything and simply test setup.')

    args = parser.parse_args()

    return args



if __name__ == '__main__':

    start_time = time.time()
    main()
    runtime = round((time.time() - start_time), 3)
    print("\n=====================    AOIDE is Finished   ====================\n")
    print("Finished in {} minutes.".format(round(runtime / 60, 3)))
