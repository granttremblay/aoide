#!/usr/bin/env python
'''
Reduce a MUSE cube | G. Tremblay
'''
import os
import sys
import glob
import time

from distutils import spawn
from shutil import copyfile

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

    print("\nRaw data directory set to {}".format(raw_data_dir))
    if os.path.exists(reduction_dir):
        print("Data products will be placed in {}".format(reduction_dir))
    else:
        os.makedirs(reduction_dir)
        print("Creating data products directory {}".format(reduction_dir))

    print("\nUsing {} processor cores for reduction.".format(args.cores))



    print("\nesorex path is {}".format(esorex_path))

    static_cal_dir = os.path.abspath(args.static_cal_dir)
    print("\nStatic Calibration Files directory is {}".format(static_cal_dir))
    print("If this is incorrect, you can specify it with --static_cal_dir.")

    if args.testsetup is True:
        sys.exit("--testsetup was set, exiting before reduction.")

    #################### CREATE SOF FILES #####################

    print("\n_________________________   Creating SOF Files    _______________________\n")
    sof.make_sof_files(raw_data_dir, reduction_dir, static_cal_dir)

    os.chdir(reduction_dir)
    print("Changed working directory to {}.".format(reduction_dir))


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


    print("=======   SCIENCE POSTPROCESSING  =======")

    # Count the number of science_scipost_N files created

    science_scipost_files = sorted(glob.glob("science_scipost_*.sof"))

    for i in range(len(science_scipost_files)):
        print("BEGGINING SCIPOST RUN #{}".format(i+1))
        os.system("OMP_NUM_THREADS={} esorex --log-file=sci_scipost_{}.log muse_scipost --filter=white,Johnson_V,Cousins_R,Cousins_I --save=cube,individual --skymodel_fraction=0.3 --skymethod=simple science_scipost_{}.sof".format(cores, i + 1, i + 1))

        os.rename("IMAGE_FOV_0001.fits", "IMAGE_FOV_0001_{}.fits".format(i+1))
        os.rename("PIXTABLE_REDUCED_0001.fits", "PIXTABLE_REDUCED_0001_{}.fits".format(i+1))
        os.rename("DATACUBE_FINAL.fits", "DATACUBE_SINGLE_FINAL_0001_{}.fits".format(i+1))


    print("=======   FINAL COMBINATION  =======")
    os.system("OMP_NUM_THREADS={} muse_exp_combine --pixfrac=0.8 --filter=white,Johnson_V,Cousins_R,Cousins_I combine.sof".format(cores))


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
    print("Your data products are in {}.".format(reduction_dir))
