#!/usr/bin/env python
'''
Aoide | Reduction & Analysis of MUSE observations
-------------------------------------------------
Dr. Grant R. Tremblay | Harvard-Smithsonian Center for Astrophysics
grant.tremblay @ cfa.harvard.edu

See the README associated with this repository for documentation & examples.
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


    print("\n=======================   Aoide | Step 1   ====================\n")
    print("AoideReduce.py: Reduce raw MUSE frames to a (mostly) final datacube.\n")

    args = parse_args()

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

    print("\n______________________  Starting MUSE Pipeline    _______________________\n")
    #################### RUN THE PIPELINE #####################

    print("=======  STARTING REDUCTION OF MUSE DATA =======")
    print("Setting number of OMP threads to {} CPU cores".format(cores))

    os.chdir(reduction_dir)
    print("Changed working directory to {}.".format(reduction_dir))


    #### Switches to skip steps, as requested #######
    skip_bias = args.skip_bias
    skip_dark = args.skip_dark
    skip_flat = args.skip_flat
    skip_arc = args.skip_arc
    skip_lsf = args.skip_lsf
    skip_twilight = args.skip_twilight
    skip_science_scibasic = args.skip_science_scibasic
    skip_std_scibasic = args.skip_std_scibasic
    skip_fluxcal = args.skip_fluxcal
    skip_scipost = args.skip_scipost
    skip_combine = args.skip_combine

    if args.skip_existing is True:
        if os.path.isfile('MASTER_BIAS.fits'):
            skip_bias = True
        if os.path.isfile('MASTER_DARK.fits'):
            skip_dark = True
        if os.path.isfile('MASTER_FLAT.fits'):
            skip_flat = True
        if os.path.isfile('WAVECAL_TABLE.fits'):
            skip_arc = True
        if os.path.isfile('LSF_PROFILE.fits'):
            skip_lsf = True
        if os.path.isfile('TWILIGHT_CUBE.fits'):
            skip_twlight = True
        if os.path.isfile('PIXTABLE_OBJECT_0001-23.fits'):
            skip_science_scibasic = True
        if os.path.isfile('PIXTABLE_STD_0001-23.fits'):
            skip_std_scibasic = True
        if os.path.isfile('PIXTABLE_STD_0001-23.fits'):
            skip_std_scibasic = True
        if os.path.isfile('STD_TELLURIC_0001.fits'):
            skip_fluxcal = True
        if os.path.isfile('DATACUBE_SINGLE_FINAL_0001_1.fits'):
            skip_scipost = True
        if os.path.isfile('DATACUBE_AOIDE_UNCLEAN.fits'):
            skip_combine = True


    if skip_bias is False:
        print("=======  CREATING MASTER BIAS =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=bias.log muse_bias --nifu=-1 --merge bias.sof".format(cores))

    if skip_dark is False:
        print("=======  CREATING MASTER DARK =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=dark.log muse_dark --nifu=-1 --merge dark.sof".format(cores))

    if skip_flat is False:
        print("=======  CREATING MASTER FLAT =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=flat.log muse_flat --nifu=-1 --merge flat.sof".format(cores))

    if skip_arc is False:
        print("=======      WAVELENGTH CALIBRATION     =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=wavecal.log muse_wavecal --nifu=-1 --resample --residuals --merge arc.sof".format(cores))

    if skip_lsf is False:
        print("=======       LINE SPREAD FUNCTION      =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=lsf.log muse_lsf --nifu=-1 --merge lsf.sof".format(cores))

    if skip_twilight is False:
        print("=======          TWILIGHT FLATS         =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=twilight.log muse_twilight twilight.sof".format(cores))

    if skip_science_scibasic is False:
        print("=======    SCIBASIC for SCI FRAMES      =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=science_scibasic.log muse_scibasic --nifu=-1 --merge science_scibasic.sof".format(cores))

    if skip_std_scibasic is False:
        print("=======    SCIBASIC for STD FRAME      =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=std_scibasic.log muse_scibasic --nifu=-1 --merge std_scibasic.sof".format(cores))

    if skip_fluxcal is False:
        print("=======        FLUX CALIBRATION         =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=fluxcal.log muse_standard --filter=white fluxcal.sof".format(cores))

    if skip_scipost is False:
        print("=======   SCIENCE POSTPROCESSING  =======")

        # Count the number of science_scipost_N files created
        science_scipost_files = sorted(glob.glob("science_scipost_*.sof"))

        # Run muse_scipost for each science exposure. We'll combine them later.
        for i in range(len(science_scipost_files)):
            print("BEGGINING SCIPOST RUN #{}".format(i+1))
            os.system("OMP_NUM_THREADS={} esorex --log-file=sci_scipost_{}.log muse_scipost --filter=white,Johnson_V,Cousins_R,Cousins_I --save=cube,individual --skymodel_fraction=0.3 --skymethod=simple science_scipost_{}.sof".format(cores, i + 1, i + 1))

            os.rename("IMAGE_FOV_0001.fits", "IMAGE_FOV_0001_{}.fits".format(i+1))
            os.rename("PIXTABLE_REDUCED_0001.fits", "PIXTABLE_REDUCED_0001_{}.fits".format(i+1))
            os.rename("DATACUBE_FINAL.fits", "DATACUBE_SINGLE_FINAL_0001_{}.fits".format(i+1))

    if skip_combine is False:
        # Combine the three reduced pixtables, trusting the default WCS solution.
        print("=======   FINAL COMBINATION  =======")
        os.system("OMP_NUM_THREADS={} esorex --log-file=combine.log muse_exp_combine --pixfrac=0.8 --filter=white,Johnson_V,Cousins_R,Cousins_I combine.sof".format(cores))
        os.rename("DATACUBE_FINAL.fits", "DATACUBE_AOIDE_UNCLEAN.fits")

def parse_args():

    parser = argparse.ArgumentParser(description="Reduce MUSE cube from raw archive files.")

    parser.add_argument('-r', '--rawdata', default='./', metavar='WORKING_DIRECTORY')

    parser.add_argument('--static_cal_dir', default='/home/grant/Software/ESO/MUSE/calib/muse-2.4.1/cal', metavar='STATIC_CALIBRATION_DIRECTORY')

    parser.add_argument('-c', '--cores', type=str, nargs='?',
                        default='6', help='Number of cores used for procesing.')

    parser.add_argument('--testsetup', action="store_true", default=False,
                        help='Flag to skip everything and simply test setup.')

    parser.add_argument('--skip_existing', action="store_true", default=False,
                        help='Flag to skip creation of existing data products.')

    parser.add_argument('--skip_bias', action="store_true", default=False,
                        help='Flag to skip master bias creation.')

    parser.add_argument('--skip_dark', action="store_true", default=False,
                        help='Flag to skip master dark creation.')

    parser.add_argument('--skip_flat', action="store_true", default=False,
                        help='Flag to skip master flat creation.')

    parser.add_argument('--skip_arc', action="store_true", default=False,
                        help='Flag to skip wavelength calibration.')

    parser.add_argument('--skip_lsf', action="store_true", default=False,
                        help='Flag to skip lsf profile creation.')

    parser.add_argument('--skip_twilight', action="store_true", default=False,
                        help='Flag to skip skyflat creation.')

    parser.add_argument('--skip_science_scibasic', action="store_true", default=False,
                        help='Flag to skip the SCIBASIC step for science frames.')

    parser.add_argument('--skip_std_scibasic', action="store_true", default=False,
                        help='Flag to skip the SCIBASIC step for standard star frames.')

    parser.add_argument('--skip_fluxcal', action="store_true", default=False,
                        help='Flag to skip flux calibration.')

    parser.add_argument('--skip_scipost', action="store_true", default=False,
                        help='Flag to skip final processing of the science observations.')

    parser.add_argument('--skip_combine', action="store_true", default=False,
                        help='Flag to skip alignment & combination of final science frames.')


    args = parser.parse_args()

    return args


if __name__ == '__main__':

    start_time = time.time()
    main()
    runtime = round((time.time() - start_time), 3)
    print("\n=====================    Aoide | Step 1 Finished   ====================\n")
    print("Finished in {} minutes.".format(round(runtime / 60, 3)))
