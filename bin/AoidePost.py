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

import time

import matplotlib.pyplot as plt

import argparse

from aoide import aoide_postprocess
from aoide import make_sky_mask


def main():
    '''
    The main Aoide function
    '''

    styleplots()

    print("\n=======================   Aoide | Step 2   ====================\n")
    print("AoidePost.py: Prepare Pipeline-reduced MUSE datacubes for Paradise\n")

    # If the user hasn't given any arguments, complain.
    if not len(sys.argv) > 1:
        print("\nNo arguments given. Type aoide -h for usage help.\n")
        sys.exit()

    args = parse_args()

    if args.dontask is True:
        interactive = False
    elif args.dontask is False:
        interactive = True


    print("_________________________   Setting Up    _______________________\n")
    print("\nAoide will be run in parallel using {} CPU cores.".format(args.cores))
    print("Working directory set to {} (change with -d)\n".format(args.working_directory))
    print("PCA fit will use {} components and {} spectra".format(args.pca_components, args.spectra))
    print("Cube 2x2 binning is set to {}".format(args.binning))
    print("A_V (foreground Galactic extinction at V-band) is set to {}.".format(args.av))

    if interactive is True:
        proceed = yes_no("Are you happy with this setup? Proceed? [yes/no]: ")

        if proceed is False:
            sys.exit("You asked to reset parameters before proceeding. Finished.")


    # Set directories & keep things clean.
    working_directory = args.working_directory
    skymask_directory = working_directory + "SKY_MASKS/"

    # Set critical filenames, based on the above directories.
    skymask_name = skymask_directory + "SKY_MASK.fits"
    pcamodel_name = skymask_directory + "PCA_SKY.fits"
    skysub_cube_name = working_directory + "DATACUBE_SKYSUB.fits"
    final_cube_name = working_directory + "DATACUBE_AOIDE_FINAL.fits"


    dirty_cube = args.input_cube
    fovimage = args.fovimage
    muse_data_extension = args.extension

    mask_sky(fovimage, skymask_directory, skymask_name,
             muse_data_extension, interactive)

    subtract_sky(dirty_cube, skysub_cube_name, pcamodel_name, skymask_name,
                 args.filter, args.spectra, args.pca_components, args.cores)

    correct_cube(skysub_cube_name, final_cube_name, args.av, args.binning)


def subtract_sky(dirty_cube, skysub_cube_name, pcamodel_name, skymask_name, filter, numspectra, components, parallel):

    print("\n________________   Fit & Subtract PCA Sky Model    _______________\n")

    print("\nPerforming PCA on sky using a median filter of {} and {} spectra.".format(
        filter, numspectra))
    aoide_postprocess.create_PCA_sky(
        dirty_cube, pcamodel_name, skymask_name, filter, numspectra, parallel)

    print("Fitting & Subtracting PCA model using {} PCA components and {} CPU cores.".format(
        components, parallel))
    print("This could take over an hour. Please be patient.")
    aoide_postprocess.subtract_PCA_sky(dirty_cube, skysub_cube_name, pcamodel_name, components, filter, parallel)

    print("\nFinally finished with PCA modeling. Created {}.".format(skysub_cube_name))


def mask_sky(fovimage, skymask_directory="SKY_MASKS/", skymask_name="SKY_MASKS/SKY_MASK.fits", muse_data_extension=1, interactive=True):

    print("______________________   Create Sky Mask    _____________________\n")

    if interactive is False:
        if os.path.exists(skymask_name):
            print("\nUsing Sky Mask {}.".format(skymask_name))
            return
        else:
            print(
                "\nTo run non-interactively, {} must already exist.".format(skymask_name))
            sys.exit(
                "ERROR: {} does not exist and --dontask is True".format(skymask_name))

    if interactive is True:
        if os.path.exists(skymask_name):
            print("\nSky Mask already exists at {}".format(skymask_name))
            edit_skymask = yes_no(
                'Would you like to edit & overwrite it? [yes/no]: ')

            if edit_skymask is False:
                print("Using existing {}, skipping this step.".format(skymask_name))
                return

    if not os.path.exists(skymask_directory):
        os.makedirs(skymask_directory)
        print("Sky Mask Directory {} not found. Creating it.".format(
            skymask_directory))

    print("Sky Mask Components will be saved in {}".format(skymask_directory))
    print("\nA window should appear.")
    print("Create boxes over BLANK SKY regions by click dragging LEFT mouse.")
    print("Make sure you draw ONLY over blank sky regions. Avoid the galaxy!")
    print("RIGHT mouse click-drag will fix mistakes.")
    print("Close the window to finish and save your Sky Mask.")

    make_mask = make_sky_mask.MaskFrame(
        fovimage, skymask_name, extension=int(muse_data_extension))
    plt.show()

    make_mask.save_mask()
    print("\nSky Mask created, saved to {}.".format(skymask_name))

def correct_cube(incube, outcube="DATACUBE_AOIDE_FINAL.fits", av, binning=False):

    print("\n__________________   Correct & Finalize Cube    _________________\n")

    aoide_postprocess.correct_cube(incube, outcube, av, binning)

    if binning is False:
        print("Finished. Corrected, unbinned cube saved to {}".format(outcube))
        print("You should copy this cube to a separate file (i.e. GALAXYNAME_DATACUBE_AOIDE_FINAL.fits) to prevent it being overwritten.")
    elif binning is True:
        print("Finished. Corrected, 2x2 binned cube saved to {}".format(outcube))
        print("You should copy this cube to a separate file (i.e. GALAXYNAME_DATACUBE_AOIDE_FINAL_BINNED.fits) to prevent it being overwritten.")

def parse_args():

    parser = argparse.ArgumentParser(
        description="PCA sky subtraction, galactic extinction correction, and preparation for MUSE cubes that are to be analyzed with Paradise")

    parser.add_argument('input_cube', metavar='CUBE_IN', type=str, nargs='?',
                        help='Input FITS datacube from which sky spectra are selected')

    parser.add_argument('output_cube', metavar='CUBE_OUT', type=str,
                        nargs='?', help='Output FITS datacube with sky residuals subtracted')

    parser.add_argument('-d', '--working_directory', metavar='WORKING_DIRECTORY', type=str,
                        default='./', nargs='?', help='Working directory in which \
                        files will be written. Default is current.')

    parser.add_argument('pca_sky', metavar='PCA', type=str,
                        nargs='?', help='Input FITS file with PCA components')

    parser.add_argument('-p', '--pca_components', type=int, nargs='?',
                        default=150, help='Number of PCA components to be used')

    parser.add_argument('-f', '--filter', type=int, nargs='?', default=30,
                        help='Size of median filter in wavelength direction to remove continuum signal before sky residual subtraction')

    parser.add_argument('-s', '--spectra', type=int, nargs='?', default=20000,
                        help='Maximum number  of spectra considered for PCA analysis. The number of selected by the MASK is always the absolute maximum and can only be reduced.')

    parser.add_argument('-i', '--fovimage', default='IMAGE_FOV_0001.fits',
                        help='FOV image for setting 2D dimesions. Default is IMG_FOV_0001.fits.')

    parser.add_argument('-c', '--cores', type=str, nargs='?',
                        default='6', help='Number of cores used for computation.')

    parser.add_argument('-a', '--av', required=True,
                        help='A_V (extinction at V-band). You can get this from NED for your source. R_V=3.1 is assumed for the galactic extinction correction.')

    parser.add_argument('-b', '--binning', action="store_true", default=False,
                        help='Flag to bin the cube 2x2.')

    parser.add_argument('-e', '--extension', default=1,
                        help='MUSE Data extension.')

    parser.add_argument('--cleanup', action="store_true", default=False,
                        help='Flag to clean up temporary working directories. Default is False.')

    parser.add_argument('--dontask', action="store_true", default=False,
                        help='Flag to suppress interactivity. SKY_MASK must be in place.')


    args = parser.parse_args()

    return args


def styleplots():
    """
    Make plots pretty and labels clear.
    """
    # plt.style.use('ggplot')

    labelsizes = 18

    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rcParams['font.size'] = labelsizes
    plt.rcParams['axes.titlesize'] = labelsizes - 2
    plt.rcParams['axes.labelsize'] = labelsizes
    plt.rcParams['xtick.labelsize'] = labelsizes
    plt.rcParams['ytick.labelsize'] = labelsizes


def yes_no(answer):
    '''
    A simple function to query the user a yes/no question.
    '''
    yes = set(['yes', 'y', 'ye', ''])
    no = set(['no', 'n'])

    while True:
        choice = input(answer).lower()
        if choice in yes:
            return True
        elif choice in no:
            return False
        else:
            print("Please respond with 'yes' or 'no'\n")


if __name__ == '__main__':

    start_time = time.time()
    main()
    runtime = round((time.time() - start_time), 3)
    print("\n=================    Aoide | Step 2 Finished   ================\n")
    print("Finished in {} minutes.".format(round(runtime / 60, 3)))
    print("You should now have a clean datacube ready for Voronoi Binning & Paradise Fitting.")
