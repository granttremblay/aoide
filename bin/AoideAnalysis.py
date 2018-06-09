#!/usr/bin/env python

import os
import sys
import glob
import time

import warnings


from aoide import make_maps
from aoide import musemovie
import imageio

import argparse


# Some things we'll be doing throw runtimewarnings that we won't care about.
warnings.filterwarnings('ignore')


def main():

    print("\n=======================   Aoide | Analysis   ====================\n")
    print("AoideAnalysis: Create Quick-look Stellar & Line (Flux/Velocity/Ratio) Maps. Also makes a cool H-alpha movie!\n")

    args = parse_args()

    current_directory = os.getcwd()
    file_save_directory = current_directory + '/quicklook/'

    if not os.path.exists(file_save_directory):
        os.makedirs(file_save_directory)
        print("Quicklook Products Directory {} not found. Creating it.".format(
            file_save_directory))
    else:
        print("Quicklook products will be saved to {}".format(file_save_directory))

    fovimage = 'IMAGE_FOV_0001.fits'
    prefix = args.prefix
    redshift = args.redshift

    kintable = '{}.kin_table.fits'.format(prefix)
    elinetable = '{}.eline_table.fits'.format(prefix)

    make_maps.make_fovimage(fovimage, name = prefix,
                            save=True, file_save_directory=file_save_directory)

    make_maps.make_stars(kintable,
                         fovimage,
                         redshift,
                         name=prefix,
                         snthresh=1,
                         velocity_thresh=None,  # in km/s
                         #cropbox=(20, 160, 65, 180),  # (x1, x2, y1, y2),
                         zoom=False,  # It will zoom on the cropbox
                         wcs=False,
                         vel_vmin=-500,
                         vel_vmax=500,
                         disp_vmin=0,
                         disp_vmax=1800,
                         colorbar_scalefactor=0.0395,  # I'm too lazy to properly implement colorbar resizing
                         project_by_median=True,
                         project_by_redshift=False,
                         save=True,
                         file_save_directory=file_save_directory)

    make_maps.make_elines(elinetable, 'Halpha',
                          fovimage,
                          redshift,
                          name=prefix,
                          snthresh=1,
                          velocity_thresh=None,  # in km/s
                          # cropbox=(20, 160, 65, 180), # (x1, x2, y1, y2),
                          zoom=False,  # It will zoom on the cropbox
                          wcs=True,
                          flux_vmin=0,
                          flux_vmax=1,
                          fluxlog=False,
                          vel_vmin=-150,
                          vel_vmax=150,
                          fwhm_vmin=0,
                          fwhm_vmax=300,
                          colorbar_scalefactor=0.047,  # I'm too lazy to properly implement colorbar resizing
                          project_by_median=True,
                          project_by_redshift=False,
                          save=True,
                          file_save_directory=file_save_directory)

    #
    #
    # make_balmer_map(Ha_flux_map, Hb_flux_map, hdr)
    # make_electron_density_map(SII6717_flux_map, SII6730_flux_map, hdr)
    #
    #
    # musemovie.make_movie(cube, redshift, center, name, thresh=thresh, frames=frames, scalefactor=scalefactor, contsub=contsub, whitebg=whitebg, linear=linear)


def parse_args():

    parser = argparse.ArgumentParser(description="Create quicklook maps, etc.")

    parser.add_argument('prefix', metavar='prefix', type=str, nargs='?',
                        help='Prefix (usually a galaxy name) that prepends all standard AOIDE-named files, including the *stellar_tab and *kin_tab files.')

    parser.add_argument('redshift', type=float,
                        nargs='?', help='Source redshift')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    start_time = time.time()
    main()
    runtime = round((time.time() - start_time), 3)
    print("\n=================    Aoide | Analysis Finished   ================\n")
    print("Finished in {} minutes.".format(round(runtime / 60, 3)))
    print("Use these cubes for Voronoi Binning & Paradise Fitting.")
