#!/usr/bin/env python
import os
import time

import numpy as np
from astropy.io import fits


from aoide import voronoi

import argparse


def main():

    print("\n=======================   Aoide | Step 3   ====================\n")
    print("AoideVoronoi.py: Voronoi Tesselation based on Stellar Continuum Signal\n")

    args = parse_args()

    input_cube = os.path.abspath(args.input_cube)
    name = args.name
    continuum_signal_band = args.band_continuum
    target_sn = args.sn
    max_area = args.max_area

    working_directory = os.path.dirname(os.path.abspath(input_cube))
    voronoi_directory = working_directory + "/voronoi/"

    if not os.path.exists(voronoi_directory):
        os.makedirs(voronoi_directory)
        print("Voronoi Products Directory {} not found. Creating it.".format(
            voronoi_directory))
    else:
        print("Voronoi products will be saved to {}".format(voronoi_directory))

    # Change to Voronoi directory
    os.chdir(voronoi_directory)
    print("Working directory set to {}".format(voronoi_directory))

    # Work on input cube
    hdu = fits.open(input_cube)
    data = hdu[0].data
    select = np.isnan(data)
    data[select] = 0
    error = np.sqrt(hdu[1].data)
    error[select] = 1e9
    hdr = hdu[0].header
    wave = np.arange(hdr['NAXIS3']) * 1.25 + hdr['CRVAL3']
    hdu.close()
    select_wave = (wave > continuum_signal_band[0]) & (
        wave < continuum_signal_band[1])
    mean_img = np.mean(data[select_wave, :, :], 0)
    err_img = np.std(data[select_wave, :, :], 0)

    # Create CONTINUUM SIGNAL and NOISE samples on which to iterate
    hdu = fits.PrimaryHDU(mean_img)
    hdu.writeto('signal_cont.fits', overwrite=True)

    hdu = fits.PrimaryHDU(err_img)
    hdu.writeto('noise_cont.fits', overwrite=True)

    voronoi.binning_fits('signal_cont.fits', voronoi_directory + 'noise_cont.fits',
                         target_sn, name, max_area=max_area, gersho=True, plot=False)

    bin_cube('{}.voronoi.pixel.dat'.format(name),
             input_cube, '{}.voronoi_rss.fits'.format(name))


def bin_cube(voronoi_pixel, input_cube, output_rss):
    hdu = fits.open(input_cube)
    data = hdu[0].data
    select = np.isnan(data)
    data[select] = 0
    error = np.sqrt(hdu[1].data)
    error[select] = 1e9
    hdr = hdu[0].header
    wave = np.arange(hdr['NAXIS3']) * 1.25 + hdr['CRVAL3']
    hdu.close()

    ascii_in = open(voronoi_pixel, 'r')
    lines = ascii_in.readlines()
    x = np.arange(len(lines) - 1, dtype=np.int16)
    y = np.arange(len(lines) - 1, dtype=np.int16)
    binNr = np.arange(len(lines) - 1, dtype=np.int16)
    for i in range(1, len(lines)):
        line = lines[i].split()
        x[i - 1] = int(line[0])
        y[i - 1] = int(line[1])
        binNr[i - 1] = int(line[2])

    rss_data = np.zeros((max(binNr), len(wave)), dtype=np.float32)
    rss_error = np.zeros((max(binNr), len(wave)), dtype=np.float32)

    for l in range(max(binNr)):
        select_bin = binNr == (l + 1)
        x_bin = x[select_bin] - 1
        y_bin = y[select_bin] - 1
        for j in range(len(x_bin)):
            rss_data[l, :] += data[:, y_bin[j], x_bin[j]]
            rss_error[l, :] += error[:, y_bin[j], x_bin[j]]**2
        rss_error[l, :] = np.sqrt(rss_error[l, :])
    rss_error[np.isnan(rss_error)] = 1e9
    rss_out = fits.HDUList([fits.PrimaryHDU(rss_data), fits.ImageHDU(
        rss_error, name='ERROR'), fits.ImageHDU(np.zeros(rss_data.shape, dtype=np.uint16), name='MASK')])
    rss_out[0].header = hdr
    rss_out[0].header['CDELT1'] = 1.25
    rss_out[0].header['CRVAL1'] = hdr['CRVAL3']
    rss_out[0].header['CRPIX1'] = 1
    rss_out.writeto(output_rss, overwrite=True, output_verify='fix')


def parse_args():

    parser = argparse.ArgumentParser(
        description="Voronoi bin a clean MUSE datacube based on S/N in the stellar continuum.")

    parser.add_argument('input_cube', metavar='CUBE_IN', type=str, nargs='?',
                        help='Input FITS MUSE cube to bin.')

    parser.add_argument('-n', '--name', metavar='GALAXY_NAME', default='NAME_OF_GALAXY', type=str,
                        nargs='?', help='Name of your galaxy.')

    parser.add_argument('-b', '--band_continuum', nargs='?', default=(6427, 6632),
                        help='Wavelength range, in Angstroms, of a line-free region of stellar continuum signal. Express as a tuple, e.g. (6427, 6632).')

    parser.add_argument('-s', '--sn', type=int, nargs='?',
                        default=50, help='Signal-to-noise ratio used for Voronoi threshold. Default is 50.')

    parser.add_argument('-m', '--max_area', type=int, nargs='?',
                        default=1000, help='Maximum area, in pixels, for any given Voronoi bin. ')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    start_time = time.time()
    main()
    runtime = round((time.time() - start_time), 3)
    print("\n=================    Aoide | Step 3 Finished   ================\n")
    print("Finished in {} minutes.".format(round(runtime / 60, 3)))
    print("Use the Voronoi-binned cube as part of your Paradise STELLAR fitting.")
