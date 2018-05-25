#!/usr/bin/env python

'''
Postprocess a MUSE cube | G. Tremblay
'''

import time

import argparse

from aoide import aoide_postprocess
from aoide import make_sky_mask


def parse_args():

	parser = argparse.ArgumentParser(description="Script to subtract sky residuals from a datacube based on a PCA spectral library")
	parser.add_argument('input_cube',metavar='CUBE_IN',type=str,nargs='?',help='Input FITS datacube from which sky spectra are selected')
	parser.add_argument('output_cube',metavar='CUBE_OUT',type=str,nargs='?',help='Output FITS datacube with sky residuals subtracted')
	parser.add_argument('pca_sky',metavar='PCA',type=str,nargs='?',help='Input FITS file with PCA components')
	parser.add_argument('-c','--components',type=int,nargs='?',default=150, help='Number of PCA components to be used')
	parser.add_argument('-f','--filter',type=int,nargs='?',default=30, help='Size of median filter in wavelength direction to remove continuum signal before sky residual subtraction')
	parser.add_argument('-p','--parallel',type=str,nargs='?',default='1', help='Number of cores used for computation.')

	args = parser.parse_args()

	return args


def main():

	args = parse_args()

if __name__ == '__main__':

    start_time = time.time()
    main()
    runtime = round((time.time() - start_time), 3)
    print("Finished in {} minutes.".format(round(runtime/60, 3)))
