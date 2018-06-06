#!/usr/bin/env python

import os

from astropy.io import fits
import numpy as np


def main():

    args = parse_args()

    name = args.name

    hdu = fits.open('{}_rss.stellar_table.fits'.format(name))
    tab = hdu[1].data
    fiber = tab.field('fiber')
    vel_fit = tab.field('vel_fit')
    vel_fit_err = tab.field('vel_fit_err')
    disp_fit = tab.field('disp_fit')
    disp_fit_err = tab.field('disp_fit_err')

    create_kin_tab('%s.voronoi.pixel.dat'.format(name), name + '_rss.kin_table.fits', fiber, vel_fit, vel_fit_err, disp_fit, disp_fit_err)


def create_kin_tab(voronoi_pixel,out_table,fiber,vel_fit,vel_fit_err,disp_fit,disp_fit_err):

  print("Creating STELLAR_KIN Table")
  ascii_in = open(voronoi_pixel,'r')
  lines = ascii_in.readlines()
  x = np.arange(len(lines)-1,dtype=np.int16)
  y = np.arange(len(lines)-1,dtype=np.int16)
  binNr = np.arange(len(lines)-1,dtype=np.int16)
  for i in range(1,len(lines)):
    line = lines[i].split()
    x[i-1]=int(line[0])
    y[i-1]=int(line[1])
    binNr[i-1]=int(line[2])
  x_cor = x
  y_cor = y
  vel_out = np.zeros(len(x),dtype=np.float32)
  vel_err_out = np.zeros(len(x),dtype=np.float32)
  disp_out = np.zeros(len(x),dtype=np.float32)
  disp_err_out = np.zeros(len(x),dtype=np.float32)
  for i in range(max(fiber)): # THIS USED TO SAY range(max(fiber))
     select = fiber[i]+1 == binNr
     vel_out[select] = vel_fit[i]
     vel_err_out[select] = vel_fit_err[i]
     disp_out[select] = disp_fit[i]
     disp_err_out[select] = disp_fit_err[i]

  columns = []
  columns.append(pyfits.Column(name='x_cor',format='J', array=x_cor))
  columns.append(pyfits.Column(name='y_cor',format='J', array=y_cor))
  columns.append(pyfits.Column(name='vel_fit',format='E', array=vel_out))
  columns.append(pyfits.Column(name='vel_fit_err',format='E', array=vel_err_out))
  columns.append(pyfits.Column(name='disp_fit',format='E', array=disp_out))
  columns.append(pyfits.Column(name='disp_fit_err',format='E', array=disp_err_out))
  new_table = pyfits.new_table(columns)
  new_table.writeto(out_table,clobber=True)


def parse_args():

    parser = argparse.ArgumentParser(
        description="PCA sky subtraction, galactic extinction correction, and preparation for MUSE cubes that are to be analyzed with Paradise")


    parser.add_argument('name', metavar='NAME', type=str, nargs='?',
                        help='Name of galaxy (used for Prefix).')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    start_time = time.time()
    main()
    runtime = round((time.time() - start_time), 3)
    print("\n=================    Aoide | Stellar Kin Table Created   ================\n")
    print("Finished in {} minutes.".format(round(runtime / 60, 3)))
    print("Use this as a partial input to your Emission Line Fitting.")
