#!/usr/bin/env python
import os

from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy.wcs import WCS
import astropy.constants as const
import astropy.units as u

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm


def make_fovimage(fovimage, name, vmin=0.001, vmax=0.5, cmap=cm.magma, nanthresh=0.0, nancolor='white', save=True, file_save_directory="./"):
    '''Show a collapsed MUSE image.'''

    data = fits.getdata(fovimage)
    wcs = WCS(fits.getheader(fovimage, 1))

    data[data < nanthresh] = np.nan

    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(projection=wcs)

    ax.coords[0].set_axislabel('Right Ascension')
    ax.coords[1].set_axislabel('Declination')

    cmap.set_bad(nancolor)  # We set NaNs to be white here.

    ax.set_title('{} \nCollapsed MUSE Image'.format(name))

    ax.tick_params(labelsize=14)

    frame = ax.imshow(data, origin='lower', # norm=LogNorm(),
                      vmin=vmin, vmax=vmax, cmap=cmap, interpolation='nearest')

    if save is True:
        # Check that the file save directory exists. If not, create it.
        if not os.path.exists(file_save_directory):
            os.makedirs(file_save_directory)
            print("Creating file save directory: {}".format(file_save_directory))
        else:
            print("Found file save directory: {}".format(file_save_directory))

        # Save the PDF figure
        pdf_file = "{}_fovimage.pdf".format(name.replace(" ", ""))
        fig.savefig(file_save_directory + pdf_file,
                    dpi=300, bbox_inches='tight')
        print("Saved FOV image to {}".format(pdf_file))


# Make Stellar Velocity & FWHM Maps

def make_stars(kintable, fovimage, redshift, name, snthresh=100, velocity_thresh=None, project_by_median=True, project_by_redshift=False, cropbox=None, zoom=False, zoombox=None, wcs=True, nancolor='white', colorbar_scalefactor=0.047, vel_vmin=-500, vel_vmax=500, disp_vmin=0, disp_vmax=300, save=True, file_save_directory="./"):

    table = fits.getdata(kintable)
    #hdr = fits.getheader(fovimage)

    x = table['x_cor']
    y = table['y_cor']

    redshift_to_vel = redshift * const.c.to(u.km / u.s).value

    # Get the 2D dimensions into which you'll paint this data
    fovdata = fits.getdata(fovimage)
    dim = fovdata.shape

    # Make empty maps of NaNs
    velmap = np.full((dim[0], dim[1]), np.nan)
    dispmap = np.full((dim[0], dim[1]), np.nan)

    # Threshold & Crop
    if velocity_thresh is not None:
        # mask = (table['vel_fit'] / table['vel_fit_err'] > snthresh) & (table['vel_fit'] < np.nanmedian(table['vel_fit'] + velocity_thresh)) & (table['vel_fit'] > np.nanmedian(table['vel_fit'] - velocity_thresh))
        mask = (table['vel_fit'] / table['vel_fit_err'] > snthresh) & (table['vel_fit'] <
                                                                       redshift_to_vel + velocity_thresh) & (table['vel_fit'] > redshift_to_vel - velocity_thresh)
    else:
        mask = ((table['vel_fit'] / table['vel_fit_err']) > snthresh)
        mask = table['vel_fit'] > 0

    # YES, y,x, in that order. I know it's confusing.
    velmap[y[mask], x[mask]] = table['vel_fit'][mask]
    dispmap[y[mask], x[mask]] = table['disp_fit'][mask]

    if cropbox is not None:
        #x1, x2, y1, y2 = cropbox
        y1, y2, x1, x2 = cropbox

        # YES, you should be confused by the below line.
        # This is *inverted* for the actual mask, but NOT for the zoom.
        # So what the user would naturally expect to be x1 is actually y1
        keep = (x > x1) & (x < x2) & (y > y1) & (y < y2)  # YES
        crop = np.logical_not(keep)

#         velmap[x[crop], y[crop]] = np.nan
        velmap[y[crop], x[crop]] = np.nan

#         dispmap[x[crop], y[crop]] = np.nan
        dispmap[y[crop], x[crop]] = np.nan

    if project_by_redshift is True:
        if project_by_median is True:
            project_by_median = False
            print("project_by_redshift=True which overrides project_by_median=True")
        velmap = velmap - redshift_to_vel
    if project_by_median is True:
        redshift = np.nanmedian(velmap) / const.c.to(u.km / u.s).value
        velmap = velmap - np.nanmedian(velmap)

    # Create the Figures

    if wcs is True:
        wcs = WCS(fits.getheader(fovimage, 1))

        velfig = plt.figure(1, figsize=(10, 10))
        velax = velfig.add_subplot(111, projection=wcs)

        dispfig = plt.figure(2, figsize=(10, 10))
        dispax = dispfig.add_subplot(111, projection=wcs)

        velax.coords[0].set_axislabel('Right Ascension')
        velax.coords[1].set_axislabel('Declination')

        dispax.coords[0].set_axislabel('Right Ascension')
        dispax.coords[1].set_axislabel('Declination')

    elif wcs is False:
        velfig = plt.figure(1, figsize=(10, 10))
        velax = velfig.add_subplot(111)

        dispfig = plt.figure(2, figsize=(10, 10))
        dispax = dispfig.add_subplot(111)

        velax.set_xlabel("X")
        velax.set_ylabel("Y")

        dispax.set_xlabel("X")
        dispax.set_ylabel("Y")

    velax.grid(False)
    dispax.grid(False)

    cmap_vel = cm.RdYlBu_r
    cmap_vel.set_bad(nancolor)

    cmap_disp = cm.plasma
    cmap_disp.set_bad(nancolor)

    if project_by_median is True or project_by_redshift is True:
        velframe = velax.imshow(velmap, origin='lower', vmin=vel_vmin,
                                vmax=vel_vmax, cmap=cmap_vel, interpolation='nearest')
        velcbar = velfig.colorbar(
            velframe, fraction=colorbar_scalefactor, pad=0.01)
        velcbar.set_label(
            r"Stellar Velocity (km s$^{{-1}}$) relative to z = {}".format(round(redshift, 4)))
    else:
        velframe = velax.imshow(velmap, origin='lower',
                                cmap=cmap_vel, interpolation='nearest')
        velcbar = velfig.colorbar(
            velframe, fraction=colorbar_scalefactor, pad=0.01)
        velcbar.set_label(r"Stellar Velocity (km s$^{-1}$)")
        vmin = None
        vmax = None

    dispframe = dispax.imshow(dispmap, origin='lower', vmin=disp_vmin,
                              vmax=disp_vmax, cmap=cmap_disp, interpolation='nearest')
    dispcbar = dispfig.colorbar(
        dispframe, fraction=colorbar_scalefactor, pad=0.01)
    dispcbar.set_label(r"Stellar Velocity Dispersion (km s$^{-1}$)")

    if zoom is True:
        if zoombox is not None:
            x1, x2, y1, y2 = zoombox
            print("Zooming to {}".format(zoombox))
        elif zoombox is None and cropbox is not None:
            print("Using cropbox as zoombox")
        elif zoombox is None and cropbox is None:
            raise Exception(
                "Zoom is TRUE but you don't have a crop or zoombox. You must specify at least one!")
        velax.set_xlim(x1, x2)
        velax.set_ylim(y1, y2)

        dispax.set_xlim(x1, x2)
        dispax.set_ylim(y1, y2)

    velax.grid(False)
    dispax.grid(False)

    # Save Everything

    if save is True:
        # Check that the file save directory exists. If not, create it.
        if not os.path.exists(file_save_directory):
            os.makedirs(file_save_directory)
            print("Creating file save directory: {}".format(file_save_directory))
        else:
            print("Found file save directory: {}".format(file_save_directory))

        # Save the PDF figure
        velfig_pdf_file = "{}_stellar_velocity.pdf".format(name.replace(" ", ""))
        velfig.savefig(file_save_directory + velfig_pdf_file, dpi=300, bbox_inches='tight')
        print("Saved Stellar Velocity Figure to {}".format(velfig_pdf_file))

        dispfig_pdf_file = "{}_stellar_dispersion.pdf".format(
            name.replace(" ", ""))
        dispfig.savefig(file_save_directory + dispfig_pdf_file,
                        dpi=300, bbox_inches='tight')
        print("Saved Stellar Velocity Dispersion Figure to {}".format(dispfig_pdf_file))

        # Save the FITS file
        vel_fits_file = "{}_stellar_velocity.fits".format(
            name.replace(" ", ""))
        hdr = WCS(fits.getheader(fovimage, 1)).to_header()
        hdu = fits.PrimaryHDU(velmap, header=hdr)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(file_save_directory + vel_fits_file,
                        overwrite=True, output_verify='silentfix')
        print("Saved Stellar Velocity Map FITS image to {}".format(vel_fits_file))

        # Save the FITS figure, along with a WCS

        disp_fits_file = "{}_stellar_dispersion.fits".format(
            name.replace(" ", ""))

        disp_fits_file = "{}_stellar_dispersion.fits".format(
            name.replace(" ", ""))
        hdu = fits.PrimaryHDU(dispmap, header=hdr)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(file_save_directory + disp_fits_file,
                        overwrite=True, output_verify='silentfix')
        print("Saved Stellar Velocity Dispersion Map FITS image to {}".format(
            disp_fits_file))

        # Make Kinemetry Table

        kinemetry_table_filename = "{}_kinemetry_table.dat".format(
            name.replace(" ", ""))

        bin_number = np.arange(1, len(table["x_cor"][mask]) + 1)
        kinemetry_table = Table([bin_number,
                                 table["x_cor"][mask],
                                 table["y_cor"][mask],
                                 table["vel_fit"][mask],
                                 table["vel_fit_err"][mask],
                                 table["disp_fit"][mask],
                                 table["disp_fit_err"][mask]],
                                names=["#", "XBIN", "YBIN", "VEL", "ER_VEL", "SIG", "ER_SIG"])
        print("Saved Kinemetry Table to {}".format(kinemetry_table_filename))
        ascii.write(kinemetry_table, file_save_directory +
                    kinemetry_table_filename, overwrite=True)


# Make Emission Line Maps

def make_elines(elinetable, eline, fovimage, redshift, name, snthresh=4, velocity_thresh=None, project_by_median=True, project_by_redshift=False, cropbox=None, zoom=False, zoombox=None, wcs=True, nancolor='white', colorbar_scalefactor=0.047, flux_vmin=0, flux_vmax=1, fluxlog=False, vel_vmin=-500, vel_vmax=500, fwhm_vmin=0, fwhm_vmax=300, save=True, file_save_directory="./"):

    table = fits.getdata(elinetable, 1)

    # Annoyingly, the tables sometimes have different naming for x_cor, y_cor.
    if 'x_cor' in table.names:
        x = table['x_cor']
        y = table['y_cor']
    elif 'x_cor_1' in table.names:
        x = table['x_cor_1']
        y = table['y_cor_1']
    else:
        sys.exit(
            "Something may be wrong with reading the table, cannot find x_cor or x_cor_1")

    # Get the 2D dimensions into which you'll paint this data
    fovdata = fits.getdata(fovimage)
    dim = fovdata.shape

    # Make empty maps of NaNs
    fluxmap = np.full((dim[0], dim[1]), np.nan)
    velmap = np.full((dim[0], dim[1]), np.nan)
    fwhmmap = np.full((dim[0], dim[1]), np.nan)

    mask = table['{}_flux'.format(eline)] / \
        table['{}_flux_err'.format(eline)] > snthresh

    fluxmap[y[mask], x[mask]] = table['{}_flux'.format(eline)][mask]
    velmap[y[mask], x[mask]] = table['{}_vel'.format(eline)][mask]
    fwhmmap[y[mask], x[mask]] = table['{}_fwhm'.format(eline)][mask]

    if cropbox is not None:
        x1, x2, y1, y2 = cropbox

        # YES, you should be confused by the below line.
        # This is *inverted* for the actual mask, but NOT for the zoom.
        # So what the user would naturally expect to be x1 is actually y1
        # keep = (x > y1) & (x < y2) & (y > x1) & (y < x2) # YES
        keep = (x > x1) & (x < x2) & (y > y1) & (y < y2)
        crop = np.logical_not(keep)

        fluxmap[x[crop], y[crop]] = np.nan
        velmap[x[crop], y[crop]] = np.nan
        fwhmmap[x[crop], y[crop]] = np.nan

    if project_by_redshift is True:
        if project_by_median is True:
            project_by_median = False
            print("project_by_redshift=True which overrides project_by_median=True")
        velmap = velmap - redshift * const.c.to(u.km / u.s).value
    if project_by_median is True:
        redshift = np.nanmedian(velmap) / const.c.to(u.km / u.s).value
        velmap = velmap - np.nanmedian(velmap)

    # Create the Figures

    if wcs is True:
        wcs = WCS(fits.getheader(fovimage, 1))

        fluxfig = plt.figure(1, figsize=(10, 10))
        fluxax = fluxfig.add_subplot(111, projection=wcs)

        velfig = plt.figure(2, figsize=(10, 10))
        velax = velfig.add_subplot(111, projection=wcs)

        fwhmfig = plt.figure(3, figsize=(10, 10))
        fwhmax = fwhmfig.add_subplot(111, projection=wcs)

        dispfig = plt.figure(4, figsize=(10, 10))
        dispax = dispfig.add_subplot(111, projection=wcs)

        fluxax.coords[0].set_axislabel('Right Ascension')
        fluxax.coords[1].set_axislabel('Declination')

        velax.coords[0].set_axislabel('Right Ascension')
        velax.coords[1].set_axislabel('Declination')

        fwhmax.coords[0].set_axislabel('Right Ascension')
        fwhmax.coords[1].set_axislabel('Declination')

        dispax.coords[0].set_axislabel('Right Ascension')
        dispax.coords[1].set_axislabel('Declination')

    elif wcs is False:
        fluxfig = plt.figure(1, figsize=(10, 10))
        fluxax = fluxfig.add_subplot(111)

        velfig = plt.figure(2, figsize=(10, 10))
        velax = velfig.add_subplot(111)

        fwhmfig = plt.figure(3, figsize=(10, 10))
        fwhmax = fwhmfig.add_subplot(111)

        dispfig = plt.figure(4, figsize=(10, 10))
        dispax = dispfig.add_subplot(111)

        fluxax.set_xlabel("X")
        fluxax.set_ylabel("Y")

        velax.set_xlabel("X")
        velax.set_ylabel("Y")

        fwhmax.set_xlabel("X")
        fwhmax.set_ylabel("Y")

        dispax.set_xlabel("X")
        dispax.set_ylabel("Y")

    fluxax.grid(False)
    velax.grid(False)
    fwhmax.grid(False)
    dispax.grid(False)

    cmap_flux = cm.magma
    cmap_flux.set_bad(nancolor)

    cmap_vel = cm.RdBu_r
    cmap_vel.set_bad(nancolor)

    cmap_fwhm = sns.cubehelix_palette(20, light=0.95, dark=0.15, as_cmap=True)
    cmap_fwhm.set_bad(nancolor)

    if fluxlog is True:
        # Ensure that vmin is positive (presumably the user won't pick something less than zero, but they might pick zero)
        flux_vmin += 0.0001
        fluxframe = fluxax.imshow(fluxmap, origin='lower', norm=LogNorm(
        ), vmin=flux_vmin, vmax=flux_vmax,  cmap=cmap_flux, interpolation='nearest')
    elif fluxlog is False:
        fluxframe = fluxax.imshow(fluxmap, origin='lower', vmin=flux_vmin,
                                  vmax=flux_vmax,  cmap=cmap_flux, interpolation='nearest')
    fluxcbar = fluxfig.colorbar(
        fluxframe, fraction=colorbar_scalefactor, pad=0.01)
    fluxcbar.set_label(
        r"{} Flux ($\times10^{{-16}}$ erg s$^{{-1}}$ cm$^{{-2}}$)".format(eline))

    if project_by_median is True or project_by_redshift is True:
        velframe = velax.imshow(velmap, origin='lower', vmin=vel_vmin,
                                vmax=vel_vmax, cmap=cmap_vel, interpolation='nearest')
        velcbar = velfig.colorbar(
            velframe, fraction=colorbar_scalefactor, pad=0.01)
        velcbar.set_label(r"{} Velocity (km s$^{{-1}}$)".format(eline))
    else:
        velframe = velax.imshow(velmap, origin='lower',
                                cmap=cmap_vel, interpolation='nearest')
        velcbar = velfig.colorbar(
            velframe, fraction=colorbar_scalefactor, pad=0.01)
        velcbar.set_label(r"{} Velocity (km s$^{{-1}}$)".format(eline))
        vmin = None
        vmax = None

    fwhmframe = fwhmax.imshow(fwhmmap, origin='lower', vmin=fwhm_vmin,
                              vmax=fwhm_vmax, cmap=cmap_fwhm, interpolation='nearest')
    fwhmcbar = fwhmfig.colorbar(
        fwhmframe, fraction=colorbar_scalefactor, pad=0.01)
    fwhmcbar.set_label(r"{} Velocity FWHM (km s$^{{-1}}$)".format(eline))

    dispframe = dispax.imshow(fwhmmap / 2.35, origin='lower', vmin=fwhm_vmin /
                              2.35, vmax=fwhm_vmax / 2.35, cmap=cmap_fwhm, interpolation='nearest')
    dispcbar = dispfig.colorbar(
        dispframe, fraction=colorbar_scalefactor, pad=0.01)
    dispcbar.set_label(r"{} Velocity Dispersion (km s$^{{-1}}$)".format(eline))

    if zoom is True:
        if zoombox is not None:
            x1, x2, y1, y2 = zoombox
            print("Zooming to {}".format(zoombox))
        elif zoombox is None and cropbox is not None:
            print("Using cropbox as zoombox")
        elif zoombox is None and cropbox is None:
            raise Exception(
                "Zoom is TRUE but you don't have a crop or zoombox. You must specify at least one!")
        fluxax.set_xlim(x1, x2)
        fluxax.set_ylim(y1, y2)

        velax.set_xlim(x1, x2)
        velax.set_ylim(y1, y2)

        fwhmax.set_xlim(x1, x2)
        fwhmax.set_ylim(y1, y2)

        dispax.set_xlim(x1, x2)
        dispax.set_ylim(y1, y2)

    fluxax.grid(False)
    velax.grid(False)
    fwhmax.grid(False)
    dispax.grid(False)

    # Save Everything

    if save is True:
        # Check that the file save directory exists. If not, create it.
        if not os.path.exists(file_save_directory):
            os.makedirs(file_save_directory)
            print("Creating file save directory: {}".format(file_save_directory))
        else:
            print("Found file save directory: {}".format(file_save_directory))

        # Save the PDF figures

        fluxfig_pdf_file = "{}_{}_flux.pdf".format(
            name.replace(" ", ""), eline)
        velfig_pdf_file = "{}_{}_vel.pdf".format(name.replace(" ", ""), eline)
        fwhmfig_pdf_file = "{}_{}_fwhm.pdf".format(
            name.replace(" ", ""), eline)
        dispfig_pdf_file = "{}_{}_disp.pdf".format(
            name.replace(" ", ""), eline)

        fluxfig.savefig(file_save_directory + fluxfig_pdf_file,
                        dpi=300, bbox_inches='tight')
        velfig.savefig(file_save_directory + velfig_pdf_file,
                       dpi=300, bbox_inches='tight')
        fwhmfig.savefig(file_save_directory + fwhmfig_pdf_file,
                        dpi=300, bbox_inches='tight')
        dispfig.savefig(file_save_directory + dispfig_pdf_file,
                        dpi=300, bbox_inches='tight')
        print("Saved {} flux map to {}".format(eline, fluxfig_pdf_file))
        print("Saved {} velocity map to {}".format(eline, velfig_pdf_file))
        print("Saved {} velocity FWHM & Dispersion maps to {} & {}".format(
            eline, fwhmfig_pdf_file, dispfig_pdf_file))

        # Save the FITS files

        flux_fits_file = "{}_{}_flux.fits".format(name.replace(" ", ""), eline)
        vel_fits_file = "{}_{}_vel.fits".format(name.replace(" ", ""), eline)
        fwhm_fits_file = "{}_{}_fwhm.fits".format(name.replace(" ", ""), eline)
        disp_fits_file = "{}_{}_disp.fits".format(name.replace(" ", ""), eline)

        hdr = WCS(fits.getheader(fovimage, 1)).to_header()

        fluxhdu = fits.PrimaryHDU(fluxmap, header=hdr)
        velhdu = fits.PrimaryHDU(velmap, header=hdr)
        fwhmhdu = fits.PrimaryHDU(fwhmmap, header=hdr)
        disphdu = fits.PrimaryHDU(fwhmmap / 2.35, header=hdr)

        fluxhdu.writeto(file_save_directory + flux_fits_file,
                        overwrite=True, output_verify='silentfix')
        velhdu.writeto(file_save_directory + vel_fits_file,
                       overwrite=True, output_verify='silentfix')
        fwhmhdu.writeto(file_save_directory + fwhm_fits_file,
                        overwrite=True, output_verify='silentfix')
        disphdu.writeto(file_save_directory + disp_fits_file,
                        overwrite=True, output_verify='silentfix')

        print("Saved {} flux map to {}".format(eline, flux_fits_file))
        print("Saved {} flux velocity to {}".format(eline, vel_fits_file))
        print("Saved {} velocity FWHM and dispersion maps to {} & {}".format(
            eline, fwhm_fits_file, disp_fits_file))


def make_electron_density_map(sii_6716_image, sii_6730_image, hdr):

    ratio = sii_6716_image / sii_6730_image

    # assuming T = 10^4 K, following eq. (3) here: https://arxiv.org/pdf/1311.5041.pdf

    log_ne_per_cm3 = 0.053 * np.tan(-3.0553 * ratio + 2.8506) + \
        6.98 - 10.6905 * ratio + \
        9.9186 * ratio**2 - 3.5442 * ratio**3

    # Mask unrealistic values
    log_ne_mask1 = log_ne_per_cm3 < 0
    log_ne_mask2 = log_ne_per_cm3 > 6

    log_ne_per_cm3[log_ne_mask1] = np.nan
    log_ne_per_cm3[log_ne_mask2] = np.nan

    hdu = fits.PrimaryHDU(log_ne_per_cm3, header=hdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto('ne.fits', overwrite=True, output_verify='silentfix')


def make_balmer_map(ha_image, hb_image, hdr):

    ratio_observed = ha_image / hb_image
    ratio_intrinsic = 2.86
    k_alpha = 2.63
    k_beta = 3.71

    ebv = (2.5 / (k_beta - k_alpha)) * \
        np.log10(ratio_observed / ratio_intrinsic)
    ebv[ebv < 0] = np.nan  # Additional masking

    av = 4.05 * ebv
    nh = 1.8e21 * av  # VERY rough, from Predehl & Schmitt, in atoms / cm2

    hdu = fits.PrimaryHDU(ebv, header=hdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto('ebv.fits', overwrite=True, output_verify='silentfix')

    hdu = fits.PrimaryHDU(av, header=hdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto('av.fits', overwrite=True, output_verify='silentfix')

    hdu = fits.PrimaryHDU(nh, header=hdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto('nh.fits', overwrite=True, output_verify='silentfix')

    return None
