from astropy.io import fits
import numpy as np

import matplotlib.pyplot as plt

import astropy.constants as const


# Read in the WCS from the data cube,
# so that we can directly paste it into the
# maps we're making

########## EDIT STUFF HERE ################


z = 0.0821  # Using the Ca II H+K line Stellar Systemic Redshift
cz = z * const.c.to('km/s')

fov = 'IMAGE_FOV_0001.fits'
inmask = 'Ha.fits'
stars = 'A2597.kin_table.fits'
gas = 'A2597.eline_table.fits'

#############################################



#### MAKE FOV

hdu = fits.open(fov)
hdr = hdu[1].header
dim = hdu[1].data.shape
hdu.close()

#### MAKE MASK

maskdata = fits.open(inmask)
thresh = maskdata[0].data

sourcemask = (thresh < 0.001)

#### MAKE STARS

# Read in the stellar table, from which we'll make stellar maps
stellar_hdu = fits.open(stars)
stellar_tab = stellar_hdu[1].data

# Read in the emission line table, from which we'll make gas maps
eline_hdu = fits.open(gas)
eline_tab = eline_hdu[1].data


stellar_columns = stellar_hdu[1].header
eline_columns = eline_hdu[1].header

stellar_x_cor = stellar_tab.field('x_cor')
stellar_y_cor = stellar_tab.field('y_cor')

# Populate the Stellar Velocity and FWHM
stellar_vel = stellar_tab.field('vel_fit')
stellar_vel_err = stellar_tab.field('vel_fit_err')

stellar_fwhm = stellar_tab.field('disp_fit')
stellar_fwhm_err = stellar_tab.field('disp_fit_err')


stellar_select = (stellar_fwhm_err < 40.) & (stellar_vel - cz.value > -200. ) & (stellar_vel - cz.value < 150.)

median_stellar_vel = np.median(stellar_vel[stellar_select])

stellar_vel_map = np.full((dim[0], dim[1]), np.nan)
stellar_vel_map[stellar_x_cor[stellar_select], stellar_y_cor[stellar_select]] = stellar_vel[stellar_select] - np.median(stellar_vel[stellar_select])



hdu = fits.PrimaryHDU(stellar_vel_map)
hdu.header = hdr
hdu.writeto('stellar_vel_map.fits',clobber=True,output_verify='fix')


stellar_fwhm_map = np.full((dim[0], dim[1]), np.nan)
stellar_fwhm_map[stellar_x_cor[stellar_select], stellar_y_cor[stellar_select]] = stellar_fwhm[stellar_select]
hdu = fits.PrimaryHDU(stellar_fwhm_map)
hdu.header = hdr
hdu.writeto('stellar_fwhm_map.fits',clobber=True,output_verify='fix')


##### MAKE GAS

eline_x_cor = eline_tab.field('x_cor')
eline_y_cor = eline_tab.field('y_cor')

Ha_flux = eline_tab.field('Halpha_flux')
Ha_flux_err = eline_tab.field('Halpha_flux_err')
Ha_vel = eline_tab.field('Halpha_vel')
Ha_vel_err = eline_tab.field('Halpha_vel_err')
Ha_fwhm = eline_tab.field('Halpha_fwhm')
Ha_fwhm_err = eline_tab.field('Halpha_fwhm_err')


Hb_flux = eline_tab.field('Hbeta_flux')
Hb_flux_err = eline_tab.field('Hbeta_flux_err')
Hb_vel = eline_tab.field('Hbeta_vel')
Hb_vel_err = eline_tab.field('Hbeta_vel_err')
Hb_fwhm = eline_tab.field('Hbeta_fwhm')
Hb_fwhm_err = eline_tab.field('Hbeta_fwhm_err')


OIII5007_flux = eline_tab.field('OIII5007_flux')
OIII5007_vel = eline_tab.field('OIII5007_vel')
OIII5007_fwhm = eline_tab.field('OIII5007_fwhm')

OI6300_flux = eline_tab.field('OI6300_flux')
OI6300_vel = eline_tab.field('OI6300_vel')
OI6300_fwhm = eline_tab.field('OI6300_fwhm')

SII6717_flux = eline_tab.field('SII6717_flux')
SII6717_flux_err = eline_tab.field('SII6717_flux')

SII6730_flux = eline_tab.field('SII6730_flux')
SII6730_flux_err = eline_tab.field('SII6730_flux')

NII6583_flux = eline_tab.field('NII6583_flux')
NII6583_flux_err = eline_tab.field('NII6583_flux_err')

gas_select = (Ha_fwhm < 1000.0)

Ha_flux_map = np.full((dim[0],dim[1]) ,np.nan)
Ha_vel_map = np.full((dim[0],dim[1]), np.nan)
Ha_fwhm_map = np.full((dim[0],dim[1]), np.nan)

Ha_flux_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = Ha_flux[gas_select]
Ha_vel_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = Ha_vel[gas_select] - cz.value
Ha_fwhm_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = Ha_fwhm[gas_select]

Ha_flux_map[sourcemask] = np.nan
Ha_vel_map[sourcemask] = np.nan
Ha_fwhm_map[sourcemask] = np.nan


hdu = fits.PrimaryHDU(Ha_flux_map)
hdu.header = hdr
hdu.writeto('Ha_flux_map.fits',clobber=True,output_verify='silentfix')

hdu = fits.PrimaryHDU(Ha_fwhm_map)
hdu.header = hdr
hdu.writeto('Ha_fwhm_map.fits',clobber=True,output_verify='silentfix')

hdu = fits.PrimaryHDU(Ha_vel_map)
hdu.header = hdr
hdu.writeto('Ha_vel_map.fits',clobber=True,output_verify='silentfix')


Hb_flux_map = np.zeros((dim[0],dim[1]))
Hb_vel_map = np.zeros((dim[0],dim[1]))
Hb_fwhm_map = np.zeros((dim[0],dim[1]))

Hb_flux_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = Hb_flux[gas_select]
Hb_vel_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = Hb_vel[gas_select] - cz.value
Hb_fwhm_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = Hb_fwhm[gas_select]

Hb_flux_map[sourcemask] = np.nan
Hb_vel_map[sourcemask] = np.nan
Hb_fwhm_map[sourcemask] = np.nan


hdu = fits.PrimaryHDU(Hb_flux_map)
hdu.header = hdr
hdu.writeto('Hb_flux_map.fits',clobber=True,output_verify='silentfix')

hdu = fits.PrimaryHDU(Hb_fwhm_map)
hdu.header = hdr
hdu.writeto('Hb_fwhm_map.fits',clobber=True,output_verify='silentfix')

hdu = fits.PrimaryHDU(Hb_vel_map)
hdu.header = hdr
hdu.writeto('Hb_vel_map.fits',clobber=True,output_verify='silentfix')



OIII5007_flux_map = np.zeros((dim[0],dim[1]))
OIII5007_vel_map = np.zeros((dim[0],dim[1]))
OIII5007_fwhm_map = np.zeros((dim[0],dim[1]))

OIII5007_flux_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = OIII5007_flux[gas_select]
OIII5007_vel_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = OIII5007_vel[gas_select] - cz.value
OIII5007_fwhm_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = OIII5007_fwhm[gas_select]

OIII5007_flux_map[sourcemask] = np.nan
OIII5007_vel_map[sourcemask] = np.nan
OIII5007_fwhm_map[sourcemask] = np.nan

hdu = fits.PrimaryHDU(OIII5007_flux_map)
hdu.header = hdr
hdu.writeto('OIII5007_flux_map.fits',clobber=True,output_verify='silentfix')

hdu = fits.PrimaryHDU(OIII5007_fwhm_map)
hdu.header = hdr
hdu.writeto('OIII5007_fwhm_map.fits',clobber=True,output_verify='silentfix')

hdu = fits.PrimaryHDU(OIII5007_vel_map)
hdu.header = hdr
hdu.writeto('OIII5007_vel_map.fits',clobber=True,output_verify='silentfix')



OI6300_flux_map = np.zeros((dim[0],dim[1]))
OI6300_vel_map = np.zeros((dim[0],dim[1]))
OI6300_fwhm_map = np.zeros((dim[0],dim[1]))

OI6300_flux_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = OI6300_flux[gas_select]
OI6300_vel_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = OI6300_vel[gas_select] - cz.value
OI6300_fwhm_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = OI6300_fwhm[gas_select]

OI6300_flux_map[sourcemask] = np.nan
OI6300_vel_map[sourcemask] = np.nan
OI6300_fwhm_map[sourcemask] = np.nan

hdu = fits.PrimaryHDU(OI6300_flux_map)
hdu.header = hdr
hdu.writeto('OI6300_flux_map.fits',clobber=True,output_verify='silentfix')

hdu = fits.PrimaryHDU(OI6300_fwhm_map)
hdu.header = hdr
hdu.writeto('OI6300_fwhm_map.fits',clobber=True,output_verify='silentfix')

hdu = fits.PrimaryHDU(OI6300_vel_map)
hdu.header = hdr
hdu.writeto('OI6300_vel_map.fits',clobber=True,output_verify='silentfix')


SII_select = (SII6730_flux / SII6730_flux_err) > 3

SII6717_flux_map = np.zeros((dim[0], dim[1]))
SII6717_flux_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = SII6717_flux[gas_select]


SII6730_flux_map = np.zeros((dim[0], dim[1]))
SII6730_flux_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = SII6730_flux[gas_select]


SII6717_flux_map[sourcemask] = np.nan
SII6730_flux_map[sourcemask] = np.nan



hdu = fits.PrimaryHDU(SII6730_flux_map)
hdu.header = hdr
hdu.writeto('SII6730_flux_map.fits',clobber=True,output_verify='silentfix')


SII_select = (NII6583_flux / NII6583_flux_err) > 3

SII6717_flux_map = np.zeros((dim[0], dim[1]))
SII6717_flux_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = SII6717_flux[gas_select]


NII6583_flux_map = np.zeros((dim[0], dim[1]))
NII6583_flux_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = NII6583_flux[gas_select]


SII6717_flux_map[sourcemask] = np.nan
NII6583_flux_map[sourcemask] = np.nan



hdu = fits.PrimaryHDU(NII6583_flux_map)
hdu.header = hdr
hdu.writeto('NII6583_flux_map.fits',clobber=True,output_verify='silentfix')



# hdu = fits.PrimaryHDU(SII6730_fwhm_map)
# hdu.header = hdr
# hdu.writeto('SII6730_fwhm_map.fits',clobber=True,output_verify='silentfix')

# hdu = fits.PrimaryHDU(SII6730_vel_map)
# hdu.header = hdr
# hdu.writeto('SII6730_vel_map.fits',clobber=True,output_verify='silentfix')




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

    ebv = (2.5 / (k_beta - k_alpha)) * np.log10(ratio_observed / ratio_intrinsic)
    ebv[ebv < 0] = np.nan # Additional masking

    av = 4.05 * ebv
    nh = 1.8e21 * av # VERY rough, from Predehl & Schmitt, in atoms / cm2


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


make_balmer_map(Ha_flux_map, Hb_flux_map, hdr)
make_electron_density_map(SII6717_flux_map, SII6730_flux_map, hdr)


