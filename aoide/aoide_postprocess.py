
from __future__ import division, print_function

import os
try:
    import pyfits
except:
    from astropy.io import fits as pyfits
import numpy
from scipy import ndimage
from scipy import interpolate
import multiprocessing
from multiprocessing import Pool
from multiprocessing import cpu_count
try:
    import copy_reg
except:
    import copyreg as copy_reg
from types import *


def subtract_sky(PIXTAB_in, SKYSPEC, PIXTAB_out, parallel=1):

    hdu_PIX = pyfits.open(PIXTAB_in)
    try:
        RVcorr = float(hdu_PIX[0].header['ESO DRS MUSE PIXTABLE RVCORR'])
    except:
        RVcorr = None

    wave = hdu_PIX['lambda'].data[:, 0]
    hdu_SKY = pyfits.open(SKYSPEC)
    tab = hdu_SKY[1].data
    hdu_SKY.close()
    wave_sky = tab.field('lambda')
    spec_sky = tab.field('data')
    var_sky = tab.field('stat')

    if RVcorr is not None:
        z = 1 + (RVcorr / 300000.0)
    else:
        z = 1
    out_sky = interpolate.interp1d(
        wave_sky * z, spec_sky, bounds_error=False, fill_value=spec_sky[0])
    out_var = interpolate.interp1d(
        wave_sky * z, var_sky, bounds_error=False, fill_value=var_sky[0])
    if parallel == 'auto':
        cpus = cpu_count()
    else:
        cpus = int(parallel)
    if cpus > 1:
        pool = Pool(processes=cpus)
        results = []
        wave_split = numpy.array_split(wave)

        def thread(out_sky, out_var, wave):
            return [out_sky(wave), out_var(wave)]

        for i in range(cpus):
            results.append(pool.apply_async(
                thread, args=(out_sky, out_var, wave_split[i])))
        pool.close()
        pool.join()
        sky_out = []
        var_out = []
        for i in range(cpus):
            sky_out.append(results[i].get()[0])
            var_out.append(results[i].get()[1])
        PIX_sky = numpy.concatenate(sky_out)
        stat_sky = numpy.concatenate(var_out)

    else:
        PIX_sky = out_sky(wave)
        stat_sky = out_var(wave)
    hdu_PIX['data'].data[:, 0] = hdu_PIX['data'].data[:, 0] - PIX_sky
    hdu_PIX['stat'].data[:, 0] = hdu_PIX['stat'].data[:, 0] + stat_sky
    hdu_PIX['xpos'].header.update('EXTNAME', 'xpos')
    hdu_PIX['ypos'].header.update('EXTNAME', 'ypos')
    hdu_PIX['lambda'].header.update('EXTNAME', 'lambda')
    hdu_PIX['data'].header.update('EXTNAME', 'data')
    hdu_PIX['dq'].header.update('EXTNAME', 'dq')
    hdu_PIX['stat'].header.update('EXTNAME', 'stat')
    hdu_PIX['origin'].header.update('EXTNAME', 'origin')
    hdu_PIX.writeto(PIXTAB_out, overwrite=True)


def interpolate_sky(SKYSPEC_in, MJD_in, SKYSPEC_out, MJD_out, header_out='', order=1, sim_err=100):
    for i in range(len(SKYSPEC_in)):
        hdu = pyfits.open(SKYSPEC_in[i])
        tab = hdu[1].data
        wave = tab.field('lambda')
        spec = tab.field('data')
        var = tab.field('stat')
        if i == 0:
            spec_in = numpy.zeros((len(SKYSPEC_in), len(spec)))
            err_in = numpy.zeros((len(SKYSPEC_in), len(spec)))
        spec_in[i, :] = spec[:]
        err_in[i, :] = numpy.sqrt(var[:])
    out_spec = numpy.zeros(len(spec))
    out_var = numpy.zeros(len(var))

    MJD = MJD_in - MJD_out
    for i in range(len(spec)):
        fit = numpy.polyfit(MJD, spec_in[:, i], order)
        out_spec[i] = fit[order]
        error_spec = numpy.zeros(sim_err)
        for s in range(sim_err):
            fit_sim = numpy.polyfit(MJD, numpy.random.normal(
                spec_in[:, i], err_in[:, i]), order)
            error_spec[s] = fit_sim[order]
        out_var[i] = numpy.std(error_spec)**2

    columns = []
    columns.append(pyfits.Column(
        name='lambda', format='D', unit='Angstrom', array=wave))
    columns.append(pyfits.Column(name='data', format='D',
                                 unit='10**(-20)*erg/s/cm**2/Angstrom', array=out_spec))
    columns.append(pyfits.Column(name='stat', format='D',
                                 unit='(10**(-20)*erg/s/cm**2/Angstrom)**2', array=out_var))
    columns.append(pyfits.Column(name='dq', format='J',
                                 array=numpy.zeros(len(out_spec))))
    hdu_tab = pyfits.new_table(columns)
    hdu_prim = pyfits.PrimaryHDU()
    hdu = pyfits.HDUList([hdu_prim, hdu_tab])
    hdu.writeto(SKYSPEC_out, overwrite=True)


def create_PCA_sky(cube_in, PCA_out, sky_mask, cont_filt=50, spectra=20000, parallel='auto'):

    hdu = pyfits.open(cube_in)
    cube = hdu[1].data
    dim = cube.shape
    hdu.close()

    hdu = pyfits.open(sky_mask)
    mask = hdu[0].data

    nan = numpy.sum(numpy.isnan(cube[5:-5, :, :]), 0)
    cube[numpy.isnan(cube)] = 0.0
    sky = cube[:, (mask == 1) & (nan == 0)].T
    indices = numpy.arange(sky.shape[0])
    numpy.random.shuffle(indices)
    smax = numpy.min([spectra, sky.shape[0]])
    out = (((sky) / 10000.0).astype(numpy.float64))[indices[:smax], :]
    if parallel == 'auto':
        cpus = cpu_count()
    else:
        cpus = int(parallel)
    if cpus > 1:
        pool = Pool(processes=cpus)
        results = []
        for i in range(smax):
            results.append(pool.apply_async(
                ndimage.filters.median_filter, args=(out[i, :], cont_filt)))
        pool.close()
        pool.join()
        smooth_out = numpy.zeros_like(out)
        for c in range(smax):
            spec = results[c].get()
            smooth_out[c, :] = spec
    else:
        smooth_out = ndimage.filters.median_filter(out, (1, cont_filt))
    out = out - smooth_out
    #hdu = pyfits.PrimaryHDU(out)
    # hdu.writeto('test.fits',clobber=True)
    M = numpy.dot(out, out.T).T
    e, EV = numpy.linalg.eigh(M)
    tmp = numpy.dot(out.T, EV).T
    V = tmp[::-1]
    # S=numpy.sqrt(e)[::-1]

    hdu = pyfits.PrimaryHDU(V[:dim[0], :])
    hdu.writeto(PCA_out, overwrite=True)


def fit_PCA(cube, idx_x, idx_y, pca_specs, cont_filt, select):
    dim = cube.shape
    clean_cube = numpy.zeros((numpy.sum(select), dim[1]))
    for m in range(len(idx_x)):
        spec = cube[:, m]
        if numpy.sum(numpy.isnan(spec)) < 10:
            spec[numpy.isnan(spec)] = 0
            smooth_spec = ndimage.filters.median_filter(spec, (cont_filt))
            out = numpy.linalg.lstsq(
                pca_specs[:, select].T, (spec - smooth_spec)[select])
            spec_sky = numpy.dot(pca_specs[:, select].T, out[0])
        else:
            spec_sky = 0.0
        clean_cube[:, m] = spec[select] - spec_sky
        # if m>20000:
        #    break
    return clean_cube


def subtract_PCA_sky(cube_in, cube_out, PCA_spec, components=150, cont_filt=40, parallel='auto'):

    hdu = pyfits.open(PCA_spec)
    pca_specs = hdu[0].data[:components, :] * 10000.0
    hdu.close()
    hdu = pyfits.open(cube_in)
    cube = hdu[1].data
    hdr = hdu[1].header
    dim = cube.shape
    wave = numpy.arange(dim[0]) * hdr['CD3_3'] + hdr['CRVAL3']
    select = ((wave > 5565) & (wave < 5585)) | ((wave > 5875) & (wave < 5905)) | ((wave > 6290) & (wave < 6310)) | ((wave > 6350) & (wave < 6370)) | ((wave > 6450) & (
        wave < 6610)) | ((wave > 6800) & (wave < 7020)) | ((wave > 7220) & (wave < 8160)) | ((wave > 8250) & (wave < 8680)) | ((wave > 8740) & (wave < 9180))

    clean_cube = numpy.zeros((numpy.sum(select), dim[1], dim[2]))
    clean_cube[:, :, :] = cube[select, :, :]
    if parallel == 'auto':
        cpus = cpu_count()
    else:
        cpus = int(parallel)
    idx = numpy.indices((dim[1], dim[2]))
    idx_x = idx[1].flatten()
    idx_y = idx[0].flatten()
    if cpus > 1:
        pool = Pool(processes=cpus)
        idx_y_split = numpy.array_split(idx_y, cpus)
        idx_x_split = numpy.array_split(idx_x, cpus)
        results = []
        for c in range(cpus):
            results.append(pool.apply_async(fit_PCA, args=(
                cube[:, idx_y_split[c], idx_x_split[c]], idx_x_split[c], idx_y_split[c], pca_specs, cont_filt, select)))
        pool.close()
        pool.join()
        for c in range(cpus):
            out = results[c].get()
            clean_cube[:, idx_y_split[c], idx_x_split[c]] = out[:, :]
    else:

        clean_rss = fit_PCA(cube[:, idx_y, idx_x], idx_y,
                            idx_x, pca_specs, cont_filt, select)
        clean_cube[:, idx_y, idx_x] = clean_rss[:, :]

    cube[select, :, :] = clean_cube
    hdu.writeto(cube_out, overwrite=True)


def gal_extinct(wave, A_V, Rv=3.1):
    m = wave / 10000.0
    x = 1 / m
    y = (x - 1.82)
    ax = 1 + (0.17699 * y) - (0.50447 * y**2) - (0.02427 * y**3) + \
        (0.72085 * y**4) + (0.01979 * y**5) - \
        (0.77530 * y**6) + (0.32999 * y**7)
    bx = (1.41338 * y) + (2.28305 * y**2) + (1.07233 * y**3) - (5.38434 *
                                                                y**4) - (0.62251 * y**5) + (5.30260 * y**6) - (2.09002 * y**7)

    Arat = ax + (bx / Rv)
    Alambda = Arat * A_V
    return Alambda


def correct_cube(incube, outcube, A_V, star_position=[], offset=(0, 0), binning=True):
    hdu = pyfits.open(incube)
    data = hdu[1].data / 1e4
    variance = hdu[2].data / (1e4**2)
    header = hdu[1].header
    hdu.close()

    dim = data.shape
    wave = numpy.arange(dim[0]) * header['CD3_3'] + header['CRVAL3']
    if binning:
        dim_out = (dim[0], (dim[1] - offset[1]) //
                   2, (dim[2] - offset[0]) // 2)
        cube = numpy.zeros(dim_out, dtype=numpy.float64)
        error = numpy.zeros(dim_out, dtype=numpy.float64)
        for x in range(dim_out[2]):
            for y in range(dim_out[1]):
                cube[:, y, x] = numpy.sum(numpy.sum(data[:, 2 * y + offset[1]:2 * (
                    y + 1) + 1 + offset[1], 2 * x + offset[0]:2 * (x + 1) + offset[0]], axis=1), axis=1)

                error[:, y, x] = numpy.sqrt(numpy.sum(numpy.sum(variance[:, 2 * y + offset[1]:2 * (
                    y + 1) + 1 + offset[1], 2 * x + offset[0]:2 * (x + 1) + offset[0]], axis=1), axis=1))
    else:
        cube = data
        error = numpy.sqrt(variance)

    nan = numpy.isnan(cube)
    cube[nan] = 0.0

    if len(star_position) > 0:
        if binning:
            star_spec = numpy.sum(numpy.sum(
                cube[:, star_position[1] // 2 - 2:star_position[1] // 2 + 2, star_position[0] // 2 - 2:star_position[0] // 2 + 2], 1), 1)
        else:
            star_spec = numpy.sum(numpy.sum(
                cube[:, star_position[1] - 4:star_position[1] + 4, star_position[0] - 4:star_position[0] + 4], 1), 1)
        telluric = [[6800, 6950], [7130, 7220], [
            7550, 7700], [8150, 8300], [9000, 9280]]
        clean_star = numpy.array(star_spec).astype(numpy.float64)
        for i in range(len(telluric)):
            blue = telluric[i][0]
            red = telluric[i][1]
            select_blue = (wave > blue - 10) & (wave < blue + 10)
            select_red = (wave > red - 10) & (wave < red + 10)
            select_telluric = (wave > blue) & (wave < red)
            flux_blue = numpy.median(star_spec[select_blue])
            flux_red = numpy.median(star_spec[select_red])
            m = (flux_red - flux_blue) / (red - blue)
            z = flux_blue - m * blue
            clean_star[select_telluric] = wave[select_telluric] * m + z
        ratio = clean_star / star_spec
        cube = cube * ratio[:, numpy.newaxis, numpy.newaxis]
        error = error * ratio[:, numpy.newaxis, numpy.newaxis]

    cube = cube / 10**(gal_extinct(wave, A_V)
                       [:, numpy.newaxis, numpy.newaxis] / -2.5)
    error = error / 10**(gal_extinct(wave, A_V)
                         [:, numpy.newaxis, numpy.newaxis] / -2.5)
    hdus = []
    hdus.append(pyfits.PrimaryHDU(cube.astype(numpy.float32)))
    hdus.append(pyfits.ImageHDU(error.astype(numpy.float32), name='ERROR'))
    hdu = pyfits.HDUList(hdus)

    if binning:
        header['CRPIX1'] = (float(header['CRPIX1']) - offset[0]) / 2.0
        header['CRPIX2'] = (float(header['CRPIX2']) - offset[1]) / 2.0
        header['CD1_1'] = float(header['CD1_1']) * 2.0
        header['CD2_2'] = float(header['CD2_2']) * 2.0

    header['CDELT3'] = header['CD3_3']
    header['BUNIT'] = '10**(-16)*erg/s/cm**2/Angstrom'
    header['EXTINCT'] = A_V
    hdu[0].header = header
    hdu.verify('silentfix')
    hdu.writeto(outcube, overwrite=True)
