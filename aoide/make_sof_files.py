#!/usr/bin/env python

from __future__ import division, print_function
import os
import glob

from astropy.io import fits


def make_sof_files(raw_data_dir, reduction_dir, static_cal_dir, science_targname):

    # Work in the reduction directory. We assume it's already createdself.

    os.chdir(reduction_dir)

    rawfiles = glob.glob(raw_data_dir + "/*.fits.fz")

    # Initialize empty lists
    science_files = []  # Science frames (i.e. OBJECT exposures go here).
    # Bias frames. We really want to only use 10 of these, often there will be more.
    bias_files = []
    dark_files = []  # Dark frames
    flat_files = []  # Lamp flats
    arc_files = []  # Wavelength calibration files (arc lamps)
    twilight_files = []  # Sky flats
    std_files = []  # STD observations

    print("I'm told your science target is {}. Make sure this is correct.".format(
        science_targname))
    print("Sorting raw data files.")

    for file in sorted(rawfiles):

        tags = fits.getval(file, "ESO DPR TYPE")
        print(file.split('/')[-1] + " is " + tags)

        if tags == 'OBJECT':
            if fits.getval(file, "OBJECT") == science_targname:
                print("Found science observation of {}, appending to science_files.".format(
                    fits.getval(file, "OBJECT")))
                science_files.append(file)
            else:
                print("Found science observation of {}, inconsistent with stated target name of {}. Skipping.".format(
                    fits.getval(file, "OBJECT"), science_targname))

        if tags == 'BIAS':
            bias_files.append(file)

        if tags == 'DARK':
            dark_files.append(file)

        if tags == 'FLAT,LAMP':
            flat_files.append(file)

        if tags == 'WAVE':
            arc_files.append(file)

        if tags == 'FLAT,SKY':
            twilight_files.append(file)

        if tags == 'STD':
            print("Found STD observation of {}".format(
                fits.getval(file, "OBJECT")))
            std_files.append(file)

    print("Finished sort.")
    print("Found {} science frames.".format(len(science_files)))
    print("Found {} bias frames.".format(len(bias_files)))
    print("Found {} dark frames.".format(len(dark_files)))
    print("Found {} lamp flats.".format(len(flat_files)))
    print("Found {} wavecal/arc exposures.".format(len(arc_files)))
    print("Found {} sky/twilight flats.".format(len(twilight_files)))
    print("Found {} fluxcal standard frames".format(len(std_files)))

    if len(science_files) == 0:
        print("WARNING: No science files found! Make sure you've properly specified your target name.")

    ## MAKE BIAS.SOF ########

    bias_sof_file = "bias.sof"

    if os.path.exists(bias_sof_file):
        print("Old {} found. Overwriting it.".format(bias_sof_file))
        os.remove(bias_sof_file)

    if len(bias_files) > 10:
        print("Only using the last 10 bias frames - you don't need more than that.")
        bias_files_filtered = bias_files[-10:]
    else:
        bias_files_filtered = bias_files

    for file in bias_files_filtered:
        sof = open(bias_sof_file, mode="a")
        sof.write(file + "    BIAS\n")
    sof.close()

    print("----------------  {} FILE CREATED ---------------\n".format(bias_sof_file.upper()))
    print(open(bias_sof_file, "r").read())
    print("-------------------------------------------------------")

    ## MAKE DARK.SOF ########

    dark_sof_file = "dark.sof"

    if os.path.exists(dark_sof_file):
        print("Old {} found. Overwriting it.".format(dark_sof_file))
        os.remove(dark_sof_file)

    for file in dark_files:
        sof = open(dark_sof_file, mode="a")
        sof.write(file + "    DARK\n")

    sof.write("MASTER_BIAS.fits    MASTER_BIAS\n")
    sof.close()

    print("----------------  {} FILE CREATED ---------------\n".format(dark_sof_file.upper()))
    print(open(dark_sof_file, "r").read())
    print("-------------------------------------------------------")

    ## MAKE FLAT.SOF ########

    flat_sof_file = "flat.sof"

    if os.path.exists(flat_sof_file):
        print("Old {} found. Overwriting it.".format(flat_sof_file))
        os.remove(flat_sof_file)

    if len(flat_files) > 10:
        print("Only using the last 10 flat frames - you don't need more than that.")
        flat_files_filtered = flat_files[-10:]
    else:
        flat_files_filtered = flat_files

    for file in flat_files_filtered:
        sof = open(flat_sof_file, mode="a")
        sof.write(file + "    FLAT\n")

    sof.write("MASTER_BIAS.fits    MASTER_BIAS\n")
    sof.write("MASTER_DARK.fits    MASTER_DARK\n")
    sof.close()

    print("----------------  {} FILE CREATED ---------------\n".format(flat_sof_file.upper()))
    print(open(flat_sof_file, "r").read())
    print("-------------------------------------------------------")

    ## MAKE WAVECAL.SOF ########

    arc_sof_file = "arc.sof"

    if os.path.exists(arc_sof_file):
        print("Old {} found. Overwriting it.".format(arc_sof_file))
        os.remove(arc_sof_file)

    if len(arc_files) > 10:
        print("Only using the last 10 arc frames - you don't need more than that.")
        arc_files_filtered = arc_files[-10:]
    else:
        arc_files_filtered = arc_files

    for file in arc_files_filtered:
        sof = open(arc_sof_file, mode="a")
        sof.write(file + "    ARC\n")

    sof.write("MASTER_BIAS.fits    MASTER_BIAS\n")
    sof.write("TRACE_TABLE.fits    TRACE_TABLE\n")
    sof.write("{}/line_catalog.fits    LINE_CATALOG\n".format(static_cal_dir))
    sof.close()

    print("----------------  {} FILE CREATED ---------------\n".format(arc_sof_file.upper()))
    print(open(arc_sof_file, "r").read())
    print("-------------------------------------------------------")

    ## MAKE LSF.SOF ########

    lsf_sof_file = "lsf.sof"

    if os.path.exists(lsf_sof_file):
        print("Old {} found. Overwriting it.".format(lsf_sof_file))
        os.remove(lsf_sof_file)

    if len(arc_files) > 10:
        print("Only using the last 10 lsf frames - you don't need more than that.")
        arc_files_filtered = arc_files[-10:]
    else:
        arc_files_filtered = arc_files

    for file in arc_files_filtered:
        sof = open(lsf_sof_file, mode="a")
        sof.write(file + "    ARC\n")

    sof.write("MASTER_BIAS.fits    MASTER_BIAS\n")
    sof.write("TRACE_TABLE.fits    TRACE_TABLE\n")
    sof.write("MASTER_FLAT.fits    MASTER_FLAT\n")
    sof.write("MASTER_DARK.fits    MASTER_DARK\n")
    sof.write("WAVECAL_TABLE.fits    WAVECAL_TABLE\n")
    sof.write("{}/line_catalog.fits    LINE_CATALOG\n".format(static_cal_dir))
    sof.close()

    print("----------------  {} FILE CREATED ---------------\n".format(lsf_sof_file.upper()))
    print(open(lsf_sof_file, "r").read())
    print("-------------------------------------------------------")

    ## MAKE TWILIGHT.SOF ########

    twilight_sof_file = "twilight.sof"

    if os.path.exists(twilight_sof_file):
        print("Old {} found. Overwriting it.".format(twilight_sof_file))
        os.remove(twilight_sof_file)

    for file in twilight_files:
        sof = open(twilight_sof_file, mode="a")
        sof.write(file + "    SKYFLAT\n")

    sof.write("MASTER_BIAS.fits    MASTER_BIAS\n")
    sof.write("TRACE_TABLE.fits    TRACE_TABLE\n")
    sof.write("MASTER_FLAT.fits    MASTER_FLAT\n")
    sof.write("MASTER_DARK.fits    MASTER_DARK\n")
    sof.write("WAVECAL_TABLE.fits    WAVECAL_TABLE\n")
    sof.write(
        "{}/geometry_table_wfm.fits    GEOMETRY_TABLE\n".format(static_cal_dir))
    sof.write("{}/vignetting_mask.fits    VIGNETTING_MASK\n".format(static_cal_dir))
    sof.close()

    print("----------------  {} FILE CREATED ---------------\n".format(twilight_sof_file.upper()))
    print(open(twilight_sof_file, "r").read())
    print("-------------------------------------------------------")

    ## MAKE STD_SCIBASIC.SOF ########

    std_sof_file = "std_scibasic.sof"

    if os.path.exists(std_sof_file):
        print("Old {} found. Overwriting it.".format(std_sof_file))
        os.remove(std_sof_file)

    if len(std_files) > 1:
        print("Found more than 1 STD observation. Check that this is correct. I'm only using the last one!")
        std_files_filtered = std_files[-1]
    else:
        std_files_filtered = std_files

    for file in std_files_filtered:
        sof = open(std_sof_file, mode="a")
        sof.write(file + "    STD\n")

    sof.write("MASTER_BIAS.fits    MASTER_BIAS\n")
    sof.write("TRACE_TABLE.fits    TRACE_TABLE\n")
    sof.write("MASTER_FLAT.fits    MASTER_FLAT\n")
    sof.write("MASTER_DARK.fits    MASTER_DARK\n")
    sof.write("TWILIGHT_CUBE.fits    TWILIGHT_CUBE\n")
    sof.write("WAVECAL_TABLE.fits    WAVECAL_TABLE\n")
    sof.write(
        "{}/geometry_table_wfm.fits    GEOMETRY_TABLE\n".format(static_cal_dir))
    sof.write("{}/badpix_table.fits    BADPIX_TABLE\n".format(static_cal_dir))
    sof.close()

    print("----------------  {} FILE CREATED ---------------\n".format(std_sof_file.upper()))
    print(open(std_sof_file, "r").read())
    print("-------------------------------------------------------")

    ## MAKE ALIGN.SOF ########

    align_sof_file = "align.sof"

    if os.path.exists(align_sof_file):
        print("Old {} found. Overwriting it.".format(align_sof_file))
        os.remove(align_sof_file)

    print("Creating {} based on counting {} science frames, assuming standard naming conventions.".format(
        align_sof_file, len(science_files)))

    for i in range(len(science_files)):
        sof = open(align_sof_file, mode="a")
        sof.write("IMAGE_FOV_0001_{}.fits".format(i + 1) + "    IMAGE_FOV\n")
    sof.close()

    print("\n----------------  {} FILE CREATED ---------------\n".format(align_sof_file.upper()))
    print(open(align_sof_file, "r").read())
    print("-------------------------------------------------------")

    ## MAKE SCIENCE_SCIBASIC.SOF ########

    science_scibasic_sof_file = "science_scibasic.sof"

    if os.path.exists(science_scibasic_sof_file):
        print("Old {} found. Overwriting it.".format(science_scibasic_sof_file))
        os.remove(science_scibasic_sof_file)

    for file in science_files:
        sof = open(science_scibasic_sof_file, mode="a")
        sof.write(file + "    OBJECT\n")

    sof.write("MASTER_BIAS.fits    MASTER_BIAS\n")
    sof.write("TRACE_TABLE.fits    TRACE_TABLE\n")
    sof.write("MASTER_FLAT.fits    MASTER_FLAT\n")
    sof.write("MASTER_DARK.fits    MASTER_DARK\n")
    sof.write("WAVECAL_TABLE.fits    WAVECAL_TABLE\n")
    sof.write("TWILIGHT_CUBE.fits    TWILIGHT_CUBE\n")
    sof.write(
        "{}/geometry_table_wfm.fits    GEOMETRY_TABLE\n".format(static_cal_dir))
    sof.write("{}/badpix_table.fits    BADPIX_TABLE\n".format(static_cal_dir))
    sof.close()

    print("----------------  {} FILE CREATED ---------------\n".format(science_scibasic_sof_file.upper()))
    print(open(science_scibasic_sof_file, "r").read())
    print("-------------------------------------------------------")

    ## MAKE SCIENCE_SCIPOST_N.SOF ########

    # scrub the files
    for i in range(len(science_files)):
        if os.path.exists("science_scipost_{}.sof".format(i + 1)):
            print("science_scipost_{}.sof found. Overwriting it.".format(i + 1))
            os.remove("science_scipost_{}.sof".format(i + 1))

    for i in range(len(science_files)):
        sof = open("science_scipost_{}.sof".format(i + 1), mode="a")
        for j in range(24):
            # you need to zero pad the first 9 numbers
            sof.write("PIXTABLE_OBJECT_000{}-{:02}".format(i +
                                                           1, j + 1) + "    PIXTABLE_OBJECT\n")
        sof.write("STD_RESPONSE_0001.fits    STD_RESPONSE\n")
        sof.write("STD_TELLURIC_0001.fits    STD_TELLURIC\n")
        sof.write("LSF_PROFILE.fits    LSF_PROFILE\n")
        sof.write(
            "{}/astrometry_wcs_wfm.fits    ASTROMETRY_WCS\n".format(static_cal_dir))
        sof.write("{}/sky_lines.fits    SKY_LINES\n".format(static_cal_dir))
        sof.write("{}/extinct_table.fits    EXTINCT_TABLE\n".format(static_cal_dir))
        sof.write("{}/filter_list.fits    FILTER_LIST\n".format(static_cal_dir))
        sof.close()

    ## MAKE COMBINE.SOF ########

    combine_sof_file = "combine.sof"

    if os.path.exists(combine_sof_file):
        print("Old {} found. Overwriting it.".format(combine_sof_file))
        os.remove(combine_sof_file)

    print("Creating {} based on counting {} science frames, assuming standard naming conventions.".format(
        combine_sof_file, len(science_files)))

    for i in range(len(science_files)):
        sof = open(combine_sof_file, mode="a")
        sof.write("PIXTABLE_REDUCED_0001_{}.fits".format(
            i + 1) + "    PIXTABLE_REDUCED\n")
    sof.write("{}/filter_list.fits    FILTER_LIST".format(static_cal_dir))
    sof.close()

    print("\n----------------  {} FILE CREATED ---------------\n".format(combine_sof_file.upper()))
    print(open(combine_sof_file, "r").read())
    print("-------------------------------------------------------")
