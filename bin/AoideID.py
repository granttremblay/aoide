#!/usr/bin/env python

import glob
from astropy.io import fits

print("\n======================    Welcome to AoideID   ====================\n")
print("Taking inventory of your MUSE raw data directory.\n")


rawfiles = glob.glob("*.fits.fz")

if len(rawfiles) == 0:
    sys.exit("ERROR. I can find no files here. Make sure you run me in a MUSE raw data directory, and that the FITS files are Fpacked with a *.fz extension.")

# Initialize empty lists
science_files = []  # Science frames (i.e. OBJECT exposures go here).
bias_files = [] # Bias frames. We really want to only use 10 of these, often there will be more.
dark_files = []  # Dark frames
flat_files = []  # Lamp flats
arc_files = []  # Wavelength calibration files (arc lamps)
twilight_files = []  # Sky flats
std_files = []  # STD observations



for file in sorted(rawfiles):

    tags = fits.getval(file, "ESO DPR TYPE")

    if tags == 'OBJECT':
            science_files.append(file)

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
        std_files.append(file)

print("Finished sort.")

print("\n==============================================================")
print("Found {} science frames.".format(len(science_files)))
print("Found {} bias frames.".format(len(bias_files)))
print("Found {} dark frames.".format(len(dark_files)))
print("Found {} lamp flats.".format(len(flat_files)))
print("Found {} wavecal/arc exposures.".format(len(arc_files)))
print("Found {} sky/twilight flats.".format(len(twilight_files)))
print("Found {} fluxcal standard frames".format(len(std_files)))

print("\n======================    SCIENCE FILES   ====================\n")

for file in science_files:
    print("{} is {} | OBJECT IS {} | EXPTIME is {}".format(file, fits.getval(file, "ESO DPR TYPE"), fits.getval(file, "OBJECT"), fits.getval(file, "EXPTIME")))

print("\n======================    STD FILES   ====================\n")

for file in std_files:
    print("{} is {} | EXPTIME is {}".format(file, fits.getval(file, "ESO DPR TYPE"), fits.getval(file, "EXPTIME")))


print("\n======================    BIAS FILES   ====================\n")

for file in bias_files:
    print("{} is {} | EXPTIME is {}".format(file, fits.getval(file, "ESO DPR TYPE"), fits.getval(file, "EXPTIME")))


print("\n======================    DARK FILES   ====================\n")

for file in dark_files:
    print("{} is {} | EXPTIME is {}".format(file, fits.getval(file, "ESO DPR TYPE"), fits.getval(file, "EXPTIME")))

print("\n======================    FLAT FILES   ====================\n")

for file in flat_files:
    print("{} is {} | EXPTIME is {}".format(file, fits.getval(file, "ESO DPR TYPE"), fits.getval(file, "EXPTIME")))


print("\n======================    ARC FILES   ====================\n")

for file in arc_files:
    print("{} is {} | EXPTIME is {}".format(file, fits.getval(file, "ESO DPR TYPE"), fits.getval(file, "EXPTIME")))

print("\n======================    TWILIGHT FILES   ====================\n")

for file in twilight_files:
    print("{} is {} | EXPTIME is {}".format(file, fits.getval(file, "ESO DPR TYPE"), fits.getval(file, "EXPTIME")))

print("\nBEFORE RUNNING AoideReduce, CHECK THAT THE ABOVE LOOKS OKAY\n")
print("\n======================    FINISHED ID RUN   ====================\n")
