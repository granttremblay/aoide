#!/usr/bin/env python

import glob
from astropy.io import fits

rawfiles = glob.glob("*.fits.fz")

for file in rawfiles:
    print("{} | DPR TYPE {} | TARGNAME is {} | EXPTIME is {}s".format(file, fits.getval(file, "ESO DPR TYPE"), fits.getval(file, "OBJECT"), fits.getval(file, "EXPTIME")))
