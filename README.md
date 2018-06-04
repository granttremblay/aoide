# Aoide

*"[Aoide](https://en.wikipedia.org/wiki/Aoide),
one of the nine daughters of Zeus and Mnemosyne..."*
<!--
<img src="misc/A2597_movie.gif" alt="MUSE is awesome" style="width: 200px;"/> -->

#### Dr. Grant R. Tremblay | Astrophysicist | Harvard-Smithsonian Center for Astrophysics


Aoide is a suite of Python tools for [MUSE](https://www.eso.org/sci/facilities/develop/instruments/muse.html),
an optical Integral Field Unit (IFU) spectrograph on ESO's Very Large Telescope.


### Step 0: Optional Convenience Aliases

In my `.bashrc` on Linux (or `.bash_profile` on macOS), I have set the following BASH aliases:
```
alias aoideid='/home/grant/Repositories/aoide/bin/AoideID.py'
alias aoidereduce='/home/grant/Repositories/aoide/bin/AoideReduce.py'
alias aoidepost='/home/grant/Repositories/aoide/bin/AoidePost.py'
```

You may want to do the same! In the below examples, I assume that Aoide scripts
such as `AoideReduce.py` are in your PATH. This will not automatically be the case, of course.
You could fix that by putting
```
export PATH=$PATH:/path/to/aoide/bin
```
in your `.bashrc`.

### Step 1: Check for Raw Data Completeness

`cd` to the `raw_data_directory` in which you have placed all `*.fits.fz` raw MUSE data downloaded from the [ESO Archive](http://archive.eso.org).

Run
```
AoideReduce.py
```
with no arguments. It will print simple information about the directory contents. Make sure it looks okay (i.e., you have the correct data for your requested observation). You can save this list by, e.g., typing
```
AoideReduce.py > contents.txt
tail -f contents.txt
```


### Step 2: Use Aoide to run the MUSE Pipeline

`AoideReduce.py` will

1. Take inventory of the contents of a MUSE raw data directory (make sure your FITS files are straight from the archive, e.g. in `*fits.fz` format.)
2. Create `.sof` ("set of frames") files, based upon requirements outlined in the [MUSE Pipeline User Manual](https://www.eso.org/sci/software/pipelines/muse/muse-pipe-recipes.html), for use with the `esorex muse_*` recipes.
3. Run the ESO MUSE Pipeline, again following all basic steps in the Pipeline Manual.  

On an Ubuntu 18.04 workstation with an Intel Xeon E5-1650 v3 (6 cores, 3.8 GHz) and 64 GB of RAM,
one run of `AoideReduce` takes ~90 minutes for a three-pointing science observation. Peak RAM useage
approaches 60 GB during the `muse_scipost` and `muse_exp_combine` steps, so beware.
If your processor supports hyperthreading (virtual cores), note that treating these
as "real" processors is useless, and will likely *decrease* performance.


### Step 2: Fit & Subtract Sky Residuals, Bin & Extinction-correct the Datacube

`AoidePost.py` will ....

```
aoidepost "DATACUBE_AOIDE_UNCLEAN.fits" -c 6 -av 0.166
```


![before_pca](misc/before_pca.png)
![after_pca](misc/after_pca.png)
