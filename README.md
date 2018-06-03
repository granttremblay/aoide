# Aoide

*"[Aoide](https://en.wikipedia.org/wiki/Aoide),
one of the nine daughters of Zeus and Mnemosyne..."*

#### Dr. Grant R. Tremblay | Astrophysicist | Harvard-Smithsonian Center for Astrophysics


Aoide is a suite of Python tools for [MUSE](https://www.eso.org/sci/facilities/develop/instruments/muse.html),
an optical Integral Field Unit (IFU) spectrograph on ESO's Very Large Telescope.


### Step 1: Check for Raw Data Completeness

`cd` to the `raw_data_directory` in which you have placed all `*.fits.fz` raw MUSE data downloaded from the [ESO Archive](http://archive.eso.org).

Run
```
AoideReduce.py > contents.txt
tail -f contents.txt
```
with no arguments. It will print simple information about the directory contents. Make sure it looks okay (i.e., you have the correct data for your requested observation). You can save this list by, e.g., typing
```
AoideReduce.py > contents.txt
tail -f contents.txt
```

### Step 2: Use Aoide to run the MUSE Pipeline




![before_pca](misc/before_pca.png)
![after_pca](misc/after_pca.png)
