# Colliding Partons at the LHC
This repository contains code used in the MSc TP Dissertation "Colliding Partons at the LHC" by Hayden
Hollenbeck at the University of Edinburgh.

## Table of Contents
* [Work Declaration](#work-declaration)
* [Dependencies](#dependencies)
* [Build](#build)
* [Running](#running)

## Work Declaration
All work in this repository is my own, except for the following.
Provided by E.R. Nocera:
* `babar/benchmark.cpp`
* `ppPI/data/filter.py` initial template was provided; substantial changes made since.
* `ppPI/src` scripts for parallel runs of unpolarized hadron code were made by myself.

Experimental data downloaded from [HepData](https://www.hepdata.net)
* `ppPI/data/\*` listed as .csv files.

## Dependencies
The following libraries are required to build and run the reweighting code:
* [LHAPDF6](https://lhapdf.hepforge.org)
* [NumPy](https://numpy.org)

Plotting requires MatplotLib:
* [MatplotLib](https://matplotlib.org)

Building the BABAR code requires:
* [APFEL](https://apfel.mi.infn.it)
* [yaml-cpp](https://github.com/jbeder/yaml-cpp)
* [LHAPDF6](https://lhapdf.hepforge.org)

Building the COMPASS code requires:
* [APFEL++](https://github.com/vbertone/apfelxx)
* [yaml-cpp](https://github.com/jbeder/yaml-cpp)
* [LHAPDF6](https://lhapdf.hepforge.org)

PDF and FF sets may be downloaded directly from [LHAPDF](https://lhapdf.hepforge.org/pdfsets.html).

## Build

To build the desired code, simply enter the directory with code to be built and run `make`.

## Running

The BABAR and COMPASS code are in standalone directories, and may be run directly from their directories.

Because of the lengthy calculations used in proton-proton collisions, the running is split into several steps.
First, run

```[bash]
ppPI/data $ python filter.py
```

To filter and reformat experimental data sets. Data sets used in reweighting are stored in their own folder
in `ppPI/data` in a CSV format. Parsing more experiments should be hard coded into `filter.py`.
Next, use the `run_all.sh` command in the `ppPI/src` directory to compute all of the cross sections.
Parameters for computing cross sections are defined in `ppPI/theory/experiment/experiment.txt`. Once
calculations have finished, `ppPI/theory/post_processing.py` should be run to generate the appropriate
conversion and formatting of the data.

Now the data is ready to be reweighted. Create an experimental set file containing a list of the experimental
data sets to be used in reweighting. In `ppPI/reweight`, run

```[bash]
ppPI/reweight $ ./make_set.sh set_name cutoff exp_set_file
```

This will compute the reweightings and create an unweighted set.
