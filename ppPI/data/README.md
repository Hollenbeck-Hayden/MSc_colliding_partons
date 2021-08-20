# Filtering data

The scope of this module is to filter the raw data (typically as available from
Hepdata) and to recast it in a common format.

This module includes:
- a separate data folder for each experiment that contains two files
  - a `metadata.yaml` file: this is written by the user according to the following structure:
    ```
    - arxiv: arxiv number
    - paper: details of the publication
    - DOI: web link to the DOI
    - hepdata: web link to the relevant hepdata entry
    - nameexp: name of the experiment in the format <name_energy_hadron>
    - ndata: number of data points
    - nsys: number of experimental systematic uncertainties
    - systype: a list of labels, for each systematic uncertainty, that could be either UNCORR or CORR
    - namesys: a list of labels, for each systematic uncertainty, that could be either ADD or MULT
    - obs: the measured bservable
    - unit: the unit of the measured observable
    - cme: the centre-of-mass energy of the pp collision
    - rap: the rapidity range [ymin,ymax] of the measurement
    - info: any additional information
    ```
  - a raw data file, in `.csv` format, typically downloaded from the Hepdata repository, see https://www.hepdata.net/
- a `filter.py` file in which the raw data in each of the seprate folders is recasted into a common format. This common format is as follows.
  ```
  nameexp
  ndata nsys
  cme
  ymin ymax
  id    pT [GeV]  obs  stat  sys1  ...  sysn
  ..    ..        ..   ..    ..    ...  ..
  ```
  where `obs` is usually E d3sigma/dp3 [mb GeV-2], `stat` is the statistical uncertainty and `sys1`, ..., `sysn` are the n systematic uncertainties.