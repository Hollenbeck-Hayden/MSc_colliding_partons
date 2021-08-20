# Computation of theoretical predictions

The scope of this module is to compute theoretical predictions corresponding to the data sets for pion productin in pp collisions stored in `../data`.

This module includes:
- a separate folder for each experiment that contains
  - a `<nameexp>.txt` file: this is written by the user according to the following structure:
    ```
    pT integration: 0 for NO, 1 for YES
    y  integration: 0 for NO, 1 for YES
    perturbative order: 0 for LO, 1 for NLO
    PDF set: 3 for NNPDF_31_nlo_as_0118
    FF set: 4 for MAPFF1.0
    choice for initial hadrons: 0 for proton, 1 for antiproton
    choice of cross section: 1 for dsigma/dy/dpt2, 2 for E*dsigma/d3p, 3 for dsigma/dy/dpt, 4 for pT3*dsigma/dy/dpT
    rescaling factor for factorisation and renormalisation scales
    ```
  - a `results` folder in which the results are stored.

The results can be produced by running hte code available in `../src`. The code can be compiled by tunning `make` and run with `./run.sh <nmeexp> <irep>` where `<nameexp>` is the name of the experiment, as specified in the `../data` folder and `<irep>` is the replica in the MAPFF1.0 FF set.