

## The **F**ramework for **A**quatic **B**iogeochemical **M**odels (**FABM**)

FABM is a Fortran 2003 programming framework for biogeochemical models of marine and freshwater systems.

FABM is published in "Bruggeman, J., Bolding, K., 2014. A general framework for aquatic biogeochemical models. Environmental Modelling & Software 61: 249–265. DOI: [10.1016/j.envsoft.2014.04.002](http://dx.doi.org/10.1016/j.envsoft.2014.04.002)"

Further information is provided on [the FABM wiki](http://fabm.net/wiki)

## This is a FABM version 0 that supports AGG model (version 1.0.0).
please go to `src/models/hzg/agg` for  [source code of AGG](https://github.com/en-pei/fabm0_agg/blob/8461eac78cb93c67235f7fedb02a2ac0a1b175eb/src/models/hzg/agg/) .

Setups for the `AGG` model see:

`src/models/hzg/agg/labsetup`   [link](https://github.com/en-pei/fabm0_agg/tree/master/src/models/hzg/agg/labsetup) .


`src/models/hzg/agg/fieldsetup`  [link](https://github.com/en-pei/fabm0_agg/tree/f8f95f43f2f9844545993d3a766b418f182806d2/src/models/hzg/agg/fieldsetup) . 



More info on the `AGG` model see README file:  
https://github.com/en-pei/fabm0_agg/blob/master/src/models/hzg/agg/README.md

[![DOI](https://zenodo.org/badge/923541153.svg)](https://doi.org/10.5281/zenodo.14802720)


## Parameter variation and post-processing 
Please see `/parallel` for scripts of parallel computing and systematic variation of parameters.

Please see `/parallel/Plot` for scripts of postprocessing and plotting scripts.


## Plotting for the settling manuscript
Please see `settling_plotting.py` and run with Python3 for reproduce the figures.

