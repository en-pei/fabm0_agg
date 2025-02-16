# AGG 1.0.0

This is `AGG`, a distribution-based aggregation-SPM model. 

It is coupled to other models and physical driver with the [FABM framework](https://github.com/fabm-model/fabm/wiki) .

The [iow_spm](https://github.com/en-pei/fabm0_agg/tree/d770ff65ddb2882371e50d3a071b0a2067c3b0d3/src/models/iow/spm) model is used to calculate resuspension. See subfolder: `/src/models/iow/spm` for source code.

## Citation

If you use this model, please cite:
- "Li, E., & Wirtz, K. (2025). _AGG model (Version 1.0.0)_. GitHub. Retrieved from https://github.com/en-pei/fabm0_agg/tree/2f265837b08193e02bdabb5f0a7f1deb51ff8614/src/models/hzg/agg"


Additionally, please cite the FABM framework:
- "Bruggeman, J., Bolding, K., 2014. A general framework for aquatic biogeochemical models. Environmental Modelling & Software 61: 249–265. DOI: 10.1016/j.envsoft.2014.04.002"

## License
Licensed under the GPL-2.0 


<!-- GETTING STARTED -->
## Getting Started
### Prerequisites
#### Compilers needed and versions
Refer to https://github.com/fabm-model/fabm/wiki/GOTM


### Installation
Used the GOTM version:

https://github.com/en-pei/gotmj_fabm0_strand.git


Then go to `extern/fabm` to clone the fabm0_agg fabm

https://github.com/en-pei/fabm0_agg.git

   ```sh
git submodule update --init –recursive
git submodule update --force --recursive --init --remote
   ```

If this does not work then use the official fabm within gotm then delete the fabm folder then get the fabm0_agg fabm (the fabm that contains  `agg.F90 `).


Note that two version of `agg.F90` are given in the folder ( `agg_lab.F90 ` and  `agg_field.F90 ` F90), by copying the lab version or the field version to  `agg.F90 ` then compile will give the right model version. Field version is for the TMZ (turbidity maximum zone) cases with do_bottom routine where erosion’s influence on size is calculated.




### cmake command: build `gotm` executable with FABM base
   ```sh
    cmake $GOTMDIR/src -DFABM_BASE=$FABMDIR/src
   ```
Successful compiling will result in a gotm executable. 





<!-- USAGE EXAMPLES -->
## Usage

### Run GOTM
   ```sh
   ./gotm
   ```
Setup files yamls are needed for running gotm. See next step.


### Set up (yaml files)
Get setup yaml and meteo file from:

 `src/models/hzg/agg/labsetup ` for a lab case.

 `src/models/hzg/agg/fieldsetup ` for a field case.

 
2.1.1 gotm.yaml: parameters of physics model, depth, layers, time… etc

2.1.2 fabm.yaml: parameters of each biogeochemical model.

2.1.3 output.yaml: specify which output to write to nc and which output frequency.




<p align="right">(<a href="#readme-top">back to top</a>)</p>




### Top contributors:
Enpei Li, Kai Wirtz


<!-- CONTACT -->
## Contact

Enpei Li - Li@bafg.de

Kai Wirtz - kai.wirtz@hereon.de

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
<!--* [](Arne...)?-->
[comment]: # (tesetttest)
* Dr. Richard Hofmeister 


<p align="right">(<a href="#readme-top">back to top</a>)</p>

