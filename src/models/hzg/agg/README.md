# AGG

This is AGG, developed as part of the FABM framework.

## Citation

If you use this model, please cite:
- Enpei Li. "AGG" [DOI or link to your publication or Zenodo archive].

Additionally, please cite the FABM framework:
- "Bruggeman, J., Bolding, K., 2014. A general framework for aquatic biogeochemical models. Environmental Modelling & Software 61: 249–265. DOI: 10.1016/j.envsoft.2014.04.002"

## License
Licensed under the Apache License, Version 2.0
(http://www.apache.org/licenses/LICENSE-2.0)


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
Successful compiling will result in a gotm executable. Run by ./gotm
Setup files yamls are needed for running gotm. See next step.

_Below is an example of how you can instruct your audience on installing and setting up your app. This template doesn't rely on any external dependencies or services._

## Run GOTM
### Set up (yaml files)

Get setup yaml and meteo file from: https://github.com/en-pei/setupyamls.git
2.1.1 gotm.yaml: parameters of physics model, depth, layers, time… etc
2.1.2 fabm.yaml: parameters of each biogeochemical model.
2.1.3 output.yaml: specify which output to write to nc and which output frequency.




<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

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

* [](Arne...)
* [](Richard ...)


<p align="right">(<a href="#readme-top">back to top</a>)</p>

