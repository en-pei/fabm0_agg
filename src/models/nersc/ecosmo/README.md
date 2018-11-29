ECOSMO
======

Switches:
---------

use_cyanos: [.true./.false.] switch cyanobacteria prognostic dynamics on/off,
            default is .true.
            in case of .false. remove initialization of variable "bg" and "bgchl"

use_chl: [.true./.false.] switch chlorophyll prognostic dynamics on/off.
         default is .true., in case of .false. remove initialization of 
         variables "diachl","flachl", and "bgchl".
         If .true. "exphy" parameter in the yaml configuration  will be
         used as linear extinction factor to the chlorophyll concentration.

couple_co2: [.true./.false.] use coupling of dic rates to external carbonate model. This is tested with the pml/carbonate model. default is .false.
            If .true. use "dic_target" and "alk_target" under "coupling" in the ECOSMO configuration.
