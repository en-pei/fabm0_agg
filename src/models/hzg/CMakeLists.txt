add_library(fabm_models_hzg OBJECT
            hzg_model_library.F90
            omexdia_p/omexdia_p.F90
            ecosmo/ecosmo.F90
            test/test.F90
            icealgea/icealgea.F90
           )

add_dependencies(fabm_models_hzg fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------
