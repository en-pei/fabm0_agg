! Copyright 2017 Helmholtz-Zentrum Geesthacht
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Richard Hofmeister

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_test --- new biogeochemical model
!
! !INTERFACE:
   module fabm_hzg_test
!
! !DESCRIPTION:
!
! This ecosystem model has to be developed
!
! !USES:
   use fabm_types
   use fabm_expressions

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_hzg_test

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hzg_test
!     Variable identifiers
      type (type_state_variable_id)         :: id_tracer
      type (type_diagnostic_variable_id)    :: id_rate
      type (type_dependency_id)             :: id_temp, id_par

!     Model parameters
      real(rk) :: param

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom
      procedure :: get_light_extinction

   end type type_hzg_test
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the test model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the test namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hzg_test),intent(inout),target  :: self
   integer,              intent(in)            :: configunit


   call self%get_parameter(self%param, 'param', '1/day', 'params long description', default=0.001_rk, scale_factor=1.0_rk/86400._rk)

   call self%register_state_variable(self%id_tracer, 'tracer', 'mg/m3', 'reactive tracer', minimum=0.0_rk)
   call self%register_diagnostic_variable(self%id_rate, 'degrad_rate', 'mg/m3/s', 'reactive tracer degradation')

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)

   return

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of test model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_hzg_test),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk) :: temp,par
   real(rk) :: tracer


   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp,temp)
   _GET_(self%id_tracer, tracer)

   _SET_ODE_(self%id_tracer, -self%param*tracer)
   _SET_DIAGNOSTIC_(self%id_rate, -self%param*tracer)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC






!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for the test model
!
! !INTERFACE:

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_hzg_test),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   real(rk) :: temp
!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,temp)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC





!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Bottom fluxes for the test model
!
! !INTERFACE:

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_hzg_test),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_

   _HORIZONTAL_LOOP_BEGIN_

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC


   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
   class (type_hzg_test), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.

   !_GET_(self%id_det, det)

   !_SET_EXTINCTION_( 1.0*det )

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction

   end module fabm_hzg_test

