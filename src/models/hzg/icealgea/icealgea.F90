! Copyright 2017 Helmholtz-Zentrum Geesthacht
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Richard Hofmeister and Deborah Benkort

#include "fabm_driver.h"

!  fabm_hzg_icealgea --- new biogeochemical model
!
   module fabm_hzg_icealgea

! !DESCRIPTION:
!
! This ecosystem model has to be developed
!

   use fabm_types
   use fabm_expressions

   implicit none

   private

   public type_hzg_icealgea

   type,extends(type_base_model) :: type_hzg_icealgea
!     Variable identifiers
      type (type_surface_state_variable_id)          :: id_phy
      type (type_horizontal_diagnostic_variable_id)  :: id_primprod
      type (type_horizontal_diagnostic_variable_id)  :: id_hbio
      type (type_dependency_id)                      :: id_temp,id_par
      type (type_horizontal_dependency_id)           :: id_h, id_c, id_icetemp
      type (type_horizontal_dependency_id)           :: id_icedh

!     Model parameters
      real(rk) :: rgr,mort

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do_surface

   end type type_hzg_icealgea

   real(rk), parameter :: one_pr_day=1.0_rk/86400.0_rk
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!
! Initialise the icealgea model
!
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the icealgea namelist is read and the variables exported
!  by the model are registered with FABM.
!
   class (type_hzg_icealgea),intent(inout),target  :: self
   integer,                  intent(in)            :: configunit
   real(rk)                                        :: ws

   call self%get_parameter(self%mort, 'mortality', '1/day', 'phytoplankton mortality', default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rgr,  'rgr', '1/day', 'relative growth rate', default=0.5_rk, scale_factor=one_pr_day)

   call self%register_surface_state_variable(self%id_phy, 'phy', 'mg/m2', 'ice phytoplankton', minimum=0.0_rk)
   
   call self%register_horizontal_diagnostic_variable(self%id_primprod, 'primprod','mg/m3/s', 'primary production')
   call self%register_horizontal_diagnostic_variable(self%id_hbio, 'hbio','m', 'thickness of biological active layer')

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_h,'icethickness','m','ice thickness')
   call self%register_dependency(self%id_c,'iceconc',     '', 'ice concentration')
   call self%register_dependency(self%id_icetemp,'icetemp',     'degC', 'ice temperature')
   call self%register_dependency(self%id_icedh,'icedh',     'm/s', 'rate of change of ice thickness')

   return

   end subroutine initialize



!
! Surface fluxes for the icealgea model

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_hzg_icealgea),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: water_temp,ice_top_temp,phy
   real(rk) :: icethickness,iceconc,ice_dh
   real(rk) :: hbio


   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp, water_temp)
   !_GET_(self%id_nit, water_nit)
   _GET_HORIZONTAL_(self%id_phy, phy)
   _GET_HORIZONTAL_(self%id_h, icethickness)
   _GET_HORIZONTAL_(self%id_c, iceconc)
   _GET_HORIZONTAL_(self%id_icetemp, ice_top_temp)
   _GET_HORIZONTAL_(self%id_icedh, ice_dh)
  
   hbio = 0.4_rk * icethickness


   _SET_SURFACE_ODE_(self%id_phy, phy*self%rgr - phy*phy*self%mort)

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_primprod, phy*self%rgr)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_hbio, hbio)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface


   end module fabm_hzg_icealgea

