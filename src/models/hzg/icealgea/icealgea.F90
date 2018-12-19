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
      type (type_surface_state_variable_id)          :: id_no3, id_nh4, id_pho, id_sil
      type (type_surface_state_variable_id)          :: id_IA, id_det
      type (type_horizontal_diagnostic_variable_id)  :: id_primprod, id_denit
      type (type_horizontal_diagnostic_variable_id)  :: id_hbio
      type (type_dependency_id)                      :: id_temp, id_par
      type (type_horizontal_dependency_id)           :: id_h, id_c, id_icetemp, id_icesalt
      type (type_horizontal_dependency_id)           :: id_icedh

!     Model parameters
      real(rk) :: muIA
      real(rk) :: TctrlIA
      real(rk) :: aaIA
      real(rk) :: M0
      real(rk) :: rM
      real(rk) :: R0
      real(rk) :: rR
      real(rk) :: Pm
      real(rk) :: rNO3
      real(rk) :: rPO4
      real(rk) :: rSi
      real(rk) :: rnit
      real(rk) :: vn

    
      contains

!     Model procedures
      procedure :: initialize
      procedure :: do_surface

   end type type_hzg_icealgea

   real(rk), parameter :: one_pr_day=1.0_rk/86400.0_rk
   real(rk), parameter :: one_pr_hour = 1.0_rk/3600.0_rk
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

   call self%get_parameter(self%muIA, 'muIA', '1/day', 'ice algae maximum growth rate', default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%TctrlIA, 'TctrlIA', '1/degC', 'temp growth dependency', default=0.001_rk)
   call self%get_parameter(self%M0, 'M0', '1/day', 'maximum mortality rate', default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rM,  'rM', '1/degC', 'mortality temp coefficient', default=0.5_rk)

   call self%get_parameter(self%R0,  'R0', '1/day', 'maximum respiration rate', default=0.5_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rR,  'rR', '1/degC', 'respiration temp coefficient', default=0.5_rk)
   call self%get_parameter(self%aaIA,  'aaIA', 'mgC/mgchla/h/Einst/m/s', 'photosynthesis ef-cy', default=0.5_rk, scale_factor=one_pr_hour)
   call self%get_parameter(self%Pm,  'Pm', 'mgC/mgchla/h', 'maximum photosynhtetic rate', default=0.5_rk, scale_factor=one_pr_hour)
   call self%get_parameter(self%rNO3,  'rNO3', 'mmol/m**3', 'NO3 half saturation', default=0.5_rk)
   call self%get_parameter(self%rPO4,  'rPO4', 'mmol/m**3', 'PO4 half saturation', default=0.5_rk)
   call self%get_parameter(self%rSi,  'rSi', 'mmol/m**3', 'SiO2 half saturation', default=0.5_rk)
   call self%get_parameter(self%rnit,  'rnit', '1/day', 'nitrification rate', default=0.5_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%vn,  'vn', 'mmolN/m**3', 'preferential uptake of nitrate half saturation', default=0.5_rk)

   call self%register_surface_state_variable(self%id_IA, 'IA', 'mg/m2', 'ice algae', minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_no3, 'no3', 'mmol/m3', 'nitrate', minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_nh4, 'nh4', 'mmol/m3', 'nitrate', minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_pho, 'pho', 'mmol/m3', 'phosphate', minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_sil, 'sil', 'mmol/m2', 'silicate', minimum=0.0_rk)
   
   call self%register_horizontal_diagnostic_variable(self%id_primprod, 'primprodIA','mg/m3/s', 'primary production')
   call self%register_horizontal_diagnostic_variable(self%id_hbio, 'hbio','m', 'thickness of biological active layer')

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_h,'icethickness','m','ice thickness')
   call self%register_dependency(self%id_c,'iceconc',     '', 'ice concentration')
   call self%register_dependency(self%id_icetemp,'icetemp',     'degC',  'ice temperature')
   call self%register_dependency(self%id_icesalt,'icesalt',     'PSU', 'ice salinity')
   call self%register_dependency(self%id_icedh,'icedh',     'm/s', 'rate of change of ice thickness')

   return

   end subroutine initialize



!
! Surface fluxes for the icealgea model

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_hzg_icealgea),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: no3, nh4, pho, sil
   real(rk) :: water_temp, water_no3, water_pho, water_sil 
   real(rk) :: ice_top_temp, ice_salt
   real(rk) :: IA
   real(rk) :: up_no3, up_pho, up_sil
   real(rk) :: pno3, rhs_nit, rhs_amm
   real(rk) :: blight, tdep, sdep, IA_prod, prod, IA_R, IA_M 
   real(rk) :: icethickness, iceconc, ice_dh
   real(rk) :: hbio
   real(rk) :: par


   _HORIZONTAL_LOOP_BEGIN_

   !_GET_(self%id_temp, water_temp)
   !_GET_(self%id_no3, water_no3)
   !_GET_(self%id_pho, water_pho)
   !_GET_(self%id_sil, water_sil)
   _GET_(self%id_par, par) 
   _GET_HORIZONTAL_(self%id_no3, no3)
   _GET_HORIZONTAL_(self%id_nh4, nh4)
   _GET_HORIZONTAL_(self%id_pho, pho)
   _GET_HORIZONTAL_(self%id_sil, sil)
   _GET_HORIZONTAL_(self%id_IA, IA)
   _GET_HORIZONTAL_(self%id_h, icethickness)
   _GET_HORIZONTAL_(self%id_c, iceconc)
   _GET_HORIZONTAL_(self%id_icetemp, ice_top_temp)
   _GET_HORIZONTAL_(self%id_icesalt, ice_salt)
   _GET_HORIZONTAL_(self%id_icedh, ice_dh)

   ! remineralisation rate

  
   ! nutrient limitation factors
   up_no3 = no3 + nh4/(self%rNO3 + (no3 + nh4))
   up_pho = pho/(self%rPO4 + pho)
   up_sil = sil/(self%rSi + sil)

   ! light limitation
   blight = max((tanh(self%aaIA *par) / self%Pm), 0.0_rk)

   ! temperature and salt dependence
   Tdep = exp(self%TctrlIA * ice_top_temp)
   
   Sdep = exp(-(2.16 - 8.30 * (10**-5) * (ice_salt ** 2.11) - 0.55 * log(ice_salt)) ** 2)  

   ! ice algae production and nutrient uptake
   IA_prod = Tdep * Sdep * min(blight, up_no3, up_pho, up_sil)
   prod = self%muIA * IA_prod * IA
 
   ! ice algae losses
   IA_R = self%R0 * self%muIA * exp(self%rR * ice_top_temp) 
   IA_M = self%M0 * exp(self%rM * ice_top_temp)

   ! nitrate and ammonium 
   !ent_nit = 
   pno3 = self%vn / (self%vn + nh4)
   rhs_nit = (-self%muIA * IA_prod * IA * pno3) + (self%rnit * nh4)
   rhs_amm = (-self%muIA * IA_prod * IA * (1-pno3)) - (self%rnit * nh4)     

   ! Thickness of biologically active layer
   hbio = 0.4_rk * icethickness


   _SET_SURFACE_ODE_(self%id_IA, (self%muIA * IA_prod - IA_R - IA_M)*IA )
   _SET_SURFACE_ODE_(self%id_no3, rhs_nit)
   _SET_SURFACE_ODE_(self%id_nh4, rhs_amm)

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_primprod, prod)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_hbio, hbio)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface


   end module fabm_hzg_icealgea

