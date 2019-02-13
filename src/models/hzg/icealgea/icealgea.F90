! Copyright 2017 Helmholtz-Zentrum Geesthacht
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Richard Hofmeister and Deborah Benkort

#include "fabm_driver.h"

!  !MODULE: fabm_hzg_icealgea --- new biogeochemical model
!
!  !INTERFACE:
   module fabm_hzg_icealgea

! !DESCRIPTION:
!
! This ecosystem model has to be developed
!

! !USES:
   use fabm_types
   use fabm_expressions

   implicit none

!  default: all is private
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_hzg_icealgea
!
! !PRIVATE DERIVED TYPES:
   type,extends(type_base_model) :: type_hzg_icealgea
!     Variable identifiers
      type (type_surface_state_variable_id)          :: id_Ino3, id_Inh4, id_Ipho, id_Isil
      type (type_surface_state_variable_id)          :: id_Ioxy, id_Ialg, id_Idet
      type (type_state_variable_id)                  :: id_Wno3, id_Wnh4, id_Wpho, id_Wsil
      type (type_state_variable_id)                  :: id_Wdia, id_Wdet, id_Wmesozoo, id_Wmicrozoo, id_Woxy
      type (type_horizontal_diagnostic_variable_id)  :: id_primprod, id_limit, id_denit, id_BAL, id_entrap, id_flush
      type (type_dependency_id)                      :: id_temp, id_par, id_U0
      type (type_horizontal_dependency_id)           :: id_h, id_hs, id_c, id_icetemp, id_icesalt, id_brsalt, id_balpar
      type (type_horizontal_dependency_id)           :: id_hbio, id_icedh

!     Model parameters
      real(rk) :: muIA
      real(rk) :: TctrlIA
      real(rk) :: aaIA
      real(rk) :: rM
      real(rk) :: rGr
      real(rk) :: Pm
      real(rk) :: rNO3
      real(rk) :: rPO4
      real(rk) :: rSi
      real(rk) :: rnit
      real(rk) :: rrem
      real(rk) :: vn
      real(rk) :: C0
      real(rk) :: Catt_ice
      real(rk) :: Catt_sno
      real(rk) :: ice_alb
      real(rk) :: sno_alb
      real(rk) :: v
      real(rk) :: Ui
      real(rk) :: Cio
      real(rk) :: D

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do_surface

   end type type_hzg_icealgea

! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk
   real(rk), parameter :: one_pr_hour = 1.0_rk/3600.0_rk
   real(rk)            :: redf(20) = 0.0_rk

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

!  Set Redfield ratios:
   redf(1)  = 6.625_rk      !C_N
   redf(2)  = 106.0_rk      !C_P
   redf(3)  = 6.625_rk      !C_SiO
   redf(4)  = 16.0_rk       !N_P
   redf(5)  = 1.0_rk        !N_SiO
   redf(6)  = 12.01_rk      !C_Cmg
   redf(7)  = 44.6608009_rk !O2mm_ml
   redf(8)  = 14.007_rk     !N_Nmg
   redf(9)  = 30.97_rk      !P_Pmg
   redf(10) = 28.09_rk     !Si_Simg    

!------------------------------------------------------------------------------------

   call self%get_parameter(self%muIA,    'muIA',    '1/day',                  'ice algae maximum growth rate',                  default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%TctrlIA, 'TctrlIA', '1/degC',                 'temperature growth dependency',                  default=0.001_rk)
   call self%get_parameter(self%rM,      'rM',      '1/day',                  'ice algae mortality rate',                       default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rGr,     'rGr',     '1/day',                  'ice algae grazing rate',                         default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%aaIA,    'aaIA',    'mgC/mgchla/h/Einst/m/s', 'photosynthesis ef-cy',                           default=0.5_rk,   scale_factor=one_pr_hour)
   call self%get_parameter(self%Pm,      'Pm',      'mgC/mgchla/h',           'maximum photosynhtetic rate',                    default=0.5_rk,   scale_factor=one_pr_hour)
   call self%get_parameter(self%rNO3,    'rNO3',    'mmolN/m**3',             'NO3 half saturation',                            default=0.5_rk,   scale_factor=redf(1)*redf(6))
   call self%get_parameter(self%rPO4,    'rPO4',    'mmolP/m**3',             'PO4 half saturation',                            default=0.5_rk,   scale_factor=redf(2)*redf(6))
   call self%get_parameter(self%rSi,     'rSi',     'mmolSi/m**3',            'SiO2 half saturation',                           default=0.5_rk,   scale_factor=redf(3)*redf(6))
   call self%get_parameter(self%rnit,    'rnit',    '1/day',                  'nitrification rate',                             default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rrem,    'rrem',    '1/day',                  'remineralisation rate',                          default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%vn,      'vn',      'mmolN/m**3',             'preferential uptake of nitrate half saturation', default=0.5_rk,   scale_factor=redf(1)*redf(6))
   call self%get_parameter(self%C0,      'C0',      '',                       'part of incoming radiation absorbed',            default=0.3_rk)
   call self%get_parameter(self%Catt_ice,'Catt_ice','1/m',                    'ice attenuation coefficent',                     default=0.5_rk)                
   call self%get_parameter(self%Catt_sno,'Catt_sno','1/m',                    'snow attenuation coefficent',                    default=0.5_rk)                
   call self%get_parameter(self%ice_alb, 'ice_alb', '',                       'albedo in the ice',                              default=0.5_rk)                
   call self%get_parameter(self%sno_alb, 'sno_alb', '',                       'albedo in the snow',                             default=0.5_rk)                
   call self%get_parameter(self%v,       'v',       'm**2/s',                 'kinematic viscosity of seawater',                default=0.0_rk)                
   call self%get_parameter(self%Ui,      'Ui',      'm/s',                    'horizontal velocity of ice',                     default=0.0_rk)                
   call self%get_parameter(self%Cio,     'Cio',     'm**2/s',                 'drag coefficient',                               default=0.0_rk)                
   call self%get_parameter(self%D,       'D',       'm**2/s',                 'molecular diffusion coefficient',                default=0.0_rk)                
 
   call self%register_surface_state_variable(self%id_Ialg,    'Ialg',   'mgC/m**3',   'ice algae',      minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Idet,    'Idet',   'mgC/m**3',   'detritus',       minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Ino3,    'Ino3',   'mgC/m**3',   'nitrate',        minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Inh4,    'Inh4',   'mgC/m**3',   'nitrate',        minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Ipho,    'Ipho',   'mgC/m**3',   'phosphate',      minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Isil,    'Isil',   'mgC/m**3',   'silicate',       minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Ioxy,    'Ioxy',   'mmolO2/m**3','oxygen',         minimum=0.0_rk)
   
   call self%register_horizontal_diagnostic_variable(self%id_primprod,  'primprodIA', 'mgC/m3/s',     'primary production')
   call self%register_horizontal_diagnostic_variable(self%id_limit,     'limitIA',    ''        ,     'limitation factor for IA')
   call self%register_horizontal_diagnostic_variable(self%id_BAL,       'BAL',        'm',            'thickness of BAL')
   call self%register_horizontal_diagnostic_variable(self%id_entrap,    'entrap',     'mmolN/m**3',   'nutrient entrap')
   call self%register_horizontal_diagnostic_variable(self%id_flush,     'flush',      'mmolN/m**3',   'nutrient release')
   call self%register_horizontal_diagnostic_variable(self%id_denit,     'denit_ice',  'mmolN/m**3/s', 'denitrification rate')

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   !call self%register_dependency(self%id_U0,  standard_variables%x_velocity)

   call self%register_dependency(self%id_h,       'icethickness', 'm',    'ice thickness')
   call self%register_dependency(self%id_hs,      'snowthickness','m',    'snow thickness')
   call self%register_dependency(self%id_hbio,    'balthickness', 'm',    'BAL thickness')
!   call self%register_dependency(self%id_c,      'iceconc',      '',     'ice concentration')
   call self%register_dependency(self%id_icetemp, 'icetemp',      'degC', 'ice temperature')
   call self%register_dependency(self%id_icesalt, 'icesalt',      'PSU',  'ice salinity')
   call self%register_dependency(self%id_brsalt,  'brsalt',       'PSU',  'ice brine salinity')
   call self%register_dependency(self%id_balpar,  'balpar',       '',     'PAR at the BAL')
   call self%register_dependency(self%id_icedh,   'icedh',        'm/s',  'rate of change of ice thickness')

   call self%register_state_dependency(self%id_Wno3, 'Wno3', 'mmolN/m**3',  'nitrate come from surface water column')
   call self%register_state_dependency(self%id_Wnh4, 'Wnh4', 'mmolN/m**3',  'ammonium come from surface water column')
   call self%register_state_dependency(self%id_Wpho, 'Wpho', 'mmolP/m**3',  'phosphate come from surface water column')
   call self%register_state_dependency(self%id_Wsil, 'Wsil', 'mmolSi/m**3', 'silicate come from surface water column')
   call self%register_state_dependency(self%id_Woxy, 'Woxy', 'mmolO2/m**3', 'oxygen come from surface water column')

   return

   end subroutine initialize

!
! Surface fluxes for the icealgea model

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_hzg_icealgea),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: Ino3, Inh4, Ipho, Isil, Ioxy
   real(rk) :: Wno3, Wnh4, Wpho, Wsil 
   real(rk) :: ice_top_temp, ice_salt, br_salt
   real(rk) :: Ialg, Idet
   real(rk) :: up_Ino3, up_Ipho, up_Isil
   real(rk) :: flux_nit, ent_Inh4, ent_Ipho, ent_Isil, entrap, flush
   real(rk) :: pno3, rhs_Init, rhs_Iamm
   real(rk) :: blight, tdep, sdep
   real(rk) :: rem, IA_prod, Prod, IA_loss, IA_gr 
   real(rk) :: icethickness, iceconc, ice_dh
   real(rk) :: snowthickness, hbio
   real(rk) :: alb, par, bio_par
   real(rk) :: hv, diff_nit, U0

      
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_Wno3, Wno3)
   _GET_(self%id_Wnh4, Wnh4)
   _GET_(self%id_Wpho, Wpho)
   _GET_(self%id_Wsil, Wsil)
   _GET_(self%id_par, par) 
   _GET_HORIZONTAL_(self%id_Ino3, Ino3)
   _GET_HORIZONTAL_(self%id_Inh4, Inh4)
   _GET_HORIZONTAL_(self%id_Ipho, Ipho)
   _GET_HORIZONTAL_(self%id_Isil, Isil)
   _GET_HORIZONTAL_(self%id_Ialg, Ialg)
   _GET_HORIZONTAL_(self%id_Idet, Idet)
   _GET_HORIZONTAL_(self%id_h, icethickness)
   _GET_HORIZONTAL_(self%id_hs, snowthickness)
   _GET_HORIZONTAL_(self%id_hbio, hbio)
   _GET_HORIZONTAL_(self%id_c, iceconc)
   _GET_HORIZONTAL_(self%id_icetemp, ice_top_temp)
   _GET_HORIZONTAL_(self%id_icesalt, ice_salt)
   _GET_HORIZONTAL_(self%id_brsalt, br_salt)
   _GET_HORIZONTAL_(self%id_icedh, ice_dh)
   _GET_HORIZONTAL_(self%id_balpar, bio_par)

   hbio = 0.05
   U0   = 0.008

   ! par on the bottom ice
   if(snowthickness > 0.0_rk) then
	alb = self%sno_alb
   else 
	alb = self%ice_alb
   end if

   bio_par = par * (1.0 - alb) * self%C0 * exp((-self%Catt_ice * icethickness) - (self%Catt_sno * snowthickness))
   
   ! loss terms
   IA_loss = max(0.0_rk, (self%rM * Ialg))
   
   ! remineralisation rate
   rem = self%rrem * Idet    
   
   ! grazing term
   IA_gr = max(0.0_rk, self%rGr * Ialg)

   ! nutrient limitation factors
   up_Ino3 = Ino3 / (self%rNO3 + Ino3) !Ino3 + Inh4/(self%rNO3 + (Ino3 + Inh4))
   up_Ipho = Ipho/(self%rPO4 + Ipho)
   up_Isil = Isil/(self%rSi + Isil)

   ! light limitation
   blight = max(1 - exp(-(self%aaIA * bio_par) / self%Pm), 0.0_rk)

   ! temperature and salt dependence
   Tdep = exp(self%TctrlIA * ice_top_temp)
   Sdep = exp(-(2.16 - 8.30 * (10**-5) * (br_salt ** 2.11) - 0.55 * log(br_salt)) ** 2)  

   ! ice algae production and nutrient uptake
   IA_prod = Tdep * Sdep * min(blight, up_Ino3) 
   Prod = self%muIA * Ialg * IA_prod
 
   ! nutrients entrapment and flushing
   entrap = max(0.0_rk,ice_dh) * max(0.0_rk,  Wno3 - Ino3) 
   flush  = min(0.0_rk,ice_dh) * Ino3

   flux_nit = entrap + flush  

   ! nutrients  diffusion
   hv = (self%v / abs(self%Ui - U0)) * (self%Cio ** (-0.5))
   diff_nit = (self%D * hv) * ((Wno3 - Ino3)/ hbio)  

   pno3 = self%vn / (self%vn + Inh4)
   rhs_Init = - Prod + rem + diff_nit + flux_nit !+ ent_Init
   !rhs_Iamm = (-self%muIA * IA_prod * IA * (1-pno3)) - (self%rnit * nh4) + ent_nh4


   _SET_SURFACE_ODE_(self%id_Ialg, Prod - IA_loss)!- (IA_R * Ialg) - (IA_M * Ialg))
   _SET_SURFACE_ODE_(self%id_Idet, IA_loss - rem )
   _SET_SURFACE_ODE_(self%id_Ino3, rhs_Init )
   _SET_SURFACE_ODE_(self%id_Inh4, rhs_Iamm)

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_primprod, Prod)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_limit, IA_prod)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_BAL, diff_nit)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_flush, flux_nit)
  

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface


   end module fabm_hzg_icealgea

