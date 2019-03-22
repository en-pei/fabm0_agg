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
      type (type_horizontal_diagnostic_variable_id)  :: id_primprod, id_l_limit, id_n_limit, id_denit, id_BAL, id_entrap, id_flush, id_sal_limit, id_temp_limit
      type (type_horizontal_diagnostic_variable_id)  :: id_remineralisation, id_mortality, id_grazing, id_lightIA, id_pno3, id_hice, id_Walg
      type (type_dependency_id)                      :: id_temp, id_U0
      type (type_horizontal_dependency_id)           :: id_air_temp, id_par
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
      real(rk) :: rgrow
      real(rk) :: rmelt
      real(rk) :: dice
      real(rk) :: dwater 
      logical  :: couple_ice

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
   redf(1)  = 6.625_rk      !molC/molN
   redf(2)  = 106.0_rk      !C_P
   redf(3)  = 6.625_rk      !C_SiO
   redf(4)  = 16.0_rk       !N_P
   redf(5)  = 1.0_rk        !N_SiO
   redf(6)  = 12.01_rk      !gC/molC
   redf(7)  = 44.6608009_rk !O2mm_ml
   redf(8)  = 14.007_rk     !N_Nmg
   redf(9)  = 30.97_rk      !P_Pmg
   redf(10) = 28.09_rk      !Si_Simg    

!------------------------------------------------------------------------------------

   call self%get_parameter(self%muIA,    'muIA',    '1/day',                  'ice algae maximum growth rate',                  default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%TctrlIA, 'TctrlIA', '1/degC',                 'temperature growth dependency',                  default=0.001_rk)
   call self%get_parameter(self%rM,      'rM',      '1/day',                  'ice algae mortality rate',                       default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rGr,     'rGr',     '1/day',                  'ice algae grazing rate',                         default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%aaIA,    'aaIA',    'mgC/mgchla/h/Einst/m/s', 'photosynthesis ef-cy',                           default=0.5_rk,   scale_factor=one_pr_hour)
   call self%get_parameter(self%Pm,      'Pm',      'mgC/mgchla/h',           'maximum photosynhtetic rate',                    default=0.5_rk,   scale_factor=one_pr_hour)
   call self%get_parameter(self%rNO3,    'rNO3',    'mmolN/m**3',             'NO3 half saturation',                            default=0.5_rk)!,   scale_factor=redf(1)*redf(6))
   call self%get_parameter(self%rPO4,    'rPO4',    'mmolP/m**3',             'PO4 half saturation',                            default=0.5_rk)!,   scale_factor=redf(2)*redf(6))
   call self%get_parameter(self%rSi,     'rSi',     'mmolSi/m**3',            'SiO2 half saturation',                           default=0.5_rk)!,   scale_factor=redf(3)*redf(6))
   call self%get_parameter(self%rnit,    'rnit',    '1/day',                  'nitrification rate',                             default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rrem,    'rrem',    '1/day',                  'remineralisation rate',                          default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%vn,      'vn',      'mmolN/m**3',             'preferential uptake of nitrate half saturation', default=0.5_rk)!,   scale_factor=redf(1)*redf(6))
   call self%get_parameter(self%C0,      'C0',      '',                       'part of incoming radiation absorbed',            default=0.3_rk)
   call self%get_parameter(self%Catt_ice,'Catt_ice','1/m',                    'ice attenuation coefficent',                     default=0.5_rk)                
   call self%get_parameter(self%Catt_sno,'Catt_sno','1/m',                    'snow attenuation coefficent',                    default=0.5_rk)                
   call self%get_parameter(self%ice_alb, 'ice_alb', '',                       'albedo in the ice',                              default=0.5_rk)                
   call self%get_parameter(self%sno_alb, 'sno_alb', '',                       'albedo in the snow',                             default=0.5_rk)                
   call self%get_parameter(self%v,       'v',       'm**2/s',                 'kinematic viscosity of seawater',                default=0.0_rk)                
   call self%get_parameter(self%Ui,      'Ui',      'm/s',                    'horizontal velocity of ice',                     default=0.0_rk)                
   call self%get_parameter(self%Cio,     'Cio',     'm**2/s',                 'drag coefficient',                               default=0.0_rk)                
   call self%get_parameter(self%D,       'D',       'm**2/s',                 'molecular diffusion coefficient',                default=0.0_rk)                
   call self%get_parameter(self%rgrow,   'rgrow',   'm/d',                    'growth rate of ice',                             default=0.0_rk)                
   call self%get_parameter(self%dice,    'dice',    'kg/m**3',                'density od sea ice',                             default=0.0_rk)
   call self%get_parameter(self%dwater,  'dwater',  'kg/m**3',                'density of the sea water',                       default=0.0_rk)  
   call self%get_parameter(self%couple_ice, 'couple_ice', '', 'switch coupling to ice', default=.false.)

   call self%register_surface_state_variable(self%id_Ialg,    'Ialg',   'mgC/m**2',   'ice algae', initial_value=0.1_rk, minimum=0.1_rk)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_Ialg, scale_factor=1.0_rk/redf(1)/redf(6))
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_Ialg, scale_factor=1.0_rk/redf(2)/redf(6))
   call self%add_to_aggregate_variable(standard_variables%total_silicate,  self%id_Ialg, scale_factor=1.0_rk/redf(3)/redf(6))

   call self%register_surface_state_variable(self%id_Idet,    'Idet',   'mgC/m**2',   'detritus',       minimum=0.0_rk)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_Idet, scale_factor=1.0_rk/redf(1)/redf(6))
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_Idet, scale_factor=1.0_rk/redf(2)/redf(6))
   call self%add_to_aggregate_variable(standard_variables%total_silicate,  self%id_Idet, scale_factor=1.0_rk/redf(3)/redf(6))

   call self%register_surface_state_variable(self%id_Ino3,    'Ino3',   'mmolN/m**2',   'nitrate',        minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Inh4,    'Inh4',   'mmolN/m**2',   'nitrate',        minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Ipho,    'Ipho',   'mmolP/m**2',   'phosphate',      minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Isil,    'Isil',   'mmolSi/m**2',   'silicate',      minimum=0.0_rk)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_Ino3)!,scale_factor=1.0_rk/redf(1)/redf(6))
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_Inh4)!,scale_factor=1.0_rk/redf(1)/redf(6))
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_Ipho)!,scale_factor=1.0_rk/redf(2)/redf(6))
   call self%add_to_aggregate_variable(standard_variables%total_silicate,  self%id_Isil)!,scale_factor=1.0_rk/redf(3)/redf(6))
   call self%register_surface_state_variable(self%id_ioxy,    'ioxy',   'mmolo2/m**2','oxygen',         minimum=0.0_rk)
   
   call self%register_horizontal_diagnostic_variable(self%id_primprod,  'primprodIA', 'mgC/m2/s',     'primary production')
   call self%register_horizontal_diagnostic_variable(self%id_remineralisation,  'remIA', 'mgC/m3/s',     'remineralisation')
   call self%register_horizontal_diagnostic_variable(self%id_mortality, 'mortalityIA','mgC/m3/s',     'mortality IA')
   call self%register_horizontal_diagnostic_variable(self%id_grazing,   'grazingIA',  'mgC/m3/s',     'mortality IA by grazing')
   call self%register_horizontal_diagnostic_variable(self%id_l_limit,   'l_limitIA',  ''        ,     'light limitation factor for ia')
   call self%register_horizontal_diagnostic_variable(self%id_sal_limit,   'sal_limitIA',  ''        ,     'salinity limitation factor for ia')
   call self%register_horizontal_diagnostic_variable(self%id_temp_limit,   'temp_limitIA',  ''        ,     'temperature limitation factor for ia')
   call self%register_horizontal_diagnostic_variable(self%id_n_limit,   'n_limitIA',  ''        ,     'nutrient limitation factor for ia')
   call self%register_horizontal_diagnostic_variable(self%id_lightIA,   'lightIA',    ''        ,     'light available for IA')
   call self%register_horizontal_diagnostic_variable(self%id_hice,      'hice',           'm',            'thickness of ice pack')
   call self%register_horizontal_diagnostic_variable(self%id_BAL,       'BAL',        'm',            'thickness of BAL')
   call self%register_horizontal_diagnostic_variable(self%id_entrap,    'entrap',     'mmolN/m**3',   'nutrient entrap')
   call self%register_horizontal_diagnostic_variable(self%id_flush,     'flush',      'mmolN/m**3',   'nutrient release')
   call self%register_horizontal_diagnostic_variable(self%id_pno3,      'pno3',       '',             '')
   call self%register_horizontal_diagnostic_variable(self%id_denit,     'denit_ice',  'mmolN/m**3/s', 'denitrification rate')
   call self%register_horizontal_diagnostic_variable(self%id_Walg,     'phyto_wat',  'mmolN/m**3/s', 'phytoplankton come from water')

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_horizontal_dependency(self%id_par, standard_variables%surface_downwelling_shortwave_flux)
   call self%register_horizontal_dependency( self%id_air_temp, standard_variables%surface_temperature)


   call self%register_dependency(self%id_h,       'icethickness', 'm',    'ice thickness')
   !call self%register_dependency(self%id_hs,      'snowthickness','m',    'snow thickness')
   !call self%register_dependency(self%id_hbio,    'balthickness', 'm',    'BAL thickness')
   !call self%register_dependency(self%id_icetemp, 'icetemp',      'degC', 'ice temperature')
   !call self%register_dependency(self%id_icesalt, 'icesalt',      'PSU',  'ice salinity')
   !call self%register_dependency(self%id_brsalt,  'brsalt',       'PSU',  'ice brine salinity')
   !call self%register_dependency(self%id_balpar,  'balpar',       '',     'PAR at the BAL')
   call self%register_dependency(self%id_icedh,   'icedh',        'm/s',  'rate of change of ice thickness')

   call self%register_state_dependency(self%id_Wno3, 'Wno3', 'mmolN/m**3',  'nitrate come from surface water column')
   call self%register_state_dependency(self%id_Wnh4, 'Wnh4', 'mmolN/m**3',  'ammonium come from surface water column')
   call self%register_state_dependency(self%id_Wpho, 'Wpho', 'mmolP/m**3',  'phosphate come from surface water column')
   call self%register_state_dependency(self%id_Wsil, 'Wsil', 'mmolSi/m**3', 'silicate come from surface water column')
   call self%register_state_dependency(self%id_Woxy, 'Woxy', 'mmolO2/m**3', 'oxygen come from surface water column')
   call self%register_state_dependency(self%id_Wdia, 'Wdia', 'mgC/m**3',    'diatom come from surface water column')
   call self%register_state_dependency(self%id_Wdet, 'Wdet', 'mgC/m**3',    'detritus come from surface water column')

   return

   end subroutine initialize

!
! Surface fluxes for the icealgea model

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_hzg_icealgea),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: Ino3, Inh4, Ipho, Isil, Ioxy
   real(rk) :: Wno3, Wnh4, Wpho, Wsil, Wdia, Wdet 
   real(rk) :: ice_top_temp, ice_salt, br_salt
   real(rk) :: Ialg, Idet
   real(rk) :: up_Ino3, up_Ipho, up_Isil
   real(rk) :: flux_nit, ent_Inh4, ent_Ipho, ent_Isil, flush
   real(rk) :: pno3, rhs_Init, rhs_Iamm, rhs_Ipho, rhs_Isil
   real(rk) :: blight, tdep, sdep
   real(rk) :: rem, rem_amn, nit, IA_prod, Prod, IA_loss, IA_gr 
   real(rk) :: icethickness, iceconc, ice_dh
   real(rk) :: snowthickness, hbio
   real(rk) :: alb, par, bio_par
   real(rk) :: hv, diff_nit, U0
   real(rk) :: remmolN,uptMolN
   real(rk) :: molN2gC, molP2gC, molS2gC
   real(rk) :: entrap_no3,flush_no3, flux_no3
   real(rk) :: entrap_nh4,flush_nh4, flux_nh4
   real(rk) :: entrap_pho,flush_pho, flux_pho
   real(rk) :: entrap_sil,flush_sil, flux_sil
   real(rk) :: flush_IA, flush_det
   real(rk) :: entrap_IA, entrap_det
   real(rk) :: air_temp   
 
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_Wno3, Wno3)
   _GET_(self%id_Wnh4, Wnh4)
   _GET_(self%id_Wpho, Wpho)
   _GET_(self%id_Wsil, Wsil)
   _GET_(self%id_Wdia, Wdia)
   _GET_(self%id_Wdet, Wdet)
   _GET_HORIZONTAL_(self%id_air_temp, air_temp)
   _GET_HORIZONTAL_(self%id_par, par) 
   _GET_HORIZONTAL_(self%id_Ino3, Ino3)
   _GET_HORIZONTAL_(self%id_Inh4, Inh4)
   _GET_HORIZONTAL_(self%id_Ipho, Ipho)
   _GET_HORIZONTAL_(self%id_Isil, Isil)
   _GET_HORIZONTAL_(self%id_Ialg, Ialg)
   _GET_HORIZONTAL_(self%id_Idet, Idet)
   _GET_HORIZONTAL_(self%id_h, icethickness)
   !_GET_HORIZONTAL_(self%id_hs, snowthickness)
   _GET_HORIZONTAL_(self%id_hbio, hbio)
   _GET_HORIZONTAL_(self%id_c, iceconc)
   _GET_HORIZONTAL_(self%id_icetemp, ice_top_temp)
   _GET_HORIZONTAL_(self%id_icesalt, ice_salt)
   _GET_HORIZONTAL_(self%id_brsalt, br_salt)
   _GET_HORIZONTAL_(self%id_icedh, ice_dh)
   _GET_HORIZONTAL_(self%id_balpar, bio_par)


   if(icethickness > 0.001_rk) then 
   snowthickness = 0.0_rk   
   hbio = 0.05
   U0   = 0.008

  !ice_dh = ice_dh  *  one_pr_day  
   

   ! molN/gC = 1/(molC/molN * gC/molC) 
   molN2gC=(1/(redf(1)*redf(6)))
   molP2gC=(1/(redf(2)*redf(6)))
   molS2gC=(1/(redf(3)*redf(6)))

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

   ! nitrification rate
   nit = self%rnit * Inh4 
   
   ! nutrient limitation factors
   up_Ino3 = Ino3/(self%rNO3 + Ino3)
   up_Ipho = Ipho/(self%rPO4 + Ipho)
   up_Isil = Isil/(self%rSi + Isil)

   ! light limitation
   blight = max(1 - exp(-(self%aaIA * bio_par) / self%Pm), 0.0_rk)

   ! temperature and salt dependence
   !Tdep = exp(self%TctrlIA * ice_top_temp)
   !Sdep = exp(-(2.16 - 8.30 * (10**-5) * (br_salt ** 2.11) - 0.55 * log(br_salt)) ** 2)  
   Tdep = 1
   Sdep = 1

   ! ice algae production and nutrient uptake
   IA_prod = Tdep * Sdep * min(blight, up_Ino3, up_Ipho, up_Isil) 
   Prod = self%muIA * Ialg * IA_prod

   ! grazing term
  ! IA_gr = max(0.0_rk, (self%rGr * prod))

   ! nutrients entrapment and flushing
   entrap_no3 = molN2gC * (max(0.0_rk, ice_dh) * max(0.0_rk, Wno3)) * 0.7 * (self%dice/self%dwater) 
   flush_no3  = (min(0.0_rk, ice_dh/icethickness) * max(0.0_rk, Ino3))
   flux_no3   = entrap_no3 + flush_no3
 
   entrap_nh4 = molN2gC * (max(0.0_rk, ice_dh) * max(0.0_rk, Wnh4)) * 0.7 * (self%dice/self%dwater)
   flush_nh4  = (min(0.0_rk, ice_dh/icethickness) * max(0.0_rk, Inh4)) 
   flux_nh4   = entrap_nh4 + flush_nh4

   entrap_sil = molS2gC * (max(0.0_rk, ice_dh) * max(0.0_rk, Wsil)) * 0.7 * (self%dice/self%dwater) 
   flush_sil  = (min(0.0_rk, ice_dh/icethickness) * max(0.0_rk, Isil)) 
   flux_sil   = entrap_sil + flush_sil

   entrap_pho = molP2gC * (max(0.0_rk, ice_dh) * max(0.0_rk, Wpho)) * 0.7 * (self%dice/self%dwater) 
   flush_pho  = (min(0.0_rk, ice_dh/icethickness) * max(0.0_rk, Ipho))
   flux_pho   = entrap_pho + flush_pho

   ! 
   flush_IA   = min(0.0_rk, ice_dh/icethickness) * max(0.0_rk, Ialg)
   entrap_IA  = (max(0.0_rk, ice_dh) * max(0.0_rk, Wdia))  * (self%dice/self%dwater) 
   flush_det  = min(0.0_rk, ice_dh/icethickness) * max(0.0_rk, Idet)
   entrap_det = (max(0.0_rk, ice_dh) * max(0.0_rk, Wdet))  * (self%dice/self%dwater) 

   ! nutrients  diffusion
   hv = (self%v / abs(self%Ui - U0)) * (self%Cio ** (-0.5))
   diff_nit = (self%D * hv) * ((Wno3 - Ino3)/ hbio)  
  
   ! nutrients change 
   pno3 = self%vn / (self%vn + Inh4)

   rhs_Init = molN2gC * (- Prod * pno3) + nit 
   rhs_Iamm = molN2gC * ((- Prod * (1-pno3)) + rem) - nit

   !rhs_Iamm = remMolN - uptMolN - nit    
   !write (*,*)'uptMolN,remMolN:',uptMolN,remMolN
   ! mmolN/m2/s = 1/(molC/molN * gC/molC) * (mgC * - ) + mgC/m2/s - mmolN/m2 
   !            = molN/gC*                *  mgC
   rhs_Ipho = (1/(redf(2)*redf(6))) * (- Prod + rem)
   rhs_Isil = (1/(redf(3)*redf(6))) * (- Prod + rem)

   _SET_SURFACE_ODE_(self%id_Ialg, Prod - IA_loss + flush_IA + entrap_IA )
   _SET_SURFACE_ODE_(self%id_Idet, IA_loss  - rem + flush_det + entrap_det)
   _SET_SURFACE_ODE_(self%id_Ino3, rhs_Init + flux_no3)
   _SET_SURFACE_ODE_(self%id_Inh4, rhs_Iamm + flux_nh4)
   _SET_SURFACE_ODE_(self%id_Ipho, rhs_Ipho + flux_pho)
   _SET_SURFACE_ODE_(self%id_Isil, rhs_Isil + flux_sil)
 
   if(self%couple_ice) then 
   _SET_SURFACE_EXCHANGE_(self%id_Wno3, -flux_no3 / molN2gC)
   _SET_SURFACE_EXCHANGE_(self%id_Wnh4, -flux_nh4 / molN2gC)
   _SET_SURFACE_EXCHANGE_(self%id_Wsil, -flux_sil / molS2gC)
   _SET_SURFACE_EXCHANGE_(self%id_Wpho, -flux_pho / molP2gC)
   _SET_SURFACE_EXCHANGE_(self%id_Wdia, -flush_IA - entrap_IA )
   _SET_SURFACE_EXCHANGE_(self%id_Wdet, -flush_det - entrap_det )
   end if 

   end if 

  !  write (*,*)'dh_ice , water no3', ice_dh, Wno3

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_primprod, Prod)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_remineralisation, nit)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_mortality, IA_loss)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_grazing, IA_gr)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_l_limit, blight)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sal_limit, Sdep)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_temp_limit, Tdep)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_n_limit, IA_prod)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_lightIA, par)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_BAL, ice_dh)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pno3, pno3)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_entrap, entrap_IA)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_flush, flush_IA)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_hice, icethickness)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Walg, Wdia)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface


   end module fabm_hzg_icealgea

