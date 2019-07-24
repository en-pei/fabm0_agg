! Copyright 2019 Helmholtz-Zentrum Geesthacht
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Richard Hofmeister

#include "fabm_driver.h"
#define AGG_WO_CHL
#define LESS_IFS

!  hzg_agg --- aggregation model for aggregation of lithogenic suspended matter
!              and detritus

   module hzg_agg

   use fabm_types
   use fabm_expressions

   private

   public type_hzg_agg, hzg_agg_init, hzg_agg_do  !
   public hzg_agg_get_vertical_movement

   type,extends(type_base_model) :: type_hzg_agg
!     Variable identifiers
      type (type_state_variable_id)      :: id_aggorg,id_agglpm
      type (type_state_variable_id)      :: id_phyn,id_detn,id_lpm
#ifndef AGG_WO_CHL
      type (type_state_variable_id)      :: id_chl,id_aggchl
#endif
      type (type_state_variable_id)      :: id_phyc,id_detc
      type (type_state_variable_id)      :: id_doc
      type (type_diagnostic_variable_id) :: id_aggvol,id_G,id_Breakup,id_ws
      type (type_diagnostic_variable_id) :: id_esd
      type (type_dependency_id) :: id_eps,id_num
      
!     Model parameters
      real(rk) :: dens_lpm
      real(rk) :: specvol_org
      real(rk) :: coagulation_rate
      real(rk) :: deg_chl
      real(rk) :: agg_porosity
      logical  :: use_lpm
      logical  :: use_phyn
      logical  :: use_phyc
#ifndef AGG_WO_CHL
      logical  :: use_chl
#endif
      logical  :: use_detn
      logical  :: use_detc
      logical  :: use_doc
      real(rk) :: dens_org
      real(rk) :: NC_agg
      real(rk) :: NP_agg
      real(rk) :: org2N
      real(rk) :: breakup_factor
      real(rk) :: doc_min
      real(rk) :: doc_mean
      real(rk) :: max_size
      integer  :: size_method
      real(rk) :: tep_remin

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: get_vertical_movement
      procedure :: sinking_velocity
      procedure :: meansize

   end type

   contains



   subroutine initialize(self,configunit)

   implicit none

   class (type_hzg_agg),intent(inout),target  :: self
   integer,             intent(in)            :: configunit

   real(rk), parameter :: secs_pr_day = 86400.d0 ! s
   real(rk), parameter :: one_pr_day = 1.0d0/86400.d0 ! 1/s
   real(rk)            :: aggorg_init=0.3        ! µg/l
   real(rk)            :: aggchl_init=0.1        ! µg/l
   real(rk)            :: agglpm_init=0.1        ! mg/l
   character(len=64)   :: phyn_variable=''
   character(len=64)   :: phyc_variable=''
   character(len=64)   :: detn_variable=''
   character(len=64)   :: detc_variable=''
   character(len=64)   :: lpm_variable=''
   character(len=64)   :: chl_variable=''
   character(len=64)   :: doc_variable=''
            

! in [g/mol]
#define _NMASS_ 14.d0
#define _CMASS_ 12.d0

   call self%get_parameter(self%dens_lpm, 'dens_lpm', 'kg/m3', 'density of lithogenic suspended matter', default=2100.0_rk)
   call self%get_parameter(self%dens_org, 'dens_org', 'kg/m3', 'density of organic suspended matter', default=1000.0_rk)
   call self%get_parameter(self%coagulation_rate, 'coagulation_rate', '', 'coagulation rate', default=0.5_rk)
   call self%get_parameter(self%agg_porosity, 'agg_porosity', '1/1', 'aggregates porosity', default=0.98_rk)
   call self%get_parameter(self%deg_chl, 'deg_chl', '1/d', 'degradation rate of chl', default=0.5_rk, scale_factor=one_pr_day)

   ! dependencies
   call self%get_parameter(lpm_variable, 'lpm_variable', '', 'variable name of lithogenic suspended matter', default='')
   self%use_lpm = lpm_variable.ne.''
   call self%get_parameter(phyn_variable, 'phyn_variable', '', 'variable name of phytoplankton', default='')
   self%use_phyn = phyn_variable.ne.''
   call self%get_parameter(phyc_variable, 'phyc_variable', '', 'variable name of phytoplankton N', default='')
   self%use_phyc = phyc_variable.ne.''
#ifndef AGG_WO_CHL
   call self%get_parameter(chl_variable, 'chl_variable', '', 'variable name of chl', default='')
   self%use_chl = self%chl_variable.ne.''
#endif
   call self%get_parameter(detn_variable, 'detn_variable', '', 'variable name of detritus N', default='')
   self%use_detn = detn_variable.ne.''
   call self%get_parameter(detc_variable, 'detc_variable', '', 'variable name of detritus C', default='')
   self%use_detc = detc_variable.ne.''
   call self%get_parameter(doc_variable, 'doc_variable', '', 'variable name of DOC', default='')
   self%use_doc = doc_variable.ne.''

   self%NC_agg = 16.0_rk/106.0_rk  ! Redfield
   self%NP_agg = 16.0_rk           ! Redfield
   self%org2N  = 16.0_rk/3550.0_rk ! g->molN

   call self%get_parameter(self%doc_min, 'doc_min', 'mmolC/m3', 'minimum doc concentration', default=1.0_rk)
   call self%get_parameter(self%doc_mean,'doc_mean', 'mmolC/m3', 'mean doc concentration', default=100.0_rk)
   ! minimim factor for coagulation rate is doc_min/opt_doc

   call self%get_parameter(self%breakup_factor, 'breakup_factor', 's**0.5/m**2', 'breakup factor', default=12068.0_rk) ! from Maerz.etal2010 [s**0.5/m**2]
   call self%get_parameter(self%max_size, 'max_size', 'm', 'maximum aggregate size', default=250.0e-6_rk)
   call self%get_parameter(self%size_method, 'size_method', '-', 'size method: 1: sigmoid, 2: Xu et al. (2008) steady-state, 3: Xu ea. (2008) Bale experiment,  4: Winterwerp (1998)', default=2)
   call self%get_parameter(self%tep_remin, 'tep_remin', '1/d', 'remineralization rate of TEP', default=0.0_rk, scale_factor=one_pr_day)

   ! Register state variables
   call self%register_state_variable(self%id_aggorg,'aggorg','mg/m**3', &
                          'concentration of biomass in aggregates',    &
                          aggorg_init,minimum=_ZERO_, &
                          no_river_dilution=.true.)

   call self%register_state_variable(self%id_agglpm,'agglpm', &
                          'mg/l','concentration of LPM in aggregates', &
                          agglpm_init,minimum=_ZERO_, &
                          no_river_dilution=.true.)

#ifndef AGG_WO_CHL
   if (self%use_chl) &
      call self%register_state_variable(self%id_aggchl,'aggchl', &
                          'mg/m**3','chlorophyll in aggregates', &
                          agglpm_init,minimum=_ZERO_, &
                          no_river_dilution=.true.)
#endif

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_G,'G','1/s','turbulent shear',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_breakup,'Breakup','1/d','breakup rate',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_aggvol,'Vol_agg','m**3/m**3','relative Volume of aggregates',time_treatment=time_treatment_last)
    call self%register_diagnostic_variable(self%id_ws,'ws','m/s','sinking velocity',time_treatment=time_treatment_last)
    call self%register_diagnostic_variable(self%id_esd,'esd','m','mean ESD',time_treatment=time_treatment_last)

   ! Register conserved quantities
   
   ! Register dependencies
   if (self%use_lpm) call self%register_state_dependency(self%id_lpm,lpm_variable)
   if (self%use_phyn) call self%register_state_dependency(self%id_phyn,phyn_variable)
   if (self%use_phyc) call self%register_state_dependency(self%id_phyn,phyc_variable)
   if (self%use_detn) call self%register_state_dependency(self%id_detn,detn_variable)
   if (self%use_detc) call self%register_state_dependency(self%id_detc,detc_variable)
#ifndef AGG_WO_CHL
   if (self%use_chl) call self%register_state_dependency(self%id_chl,chl_variable)
#endif
   if (self%use_doc) call self%register_state_dependency(self%id_doc,doc_variable)

   call self%register_dependency(self%id_eps, standard_variables%turbulence_dissipation)
   call self%register_dependency(self%id_num, standard_variables%eddy_viscosity)
   
   return

   end subroutine initialize


   subroutine do(self,_ARGUMENTS_DO_)
   implicit none

   class (type_hzg_agg),       intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

   real(rk)                   :: lpm,phyn,phyc,detn,detc
   real(rk)                   :: doc, tep
   real(rk)                   :: G,eps,num_turb
   real(rk)                   :: coagulation, decomposition, breakup
   real(rk)                   :: A1_lpm,A2_lpm, Loss_lpm
   real(rk)                   :: A1_phyn,A2_phyn
   real(rk)                   :: A1_detn,A2_detn
   real(rk)                   :: sms,aggorg,agglpm
#ifndef AGG_WO_CHL
   real(rk)                   :: aggchl,phychl
#endif
   real(rk)                   :: Vol_agg
   real(rk)                   :: num_water=1.1d-3/1025_rk

   _LOOP_BEGIN_
   
   _GET_STATE_(self%id_agglpm,agglpm)
   _GET_STATE_(self%id_aggorg,aggorg)
#ifndef AGG_WO_CHL
   if (self%use_chl) &
     _GET_STATE_(self%id_aggchl,aggchl)
#endif

   ! total volume
   Vol_agg=1.d-3 * (agglpm/self%dens_lpm + 1.d-3*aggorg/self%dens_org)/(_ONE_-self%agg_porosity) ! [m**3/m**3]

   ! Retrieve dependencies
   _GET_DEPENDENCY_(self%id_eps, eps) ! dissipation [m**2/s**3]
   _GET_DEPENDENCY_(self%id_num, num_turb) ! kinematic (turbulent) viscosity [m**2/s]
   G = sqrt(eps/(num_turb + num_water)) ! turbulent shear

   breakup = self%breakup_factor * G**1.5d0 * self%meansize(agglpm+1.d-3*aggorg,G)**2
   !breakup = self%breakup_factor * G**1.5d0 * self%max_size**2 ! [1/s]
   !breakup_factor is in Xu.etal2008: efficiency*sqrt(viscosity/yield_strength)*2 with
   !yield_strength=10.e-10 and efficiency=?

   !_GET_STATE_(self%id_doc,doc)
   doc=self%doc_min
   tep=0.64d0*doc ! 64% of DOC formed TEP
   coagulation = self%coagulation_rate * G * max(self%doc_min,doc)/(doc+self%doc_mean)

   ! bacterial decomposition of TEP binding
   ! so far, organic content is not hydrolysed
   ! agg_org -> DetN agglpm -> lpm
   decomposition = self%tep_remin ! [1/s]

#ifndef LESS_IFS
   if (self%use_phyn) then
#endif
      _GET_STATE_(self%id_phyn,phyn) !
      A1_phyn=coagulation/self%dens_org * self%org2N *1.d-6 * phyn**2
      A2_phyn=coagulation * Vol_agg * phyn

      sms=A1_phyn + A2_phyn ![mmol/m**3/s]
      _SET_ODE_(self%id_phyn,-sms)
      _SET_ODE_(self%id_aggorg,sms/self%org2N)
      if (self%use_phyc) then
         _SET_ODE_(self%id_phyc,-sms/self%NC_agg)
         ! possibly _GET_STATE_(self%id_phyc,phyc)
         ! possibly better: _SET_ODE_(self%id_phyc,-coagulation*Vol_agg*phyc)
         !_SET_ODE_(self%id_aggorg,_CMASS_/self%NC_agg*sms)
      endif

#ifndef AGG_WO_CHL      
      if (self%use_chl) then
         _GET_STATE_(self%id_chl,phychl)
         _SET_ODE_(self%id_chl,phychl/phyn*(-sms))
         _SET_ODE_(self%id_aggchl,phychl/phyn*sms-self%deg_chl*aggchl)
      endif
#endif
#ifndef LESS_IFS
   end if

   if (self%use_detn) then
#endif
      _GET_STATE_(self%id_detn,detn) !
      A1_detn=coagulation/self%dens_org *self%org2N * 1.d-6 * detn**2
      A2_detn=coagulation * Vol_agg * detn

      sms= A1_detn + A2_detn - (decomposition + breakup)*aggorg*self%org2N ![mmolN/m3/s]
      _SET_ODE_(self%id_detn,-sms)
      _SET_ODE_(self%id_aggorg,sms/self%org2N)
      if (self%use_detc) then
         _SET_ODE_(self%id_detc,_ONE_/self%NC_agg*(-sms))
      endif
#ifndef LESS_IFS
   end if
#endif

   if (self%use_lpm) then
      _GET_STATE_(self%id_lpm,lpm) !
      A1_lpm=coagulation/self%dens_lpm * 1.d-3 * lpm**2
      A2_lpm=coagulation * Vol_agg * lpm
      Loss_lpm =  (decomposition + breakup)*agglpm ![g/m**3/s]

      _SET_ODE_(self%id_lpm,Loss_lpm - A1_lpm - A2_lpm)
      _SET_ODE_(self%id_agglpm,A1_lpm + A2_lpm - Loss_lpm)
   endif

   _SET_DIAGNOSTIC_(self%id_G,G)
   _SET_DIAGNOSTIC_(self%id_breakup,breakup)
   _SET_DIAGNOSTIC_(self%id_aggvol,Vol_agg)
   _SET_DIAGNOSTIC_(self%id_ws,self%sinking_velocity(aggorg,agglpm,G))
   _SET_DIAGNOSTIC_(self%id_esd,self%meansize(aggorg*1.d-3+agglpm,G))

   _LOOP_END_

   end subroutine do




   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)

   implicit none

   class (type_hzg_agg), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

   real(rk)                :: aggorg, agglpm, agg_mass, Vol_agg
   real(rk)                :: rho_part,rho_water,visc
   real(rk)                :: G,eps,num_turb
   real(rk)                :: ws
   real(rk)                :: num_water=1.1d-3/1025_rk
   
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_agglpm,agglpm)
   _GET_STATE_(self%id_aggorg,aggorg)

   ! Retrieve dependencies
   _GET_DEPENDENCY_(self%id_eps, eps) ! dissipation [m**2/s**3]
   _GET_DEPENDENCY_(self%id_num, num_turb) ! kinematic (turbulent) viscosity [m**2/s]
   G = sqrt(eps/(num_turb + num_water)) ! turbulent shear

   ws = self%sinking_velocity(aggorg,agglpm,G)

   _SET_VERTICAL_MOVEMENT_(self%id_agglpm,ws)
   _SET_VERTICAL_MOVEMENT_(self%id_aggorg,ws)
#ifndef AGG_WO_CHL
   if (self%use_chl) &
     _SET_VERTICAL_MOVEMENT_(self%id_aggchl,ws)
#endif
   _SET_DIAGNOSTIC_(self%id_ws,ws)

   _LOOP_END_

   end subroutine get_vertical_movement


   real(rk) function sinking_velocity(self,aggorg,agglpm,G)
   implicit none
   class(type_hzg_agg)        :: self
   real(rk), intent(in)       :: aggorg
   real(rk), intent(in)       :: agglpm
   real(rk), intent(in)       :: G
   real(rk)                   :: rho_part,rho_water,visc
   real(rk)                   :: agg_mass,Vol_agg

   rho_water = 1025.d0 ! [kg/m**3]
   visc = 1.1d-3 ! dynamic viscosity for about 17 degC water [kg/(m*s)]
   Vol_agg=1.d-3*(agglpm/self%dens_lpm + 1.d-3*aggorg/self%dens_org)/(_ONE_-self%agg_porosity)
   agg_mass = 1.d-3*aggorg+agglpm
   rho_part = (_ONE_-self%agg_porosity)*1.d-3*agg_mass/Vol_agg + self%agg_porosity*rho_water

   !Stokes law:
   sinking_velocity = 2.d0*(rho_part - rho_water)*9.81d0/(9.d0*visc) * \
         (self%meansize(agg_mass,G)/2.d0)**2

   end function sinking_velocity


   real(rk) function meansize(self,agg_mass,G)
   implicit none
   class(type_hzg_agg)        :: self
   real(rk), intent(in)       :: agg_mass
   real(rk), intent(in)       :: G
   real(rk)                   :: modesize,sigma
   real(rk), parameter        :: minsize= 50.d-6 ! [m] minimum size of arregates
   real(rk), parameter        :: k_size = 50.d0  ! [g/m**3] half-saturation constant for size distribution

   if (self%size_method == 1) then
     ! Get mean size from log-normal distribution:
     !   most probable (mode) size is approaching self%max_size for high SPM
     modesize = minsize+(self%max_size - minsize)*agg_mass/(agg_mass+k_size)
   else if (self%size_method == 2) then
     ! from 1d experiments in Xu.etal2008
     ! equilibrium D50 value depending on agg_mass [g/l] and G
     ! resulting from steady-state experiment
     modesize=0.0001_rk + 0.0004_rk * 1.d-3*agg_mass/sqrt(G)
   else if (self%size_method == 3) then
     ! from 1d experiments in Xu.etal2008
     ! equilibrium D50 value depending on agg_mass [g/l] and G
     ! resulting from the tank experiment
     modesize=0.0001_rk + 0.0003_rk*1.d-3*agg_mass/G
   else if (self%size_method == 4) then
     ! Winterwerp et al (1998)
     modesize = 4.e-6 + 4.e-3*1.e-3*agg_mass/sqrt(G)
   end if
   meansize = min(modesize,self%max_size)
#if 0
   modesize = min(modesize,self%max_size)
   !   width of distribution increases with SPM up to sigma=0.75
   sigma = 0.75d0*agg_mass/(agg_mass+k_size)
   meansize = modesize*exp(1.5d0 * sigma**2) 
#endif

   end function meansize

   end module hzg_agg

