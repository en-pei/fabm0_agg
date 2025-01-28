! Copyright 2025 Helmholtz-Zentrum Hereon
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Enpei Li, Kai Wirtz, Richard Hofmeister

#include "fabm_driver.h"
#define AGG_WO_CHL
#define LESS_IFS

!  hzg_agg --- aggregation model without bottom layer. With options of dynamical size and diagnostic size calculation.

   module hzg_agg

   use fabm_types
   use fabm_expressions
   use fabm_standard_variables, only : type_standard_variable_set

   private

   public type_hzg_agg, hzg_agg_init, hzg_agg_do  !
   public hzg_agg_get_vertical_movement

   type,extends(type_base_model) :: type_hzg_agg
!     Variable identifiers
      type (type_state_variable_id)      :: id_aggorg,id_agglpm,id_aggtep
      type (type_state_variable_id)      :: id_phyn,id_detn,id_lpm
#ifndef AGG_WO_CHL
      type (type_state_variable_id)      :: id_chl,id_aggchl

#endif
      type (type_state_variable_id)      :: id_phyc,id_detc
      type (type_state_variable_id)      :: id_doc
      type (type_diagnostic_variable_id) :: id_esd,id_rho_part,id_fractal_dimension
      type (type_diagnostic_variable_id) :: id_aggvol,id_G,id_Breakup,id_ws,id_aggmass,id_coagulationphy,id_coagulationlpm,id_coagulationdet,id_diffsetlpm,id_sinkinglpm,id_resuspensionlpm,id_Dsize 
      type (type_dependency_id) :: id_eps,id_num
      type (type_state_variable_id)      :: id_dD,id_Xsize 
      
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
      real(rk) :: onoff
      real(rk) :: fractal_dimension1
      real(rk) :: min_size
      real(rk) :: const_ws
      real(rk) :: ks
      real(rk) :: kd
      real(rk) :: pc, pb, ps, pd, pws, kws, pgc, pgb

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
   real(rk)            :: aggorg_init=0.3*0.001  ! g m-3      ! µg/l
   real(rk)            :: aggtep_init=1  ! g m-3    
   real(rk)            :: aggchl_init=0.1        ! µg/l !check
   real(rk)            :: agglpm_init=0.1        ! g m-3		!mg/l
   real(rk)            :: dD_init=0.00001      	! m/s !added
   real(rk)            :: Xsize_init=0.00001      	! m !added
   character(len=64)   :: phyn_variable=''
   character(len=64)   :: phyc_variable=''
   character(len=64)   :: detn_variable=''
   character(len=64)   :: detc_variable=''
   character(len=64)   :: lpm_variable=''
   character(len=64)   :: chl_variable=''
   character(len=64)   :: doc_variable=''
!   character(len=64)   :: aggmass_variable=''	 !added

            
   type(type_standard_variable_set) :: standard_variable_set

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

   self%NC_agg = 1000*16.0_rk/106.0_rk  ! Redfield
   self%NP_agg = 1000*16.0_rk           ! Redfield
   self%org2N  = 1000*16.0_rk/3550.0_rk ! g->mmolN 
!   self%min_size=1.d-5 

   call self%get_parameter(self%doc_min, 'doc_min', 'mmolC/m3', 'minimum doc concentration', default=1.0_rk)
   call self%get_parameter(self%doc_mean,'doc_mean', 'mmolC/m3', 'mean doc concentration', default=100.0_rk)
   ! minimim factor for coagulation rate is doc_min/opt_doc

   call self%get_parameter(self%breakup_factor, 'breakup_factor', 's**0.5/m**2', 'breakup factor', default=12068.0_rk) ! from Maerz.etal2010 [s**0.5/m**2]
   call self%get_parameter(self%max_size, 'max_size', 'm', 'maximum aggregate size', default=400.0e-6_rk)!250.0e-6_rk)
   call self%get_parameter(self%size_method, 'size_method', '-', 'size method: 1: sigmoid, 2: Xu et al. (2008) steady-state, 3: Xu ea. (2008) Bale experiment,  4: Winterwerp (1998)', default=2)
   call self%get_parameter(self%tep_remin, 'tep_remin', '1/d', 'remineralization rate of TEP', default=0.0_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%onoff, 'onoff', 'dimensionless', 'on or off of dD/dt term', default=0.0_rk)
   call self%get_parameter(self%fractal_dimension1, 'fractal_dimension1', 'dimensionless', 'fractal dimension from input', default=2.2_rk)
   call self%get_parameter(self%min_size, 'min_size', 'm', 'maximum aggregate size', default=1.0e-6_rk)
   call self%get_parameter(self%const_ws, 'const_ws', 'm/s', 'constant sinking velocity', default=0.0_rk)
   call self%get_parameter(self%ks, 'ks', '/', 'constant for dD_s sinking', default=1.0_rk)
   call self%get_parameter(self%kd, 'kd', '/', 'constant for dD_d differential settling', default=1.0_rk)

   call self%get_parameter(self%pc, 'pc', '/', 'power for dD_c coagulation', default=4.0_rk)
   call self%get_parameter(self%pb, 'pb', '/', 'power for dD_b breakup', default=2.0_rk)
   call self%get_parameter(self%ps, 'ps', '/', 'power for dD_s sinking', default=1.54_rk)
   call self%get_parameter(self%pd, 'pd', '/', 'power for dD_d differential settling', default=5.5_rk)
   call self%get_parameter(self%pws, 'pws', '/', 'power for sinking velocity', default=1.54_rk)
   call self%get_parameter(self%kws, 'kws', '/', 'constant for sinking velocity', default=347.5602_rk)
   call self%get_parameter(self%pgc, 'pgc', '/', 'power for G in coagulation', default=1.0_rk)
   call self%get_parameter(self%pgb, 'pgb', '/', 'power for G in breakup', default=1.5_rk)

   ! Register state variables
   call self%register_state_variable(self%id_aggorg,'aggorg','g/m**3', & !changed from mg/m**3 to g/m**3
                          'concentration of biomass in aggregates',    &
                          aggorg_init,minimum=_ZERO_, &
                          no_river_dilution=.true.)

   call self%register_state_variable(self%id_aggtep,'aggtep','g/m**3', & !changed from mg/m**3 to g/m**3
                          'tep in aggregates',    &
                          aggtep_init,minimum=_ZERO_, &
                          no_river_dilution=.true.)

   call self%register_state_variable(self%id_agglpm,'agglpm', &    	  !mg/L
                          'g/m**3','concentration of LPM in aggregates', &
                          agglpm_init,minimum=_ZERO_, &
                          no_river_dilution=.true.)

   call self%register_state_variable(self%id_dD,'dD', &
                          'm/s','derivative of D over time', &
                          dD_init,minimum=_ZERO_, &
                          no_river_dilution=.true.) ! for dynamical size

   call self%register_state_variable(self%id_Xsize,'Xsize', &
                          'm*kg/m-3','Diameter*spmc', &
                          Xsize_init,minimum=_ZERO_, &
                          no_river_dilution=.true.) !compound variable


#ifndef AGG_WO_CHL
   if (self%use_chl) &
      call self%register_state_variable(self%id_aggchl,'aggchl', &	!check
                          'mg/m**3','chlorophyll in aggregates', &
                          agglpm_init,minimum=_ZERO_, &
                          no_river_dilution=.true.)
#endif

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_G,'G','1/s','turbulent shear',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_Breakup,'Breakup','1/d','breakup rate',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_aggvol,'Vol_agg','m**3/m**3','relative Volume of aggregates',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_ws,'ws','m/s','sinking velocity',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_esd,'esd','m','mean ESD',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_aggmass,'aggmass','g m-3','mass concentration of aggregates',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_rho_part,'rho_part','kg m-3','density of aggregates',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_coagulationphy,'coagulationphy','1/s','coagulaiton phy',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_coagulationdet,'coagulationdet','1/s','coagulaiton det',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_coagulationlpm,'coagulationlpm','1/s','coagulaiton lpm',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_diffsetlpm,'diffsetlpm','1/s','differential settling caused size change lpm',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_sinkinglpm,'sinkinglpm','1/s','sinking caused size change lpm',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_resuspensionlpm,'resuspensionlpm','1/s','resuspension caused size change lpm',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_fractal_dimension,'fractal_dimension','/','fractal_dimension',time_treatment=time_treatment_last)
   call self%register_diagnostic_variable(self%id_Dsize,'Dsize','m','Dsize',time_treatment=time_treatment_last)


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

   ! The dissipation of the turbulent kinetic energy is usually abbreviated as
   ! eps with greek symbol notation $\epsilon$. Its unit is W kg-1, or,
   ! equivalently m2 s-3.
   if (standard_variable_set%contains('turbulent_kinetic_energy_dissipation')) then 
     !call self%register_dependency(self%id_eps, standard_variables%turbulent_kinetic_energy_dissipation)
     call self%register_dependency(self%id_eps, &
       type_bulk_standard_variable(name='turbulent_kinetic_energy_dissipation'))
   else
     call self%register_dependency(self%id_eps, &
       type_bulk_standard_variable(name='turbulent_kinetic_energy_dissipation',units='W kg-1', &
       cf_names='specific_turbulent_kinetic_energy_dissipation_in_sea_water'))
   endif

   ! The vertical eddy viscosity or momentum diffusivity is usually abbreviated as num with greek
   ! symbol $\nu_m$.  Its unit is m2 s-1 and it was renamed eddy_viscosity => momentum_diffusivity 
   if (standard_variable_set%contains('momentum_diffusivity')) then 
     !call self%register_dependency(self%id_num, standard_variables%momentum_diffusivity)
     call self%register_dependency(self%id_num, &
       type_bulk_standard_variable(name='momentum_diffusivity'))
   else
     call self%register_dependency(self%id_num, &
       type_bulk_standard_variable(name='momentum_diffusivity',units='m2 s-1', &
       cf_names='ocean_vertical_momentum_diffusivity'))
   endif
   
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
   real(rk)                   :: A1_lpm,A2_lpm, loss_lpm, loss_detn, loss_phyn
   real(rk)                   :: A1_phyn,A2_phyn
   real(rk)                   :: A1_detn,A2_detn
   real(rk)                   :: aggorg,agglpm, coagulation_detn, coagulation_phyn, coagulation_lpm, sms, Dsize,aggtep
#ifndef AGG_WO_CHL
   real(rk)                   :: aggchl,phychl
#endif
   real(rk)                   :: Vol_agg, aggmass, xx=1 !10 xx is used to scale the size all together
   real(rk)                   :: num_water=1.1d-3/1025_rk, rho_part
   real(rk)                   :: dD, Xsize, W, diffset_lpm, sinking_lpm, resuspension_lpm, fractal_dimension
   real(rk)	 	      :: Pi=4.D0*DATAN(1.D0), rho_water = 1025.d0, visc = 1.1d-3 ! dynamic viscosity for about 17 degC water [kg/(m*s)] 



   _LOOP_BEGIN_
   
   _GET_STATE_(self%id_agglpm,agglpm)
   _GET_STATE_(self%id_aggorg,aggorg)
   _GET_STATE_(self%id_aggtep,aggtep)
   _GET_STATE_(self%id_dD,dD)
   _GET_STATE_(self%id_Xsize,Xsize)
   _GET_STATE_(self%id_lpm,lpm)
   Dsize=Xsize/(lpm+1e-6)  ! add a small value in case lpm gets too small
   _SET_DIAGNOSTIC_(self%id_Dsize,Dsize)
!   _GET_STATE_(self%id_aggmass,aggmass) !added
#ifndef AGG_WO_CHL
   if (self%use_chl) &
     _GET_STATE_(self%id_aggchl,aggchl)
   endif ! bug?????? added endif
#endif

   ! total volume
   Vol_agg=1.d-3 * (agglpm/self%dens_lpm + aggorg/self%dens_org)/(_ONE_-self%agg_porosity) ! [m**3/m**3]
   aggmass= aggorg + agglpm !added here

   ! Retrieve dependencies
   _GET_DEPENDENCY_(self%id_eps, eps) ! dissipation [m**2/s**3]
   _GET_DEPENDENCY_(self%id_num, num_turb) ! kinematic (turbulent) viscosity [m**2/s]
   G = sqrt(eps/(num_turb + num_water)) ! turbulent shear 
   sms= (G+1e-6)/7.16755 !7.167541 !varying nf because of G

   fractal_dimension=self%fractal_dimension1

   if (self%onoff==0) then
   	breakup= self%breakup_factor*G**self%pgb*(Dsize*xx-self%min_size)**(3-fractal_dimension)*(Dsize*xx)**self%pb !varying power of Dsize

   else if (self%onoff==1) then   !negative sign added in loss_lpm
   	breakup = self%breakup_factor * G**1.5d0 * (self%meansize(aggmass,G,doc, lpm, agglpm,aggorg))**2! self%meansize(agglpm+aggorg,G)**2      
   end if


   !_GET_STATE_(self%id_doc,doc)
   doc=self%doc_min
   tep=0.64d0*doc ! 64% of DOC formed TEP
   coagulation = self%coagulation_rate * max(self%doc_min,doc)/(doc+self%doc_mean)

   ! bacterial decomposition of TEP binding
   decomposition = self%tep_remin 

#ifndef LESS_IFS
   if (self%use_phyn) then
#endif
      _GET_STATE_(self%id_phyn,phyn) !

      coagulation_phyn	= coagulation * (phyn *1.d-3/(self%org2N *self%dens_org)+Vol_agg) * phyn *G  !mmol m**-3 s-1
      loss_phyn	=	(decomposition + breakup) * aggorg * self%org2N *0.9	!90%back		!mmolN m**-3 s-1

      _SET_ODE_(self%id_phyn,-coagulation_phyn + loss_phyn) !-sms !should be in mol m**-3 s-1  !*1.d-3
      _SET_ODE_(self%id_aggorg,(coagulation_phyn - loss_phyn)/self%org2N) !sms/self%org2N  !smaller aggregation rate for phyn ! g m-3 s-1 


       _SET_DIAGNOSTIC_(self%id_coagulationphy,coagulation_phyn) !mmol m-3 s-1
      if (self%use_phyc) then			!check
         _SET_ODE_(self%id_phyc,-sms/self%NC_agg)
      endif

#ifndef AGG_WO_CHL      			
      if (self%use_chl) then			!check
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
      coagulation_detn	=	coagulation*G * (detn *1.d-3/(self%org2N *self%dens_org)+Vol_agg) * detn !*10	!mmolN m**-3 s-1
          _SET_DIAGNOSTIC_(self%id_coagulationdet,coagulation_detn)
      loss_detn	=	(decomposition + breakup) * aggorg * self%org2N	* 0.1	!10%		!mmolN m**-3 s-1

      _SET_ODE_(self%id_detn, - coagulation_detn + loss_detn)		!-sms !mmolN m**-3 s-1
      _SET_ODE_(self%id_aggorg, (coagulation_detn - loss_detn)/self%org2N)	!sms/self%org2N !g m-3 s-1
      if (self%use_detc) then							!check
         _SET_ODE_(self%id_detc,_ONE_/self%NC_agg*(-sms))
      endif
#ifndef LESS_IFS
   end if
#endif

   if (self%use_lpm) then
     _GET_STATE_(self%id_lpm,lpm) !
     _GET_STATE_(self%id_aggtep,aggtep) !
     loss_lpm =  (decomposition + breakup)*agglpm ![g/m**3/s]
         if (self%onoff==0) then  
            coagulation_lpm = coagulation * aggtep * (G**self%pgc) * 1.d-3 * (lpm**1.0) * (Dsize * xx)**(self%pc - fractal_dimension) !exponential decay change of coagulation

         else if (self%onoff==1) then  
            coagulation_lpm = coagulation * (lpm**2*1.d-3/self%dens_lpm + Vol_agg* lpm) *G			!g m**-3 s-1 
            !coagulation_lpm = coagulation*G*(Vol_agg*lpm) 			!g m**-3 s-1    !only coagulates with existing aggregates 

         end if
      _SET_ODE_(self%id_lpm, (loss_lpm - coagulation_lpm)*self%onoff) !- A1_lpm - A2_lpm
      _SET_ODE_(self%id_agglpm, (-loss_lpm + coagulation_lpm)*self%onoff) !A1_lpm + A2_lpm 
      _SET_ODE_(self%id_aggtep, - decomposition * aggtep) !exponential decay

   
   endif

   if (self%onoff==0) then
      _GET_STATE_(self%id_lpm,lpm) !added
      W=- self%sinking_velocity(Dsize*xx,aggorg,agglpm)/((Dsize*xx)**self%ps) 
      diffset_lpm=self%kd* coagulation*2*Pi*(exp(1.)-1)/(exp(3.)*fractal_dimension*Pi/6*self%dens_lpm)*(Dsize*xx)**(self%pd-fractal_dimension)*(lpm*1.d-3)*W !varying power for Dsize
      _SET_DIAGNOSTIC_(self%id_diffsetlpm,diffset_lpm)  
      sinking_lpm = self%ks*self%sinking_velocity(Dsize*xx,aggorg,agglpm)/0.014*(Dsize*xx)**self%ps*(fractal_dimension-1)/(fractal_dimension+1)  !varying power for Dsize

      _SET_DIAGNOSTIC_(self%id_sinkinglpm,sinking_lpm)


      resuspension_lpm=coagulation_lpm - breakup + diffset_lpm + sinking_lpm !dD sum this is the sum all changing term.!
      _SET_DIAGNOSTIC_(self%id_resuspensionlpm,resuspension_lpm)

! size change caused by all the subprocesses as change of counpound variable XSize!
      _SET_ODE_(self%id_Xsize,(coagulation_lpm-breakup+diffset_lpm+sinking_lpm)*(lpm+1e-6)/xx) !weighted by spmc lpm


      _SET_DIAGNOSTIC_(self%id_coagulationphy,coagulation_phyn) 
      _SET_DIAGNOSTIC_(self%id_coagulationlpm,coagulation_lpm)

   else if (self%onoff==1) then
      rho_part =   1.d-3*aggmass/(_ONE_-self%agg_porosity)/Vol_agg
      _SET_DIAGNOSTIC_(self%id_rho_part,rho_part) 
   end if
   



 
   _SET_DIAGNOSTIC_(self%id_G,G)
   _SET_DIAGNOSTIC_(self%id_Breakup,breakup) 
   _SET_DIAGNOSTIC_(self%id_aggvol,Vol_agg)
   _SET_DIAGNOSTIC_(self%id_esd,self%meansize(aggmass,G,doc,lpm, agglpm,aggorg)) 
   _SET_DIAGNOSTIC_(self%id_aggmass,aggmass)
   _SET_DIAGNOSTIC_(self%id_fractal_dimension,fractal_dimension)


   

   _LOOP_END_

   end subroutine do




   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)

   implicit none

   class (type_hzg_agg), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

!   real(rk), intent(in)    :: !
   real(rk)                :: ws,aggorg,agglpm,rho_part,rho_water,aggmass,fractal_dimension,Dsize,G!,visc
   real(rk)                :: Xsize,lpm,eps,num_turb
   real(rk)                :: num_water=1.1d-3/1025_rk, xx=1 !10
   
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_agglpm,agglpm)
   _GET_STATE_(self%id_aggorg,aggorg)
   _GET_STATE_(self%id_lpm,lpm) !added
   _GET_STATE_(self%id_Xsize,Xsize) !added

!   _GET_STATE_(self%id_aggmass,aggmass)

   ! Retrieve dependencies
   _GET_DEPENDENCY_(self%id_eps, eps) ! dissipation [m**2/s**3]
   _GET_DEPENDENCY_(self%id_num, num_turb) ! kinematic (turbulent) viscosity [m**2/s]
   G = sqrt(eps/(num_turb + num_water)) ! turbulent shear  !num_water~mu
   Dsize=Xsize/(lpm+1e-6) !!calculated locally second time
   fractal_dimension =self%fractal_dimension1 !3* (Dsize*1d6)**(-0.076527) !2.2  !!second time
   ws = self%sinking_velocity(Dsize*xx,aggorg,agglpm)
   if (self%onoff==1) then  
      _SET_VERTICAL_MOVEMENT_(self%id_agglpm,ws)
      _SET_VERTICAL_MOVEMENT_(self%id_aggorg,ws)
   else if (self%onoff==0) then
      _SET_VERTICAL_MOVEMENT_(self%id_lpm,ws)
      _SET_VERTICAL_MOVEMENT_(self%id_Xsize,ws)  !added ws for Xsize
   else 
      print*, "unknown size method for vertical movement"

   end if
   _SET_DIAGNOSTIC_(self%id_ws,ws) 

#ifndef AGG_WO_CHL
   if (self%use_chl) &
     _SET_VERTICAL_MOVEMENT_(self%id_aggchl,ws)
#endif


   _LOOP_END_

   end subroutine get_vertical_movement


   real(rk) function sinking_velocity(self,fDsize,aggorg,agglpm)
   implicit none
   class(type_hzg_agg)        :: self
   real(rk), intent(in)       :: aggorg, agglpm
   real(rk), intent(in)       :: fDsize !,G, fractal_dimension
   real(rk)                   :: rho_part,rho_water,visc
   real(rk)                   :: Vol_agg, aggmass,doc,lpm
   real(rk)                   :: dD

   !Stokes law:
   if (self%onoff==0) then !dynamical
   sinking_velocity = - self%kws*(2*fDsize)**self%pws
   else if (self%onoff==1) then !diagnostic
   rho_water = 1025.d0 ! [kg/m**3]
   visc = 1.1d-3 ! dynamic viscosity for about 17 degC water [kg/(m*s)]
   Vol_agg=1.d-3*(agglpm/self%dens_lpm + aggorg/self%dens_org)/(_ONE_-self%agg_porosity)
   aggmass = aggorg + agglpm	!1.d-3*aggorg+agglpm
   rho_part =   1.d-3*aggmass/(_ONE_-self%agg_porosity)/Vol_agg
   sinking_velocity = -0.012 !!!turned off for simple sinking velocity funciton for Dsize !!!!!!!
   !sinking_velocity= -1.d0/18*(rho_part - rho_water)*9.81d0/visc*(self%meansize(aggmass,G,doc,lpm,agglpm,aggorg))**1.54*0.01 !*0.007 ! !!working function
   ! to turn back on the sinking_velocity function here instead of a constant!!!! muted for dynamical size.
   else 
      print*, "unknown size method for calculating sinking velocity"
   end if 
   end function sinking_velocity


   real(rk) function meansize(self,aggmass,G,doc,lpm,agglpm,aggorg)
   implicit none
   class(type_hzg_agg)        :: self
   real(rk), intent(in)       :: aggmass
   real(rk), intent(in)       :: G
   real(rk)                   :: modesize,sigma, doc, lpm, agglpm,aggorg,sms,k
!   real(rk), parameter        :: minsize= 0.0001 !changed to size of pp, was 50.d-6 ! [m] minimum size of arregates
   real(rk), parameter        :: k_size = 50.d0  ! [g/m**3] half-saturation constant for size distribution !need to be checked

   if (self%size_method == 1) then
     ! Get mean size from log-normal distribution:
     !   most probable (mode) size is approaching self%max_size for high SPM
     modesize = self%min_size+(self%max_size - self%min_size)*aggmass/(aggmass+k_size)
   else if (self%size_method == 2) then
     ! from 1d experiments in Xu.etal2008
     ! equilibrium D50 value depending on aggmass [g/l] and G, aggmass[g/m-3]
     ! resulting from steady-state experiment
     modesize=0.0001_rk + 0.0004_rk * 1.d-3*aggmass/sqrt(G)
   else if (self%size_method == 3) then !modified
     ! from 1d experiments in Xu.etal2008
     ! equilibrium D50 value depending on aggmass [g/l] and G
     ! resulting from the tank experiment
     modesize=0.0001_rk + 0.0003_rk*1.d-3*aggmass/G !0.00005_rk + 0.0001_rk*1.d-3*aggmass/G
   else if (self%size_method == 4) then  !current size method
     ! Winterwerp et al (1998)
     !modesize = 70.e-6+4.e-6 + 4.e-3*1.e-3*(agglpm+aggorg)/sqrt(G)  !previously used for spmc 4000
     !modesize = 1.e-6 + 30*4.e-3*1.e-3*(agglpm+aggorg)/sqrt(G) !30
     modesize = 1.e-6 + self%coagulation_rate/self%breakup_factor*1.e-3*(agglpm+aggorg)/sqrt(G) !30
   else if (self%size_method == 5) then
     ! equilibrium of aggregation in this code when only lpm
     modesize = sqrt(0.001* self%coagulation_rate*max(self%doc_min,doc)/(doc+self%doc_mean)/(self%dens_lpm*(_ONE_-self%agg_porosity)*self%breakup_factor*sqrt(G))*(1+(_ONE_-self%agg_porosity)*0.001*lpm/0.001*agglpm*0.001*lpm)) 
   else if (self%size_method == 6) then
     ! equilibrium of aggregation in this code when only aggregates forming aggregates
     modesize = sqrt((1.e-3*aggorg/self%dens_org + 1.e-3*agglpm/self%dens_lpm)*self%coagulation_rate*max(self%doc_min,doc)/(doc+self%doc_mean)*1/((_ONE_-self%agg_porosity)*self%breakup_factor*sqrt(G))) !1.3* 
   else if (self%size_method == 7) then !1e-4*method6
    modesize = 1.e-4*sqrt((1.e-3*aggorg/self%dens_org + 1.e-3*agglpm/self%dens_lpm)*self%coagulation_rate*max(self%doc_min,doc)/(doc+self%doc_mean)*1/(_ONE_-self%agg_porosity)*self%breakup_factor*sqrt(G)) 
   else if (self%size_method == 8) then
     k=2
     modesize = sqrt(self%coagulation_rate*max(self%doc_min,doc)/((doc+self%doc_mean)*self%breakup_factor*sqrt(G))*(1/(_ONE_-self%agg_porosity)+k-1)*(k-1)*(1.e-3*agglpm/self%dens_lpm))
   else if (self%size_method == 9) then !agglpm+lpm=k*agglpm
     k=1+lpm/agglpm 
     modesize = sqrt(self%coagulation_rate*max(self%doc_min,doc)/((doc+self%doc_mean)*self%breakup_factor*sqrt(G))*(1/(_ONE_-self%agg_porosity)+k-1)*(k-1)*(1.e-3*agglpm/self%dens_lpm))
!   else if (self%onoff == 0) then 
!     modesize = Dsize  !????does it work?
   end if
   meansize = min(modesize,self%max_size) !test the min equation
#if 0
   modesize = min(modesize,self%max_size)
   !   width of distribution increases with SPM up to sigma=0.75
   sigma = 0.75d0*aggmass/(aggmass+k_size)
   meansize = modesize*exp(1.5d0 * sigma**2) 
#endif

   end function meansize

   end module hzg_agg

