! Copyright 2017 Helmholtz-Zentrum Geesthacht
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Richard Hofmeister

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_pops --- new POP chemistry model
!
! !INTERFACE:
   module fabm_hzg_pops
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
   public type_hzg_pops

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hzg_pops
!     Model variables

      ! Variable identifiers
      type (type_state_variable_id)         :: id_Cx, id_Cs, id_Cds, id_Cm, id_Cphoto, id_Ca! Concentrations
      type (type_state_variable_id)         :: id_test
      type (type_state_variable_id)         :: id_Fwa  ! Flux atmospheric
      type (type_surface_state_variable_id) :: id_Catm0 ! Concentration atm
      type (type_bottom_state_variable_id)  :: id_sedPCB ! Concentration in sediment
      ! Set dependencies
      type (type_dependency_id)             :: id_temp, id_par, id_PF, id_kphoto
      type (type_dependency_id)             :: id_Ca_D, id_Cda_D, id_Cds_D, id_C_exD
      type (type_dependency_id)             :: id_dz   ! Cell thickness
      type (type_global_dependency_id)      :: id_d    ! Days in the year
      type (type_horizontal_dependency_id)  :: id_swr0, id_u10 ! Surface swr and wind speed
      type (type_dependency_id)             :: id_dom, id_no3, id_det, id_phy, id_oxy !ecosmo parameters
      type (type_horizontal_dependency_id)  :: id_sed1, id_tbs ! sediment parameters
      type (type_horizontal_dependency_id)  :: id_frice !coupling with ice
      ! Identifiers for diagnostic variables [model outputs]
      type (type_diagnostic_variable_id)    :: id_diag !, id_C_ex
      type (type_diagnostic_variable_id)    :: id_Ct  ! total concentrations
      type (type_diagnostic_variable_id)    :: id_uv_tot, id_p_tot, id_PF_28, id_kph  ! light variables
      type (type_horizontal_diagnostic_variable_id):: id_dc, id_Cs_flux_bot, id_Cx_flux_bot !fluxes of PCB
     
!     Model parameters
      real(rk) :: diag
      real(rk) :: test
      real(rk) :: kdeso, kdd, dep
      real(rk) :: Kd, Kdom, Kbio
      real(rk) :: Q, Catm

!     Light parameters
      real, dimension(23)      :: eps_dom = (/ 0.01, 0.01, 0.01, 0.01, &
                   0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, &
                   0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, &
                   0.01, 0.01, 0.01/) ! Extinction coefficients for DOM

      real, dimension(23)      :: eps_no3 = (/ 0.00, 0.00, 0.00, 0.00, &
                   0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                   0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
                   0.00, 0.00, 0.00/) ! Extinction coefficients for NO3
! no3 absorbs UV irradiation below 250 nm, my calculations start from 290 nm

      real, dimension(23)      :: eps_28 = (/ 0.005, 0.005, 0.005, 0.005, &
                   0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, &
                   0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, &
                   0.005, 0.005, 0.005/) ! Extinction coefficients for PCB 28


      contains
!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom
      procedure :: get_light_extinction

   end type type_hzg_pops
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
   class (type_hzg_pops),intent(inout),target  :: self
   integer,              intent(in)            :: configunit

! !LOCAL VARIABLES:


  ! 1. Register model parameters
  ! 2. Register model variables
  ! 3. Register diagnostic variables
  ! 4. Register dependencies

  ! Register model parameters from yaml
   call self%get_parameter(self%Kd,    'Kd', '', 'Kd-DET for PCB153', scale_factor=1.0_rk)
   call self%get_parameter(self%Kbio,  'Kbio', '', 'Kd-phytoplankton for PCB153', scale_factor=1.0_rk)
   call self%get_parameter(self%kdeso, 'kdeso', '', 'rate of desorption from DET', scale_factor=1.0_rk)
   call self%get_parameter(self%Kdom,  'Kdom', '', 'Kd- DOM for PCB153', scale_factor=1.0_rk)
   call self%get_parameter(self%kdd,   'kdd', '', 'rate of desorption from DOM', scale_factor=1.0_rk)
   call self%get_parameter(self%Q,     'Q', '', 'quantum yield', default=0.005_rk, scale_factor=1.0_rk)
   call self%get_parameter(self%dep,   'dep', 'pg/m²/s','atmospheric deposition', scale_factor=1.0_rk)
   call self%get_parameter(self%Catm,  'Catm', 'pg/m³', 'initial atmospheric concentrations', scale_factor=1.0_rk)


  ! Register main model variables
   call self%register_state_variable(self%id_Cx,       'Cx' , 'pg/m3', 'dissolved PCB153 without atm', minimum= 0.0_rk)
   call self%register_state_variable(self%id_Cs,       'Cs' , 'pg/m3', 'sorbed on PM PCB153', vertical_movement = -0.00005787_rk, minimum = 0.0_rk)
   call self%register_state_variable(self%id_Cds,      'Cds' , 'pg/m3', 'sorbed on DOM PCB153', minimum=0.0_rk)
   call self%register_state_variable(self%id_Cphoto,   'Cphoto' , 'pg/m3', 'Photodegrated PCB153', minimum=0.0_rk)
   call self%register_state_variable(self%id_Cm,       'Cm' , 'pg/m3', 'Biodegrated PCB153', minimum=0.0_rk)
   call self%register_state_variable(self%id_test,     'test', '', 'test', minimum=-10000000.0_RK)

  ! Photolysis

   call self%register_diagnostic_variable(self%id_p_tot, 'photon_flux_tot','µE','photon flux', source=source_do_column)
   call self%register_diagnostic_variable(self%id_uv_tot,'radiation_tot','W/m2','UV irradiation', source=source_do_column)
   call self%register_diagnostic_variable(self%id_PF_28, 'PF_28','�E', 'Photon flux for PCB153')
   call self%register_diagnostic_variable(self%id_kph,   'kph','', 'photodegradation constant for PCB153')

  ! Atmosphere

   call self%register_state_variable(self%id_Ca,                'Ca' , 'pg/m3', 'Atmospheric input PCB153', minimum=0.0_rk)
   call self%register_state_variable(self%id_Fwa,               'Fwa', 'pg/m2/s', 'atm flux of PCB153', minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Catm0,     'Catm0', 'pg/m3', 'PCB test atm concentrations', initial_value=200._rk, minimum=0.0_rk)
   call self%register_horizontal_diagnostic_variable(self%id_dc,'dc','pg/m2/s', 'Atmospheric flux of PCB153')

  ! General

   call self%register_diagnostic_variable(self%id_Ct,   'Ct','pg/m3', 'total concentrtion of PCB153')
   call self%register_diagnostic_variable(self%id_diag, 'diag','mg/m3/s', 'pops without det')
 
 
  ! General standard variables

   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_frice,standard_variables%ice_fraction) 
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_swr0,standard_variables%surface_downwelling_shortwave_flux)
   call self%register_dependency(self%id_u10,standard_variables%wind_speed)
   call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
   call self%register_dependency(self%id_d,standard_variables%number_of_days_since_start_of_the_year)

  ! From ECOSMO general

   call self%register_dependency(self%id_det,'det','mgC/m3','detritus')
   call self%register_dependency(self%id_no3,'no3','mgC/m3','nitrate')
   call self%register_dependency(self%id_dom,'dom','mgC/m3','dissolved organic matter')
   call self%register_dependency(self%id_phy,'phy','mgC/m3','phytoplankton')
   call self%register_dependency(self%id_oxy,'oxy','mmolO/m3','oxygen')

  ! Sediment 

   call self%register_dependency(self%id_sed1,'sed1','mgC/m3','pool of sed det')
   call self%register_dependency(self%id_tbs,standard_variables%bottom_stress)
   call self%register_diagnostic_variable(self%id_Cs_flux_bot,'bottom_Cs_flux','mgC/m**2/s', &
         'bottom Cs exchange flux')
   call self%register_diagnostic_variable(self%id_Cx_flux_bot,'bottom_Cx_flux','mgC/m²/s', &
         'bottom Cx exchange flux diffusion')
   call self%register_bottom_state_variable(self%id_sedPCB, 'sedPCB', 'pg/m2', 'Sediment PCB153', initial_value=50._rk)

  ! Inside model dependencies (From POPs)

   call self%register_dependency(self%id_PF,'PF_28','','photon flux of pcb')
   call self%register_dependency(self%id_kphoto,'kph','','conc photo')
   call self%register_dependency(self%id_Ca_D,'Ct','','atmosph conc')
   call self%register_dependency(self%id_Cda_D,'Ct','','Cda of pcb')
   call self%register_dependency(self%id_C_exD,'Ct','','Cds of pcb')

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
   class (type_hzg_pops),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk) :: Ca, Cx, Cs, Cm, Cds, Cphoto
   real(rk) :: par, det, dom, phy, temp 
   real(rk) :: diag, test, Ct
   real(rk) :: Kd, Kdom, Kbio
   real(rk) :: kdeso, kdd
   real(rk) :: kph, PF_28
   real(rk) :: dom_rem, det_rem, ksb, ksp, ksd, POM, DiOM
   real(rk) :: kbact, kds, ksorp, BIOd
   real(rk) :: cdet, clab
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_test,test)
   _GET_(self%id_PF, PF_28)
   _GET_(self%id_kphoto, kph)


   !2. Partitioning calculations
   !2.1 At particulate organic matter (det)
   _GET_(self%id_Cx, Cx)
   _GET_(self%id_Cs, Cs)
   _GET_(self%id_Cds, Cds)
   _GET_(self%id_det, det)
   _GET_(self%id_dom, dom)
   _GET_(self%id_phy, phy)
   _GET_(self%id_temp,temp)

    
   POM = det!*0.001_rk  !mgC/m³ to gC/m³
   DiOM = dom!*0.001_rk

   cdet = Cs/POM !transfer coeffitients
   clab = Cds/DiOM


   if (det .GT. 1000._rk) then
      print*, 'OMG, again'
   endif
  
   ! 2.1 At particulate organic matter (det)

   det_rem =  0.00000002_rk * (1._rk+20._rk*(temp**2/(13._rk**2+temp**2))) !remineralization of det (from ecosmo)
   ksorp = (self%kdeso*self%Kd)/POM !kinetic sorption calculated from Schw & Böhm2016
   
   ! 2.2 At dissolved organic matter (dom)

   dom_rem = 10._rk * 0.00000002_rk * (1._rk+20._rk*(temp**2/(13._rk**2+temp**2)))*0.6_rk !remineralization of dom
   kds = (self%kdd*self%Kdom)/DiOM
   
! ksorp - sorption of det
! kdeso - desorption from det
! kds - sorption of dom
! kdd - desorption from dom
! kph - photolytic degradation (calculated in extinction subroutine)

  ! 2.3 Biodegradation
 
   _GET_(self%id_Cphoto, Cphoto)
   _GET_(self%id_Cm, Cm)
   test =0._rk  
   BIOd =0.001_rk* (dom_rem*dom + det_rem*det)
   kbact = 0.0000001_rk

!   ! 2.3 Calculate a partitioning for atmospheric concentrations

   if ((Cs/Cx) .LT. self%Kd) then          ! implementation of Kd conditions (to
                                           !       adsorb not more then its allowed by def)
       ksorp = (self%kdeso*self%Kd)/POM !kinetic sorption calculated from Schw & Böhm2016
   else
       ksorp = 0._rk 
   endif

   if ((Cds/Cx) .LT. self%Kdom) then
       kds = (self%kdd*self%Kdom)/DiOM
   else
       kds = 0._rk
   endif


   _SET_ODE_(self%id_Cx, -0.001_rk*(ksorp*Cx*POM + kds*Cx*DiOM) + kdeso*cdet*POM + &
                kdd*clab*DiOM - kph*(Cx) + det_rem*cdet*POM + 0.6_rk*dom_rem*clab*DiOM - &
                kbact*BIOd*(Cx))
    _SET_ODE_(self%id_Cds, 0.001_rk*kds*(Cx)*DiOM - kdd*clab*DiOM - &
                0.6_rk*dom_rem*clab*DiOM - 0.6_rk*kph*clab*DiOM)
    _SET_ODE_(self%id_Cs, 0.001_rk*ksorp*(Cx)*POM - kdeso*cdet*POM - det_rem*cdet*POM)
    _SET_ODE_ (self%id_Cm,  kbact*BIOd*(Cx))
    _SET_ODE_ (self%id_Cphoto,  kph*(Cx+ 0.6_rk*clab*DiOM))


    _SET_ODE_ (self%id_test,  ksorp)


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
   class (type_hzg_pops),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk)  :: HEN, Catm, u10, temp, Ctot, Ca, Cx, Cs, Cds, Catm0, dep, frice
   real(rk)  :: kw, ka, kol, dc, Fwa, d

!EOP
!-----------------------------------------------------------------------
!BOC

   _GET_(self%id_temp,temp)
   _GET_(self%id_Ca, Ca)
!   _GET_HORIZONTAL_(self%id_Catm, Catm)
   _GET_HORIZONTAL_(self%id_u10, u10)
   _GET_(self%id_Fwa, Fwa)
   _GET_GLOBAL_(self%id_d, d)
   _GET_HORIZONTAL_(self%id_Catm0, Catm0)
   _GET_HORIZONTAL_(self%id_frice, frice)

    _HORIZONTAL_LOOP_BEGIN_

         Ca = 0._rk !Cx

         HEN = EXP (-7950.0027_rk/ (temp + 273.15_rk) + 22.8517_rk )
   !henrys constant in atm m3/mol. 2.03 in Pa m3/mol

         kw = ((0.245_rk * (u10**2) + 0.061 * u10)/3600._rk) * 0.29972323_rk
         ka = (0.2_rk * u10 + 0.3_rk) * 0.1658317_rk
   !calculated by me from 3 articles (look at 1st notepad)


         dc = (self%Catm/HEN - Cx)  !concept from Odabasi et al 2008. 
         kol = 1._rk/(1._rk/(kw) + 1._rk/(ka*HEN))

         Fwa = (dc*kol/86400._rk) + self%dep    !last coefficient is a deposition of
                                                !PCB from EMEP 2018 pg/m²/s

   if (frice .LT. 1._rk) then
           Fwa = (1._rk-frice)*Fwa
   else
           Fwa = 0._rk
   endif

   _SET_SURFACE_EXCHANGE_(self%id_Cx, Fwa)

   _SET_HORIZONTAL_DIAGNOSTIC_ (self%id_dc, Fwa)

   


   _HORIZONTAL_LOOP_END_

   end subroutine do_surface

!----------------------------------------------------------------------
! !IROUTINE: Bottom fluxes for the test model
!
! !INTERFACE:

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_hzg_pops),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_

   real(rk)  :: Cs, Cx, det!, detfast!, Cs_flux_bot
   real(rk)  :: sedPCB, tbs, sed1, oxy
   real(rk)  :: Rsd, Rds, Rsa, Rsdenit, temp
   real(rk)  :: ksed, kb, ksw, kws, ksd, sid
   integer   :: d


   _HORIZONTAL_LOOP_BEGIN_


     _HORIZONTAL_LOOP_BEGIN_

     _GET_(self%id_Cs, Cs)
 !    _GET_(self%id_det, det)
!     _GET_GLOBAL_(self%id_d, d)
!     _GET_(self%id_sed1, sed1)

     _SET_BOTTOM_EXCHANGE_(self%id_Cs, -0.00005787_rk*Cs) !sinking rate = BioC(23)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Cs_flux_bot, -0.00005787_rk*Cs)


     _HORIZONTAL_LOOP_END_
     
     _HORIZONTAL_LOOP_BEGIN_


     _GET_(self%id_det,det)
     _GET_(self%id_Cs, Cs)
     _GET_(self%id_Cx, Cx)
     _GET_HORIZONTAL_(self%id_sedPCB, sedPCB)
     _GET_HORIZONTAL_(self%id_sed1, sed1)
     _GET_HORIZONTAL_(self%id_tbs,tbs)
     _GET_(self%id_oxy,oxy)
     _GET_(self%id_temp,temp)


!----citical bottom shear stress
        if (tbs .GE. 0.007_rk) then  !0.007 = BioC(34)
          Rsd=0.000025_rk !28935_rk
          Rds=0.0_rk
        else if (tbs .LT. 0.007_rk) then
          Rsd=0.0_rk
          Rds=0.0000035_rk !4051_rk
        end if


!----mineralization dependency
        if (oxy .gt. 0._rk) then  !0.007 = BioC(34)
          Rsa=0.00000001_rk*exp(0.15*temp)*1.0_rk
          Rsdenit=0.0_rk
        else if (tbs .le. 0._rk) then
          Rsa=0.0_rk
          Rsdenit=0.00000001_rk*exp(0.15*temp)*2.0_rk
        end if

!       sid = 86400._rk !seconds in day      
 
       ksd = 3.9e-10 !0.000034_rk/sid !degradation rate in sediment (Jay A. Davis, 2009 SETAC)
       kws = 4.05e-10 !0.000035_rk/sid !water-sediment diffusion (-II-)
       ksw = 1.29e-11 !0.0000012_rk/sid ! sed-water diffution
       kb = 1.1e-10 !0.0000000011575_rk !burial of PCB (constant from ECOSMO or Johannes)
       
       ksed = Rds*Cs + kws*Cx - Rsd*sedPCB - ksw*sedPCB - kb*sedPCB - ksd*sedPCB

       _SET_BOTTOM_ODE_ (self%id_sedPCB, Rds*Cs - Rsd*sedPCB - kb*sedPCB - ksd*sedPCB - 2._rk*Rsa*sedPCB - Rsdenit*sedPCB + kws*Cx - ksw*sedPCB)
       _SET_BOTTOM_EXCHANGE_(self%id_Cs, Rsd*sedPCB - Rds*Cs)
       _SET_BOTTOM_EXCHANGE_(self%id_Cx, ksw*sedPCB - kws*Cx)
       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Cs_flux_bot, Rsd*sedPCB - Rds*Cs)
       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Cx_flux_bot, ksw*sedPCB - kws*Cx)


     _HORIZONTAL_LOOP_END_

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC


!-----------------------------------------------------------------------
!BOP

   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
   class (type_hzg_pops), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_

   real(rk) :: ext1, ext2, ext3 = 0.01, ext4 = 0.01
 !1 - extinction coefficient for uvA (Kd), 2 - uvB, 3 - det, 4 - phypl
 ! (pure seawater) Hargreaves
   real(rk) :: ex, swr0, dz, z, tot_ex, PF_28, kph, frice
!   real(rk) ::swr1 
   real(rk) :: eps_dom(23), eps_28(23), exti(23), irrad_fr(23), eps_no3(23)
   real(rk) :: uv(23), p(23), p_28(23), det, dom, no3, phy, uv_tot, p_tot, test, lambda(23), Ct, Cx, DiOM, POM
!   real(rk) :: detfast
   integer  :: i

   ! Enter spatial loops (if any)
       _GET_HORIZONTAL_(self%id_swr0,swr0)
!       _GET_HORIZONTAL_(self%id_swr1,swr1)
   
 
        irrad_fr(1) = 0.0000015
        irrad_fr(2) = 0.0001109
        irrad_fr(3) = 0.0011055
        irrad_fr(4) = 0.0072166
        irrad_fr(5) = 0.0126972
        irrad_fr(6) = 0.0223206
        irrad_fr(7) = 0.0297674
        irrad_fr(9) = 0.0496409
        irrad_fr(10) = 0.0441531
        irrad_fr(11) = 0.0506105
        irrad_fr(12) = 0.0486977
        irrad_fr(13) = 0.0527212
        irrad_fr(14) = 0.0549309
        irrad_fr(15) = 0.0524046
        irrad_fr(16) = 0.0623909
        irrad_fr(17) = 0.0658141
        irrad_fr(18) = 0.0613157
        irrad_fr(19) = 0.0640530
        irrad_fr(20) = 0.0555575
        irrad_fr(21) = 0.0650886
        irrad_fr(22) = 0.0625821
        irrad_fr(23) = 0.0938600

       z       = 0._rk
       ex = 0._rk  

     _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_dz,dz)       ! Layer height (m)
         _GET_(self%id_det, det)
!         _GET_(self%id_detfast, detfast)
         _GET_(self%id_dom, dom)
         _GET_(self%id_no3, no3)
         _GET_(self%id_phy, phy)
         _GET_(self%id_Cda_D, Ct)
         _GET_(self%id_Cx, Cx)
         _GET_HORIZONTAL_(self%id_frice,frice)

    ! Absorbance by dom, no3 and pcb
    ! Wavelength dependent extinctions
    ! Set depth to centre of layer

    ! ext1=0.027, ext2 = 0.094


      POM = det*0.001_rk
      DiOM = dom*0.001_rk

      tot_ex = 0._rk
      do i = 1, 23
          if (i .LE. 7) then
             ext1 = 0.027_rk
          else
             ext1 = 0.094_rk
          endif
          exti(i) = (ext1 + 0.01*(DiOM) + 0.008*(POM) + 0.004*(phy*0.00008_rk))
          tot_ex = tot_ex + exti(i)
      end do

      z = z + dz/2._rk
      ex = ex + tot_ex*dz/2._rk

    ! Calculate UV radiation
     do i = 1, 23
          lambda(i) = 285._rk + i * 5._rk
          uv(i) = (1-frice)*swr0*1000._rk/1000._rk * 0.04_rk* irrad_fr(i) * EXP(-z/5._rk - ex)
          p(i) = uv(i) * (lambda(i)/10**9) * 0.00836_rk  !photon flux converted from (wl in m!)
          p_28(i) = p(i) * (self%eps_28(i)/exti(i)) !* Cx / tot_ex)
     enddo

     uv_tot = 0._rk
     p_tot = 0._rk
     PF_28 = 0._rk

     do i =1,23 
          uv_tot = uv_tot + uv(i) 
          p_tot = p_tot + p(i)
          PF_28 = PF_28 + p_28(i) 
     enddo


     ! Move to bottom of layer
      z = z + dz/2._rk
      ex = ex + tot_ex*dz/2._rk

       !Here we write to the diagnosticvariable the irradiation at middle of
       !cell



      _SET_DIAGNOSTIC_(self%id_uv_tot, uv_tot)
      _SET_DIAGNOSTIC_(self%id_p_tot, p_tot)
      _SET_DIAGNOSTIC_(self%id_PF_28, PF_28)

      kph = PF_28*self%Q
      
      _SET_DIAGNOSTIC_(self%id_kph, kph)


     _VERTICAL_LOOP_END_


   ! Retrieve current (local) state variable values.
   ! Leave spatial loops (if any)


   end subroutine get_light_extinction
   

   end module fabm_hzg_pops
