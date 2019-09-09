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
      type (type_state_variable_id)         :: id_Cx_28, id_Cs_28, id_Cda_28, id_Cds_28, id_Cp_28, id_Ca ! Concentrations
      type (type_state_variable_id)         :: id_wtf, id_test
      type (type_state_variable_id)         :: id_ksorp_28, id_kds_28 ! Sorption
      type (type_state_variable_id)         :: id_Fwa_28  ! Flux from water to atm
      ! Set dependencies
      type (type_horizontal_dependency_id)  :: id_swr0, id_u10 ! Surface swr and wind speed
      type (type_dependency_id)             :: id_temp, id_par, id_PF
      type (type_dependency_id)             :: id_Ca_28D, id_Cda_28D, id_Cds_28D, id_C_exD
      type (type_dependency_id)             :: id_dz   ! Cell thickness
      type (type_dependency_id)             :: id_dom, id_no3, id_det, id_phy !ecosmo
      type (type_dependency_id)             :: id_Catm_28 !coupling with atmosphere
      ! Identifiers for diagnostic variables [model outputs]
      type (type_diagnostic_variable_id)    :: id_diag, id_fuck1, id_conc!, id_C_ex 
      type (type_diagnostic_variable_id)    :: id_Ctot_28, id_Ct_28  ! total concentrations
      type (type_diagnostic_variable_id)    :: id_uv_tot, id_p_tot, id_PF_28, id_kph_28  ! light variables
      type (type_horizontal_diagnostic_variable_id):: id_dc !concentr difference between water and atm


!     Model parameters
      real(rk) :: diag
      real(rk) :: test 
      real(rk) :: kdeso_28, kdd_28
      real(rk) :: Kd_28, Kdom_28 
      real(rk) :: Q_28 

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

      real, dimension(23)      :: eps_28 = (/ 0.001, 0.001, 0.001, 0.001, &
                   0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, & 
                   0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, &
                   0.001, 0.001, 0.001/) ! Extinction coefficients for PCB 28

     
      contains
!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom
      procedure :: get_light_extinction
      procedure :: check_state

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




  ! Register model parameters
   call self%get_parameter(self%Kd_28, 'Kd_28', '', 'Kd for PCB 28', scale_factor=1.0_rk)
   call self%get_parameter(self%kdeso_28, 'kdeso_28', '', 'constant rate for desorption for 28', scale_factor=1.0_rk)
   call self%get_parameter(self%Kdom_28, 'Kdom_28', '', 'Kd for PCB 28', scale_factor=1.0_rk)
   call self%get_parameter(self%kdd_28, 'kdd_28', '', 'constant rate for desorption for 28', scale_factor=1.0_rk)
   call self%get_parameter(self%Q_28, 'Q_28', '', 'quantum yield for PCB28', default=0.005_rk, scale_factor=1.0_rk)


  ! Register model variables
   call self%register_state_variable(self%id_Cx_28,    'Cx_28' , 'pg/L', 'PCB28 dissolved', minimum=0.0_rk)
   call self%register_state_variable(self%id_Cs_28,    'Cs_28' , 'pg/g', 'PCB28 sorbed', minimum=0.0_rk)
   call self%register_state_variable(self%id_Cda_28,    'Cda_28' , 'pg/L', 'PCB28 dissolved', initial_value=1.0_rk,  minimum=0.0_rk)
   call self%register_state_variable(self%id_Cds_28,    'Cds_28' , 'pg/g', 'PCB28 sorbed', minimum=0.0_rk)
   call self%register_state_variable(self%id_Cp_28,    'Cp_28' , 'pg/g', 'PCB28 sorbed', minimum=0.0_rk)
   call self%register_state_variable(self%id_Ca,    'Ca' , 'pg/g', 'PCB28 sorbed', minimum=0.0_rk)

   call self%register_state_variable(self%id_wtf,      'wtf' , '', 'what the fuck value', minimum=0.0_rk)
   call self%register_state_variable(self%id_test,     'test', '', 'WTF',minimum=-1000.0_RK)
   call self%register_state_variable(self%id_ksorp_28, 'ksorp_28', '', 'rate of sorption for 28', minimum=0.0_rk)
   call self%register_state_variable(self%id_kds_28, 'kds_28', '', 'rate of sorption at DOM  for 28', minimum=0.0_rk)
   call self%register_state_variable(self%id_Fwa_28, 'Fwa_28', '', 'pops flux wa', minimum=0.0_rk)

   
  ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_p_tot,'photon_flux_tot','???','photon flux', source=source_do_column)
   call self%register_diagnostic_variable(self%id_uv_tot,'radiation_tot','???','UV radiation', source=source_do_column)

   call self%register_diagnostic_variable(self%id_diag, 'diag','mg/m3/s', 'reactive tracer degradation')
   call self%register_diagnostic_variable(self%id_fuck1, 'fuck_it_up','mg/m3/s', 'reactive tracer degradation')
   call self%register_diagnostic_variable(self%id_conc, 'conc','pg/L', 'Final dissolved concentration of PCB 28')
   call self%register_diagnostic_variable(self%id_Ctot_28, 'Ctot_28','mmol/l', 'total concentrtion of 28 PCB')
   call self%register_diagnostic_variable(self%id_Ct_28, 'Ct_28','mmol/l', 'total concentrtion of 28 PCB')
   call self%register_diagnostic_variable(self%id_PF_28, 'PF_28','blabla', 'photon flux for 28 PCB')
   call self%register_diagnostic_variable(self%id_kph_28, 'kph_28','blabla', 'photon flux for 28 PCB')
   call self%register_horizontal_diagnostic_variable(self%id_dc, 'dc','blabla', 'concentr 28')


  ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)  
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_swr0,standard_variables%surface_downwelling_shortwave_flux)
   call self%register_dependency(self%id_u10,standard_variables%wind_speed)
   call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
   call self%register_dependency(self%id_det,'det','mgC/m3','detritus')
   call self%register_dependency(self%id_no3,'no3','mgC/m3','nitrate')
   call self%register_dependency(self%id_dom,'dom','mgC/m3','dissolved organic matter')
   call self%register_dependency(self%id_phy,'phy','mgC/m3','phytoplankton')
   call self%register_dependency(self%id_Catm_28,'Catm_28','pg/m3','atmospheric concentrations of pcb')
   call self%register_dependency(self%id_PF,'PF_28','','photon flux of pcb')
   call self%register_dependency(self%id_Ca_28D,'Ct_28','','Ct pcb')
   call self%register_dependency(self%id_Cda_28D,'Cx_28','','Cda of pcb')
   call self%register_dependency(self%id_Cds_28D,'Cp_28','','Cds of pcb')
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
!   real(rk) :: bde209, furan, bde47
   real(rk) :: Cx_28, Cs_28, par, det, wtf, Cda_28, Cds_28, dom
   real(rk) :: diag, test, Ctot_28, Ct_28, fuck1, conc

   real(rk) :: ksorp_28, kdeso_28, Kd_28, kds_28, kdd_28, Kdom_28
   real(rk) :: kph_28, Cp_28, PF_28, Ca

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_par, par)
   _GET_(self%id_det, det)
   _GET_(self%id_Cx_28, Cx_28)
   _GET_(self%id_Cs_28, Cs_28)
   _GET_(self%id_dom, dom)
   _GET_(self%id_wtf,wtf)
   _GET_(self%id_test,test)
   _GET_(self%id_PF, PF_28)


   !1. Calculation of total pollutant amount (dissolved + adsorbed)
!   det = MAX(det,1.0)
  
   Ctot_28 = Cs_28 + Cx_28 
   _SET_DIAGNOSTIC_(self%id_Ctot_28, Ctot_28)

   !2. Partitioning calculations 
   !2.1 At particulate organic matter (det)
  
   ksorp_28 = self%kdeso_28/self%Kd_28 

   _SET_ODE_(self%id_Cx_28, (kdeso_28*(Ctot_28)) - ((ksorp_28*(det/(10**3))+kdeso_28)*Cx_28))
   _SET_ODE_(self%id_Cs_28, (ksorp_28*(Ctot_28)*(det/(10**3))) - ((ksorp_28*(det/(10**3)) + kdeso_28)*Cs_28))

   !2.2 At dissolved organic matter (dom)


   _GET_(self%id_Ca, Ca)
    Ct_28 = Cx_28 + Ca
   _SET_DIAGNOSTIC_ (self%id_Ct_28, Ct_28 )
   
   _GET_(self%id_Ca_28D, Ct_28)
   _GET_(self%id_Cda_28, Cda_28)
   _GET_(self%id_Cds_28, Cds_28)

   kds_28 = self%kdd_28/self%Kdom_28 
   
   _SET_ODE_(self%id_Cda_28,  (kdd_28*(Ct_28)) - ((kds_28*(dom/(10._rk))+kdd_28)*Cda_28))
   _SET_ODE_(self%id_Cds_28, (kds_28*(Ct_28)*(dom/(10._rk))) - ((kds_28*(dom/10._rk) + kdd_28)*Cds_28))

!   write(0,*) 'temp: ',temp

   !3. Direct photolysis calculations 

   kph_28 = PF_28*self%Q_28*100000._rk*Cda_28
   
   _GET_ (self%id_Cp_28, Cp_28)
   Cp_28 = Cda_28   

   _SET_ODE_(self%id_Cp_28, -kph_28)
   _GET_(self%id_Cds_28D, Cp_28)

   conc = Cda_28 +  (1._rk - Cp_28) 

   _SET_DIAGNOSTIC_(self%id_conc, conc)
   _SET_DIAGNOSTIC_(self%id_fuck1, kph_28)


   

!   write(0,*) 'temp: ',temp
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

   real(rk)  :: H_28, Catm_28, u10, temp, Ctot_28, Ca, Cx_28!, C_ex
   real(rk)  :: kw, ka, kol, dc, Fwa_28, bla


!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,temp)
   _GET_(self%id_Ca, Ca)
   _GET_(self%id_Catm_28, Catm_28)
   _GET_HORIZONTAL_(self%id_u10, u10)
   _GET_(self%id_Fwa_28, Fwa_28)

   Ca = Ctot_28

   H_28 = EXP (-3908.8515_rk/ (temp + 273.15_rk) + 8.9002_rk )
   !henrys constant in atm m3/mol. 2.03 in Pa m3/mol

   kw = ((0.245_rk * (u10**2) + 0.061 * u10)/3600._rk) * 0.34375_rk 
   !calculated by me from 3 articles (look at 1st notepad)
   ka = (0.2_rk * u10 + 0.3_rk) * 0.213923_rk 


   dc = 210._rk - (Ctot_28/H_28) !concept from Odabasi et al 2008. Instead of
                                    !210 is Catm
   kol = 1._rk/(kw*H_28) + 1._rk/ka

   Fwa_28 = dc/(kol*1000._rk)

   _SET_SURFACE_EXCHANGE_(self%id_Ca, Fwa_28) 

!   bla = Ctot_28 + Ca

   _SET_HORIZONTAL_DIAGNOSTIC_ (self%id_dc, Fwa_28)

!   diag = Ca
!   _SET_DIAGNOSTIC_ (self%id_diag, diag)


!   _GET_ (self%id_C_ex,C_ex)
   
!   Ctot_28 = Cx_28 + Ca 


!   _SET_DIAGNOSTIC_(self%id_C_ex, Ctot_28)


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
   class (type_hzg_pops),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_

   _HORIZONTAL_LOOP_BEGIN_

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC


   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
   class (type_hzg_pops), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_

   real(rk) :: ext1 = 0.027, ext2 = 0.094, ext3 = 0.01, ext4 = 0.01 
 !1 - extinction coefficient for uvA (Kd), 2 - uvB, 3 - det, 4 - phypl
 ! (pure seawater) Hargreaves
   real(rk) :: swr0, dz, z, tot_ex, PF_28
   real(rk) :: eps_dom(23), eps_28(23), exti(23), irrad_fr(23), eps_no3(23)
   real(rk) :: uv(23), p(23), p_28(23), det, dom, no3, phy, uv_tot, p_tot, test, lambda(23), Cx_28
   integer  :: i

   ! Enter spatial loops (if any)
       _GET_HORIZONTAL_(self%id_swr0,swr0)

        irrad_fr(1) = 0.0000015
        irrad_fr(2) = 0.0001109
        irrad_fr(3) = 0.0011055
        irrad_fr(4) = 0.0072166
        irrad_fr(5) = 0.0126972
        irrad_fr(6) = 0.0223206
        irrad_fr(7) = 0.0297674
        irrad_fr(8) = 0.0429593
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

       p_tot   = 0._rk
       z       = 0._rk
       tot_ex = 0._rk
       PF_28 = 0._rk
     _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_dz,dz)       ! Layer height (m)
         _GET_(self%id_det, det)
         _GET_(self%id_dom, dom)
         _GET_(self%id_no3, no3)
         _GET_(self%id_phy, phy)
         _GET_(self%id_Cda_28D, Cx_28)
    ! Absorbance by dom, no3 and pcb
    ! Wavelength dependent extinctions               
    ! Set depth to centre of layer

      z = z + dz/2._rk
      do i = 1, 23
!      print*, self%ext1, self%ext2, self%ext3, self%ext4, self%eps_no3(i), self%eps_dom(i)

          if (i .LE. 7) then
             exti(i) = (ext1 + self%eps_dom(i)*(dom/12011._rk) + ext3*(det/12011._rk) &
             + ext4*(phy/12011._rk) + self%eps_no3(i)*(no3/14007._rk))*dz/2._rk
          else
             exti(i) = (ext2 + self%eps_dom(i)*(dom/12011._rk) + ext3*(det/12011._rk) &
             + ext4*(phy/12011._rk) + self%eps_no3(i)*(no3/14007._rk))*dz/2._rk
          endif
           tot_ex = exti(1) + exti(2) + exti(3) + exti(4) + exti(5)&
                 + exti(6) + exti(7) + exti(8) + exti(9) + exti(10)&
                 + exti(11) + exti(12) + exti(13) + exti(14) + exti(15)&
                 + exti(16) + exti(17) + exti(18) + exti(19) + exti(20)&
                 + exti(21) + exti(22) + exti(23)
     end do
     
    ! Calculate UV radiation
      do i = 1, 23
          lambda(i) = 285._rk + i * 5._rk
          uv(i) = swr0 *0.04_rk * irrad_fr(i) * EXP(-z/5._rk - exti(i))
          p(i) = uv(i) * (lambda(i)/10**9) * 0.00836_rk  !photon flux converted from (wl in m!)
          p_28(i) = p(i) * (self%eps_28(i) * Cx_28 / exti(i))
      end do
     

       uv_tot = uv(1) + uv(2) + uv(3) + uv(4) + uv(5) + uv(6) + uv(7)&
         + uv(8) + uv(9) + uv(10) + uv(11) + uv(12) + uv(13) + uv(14)&
         + uv(15) + uv(16) + uv(17) + uv(18) + uv(19) + uv(20) + uv(21)&
         + uv(22) + uv(23)

 
       p_tot = p(1) + p(2) + p(3) + p(4) + p(5) + p(6) + p(7) + p(8)&
         + p(9) + p(10) + p(11) + p(12) + p(13) + p(14) + p(15) + p(16)&
         + p(17) + p(18) + p(19) + p(20) + p(21) + p(22) + p(23)
      
 
       PF_28 = p_28(1) + p_28(2) + p_28(3) + p_28(4) + p_28(5) + p_28(6)&
         + p_28(7) + p_28(8) + p_28(9) + p_28(10) + p_28(11) + p_28(12)&
         + p_28(13) + p_28(14) + p_28(15) + p_28(16) + p_28(17) + p_28(18)&
         + p_28(19) + p_28(20) + p_28(21) + p_28(22) + p_28(23)

    
        ! Move to bottom of layer
      z = z + dz/2._rk
      do i = 1, 23
          if (i .LE. 7) then
             exti(i) = (ext1 + self%eps_dom(i)*(dom/12011._rk) + ext3*(det/12011._rk) &
             + ext4*(phy/12011._rk) + self%eps_no3(i)*(no3/14007._rk))*dz/2._rk
          else
             exti(i) = (ext2 + self%eps_dom(i)*(dom/12011._rk) + ext3*(det/12011._rk) &
             + ext4*(phy/12011._rk) + self%eps_no3(i)*(no3/14007._rk))*dz/2._rk
          end if
           tot_ex = tot_ex + exti(i)
       end do

       !Here we write to the diagnosticvariable the irradiation at middle of
       !cell

      _SET_DIAGNOSTIC_(self%id_uv_tot, uv_tot)
      _SET_DIAGNOSTIC_(self%id_p_tot, p_tot)
      _SET_DIAGNOSTIC_(self%id_PF_28, PF_28)


     _VERTICAL_LOOP_END_


   ! Retrieve current (local) state variable values.

   !_GET_(self%id_det, det)

   !_SET_EXTINCTION_( 1.0*det )

   ! Leave spatial loops (if any)
  

   end subroutine get_light_extinction
   


   subroutine check_state(self, _ARGUMENTS_CHECK_STATE_)
   class (type_hzg_pops), intent(in) :: self
   _DECLARE_ARGUMENTS_CHECK_STATE_
   real(rk) :: Ct_28
   _LOOP_BEGIN_

   _GET_(self%id_C_exD, Ct_28)
   _SET_(self%id_Cda_28, 0.99*Ct_28)
   _SET_(self%id_Cds_28, 0.01*Ct_28)

   _LOOP_END_

   end subroutine check_state

   end module fabm_hzg_pops



