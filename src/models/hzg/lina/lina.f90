#include "fabm_driver.h"
!#ifdef _FABM_F2003_
!#include fabm_version.h
!#if _FABM_API_VERSION_ < 1
!#  error You need FABM 1.0 or later
!#endif

!---------------------------------------------------------------------------
!NPD-Type model that includes the LIght-Nutrient-Aggreation Feedback
!
!Idea: Enpei Li, Ovidio Garcia, Johannes Timm 2020
!Equations: OG 2020
!Code: JT 2020
!
!
!---------------------------------------------------------------------------

module hzg_lina
  !
  use fabm_types
  use fabm_driver
  implicit none

  public aggregation_kinetics

  private
   type type_lina_state_var
     real(rk):: Biomass_Phytoplankton,Biomass_Aggregates,Mass_Lithogenous,Biomass_Detritus
     real(rk):: DIN,DIP,TEP
     real(rk):: Phyto_N,Phyto_P, Agg_N,Agg_P,Det_N,Det_P
     real(rk):: active_Light

     !State Variables of Lina
     real(rk):: lina_X
     real(rk)::lina_A
     real(rk)::lina_E
     real(rk)::lina_L
     real(rk)::lina_D
     real(rk)::lina_N
     real(rk)::lina_P
   end type type_lina_state_var

  type type_lina_rhs
     real(rk):: Biomass_Phytoplankton,Biomass_Aggregates,Mass_Lithogenous,Biomass_Detritus
     real(rk):: DIN,DIP,TEP
     real(rk):: Phyto_N,Phyto_P, Agg_N,Agg_P,Det_N,Det_P
     real(rk):: active_Light
   
     !State Variables need to have a RHS solution
     real(rk):: lina_X
     real(rk)::lina_A
     real(rk)::lina_E
     real(rk)::lina_L
     real(rk)::lina_D
     real(rk)::lina_N
     real(rk)::lina_P

  end type type_lina_rhs

 type,extends(type_base_model),public:: type_hzg_lina

  type(type_state_variable_id) :: id_Biomass_Phytoplankton,id_Biomass_Aggregates,id_Mass_Lithogenous,id_Biomass_Detritus,id_DIN,id_DIP,id_Phyto_N,id_Phyto_P, id_Agg_N,id_Agg_P,id_Det_N,id_Det_P,id_active_Light,TEP

  type (type_dependency_id)   :: id_Temperature, id_PAR,id_Velocity,id_Tide
 
  type (type_diagnostic_variable_id) :: id_Agg_Lith_ratio,id_Agg_size,id_TEP_rate !.......many more to follow
  
! Internal states of the LINA plankton compartment
  type (type_diagnostic_variable_id) ::id_lina_QXN
  type (type_diagnostic_variable_id) ::id_lina_QXP
  type (type_diagnostic_variable_id) ::id_lina_Q_starN
  type (type_diagnostic_variable_id) ::id_lina_Q_starP
  type (type_diagnostic_variable_id) ::id_lina_muX
  type (type_diagnostic_variable_id) ::id_lina_mx
  type (type_diagnostic_variable_id) ::id_lina_wx
  type (type_diagnostic_variable_id) ::id_lina_CX
  type (type_diagnostic_variable_id) ::id_lina_gammaN
  type (type_diagnostic_variable_id) ::id_lina_gammaP
  type (type_diagnostic_variable_id) ::id_lina_cI
  type (type_diagnostic_variable_id) ::id_lina_c
  type (type_diagnostic_variable_id) ::id_lina_n
  type (type_diagnostic_variable_id) ::id_lina_alpha
  type (type_diagnostic_variable_id) ::id_lina_mu_max
  type (type_diagnostic_variable_id) ::id_lina_E
  type (type_diagnostic_variable_id) ::id_lina_E_min
  type (type_diagnostic_variable_id) ::id_lina_E_max
  type (type_diagnostic_variable_id) ::id_lina_C_dot
  type (type_diagnostic_variable_id) ::id_lina_B
  type (type_diagnostic_variable_id) ::id_linaB_star
  type (type_diagnostic_variable_id) ::id_lina_n_star
! Internal state of the LINA aggregate compartment
  type (type_diagnostic_variable_id) ::id_lina_QAN
  type (type_diagnostic_variable_id) ::id_lina_QAP


  logical::dummy_logical
  real(rk):: dummy_inital_value,dummy_parameter

 contains
   procedure :: initialize
   procedure :: do

end type type_hzg_lina

real(rk), parameter :: secs_pr_day = 86400.0_rk

contains


!This subroutine initailizes the model and reads the inital values via FABM from the file/coupling specifried in fabm.yaml and input.yaml
 subroutine initialize(self,configunit)
   class(type_hzg_lina),intent(inout),target :: self
   integer , intent(in) ::configunit
!local variables
  real(rk):: dummy_inital_value
  logical:: dummy_model_switch
  real(rk)::dummy_parameter

!local internal LINA variables
 real(rk) ::lina_QXN
  real(rk) ::lina_QXP
  real(rk) ::lina_Q_starN
  real(rk) ::lina_Q_starP
  real(rk) ::lina_muX
  real(rk) ::lina_mx
  real(rk) ::lina_wx
  real(rk) ::lina_CX
  real(rk) ::lina_gammaN
  real(rk) ::lina_gammaP
  real(rk) ::lina_cI
  real(rk) ::lina_c
  real(rk) ::lina_n
  real(rk) ::lina_alpha
  real(rk) ::lina_mu_max
  real(rk) ::lina_E
  real(rk) ::lina_E_min
  real(rk) ::lina_E_max
  real(rk) ::lina_C_dot
  real(rk) ::lina_B
  real(rk) ::linaB_star
  real(rk) ::lina_n_star
! Internal state of the LINA aggregate compartment
  real(rk) :: lina_QAN
  real(rk) ::lina_QAP




  dummy_inital_value=-999.99_rk
  dummy_model_switch=.false.

  write(0,*) 'reading parameters'

  call self%get_parameter(self%dummy_parameter,'Dummy_tunable_Parameter', default=-999.99_rk)
  call self%get_parameter(self%dummy_logical,'Dummy_logical_Parameter', default=.true.)


!Register LINA interals as diagnostic variables, allowing debugging
call self%register_diagnostic_variable(self%id_lina_QXN, 'QXN','mol-N mol-C-1', 'Phytoplankton Nitrogen to Carbon ratio', &
         output=output_instantaneous)
call self%register_diagnostic_variable(self%id_lina_QXP, 'QXP','mol-P mol-C-1', 'Phytoplankton Phosphorous to Carbon ratio', &
         output=output_instantaneous)
call self%register_diagnostic_variable(self%id_lina_Q_starN, 'QstarN','mol-N mol-C-1', 'Maximum Nitrogen to Carbon ratio', &
         output=output_instantaneous)
call self%register_diagnostic_variable(self%id_lina_QXP, 'QstarP','mol-P mol-C-1', 'Maximum Phosphorous to Carbon ratio', &
         output=output_instantaneous)
call self%register_diagnostic_variable(self%id_lina_muX, 'muX','1/d', 'Phytoplankton Growth Rate', &
         output=output_instantaneous)  
call self%register_diagnostic_variable(self%id_lina_mX, 'mX','1/d', 'Phytoplankton Mortalility Rate', &
         output=output_instantaneous)  
call self%register_diagnostic_variable(self%id_lina_wX, 'wX','1/d', 'Phytoplankton Sinking Rate', &
         output=output_instantaneous)   
call self%register_diagnostic_variable(self%id_lina_CX, 'CX','1/d', 'Phytoplankton Coagulation Rate', &
         output=output_instantaneous)   
call self%register_diagnostic_variable(self%id_lina_gammaN, 'gammaN','mol-N mol-C-1 d-1', 'Nitrogen Uptake Rate', &
         output=output_instantaneous)  
call self%register_diagnostic_variable(self%id_lina_gammaP, 'gammaP','mol-P mol-C-1 d-1', 'Phosphorus Uptake Rate', &
         output=output_instantaneous)  
call self%register_diagnostic_variable(self%id_lina_cI, 'cI','1/1', 'Light dependence in growth rate', &
         output=output_instantaneous) 
call self%register_diagnostic_variable(self%id_lina_c, 'c','1/1', 'Physiology colimitation in growth rate', &
         output=output_instantaneous) 
call self%register_diagnostic_variable(self%id_lina_n, 'n','1/1', 'Metabolic independence', &
         output=output_instantaneous)   
 call self%register_diagnostic_variable(self%id_lina_alpha, 'alpha','1/1', 'Photosysnthesis slope', &
         output=output_instantaneous)   
call self%register_diagnostic_variable(self%id_lina_mu_max, 'mu_max','1/d', 'Maximum Growth Rate', &
         output=output_instantaneous)   
call self%register_diagnostic_variable(self%id_lina_E, 'E','mol-C m-3', 'EPS production', &
         output=output_instantaneous)   
call self%register_diagnostic_variable(self%id_lina_E_min, 'E_min','mol-C m-3', 'minimum EPS production', &
         output=output_instantaneous)   
call self%register_diagnostic_variable(self%id_lina_E_max, 'E_max','mol-C m-3', 'maximum EPS production', &
         output=output_instantaneous)   
call self%register_diagnostic_variable(self%id_lina_C_dot, 'C_dot','1/d', 'Physiology stress for nutrient limitaion', &
         output=output_instantaneous)   
call self%register_diagnostic_variable(self%id_lina_B, 'B','1/1', 'EPS production parameter', &
         output=output_instantaneous)   
call self%register_diagnostic_variable(self%id_lina_B_star, 'B_star','d', 'Other EPS production parameter', &
         output=output_instantaneous) 
call self%register_diagnostic_variable(self%id_lina_n_star, 'n_star','1/1', 'Intrinsic Metabolic independence', &
         output=output_instantaneous) 

! Internal state of the LINA aggregate compartment

call self%register_diagnostic_variable(self%id_lina_QAN, 'QAN','mol-N mol-C-1', 'Aggregate Nitrogen to Carbon ratio', &
         output=output_instantaneous)
call self%register_diagnostic_variable(self%id_lina_QAP, 'QAP','mol-P mol-C-1', 'Aggregate Phosphorous to Carbon ratio', &
         output=output_instantaneous)


!Register LINA State variables

 call self%register_state_variable(self%id_Biomass_Phytoplankton,  'Biomass_Phytoplankton','µg-C/L','Biomass of free living phytoplankton in µg-carbon per liter', &
         , minimum=_ZERO_, no_river_dilution=.true. )



