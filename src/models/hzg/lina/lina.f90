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
!Equaitons: OG 2020
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
   end type type_lina_state_var

  type type_lina_rhs
  real(rk):: Biomass_Phytoplankton,Biomass_Aggregates,Mass_Lithogenous,Biomass_Detritus
     real(rk):: DIN,DIP,TEP
     real(rk):: Phyto_N,Phyto_P, Agg_N,Agg_P,Det_N,Det_P
   real(rk):: active_Light
  end type type_lina_rhs

 type,extends(type_base_model),public:: type_hzg_lina

  type(type_state_variable_id) :: id_Biomass_Phytoplankton,id_Biomass_Aggregates,id_Mass_Lithogenous,id_Biomass_Detritus,id_DIN,id_DIP,id_Phyto_N,id_Phyto_P, id_Agg_N,id_Agg_P,id_Det_N,id_Det_P,id_active_Light,TEP

  type (type_dependency_id)   :: id_Temperature, id_PAR,id_Velocity,id_Tide
 
  type (type_diagnostic_variable_id) :: id_Agg_Lith_ratio,id_Agg_size,id_TEP_rate !.......many more to follow

  logical::dummy_logical
  real(rk):: dummy_inital_value,dummy_parameter

 contains
   procedure :: initialize
   procedure :: do

end type type_hzg_lina




