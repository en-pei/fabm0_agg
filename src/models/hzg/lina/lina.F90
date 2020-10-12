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
     real(rk)::lina_X
     real(rk)::lina_A
     real(rk)::lina_E
     real(rk)::lina_L
     real(rk)::lina_D
     real(rk)::lina_N
     real(rk)::lina_P

     real(rk)::lina_QAP
     real(rk)::lina_QAN
     real(rk)::lina_QDP
     real(rk)::lina_QDN
     real(rk)::lina_QXP
     real(rk)::lina_QXN 

     !External Input
     real(rk)::Temperature

   end type type_lina_state_var

  type type_lina_rhs
     real(rk):: Biomass_Phytoplankton,Biomass_Aggregates,Mass_Lithogenous,Biomass_Detritus
     real(rk):: DIN,DIP,TEP
     real(rk):: Phyto_N,Phyto_P, Agg_N,Agg_P,Det_N,Det_P
     real(rk):: active_Light
   
     !State Variables need to have a RHS solution
     real(rk)::lina_X
     real(rk)::lina_A
     real(rk)::lina_E
     real(rk)::lina_L
     real(rk)::lina_D
     real(rk)::lina_N
     real(rk)::lina_P
     real(rk)::lina_QAP
     real(rk)::lina_QAN
     real(rk)::lina_QDP
     real(rk)::lina_QDN
     real(rk)::lina_QXP
     real(rk)::lina_QXN

  end type type_lina_rhs

  type type_lina_param
     !Type for holding tunable (non changing) parameters
     real(rk)::lina_Q_starN
     real(rk)::lina_Q_starP
     real(rk)::lina_Q0N
     real(rk)::lina_Q0P
     real(rk)::lina_rhol
     real(rk)::lina_rhow
     real(rk)::lina_alpha
     real(rk)::lina_mu_max
     real(rk)::lina_E_min
     real(rk)::lina_E_max
     real(rk)::lina_zeta
     real(rk)::lina_MI_star
     real(rk)::lina_eN
     real(rk)::lina_eP

  end type type_lina_param

 type,extends(type_base_model),public:: type_hzg_lina

  type(type_state_variable_id) :: id_Biomass_Phytoplankton,id_Biomass_Aggregates,id_Mass_Lithogenous,id_Biomass_Detritus,id_DIN,id_DIP,id_Phyto_N,id_Phyto_P, id_Agg_N,id_Agg_P,id_Det_N,id_Det_P,id_active_Light,TEP

  type (type_dependency_id)   :: id_Temperature, id_PAR,id_Velocity,id_Tide
 
  type (type_diagnostic_variable_id) :: id_Agg_Lith_ratio,id_Agg_size,id_TEP_rate !.......many more to follow
  

  !LINA State Variables
  type(type_state_variable_id)::id_lina_X
  type(type_state_variable_id)::id_lina_A
  type(type_state_variable_id)::id_lina_E
  type(type_state_variable_id)::id_lina_L
  type(type_state_variable_id)::id_lina_D
  type(type_state_variable_id)::id_lina_N
  type(type_state_variable_id)::id_lina_P

  type (type_state_variable_id) ::id_lina_QXN
  type (type_state_variable_id) ::id_lina_QXP
  type (type_state_variable_id) ::id_lina_QDN
  type (type_state_variable_id) ::id_lina_QDP
  type (type_state_variable_id) ::id_lina_QAN
  type (type_state_variable_id) ::id_lina_QAP



  ! Internal states of the LINA plankton compartment

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
  type (type_diagnostic_variable_id) ::id_lina_MI
  type (type_diagnostic_variable_id) ::id_lina_alpha
  type (type_diagnostic_variable_id) ::id_lina_mu_max
  type (type_diagnostic_variable_id) ::id_lina_E
  type (type_diagnostic_variable_id) ::id_lina_E_min
  type (type_diagnostic_variable_id) ::id_lina_E_max
  type (type_diagnostic_variable_id) ::id_lina_C_dot
  type (type_diagnostic_variable_id) ::id_lina_B
  type (type_diagnostic_variable_id) ::id_lina_B_star
  type (type_diagnostic_variable_id) ::id_lina_MI_star
  ! Internal state of the LINA aggregate compartment


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
   real(rk) ::lina_MI
   real(rk) ::lina_alpha
   real(rk) ::lina_mu_max
   real(rk) ::lina_E
   real(rk) ::lina_E_min
   real(rk) ::lina_E_max
   real(rk) ::lina_C_dot
   real(rk) ::lina_B
   real(rk) ::lina_B_star
   real(rk) ::lina_MI_star
   ! Internal state of the LINA aggregate compartment
   real(rk) ::lina_QAN
   real(rk) ::lina_QAP

   dummy_inital_value=-999.99_rk
   dummy_model_switch=.false.

   write(0,*) 'reading parameters'

   call self%get_parameter(self%dummy_parameter,'Dummy_tunable_Parameter', default=-999.99_rk)
   call self%get_parameter(self%dummy_logical,'Dummy_logical_Parameter', default=.true.)

   
   !Register LINA parameters
   call self%get_parameter(self%lina_eN,'external Nitrogen input', default=-999.99_rk)
   call self%get_parameter(self%lina_eP,'external Phosphorus input', default=-999.99_rk)
   !Register LINA interals as diagnostic variables, allowing debugging
    call self%register_diagnostic_variable(self%id_lina_Q_starN, 'QstarN','mol-N mol-C-1', 'Maximum Nitrogen to Carbon ratio',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_lina_QXP, 'QstarP','mol-P mol-C-1', 'Maximum Phosphorous to Carbon ratio', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_lina_muX, 'muX','1/d', 'Phytoplankton Growth Rate', output=output_instantaneous)  
   call self%register_diagnostic_variable(self%id_lina_mX, 'mX','1/d', 'Phytoplankton Mortalility Rate', output=output_instantaneous)  
   call self%register_diagnostic_variable(self%id_lina_wX, 'wX','1/d', 'Phytoplankton Sinking Rate', output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_lina_CX, 'CX','1/d', 'Phytoplankton Coagulation Rate',  output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_lina_gammaN, 'gammaN','mol-N mol-C-1 d-1', 'Nitrogen Uptake Rate', output=output_instantaneous)  
   call self%register_diagnostic_variable(self%id_lina_gammaP, 'gammaP','mol-P mol-C-1 d-1', 'Phosphorus Uptake Rate', output=output_instantaneous)  
   call self%register_diagnostic_variable(self%id_lina_cI, 'cI','1/1', 'Light dependence in growth rate', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_lina_c, 'c','1/1', 'Physiology colimitation in growth rate', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_lina_MI, 'n','1/1', 'Metabolic independence', output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_lina_alpha, 'alpha','1/1', 'Photosysnthesis slope', output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_lina_mu_max, 'mu_max','1/d', 'Maximum Growth Rate', output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_lina_E_min, 'E_min','mol-C m-3', 'minimum EPS production', output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_lina_E_max, 'E_max','mol-C m-3', 'maximum EPS production', output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_lina_C_dot, 'C_dot','1/d', 'Physiology stress for nutrient limitaion', output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_lina_B, 'B','1/1', 'EPS production parameter', output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_lina_B_star, 'B_star','d', 'Other EPS production parameter', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_lina_MI_star, 'n_star','1/1', 'Intrinsic Metabolic independence', output=output_instantaneous) 

   ! Internal state of the LINA aggregate compartment

  

   !Register LINA State variables

   !Inital Values are set via the FABM.yaml, we allow changes by the physical host (no river dilution). 

   call self%register_state_variable(self%id_lina_X,  'Biomass_Phytoplankton','mol-C m-3','Phytoplankton concentration', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_A,  'Biomass_Aggregates','mol-C m-3','Aggregate concentration', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_E,  'Free EPS','mol-C m-3','EPS concentration', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_L,  'Primary lithogenic particles','mol-C m-3','lithogenic concentration', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_D,  'Detritus','mol-C m-3','Detritus concentration', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_N,  'Nitrogen','mol-N m-3','Nitrogen concentration', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_P,  'Phosphorus','mol-P m-3','Phosphorus concentration', minimum=_ZERO_, no_river_dilution=.false. )
   
 
   
   call self%register_state_variable(self%id_lina_QAN, 'QAN','mol-N mol-C-1', 'Aggregate Nitrogen to Carbon ratio', output=output_instantaneous)
   call self%register_state_variable(self%id_lina_QAP, 'QAP','mol-P mol-C-1', 'Aggregate Phosphorous to Carbon ratio', output=output_instantaneous)
   call self%register_state_variable(self%id_lina_QDN, 'QDN','mol-N mol-C-1', 'Detritus Nitrogen to Carbon ratio', output=output_instantaneous)
   call self%register_state_variable(self%id_lina_QDP, 'QDP','mol-P mol-C-1', 'Detritus Phosphorous to Carbon ratio', output=output_instantaneous)
   call self%register_state_variable(self%id_lina_QXN, 'QXN','mol-N mol-C-1', 'Phytoplankton Nitrogen to Carbon ratio', output=output_instantaneous)
   call self%register_state_variable(self%id_lina_QXP, 'QXP','mol-P mol-C-1', 'Phytoplankton Phosphorous to Carbon ratio', output=output_instantaneous)
  
   !Register external dependencies
   call self%register_dependency(self%id_Temperature,standard_variables%temperature)
   !call self%register_dependency(self%id_Biomass_Phytoplankton,'phy','mol-C/m3','mole_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water')


   !..........................................................................................................................................................
   return
  end subroutine initialize




!------------------------------------------------------------------------------------------------------------------------------------------------------------------

subroutine do(self,_ARGUMENTS_DO_)
    !
    ! !INPUT PARAMETERS:
    class (type_hzg_lina),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_
    !
    !LOCAL VARIABLES:
  
   logical:: Debugout=.TRUE.
    type (type_lina_state_var), dimension(1)       :: var
    type (type_lina_rhs)       :: rhsv
    type (type_lina_param) ::param


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
   real(rk) ::lina_MI
   real(rk) ::lina_alpha
   real(rk) ::lina_mu_max
   real(rk) ::lina_E
   real(rk) ::lina_E_min
   real(rk) ::lina_E_max
   real(rk) ::lina_C_dot
   real(rk) ::lina_B
   real(rk) ::lina_B_star
   real(rk) ::lina_MI_star
   ! Internal state of the LINA aggregate compartment
   real(rk) ::lina_QAN
   real(rk) ::lina_QAP



!---------------------------------

   if (Debugout) then
     write(*,*) 'Starting DO routine'
   else
     write(*,*) 'debugging turned off'
   endif

!---------------------------------


!-----------------------------------
!Enter potential spatial loops

   !_FABM_LOOP_BEGIN_
   _ARGUMENTS_DO_
   _LOOP_BEGIN_

!----------------------------------
!get inital values for external and state variables


   _GET_(self%id_lina_X, var%lina_X)
   _GET_(self%id_lina_A, var%lina_A)
   _GET_(self%id_lina_E, var%lina_E)
   _GET_(self%id_lina_L, var%lina_L)
   _GET_(self%id_lina_D, var%lina_D)
   _GET_(self%id_lina_N, var%lina_N)
   _GET_(self%id_lina_P, var%lina_P)


   _GET_(self%id_lina_QAN,var%lina_QAN)
   _GET_(self%id_lina_QAP,var%lina_QAP)
   _GET_(self%id_lina_QDN,var%lina_QDN)
   _GET_(self%id_lina_QDP,var%lina_QDP)
   _GET_(self%id_lina_QXN,var%lina_QXN)
   _GET_(self%id_lina_QXP,var%lina_QXN)

   _GET_(self%id_Temperature, var%Temperature)  !
!-------------------------------------------------

!Copy Variables into private copy
!lina_A=var%lina_A

!---------------------------------------------
!intermediate results......
!
!
   qN=(self%lina_QXN-param%lina_Q0N)/(param%lina_Q_starN-param%lina_Q0N)              !eqn 10
   qP=(self%lina_QXP-param%lina_Q0P)/(param%lina_Q_starP-param%lina_Q0P)  !eqn 10
   gammaN=param%gamma_starN*((var%lina_N*var%lina_AN)/(var%gamma_starN + (var%lina_N*var%lina_AN)))*(1.0_rk-qN) !eqn 11
   gammaP=param%gamma_starP*((var%lina_P*var%lina_AP)/(var%gamma_starP + (var%lina_P*var%lina_AP)))*(1.0_rk-qP) !eqn 12
   R=zeta&gammaN !eqn 14
   cI=1-exp(-param%alpha*var%PAR/param%mu_max) !eqn 15
   muX=cI*c*var%mu_max-R !eqn 13
   MI=param%lina_MI_star*(1+qN) ! eqn 21
   c=lina_colimitation(qN,qP,c,cn) !colimitation implemented as function !eqn 17,18
    cdot= !eqn 20
   eta=param%lina_E_min+(param%lina_E_max-param%lina_E_in)*(1.0_rk+ 0.5_rk * tanh(param%lina_B_star*cdot-param%lina_B)) !eqn 16
 
   rhox =var%lina_rho_starx*var%lina_dx^(-var%lina_arho) * (1.0_rk-(1.0_rk- var%lina_dx^(-var%lina_arho))*cI*c)!eqn 22
   wx=1.0_rk/(18.0_rk*var%lina_muw)*(var%lina_rhox-param%lina_rhow)*var%lina_g*var%lina_dx^2.0_rk !eqn 23



!-------------------------------------------------
!Calculate RHS of the State equations these rates of change are defined by FABM to be per second

   rhsv%lina_X= (self%lina_muX - self%lina_mx -self%lina_wx -self%lina_CX) * var%lina_X !eqn 1
   rhsv%lina_A= lina_C_tot - (self%lina_wa + self%lina_ma + self%lina_kB)m * var%lina_A !eqn 2
   !(E and L ,D  are solved by the external AGG model)
   !rhsv%lina_E= lina_eta * var%lina_X - lina_h * var%lina_A !eqn 3
   !rhsv%lina_L= lina_rL - (self%lina_wl +self%lina_CL)*var%lina_L + self%lina_kB*lina_psi * var%lina_A !eqn 4 
   !rhsv%lina_D= self%lina_mx * var%lina_X - (self%lina_mD + self%lina_wD +self%lina_CD) * var%lina_D + self%lina_kB (1- lina_psi) * var%lina_A !eqn 5
   rhsv%lina_N= - self%lina_gammaN * var%lina_X + self%lina_mD * var%lina_D * var%lina_QDN + self%lina_ma * var%lina_A * var%lina_QAN +param%lina_eN !eqn 6
   rhsv%lina_P= - self%lina_gammaP * var%lina_X + self%lina_mD * var%lina_D * var%lina_QDP + self%lina_ma * var%lina_A * var%lina_QAP +param%lina_eP !eqn 7

   rhsv%lina_QXN= lina_gammaN- lina_muX * var%lina_QXN !eqn 8
   rhsv%lina_QXP= lina_gammaP- lina_muX * var%lina_QXP !eqn 9

   rhsv%lina_QDN= self%lina_mx *var%lina_X * var%lina_QXN - ( self%lina_mD + self%lina_wD + self%lina_CD) * var%lina_D * var%lina_QDN + self%lina_kB (1- lina_psi) * var%lina_A * var%lina_QAN !eqn 39

   rhsv%lina_QDP= self%lina_mx *var%lina_X * var%lina_QXP - ( self%lina_mD + self%lina_wD + self%lina_CD) * var%lina_D * var%lina_QDP + self%lina_kB (1- lina_psi) * var%lina_A * var%lina_QAP !eqn 40

   rhsv%lina_psi= var%lina_CL*var%lina_L/var%lina_A-var%lina_kB*var%lina_psi-var%lina_psi/var%lina_A* rhsv%lina_A !eqn 24
   rhsv%lina_QAN= self%lina_Cx *var%lina_X * var%lina_QXN + self%lina_CD* var%lina_D * var%lina_QDN - (self%lina_wa + self%lina_ma + self%lina_kB) * var%lina_A * var%lina_QAN !eqn 25
    rhsv%lina_QAP= self%lina_Cx *var%lina_X * var%lina_QXP + self%lina_CD* var%lina_D * var%lina_QDP - (self%lina_wa + self%lina_ma + self%lina_kB) * var%lina_A * var%lina_QAP !eqn 26

!---------------------------------------------------

!-----------------------------------------------output of diagnostic
   if Debugout then
    write(*,*) 'diagnostic variables'
    _SET_DIAGNOSTIC_(self%id_lina_QXN, lina_QXN)
    _SET_DIAGNOSTIC_(self%id_lina_QXP, lina_QXP)
    _SET_DIAGNOSTIC_(self%id_lina_Q_starN,lina_Q_starN)    
    _SET_DIAGNOSTIC_(self%id_lina_Q_starP,lina_Q_starP) 
    _SET_DIAGNOSTIC_(self%id_lina_muX, lina_muX)                   !Phytoplankton Growth rate
    _SET_DIAGNOSTIC_(self%id_lina_mX, lina_mX)
    _SET_DIAGNOSTIC_(self%id_lina_wX, lina_wX)
    _SET_DIAGNOSTIC_(self%id_lina_CX,lina_CX)
    _SET_DIAGNOSTIC_(self%id_lina_gammaN,lina_gammaN)
    _SET_DIAGNOSTIC_(self%id_lina_gammaP,lina_gammaP)
    _SET_DIAGNOSTIC_(self%id_lina_cI,lina_cI)
    _SET_DIAGNOSTIC_(self%id_lina_c,lina_c)
    _SET_DIAGNOSTIC_(self%id_lina_MI,lina_MI)
    _SET_DIAGNOSTIC_(self%id_lina_alpha,lina_alpha)
    _SET_DIAGNOSTIC_(self%id_lina_mu_max,lina_mu_max)
    _SET_DIAGNOSTIC_(self%id_lina_E_min,lina_E_min)
    _SET_DIAGNOSTIC_(self%id_lina_E_max,lina_E_max)
    _SET_DIAGNOSTIC_(self%id_lina_C_dot,lina_C_dot)
    _SET_DIAGNOSTIC_(self%id_lina_B,lina_B)
    _SET_DIAGNOSTIC_(self%id_lina_B_star,lina_B_star)
    _SET_DIAGNOSTIC_(self%id_lina_MI_star,lina_MI_star)
    
    _SET_DIAGNOSTIC_(self%id_lina_QAN,lina_QAN)
    _SET_DIAGNOSTIC_(self%id_lina_QAP,lina_QAP)
   endif

!-----------------------------------------calculate the ODE
   if Debugout then write(*,*) 'calculate ODE'
   _ADD_SOURCE_(self%id_lina_X,rhsv%lina_X)
   !_SET_ODE_(self%id_lina_X, rhsv%lina_X)
!--------------------------------------------finalize FABM loops
   !_FABM_LOOP_END_
   _LOOP_END_
  end subroutine do


  function lina_colimitation(qi,qj,n,cn) result(c)
     real(rk), intent(in)::qi,qj,n,cn
     real(rk)  ::c

     c=qi*lina_gn(qj/qi,n)*(1+qi*qj*n+cn)
  end function
  
  function lina_gn(r,n) result(g)
    real(rk),intent(in):: r,n
    real(rk) ::g
    
    g=(r-r^(1+n))/(1-r^(1+n))
  end function

