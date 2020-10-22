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



  private
   type type_lina_state_var
 
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
     real(rk)::lina_psi

     !External Input
     real(rk)::ma
     real(rk)::mx
     real(rk)::kB
     real(rk)::wa
     real(rk)::dx
     real(rk)::Cx
     real(rk)::md,wd,CD,CL,C_tot
     real(rk)::Temperature
     real(rk)::PAR
   end type type_lina_state_var

   type type_lina_rhs
    
   
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
     real(rk)::lina_psi
   end type type_lina_rhs

   type type_lina_param
     !Type for holding tunable (non changing) parameters


   end type type_lina_param

 type,extends(type_base_model),public:: type_hzg_lina


  type (type_dependency_id)   :: id_Temperature, id_PAR!,id_Velocity,id_Tide
 
  
  type(type_dependency_id)::id_mx
  type(type_dependency_id)::id_ma
  type(type_dependency_id)::id_wa
  type(type_dependency_id)::id_kB
  type(type_dependency_id)::id_dx
  type(type_dependency_id)::id_Cx
  type(type_dependency_id)::id_md
  type(type_dependency_id)::id_wD
  type(type_dependency_id)::id_CD
  type(type_dependency_id)::id_CL
  type(type_dependency_id)::id_C_tot          
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
  type (type_state_variable_id) ::id_lina_psi


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

  type (type_diagnostic_variable_id) ::id_lina_E_min
  type (type_diagnostic_variable_id) ::id_lina_E_max
  type (type_diagnostic_variable_id) ::id_lina_C_dot
!  type (type_diagnostic_variable_id) ::id_lina_B
!  type (type_diagnostic_variable_id) ::id_lina_B_star
  type (type_diagnostic_variable_id) ::id_lina_MI_star
  type (type_diagnostic_variable_id) ::id_lina_QN 
  type (type_diagnostic_variable_id) ::id_lina_QP
  type (type_diagnostic_variable_id) ::id_lina_eta
  type (type_diagnostic_variable_id) ::id_lina_R
  type (type_diagnostic_variable_id) ::id_lina_rhox
 

  ! parameters
     real(rk)::lina_Q_starN
     real(rk)::lina_Q_starP
     real(rk)::lina_Q0N
     real(rk)::lina_Q0P
     real(rk)::lina_rhol
     real(rk)::lina_rhow
     real(rk)::lina_rho_starx
     real(rk)::lina_alpha
     real(rk)::lina_mu_max
     real(rk)::lina_muw
     real(rk)::lina_E_min
     real(rk)::lina_E_max
     real(rk)::lina_zeta
     real(rk)::lina_MI_star
     real(rk)::lina_eN
     real(rk)::lina_eP  
     real(rk)::lina_gamma_starN
     real(rk)::lina_gamma_starP
     real(rk)::lina_AN
     real(rk)::lina_AP
     real(rk)::lina_B
     real(rk)::lina_B_star
     real(rk)::lina_g
     real(rk)::lina_arho
     real(rk)::lina_cn

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
   !real(rk)::dummy_parameter

   !local internal LINA variables


   dummy_inital_value=-999.99_rk
   dummy_model_switch=.false.

   write(0,*) 'reading parameters'

   call self%get_parameter(self%dummy_parameter,'Dummy_tunable_Parameter', default=-999.99_rk)
   call self%get_parameter(self%dummy_logical,'Dummy_logical_Parameter', default=.true.)

   
   !Register LINA parameters
     
   call self%get_parameter(self%lina_Q_starN, 'QstarN','mol-N mol-C-1', 'Maximum Nitrogen to Carbon ratio',  default=-999.99_rk)
   call self%get_parameter(self%lina_Q_starP, 'QstarP','mol-P mol-C-1', 'Maximum Phosphorous to Carbon ratio',default=-999.99_rk)
   call self%get_parameter(self%lina_Q0N, 'Q0N','mol-N mol-C-1', 'Baseline Nitrogen to Carbon ratio',  default=-999.99_rk)
   call self%get_parameter(self%lina_Q0P, 'Q0P','mol-P mol-C-1', 'Baseline Phosphorous to Carbon ratio',default=-999.99_rk)
   call self%get_parameter(self%lina_rhol, 'rohl','Density of Lithogenous',default=-999.99_rk)
   call self%get_parameter(self%lina_rhow, 'rohw','Density of Water',default=1035.00_rk)
   call self%get_parameter(self%lina_rho_starx, 'roh_starx','Base Density of Phytoplankton',default=-999.99_rk)
   call self%get_parameter(self%lina_alpha, 'alpha','1/1', 'Photosysnthesis slope', default=-999.99_rk)   
   call self%get_parameter(self%lina_mu_max, 'mu_max','1/d', 'Maximum Growth Rate', default=-999.99_rk)  
   call self%get_parameter(self%lina_muw, 'muw','', 'water resistance', default=-999.99_rk)   
   call self%get_parameter(self%lina_E_min, 'E_min','mol-C m-3', 'minimum EPS production',default=-999.99_rk)
   call self%get_parameter(self%lina_E_max, 'E_max','mol-C m-3', 'maximum EPS production', default=-999.99_rk)
   call self%get_parameter(self%lina_zeta, 'zeta','1/1', 'Baseline Respiration rate', default=-999.99_rk)  
   call self%get_parameter(self%lina_MI_star, 'MI_star','1/1', 'Intrinsic Metabolic independence', default=-999.99_rk)   
   call self%get_parameter(self%lina_eN,'external Nitrogen input', default=-999.99_rk)
   call self%get_parameter(self%lina_eP,'external Phosphorus input', default=-999.99_rk)   
   call self%get_parameter(self%lina_gamma_starN,'Maximum Nitrogen uptake rate', default=-999.99_rk)  
   call self%get_parameter(self%lina_gamma_starP,'Maximum Phosphorous uptake rate', default=-999.99_rk)  
   call self%get_parameter(self%lina_AN,'Agg Nitrogen? ', default=-999.99_rk)  
   call self%get_parameter(self%lina_AP,'Agg Phosphorous ?', default=-999.99_rk) 
   call self%get_parameter(self%lina_B, 'B','d', 'EPS production parameter', default=-999.99_rk) 
   call self%get_parameter(self%lina_g, 'g','m2s2', 'gravitational acceleration', default=9.81_rk) 
   call self%get_parameter(self%lina_B_star, 'B_star','d', 'Other EPS production parameter', default=-999.99_rk) 
   call self%get_parameter(self%lina_arho, 'aroh','1/1', 'Density Parameter', default=-999.99_rk) 
   call self%get_parameter(self%lina_cn, 'cn','1/1', 'Colimitation Parameter', default=-999.99_rk) 
   !Register LINA interals as diagnostic variables, allowing debugging

   call self%register_diagnostic_variable(self%id_lina_muX, 'muX','1/d', 'Phytoplankton Growth Rate', output=output_instantaneous)  
   call self%register_diagnostic_variable(self%id_lina_mX, 'mX','1/d', 'Phytoplankton Mortalility Rate', output=output_instantaneous)  
   call self%register_diagnostic_variable(self%id_lina_wX, 'wX','1/d', 'Phytoplankton Sinking Rate', output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_lina_CX, 'CX','1/d', 'Phytoplankton Coagulation Rate',  output=output_instantaneous)   
   call self%register_diagnostic_variable(self%id_lina_gammaN, 'gammaN','mol-N mol-C-1 d-1', 'Nitrogen Uptake Rate', output=output_instantaneous)  
   call self%register_diagnostic_variable(self%id_lina_gammaP, 'gammaP','mol-P mol-C-1 d-1', 'Phosphorus Uptake Rate', output=output_instantaneous)  
   call self%register_diagnostic_variable(self%id_lina_cI, 'cI','1/1', 'Light dependence in growth rate', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_lina_c, 'c','1/1', 'Physiology colimitation in growth rate', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_lina_MI, 'MI','1/1', 'Metabolic independence', output=output_instantaneous)   


 
   call self%register_diagnostic_variable(self%id_lina_C_dot, 'C_dot','1/d', 'Physiology stress for nutrient limitaion', output=output_instantaneous)   

  
  

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
   call self%register_state_variable(self%id_lina_QAN, 'QAN','mol-N mol-C-1', 'Aggregate Nitrogen to Carbon ratio', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_QAP, 'QAP','mol-P mol-C-1', 'Aggregate Phosphorous to Carbon ratio', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_QDN, 'QDN','mol-N mol-C-1', 'Detritus Nitrogen to Carbon ratio', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_QDP, 'QDP','mol-P mol-C-1', 'Detritus Phosphorous to Carbon ratio', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_QXN, 'QXN','mol-N mol-C-1', 'Phytoplankton Nitrogen to Carbon ratio', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_QXP, 'QXP','mol-P mol-C-1', 'Phytoplankton Phosphorous to Carbon ratio', minimum=_ZERO_, no_river_dilution=.false. )
   call self%register_state_variable(self%id_lina_psi, 'psi','mol-P mol-C-1', 'psi', minimum=_ZERO_, no_river_dilution=.false. )
  
   !Register external dependencies
!depending on NPZD
   call self%register_dependency(self%id_mx,'mx','d-1','NPZD phytoplankton mortality rate')
!depending on AGG
   call self%register_dependency(self%id_Cx,'Cx','d-1','Phytoplankton Coagulation rate')
   call self%register_dependency(self%id_kB,'kB','1/1','AGG kB')
   call self%register_dependency(self%id_wa,'wa','ms-1','AGG Aggregate sinking speed')
   call self%register_dependency(self%id_ma,'ma','ms-1','AGG Aggregate mortality')
   call self%register_dependency(self%id_dx,'dx','1/1','AGG Phytoplankton Size')
   call self%register_dependency(self%id_md,'md','ms-1','AGG detritus mortality')
   call self%register_dependency(self%id_wD,'wd','ms-1','AGG Detritus sinking speed')
   call self%register_dependency(self%id_CD,'CD','1/1','AGG Detritus Coagulation ')
   call self%register_dependency(self%id_CL,'CL','1/1','AGG Lithogenous Coagulation ')
   call self%register_dependency(self%id_C_tot,'C_tot','1/1','AGG total Coagulation ')
   
!depending on phyics
   call self%register_dependency(self%id_Temperature,standard_variables%temperature)
   call self%register_dependency(self%id_PAR,standard_variables%downwelling_photosynthetic_radiative_flux)
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
    type (type_lina_state_var)      :: var
    type (type_lina_rhs)       :: rhsv
  

   real(rk) ::lina_qN
   real(rk) ::lina_qP
   real(rk) ::lina_wX
   real(rk) ::lina_rhox
   real(rk)::lina_gammaN
   real(rk)::lina_gammaP
   real(rk)::lina_R
   real(rk)::lina_cI
   real(rk)::lina_muX
   real(rk)::lina_MI
   real(rk)::lina_c
   real(rk)::lina_c_dot
   real(rk)::lina_eta

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
!   _ARGUMENTS_DO_
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
 
!Get externals
   _GET_(self%id_Temperature, var%Temperature)  !
   _GET_(self%id_PAR, var%PAR)  !
   _GET_(self%id_kB, var%kB)
   _GET_(self%id_dx, var%dx)
   _GET_(self%id_mx, var%mx)
   _GET_(self%id_ma, var%ma)
   _GET_(self%id_Cx, var%Cx)
   _GET_(self%id_md, var%md)
   _GET_(self%id_CD, var%CD)
   _GET_(self%id_CL, var%CL)
   _GET_(self%id_wd, var%wd)  
   _GET_(self%id_c_tot,var%C_tot)      
!-------------------------------------------------

!Copy Variables into private copy
!lina_A=var%lina_A

!---------------------------------------------
!intermediate results......
!
!
   lina_qN=(var%lina_QXN-self%lina_Q0N)/(self%lina_Q_starN-self%lina_Q0N)              !eqn 10
   lina_qP=(var%lina_QXP-self%lina_Q0P)/(self%lina_Q_starP-self%lina_Q0P)  !eqn 10
   lina_gammaN=self%lina_gamma_starN*((var%lina_N*self%lina_AN)/(self%lina_gamma_starN + (var%lina_N*self%lina_AN)))*(1.0_rk-lina_gn(lina_qN,lina_MI)) !eqn 11
   lina_gammaP=self%lina_gamma_starP*((var%lina_P*self%lina_AP)/(self%lina_gamma_starP + (var%lina_P*self%lina_AP)))*(1.0_rk-lina_gn(lina_qP,lina_MI)) !eqn 12
   lina_R=self%lina_zeta*lina_gammaN !eqn 14
   lina_cI=1-exp(-self%lina_alpha*var%PAR/self%lina_mu_max) !eqn 15
   lina_MI=self%lina_MI_star*(1+lina_qN) ! eqn 21
   lina_c=lina_colimitation(lina_qN,lina_qP,lina_MI,self%lina_cn) !colimitation implemented as function !eqn 17,18
   lina_muX=lina_cI*lina_c*self%lina_mu_max-lina_R !eqn 13
  !   lina_c_dot= !eqn 20
   lina_c_dot=-9999.99_rk
   lina_eta=lina_EPS_production(self%lina_E_min,self%lina_E_max,lina_muX,lina_MI,var%lina_N,var%lina_P,lina_qN,lina_qP,lina_gammaN,lina_gammaP)
   lina_rhox= self%lina_rho_starx * (var%dx ** (-self%lina_arho)) *  (1.0_rk - (1.0_rk - (var%dx ** (-self%lina_arho))) * lina_cI * lina_c ) !eqn 22
   lina_wx= 1.0_rk / (18.0_rk * self%lina_muw) * (lina_rhox-self%lina_rhow) * self%lina_g * var%dx ** 2.0_rk !eqn 23
!-------------------------------------------------
!Calculate RHS of the State equations these rates of change are defined by FABM to be per second

   rhsv%lina_X= (lina_muX - var%mx -lina_wx -var%Cx) * var%lina_X !eqn 1
   rhsv%lina_A= var%C_tot - (var%wa + var%ma + var%kB) * var%lina_A !eqn 2
   !(E and L ,D  are solved by the external AGG model)
   !rhsv%lina_E= lina_eta * var%lina_X - lina_h * var%lina_A !eqn 3
   !rhsv%lina_L= lina_rL - (self%lina_wl +var%CL)*var%lina_L + self%lina_kB*lina_psi * var%lina_A !eqn 4 
   !rhsv%lina_D= var%mx * var%lina_X - (var%mD + var%wD +var%CD) * var%lina_D + self%lina_kB (1- lina_psi) * var%lina_A !eqn 5
   rhsv%lina_N= - lina_gammaN * var%lina_X + var%mD * var%lina_D * var%lina_QDN + var%ma * var%lina_A * var%lina_QAN +self%lina_eN !eqn 6
   rhsv%lina_P= - lina_gammaP * var%lina_X + var%mD * var%lina_D * var%lina_QDP + var%ma * var%lina_A * var%lina_QAP +self%lina_eP !eqn 7

   rhsv%lina_QXN= lina_gammaN- lina_muX * var%lina_QXN !eqn 8
   rhsv%lina_QXP= lina_gammaP- lina_muX * var%lina_QXP !eqn 9

   rhsv%lina_QDN= var%mx *var%lina_X * var%lina_QXN - ( var%mD + var%wD + var%CD) * var%lina_D * var%lina_QDN + var%kB *(1- var%lina_psi) * var%lina_A * var%lina_QAN !eqn 39

   rhsv%lina_QDP= var%mx *var%lina_X * var%lina_QXP - ( var%mD + var%wD + var%CD) * var%lina_D * var%lina_QDP + var%kB *(1- var%lina_psi) * var%lina_A * var%lina_QAP !eqn 40

   rhsv%lina_psi= var%CL*var%lina_L/var%lina_A-var%kB*var%lina_psi-var%lina_psi/var%lina_A* rhsv%lina_A !eqn 24
   rhsv%lina_QAN= var%Cx *var%lina_X * var%lina_QXN + var%CD* var%lina_D * var%lina_QDN - (var%wa + var%ma + var%kB) * var%lina_A * var%lina_QAN !eqn 25
   rhsv%lina_QAP= var%Cx *var%lina_X * var%lina_QXP + var%CD* var%lina_D * var%lina_QDP - (var%wa + var%ma + var%kB) * var%lina_A * var%lina_QAP !eqn 26

!---------------------------------------------------

!-----------------------------------------------output of diagnostic
   if (Debugout) then
    write(*,*) 'diagnostic variables'
   endif
    _SET_DIAGNOSTIC_(self%id_lina_QN,lina_QN)    
    _SET_DIAGNOSTIC_(self%id_lina_QP,lina_QP) 
    _SET_DIAGNOSTIC_(self%id_lina_muX, lina_muX)                   !Phytoplankton Growth rate
    _SET_DIAGNOSTIC_(self%id_lina_wX, lina_wX)
    _SET_DIAGNOSTIC_(self%id_lina_gammaN,lina_gammaN)
    _SET_DIAGNOSTIC_(self%id_lina_gammaP,lina_gammaP)
    _SET_DIAGNOSTIC_(self%id_lina_cI,lina_cI)
    _SET_DIAGNOSTIC_(self%id_lina_c,lina_c)
    _SET_DIAGNOSTIC_(self%id_lina_MI,lina_MI)
!!    _SET_DIAGNOSTIC_(self%id_lina_C_dot,lina_C_dot)
    _SET_DIAGNOSTIC_(self%id_lina_eta,lina_eta)
    _SET_DIAGNOSTIC_(self%id_lina_R,lina_R)
    _SET_DIAGNOSTIC_(self%id_lina_rhox,lina_rhox)
 

  
  

!-----------------------------------------calculate the ODE
 !  if (Debugout) then write(*,*) 'calculate ODE' endif
!    _ADD_SOURCE_ (self%id_lina_X,rhsv%lina_X)
!   __ADD_SOURCE__(self%id_lina_X, rhsv%lina_X)
!   __ADD_SOURCE__(self%id_lina_A, rhsv%lina_A)
!   __ADD_SOURCE__(self%id_lina_E, rhsv%lina_E)
!   __ADD_SOURCE__(self%id_lina_L, rhsv%lina_L)
!   __ADD_SOURCE__(self%id_lina_D, rhsv%lina_D)
!   __ADD_SOURCE__(self%id_lina_N, rhsv%lina_N)
!   __ADD_SOURCE__(self%id_lina_P, rhsv%lina_P)
!   __ADD_SOURCE__(self%id_lina_QAN,rhsv%lina_QAN)
!   __ADD_SOURCE__(self%id_lina_QAP,rhsv%lina_QAP)
!   __ADD_SOURCE__(self%id_lina_QDN,rhsv%lina_QDN)
!   __ADD_SOURCE__(self%id_lina_QDP,rhsv%lina_QDP)
!   __ADD_SOURCE__(self%id_lina_QXN,rhsv%lina_QXN)
!   __ADD_SOURCE__(self%id_lina_QXP,rhsv%lina_QXN)
!   __ADD_SOURCE__(self%id_lina_psi,rhsv%lina_psi)
 
   
   _SET_ODE_(self%id_lina_X, rhsv%lina_X)
   _SET_ODE_(self%id_lina_A, rhsv%lina_A)
   _SET_ODE_(self%id_lina_E, rhsv%lina_E)
   _SET_ODE_(self%id_lina_L, rhsv%lina_L)
   _SET_ODE_(self%id_lina_D, rhsv%lina_D)
   _SET_ODE_(self%id_lina_N, rhsv%lina_N)
   _SET_ODE_(self%id_lina_P, rhsv%lina_P)
   _SET_ODE_(self%id_lina_QAN,rhsv%lina_QAN)
   _SET_ODE_(self%id_lina_QAP,rhsv%lina_QAP)
   _SET_ODE_(self%id_lina_QDN,rhsv%lina_QDN)
   _SET_ODE_(self%id_lina_QDP,rhsv%lina_QDP)
   _SET_ODE_(self%id_lina_QXN,rhsv%lina_QXN)
   _SET_ODE_(self%id_lina_QXP,rhsv%lina_QXN)
   _SET_ODE_(self%id_lina_psi,rhsv%lina_psi)

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
  
  function lina_gn(gn_r,gn_n) result(gn_g)
    real(rk),intent(in):: gn_r,gn_n
    real(rk) ::gn_g
    
    gn_g= (gn_r - ( gn_r**(1+gn_n))) / (1-(gn_r**(1+gn_n)))
  end function

function lina_EPS_production(E_min,E_max,muX,MI,N,P,qN,qP,gammaN,gammaP) result(eta)
    real(rk),intent(in):: E_min,E_max,muX,MI,N,P,qN,qP,gammaN,gammaP
    real(rk) :: eta
    
    if (qN .le. qP) then
        eta=E_min+(E_max-E_min)*(muX-gammaN/(1.-lina_gn(qN,MI))*(1.-(1.+MI)*qN**MI+MI*qN**(1.+MI))/(1.-qN**(1.+MI))**2)*(1-lina_gn(qN,MI)) !eqn 16
    else
        eta=E_min+(E_max-E_min)*(muX-gammaP/(1.-lina_gn(qP,MI))*(1.-(1.+MI)*qP**MI+MI*qP**(1.+MI))/(1.-qP**(1.+MI))**2)*(1-lina_gn(qP,MI)) !eqn 16
    end if

  end function
end module hzg_lina
