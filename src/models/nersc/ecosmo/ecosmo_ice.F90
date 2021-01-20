#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_nersc_ecosmo_ice --- ECOSMO biogeochemical model + ice
! biogeochemistry
!
! !INTERFACE:
   module fabm_nersc_ecosmo_ice
!
! !DESCRIPTION:
!
! The ECOSMO model is based on Daewel & Schrum (JMS,2013)
!
! !USES:
   use fabm_types
   use fabm_expressions

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_nersc_ecosmo_ice
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: sedy0 = 86400.0_rk
   real(rk), parameter :: mmolm3_in_mll = 44.6608009_rk
   real(rk)            :: redf(20)=0.0_rk
   real(rk)            :: BioC(45)=0.0_rk
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_nersc_ecosmo_ice
!     Variable identifiers
      type (type_state_variable_id)         :: id_no3, id_nh4, id_pho, id_sil
      type (type_state_variable_id)         :: id_opa, id_det1, id_det2, id_dia, id_fla
      type (type_state_variable_id)         :: id_diachl, id_flachl, id_bgchl
      type (type_state_variable_id)         :: id_mesozoo, id_microzoo, id_bg,id_dom, id_oxy
      type (type_state_variable_id)         :: id_dic, id_alk
      type (type_bottom_state_variable_id)  :: id_sed1, id_sed2, id_sed3
      type (type_dependency_id)             :: id_temp, id_salt, id_par
      type (type_dependency_id)             :: id_parmean
      type (type_horizontal_dependency_id)  :: id_tbs
      type (type_horizontal_dependency_id)  :: id_sfpar, id_meansfpar
      type (type_diagnostic_variable_id)    :: id_denit, id_primprod, id_secprod
      type (type_diagnostic_variable_id)    :: id_parmean_diag
      type (type_diagnostic_variable_id)    :: id_c2chl_fla, id_c2chl_dia,id_c2chl_bg
      type (type_horizontal_diagnostic_variable_id)    :: id_tbsout
      
      type (type_state_variable_id)         :: id_fish1, id_fish2
      type (type_bottom_state_variable_id)  :: id_mb1

      ! community sinking id's
      type (type_state_variable_id)         :: id_dsnk, id_dsnk2 ! variable that advects detritus sinking speed
      type (type_diagnostic_variable_id)    :: id_snkspd, id_snkspd2 ! average calculated sinking speed

!id diagnostics for fish calculation - integration
      type (type_horizontal_diagnostic_variable_id)    :: id_meso_integr,id_micro_integr,id_det1_integr,id_det2_integr,id_fish_integr,id_fish2_integr
      type (type_horizontal_dependency_id)             :: id_mesoint,id_microint, id_det1int,id_det2int,id_fishint,id_fish2int
      type (type_horizontal_diagnostic_variable_id)    :: id_dF1onZl,id_dF1onZs,id_dF1onDet1,id_dF1onDet2
      type (type_horizontal_diagnostic_variable_id)    :: id_dF2onZl,id_dF2onF1,id_dF2onDet1,id_dF2onDet2
      type (type_horizontal_dependency_id)             :: id_f1zl,id_f1zs, id_f1det1, id_f1det2, id_f1mb1
      type (type_horizontal_dependency_id)             :: id_f2zl,id_f2f1, id_f2det1, id_f2det2, id_f2mb1

!     Model parameters
      real(rk) :: BioC(45)
      real(rk) :: zpr, frr
      real(rk) :: prefZsPs
      real(rk) :: prefZsPl
      real(rk) :: prefZsBG
      real(rk) :: prefZsD
      real(rk) :: prefZlPs
      real(rk) :: prefZlPl
      real(rk) :: prefZlBG
      real(rk) :: prefZlD
      real(rk) :: prefZlZs
      real(rk) :: surface_deposition_no3
      real(rk) :: surface_deposition_nh4
      real(rk) :: surface_deposition_pho
      real(rk) :: surface_deposition_sil
      real(rk) :: nfixation_minimum_daily_par
      real(rk) :: bg_growth_minimum_daily_rad
      real(rk) :: MAXchl2nPs, MINchl2nPs 
      real(rk) :: MAXchl2nPl, MINchl2nPl 
      real(rk) :: MAXchl2nBG, MINchl2nBG 
      real(rk) :: alfaPl, alfaPs, alfaBG
! Déborah IA parameter declaration
      real(rk) :: sinkIdet
      real(rk) :: T_lim_phy
! Ute MB parameter declaration   
      real(rk) :: mMB1,excMB1,rMB1,tempcMB1
      real(rk) :: GrMB1Z,GrMB1P,GrMB1Det,GrMB1Sed 
      real(rk) :: prefMB1Zs
      real(rk) :: prefMB1Zl
      real(rk) :: prefMB1Sed
      real(rk) :: prefMB1Det
      real(rk) :: prefMB1P
      real(rk) :: prefMB1Dom
! Ute Fish parameter declaration   
      real(rk) :: mF1,excF1,rF1,asefF1,tempcF1,rF1MB
      real(rk) :: mF2,excF2,rF2,asefF2,tempcF2,rF2MB
      real(rk) :: GrF1Zl,GrF1Zs,GrF1Det,GrF1MB1 
      real(rk) :: GrF2Zl,GrF2F1,GrF2Det,GrF2MB1 
      real(rk) :: prefF1Zl,prefF2Zl
      real(rk) :: prefF1Zs,prefF2F1
      real(rk) :: prefF1Det,prefF2Det
      real(rk) :: prefF1MB1,prefF2MB1
! community dependent sinking parameters
      real(rk) :: sinkDiaD,sinkFlaD,sinkMicD,sinkMesD,sinkBgD,sinkFish1D,sinkFish2D
 
      logical  :: use_chl,use_cyanos,couple_co2, use_fish, community_export
      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom
      procedure :: get_light_extinction
      procedure :: get_vertical_movement 

   end type type_nersc_ecosmo_ice
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the ECOSMO model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
   class (type_nersc_ecosmo_ice),intent(inout),target  :: self
   integer,                intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
!  Caglar Yumruktepe:
!  Added fabm.yaml support: parameters from yaml file are copied to
!                           BioC array. Eventually, BioC array will be removed
!                           from the model where parameter names from the yaml
!                           file will be used.
!  Added dynamic chlorophyll-a from Geider etal., 1997
!
! !LOCAL VARIABLES:
! Everything else taken from yaml file
!
   integer :: i
   ! set Redfield ratios:
   redf(1) = 6.625_rk      !C_N
   redf(2) = 106.0_rk      !C_P
   redf(3) = 6.625_rk      !C_SiO
   redf(4) = 16.0_rk       !N_P
   redf(5) = 1.0_rk        !N_SiO
   redf(6) = 12.01_rk      !C_Cmg
   redf(7) = 44.6608009_rk !O2mm_ml
   redf(8) = 14.007_rk     !N_Nmg
   redf(9) = 30.97_rk      !P_Pmg
   redf(10) = 28.09_rk     !Si_Simg
   do i=1,10
     redf(i+10) = 1._rk/redf(i)
   end do


!EOP
!-----------------------------------------------------------------------
!BOC

   call self%get_parameter(self%zpr, 'zpr', '1/day', 'zpr_long_name_needed', default=0.001_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter(self%frr, 'frr', '-', 'fraction of dissolved from det.', default=0.4_rk)
   call self%get_parameter(self%nfixation_minimum_daily_par, 'nfixation_minimum_daily_par', 'nfixation minimum daily par', default=40.0_rk)
   call self%get_parameter(self%bg_growth_minimum_daily_rad, 'bg_growth_minimum_daily_rad', 'bg growth minimum daily rad', default=120.0_rk)
   ! set surface fluxes in [mgC/m2/s]
   call self%get_parameter( self%surface_deposition_no3, 'surface_deposition_no3', 'mmolN/m**2 d', 'surface deposition no3', default=0.0_rk, scale_factor=redf(1)*redf(6)/sedy0 )
   call self%get_parameter( self%surface_deposition_nh4, 'surface_deposition_nh4', 'mmolN/m**2 d', 'surface deposition nh4', default=0.0_rk, scale_factor=redf(1)*redf(6)/sedy0 )
   call self%get_parameter( self%surface_deposition_pho, 'surface_deposition_pho', 'mmolN/m**2 d', 'surface deposition pho', default=0.0_rk, scale_factor=redf(2)*redf(6)/sedy0 )
   call self%get_parameter( self%surface_deposition_sil, 'surface_deposition_sil', 'mmolN/m**2 d', 'surface deposition sil', default=0.0_rk, scale_factor=redf(3)*redf(6)/sedy0 )
   !  change units 1/day to 1/sec and mmolN,P,Si to mmolC
   call self%get_parameter( self%BioC(1) , 'muPl',        '1/day',      'max growth rate for Pl',          default=1.30_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(2) , 'muPs',        '1/day',      'max growth rate for Ps',          default=1.10_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(3) , 'aa',          'm**2/W',     'photosynthesis ef-cy',            default=0.04_rk)
   call self%get_parameter( self%BioC(4) , 'EXw',         '1/m',        'light extinction',                default=0.041_rk)
   if (self%use_chl) then
      call self%get_parameter( self%BioC(5) , 'Exphy',       'm**2/mgCHL', 'phyto self-shading',              default=0.04_rk )
   else
      call self%get_parameter( self%BioC(5) , 'Exphy',       'm**2/mmolN', 'phyto self-shading',              default=0.04_rk, scale_factor=1.0_rk/(redf(1)*redf(6)) )
   end if
      call self%get_parameter( self%extdet , 'Exdet',       'm**2/molC', 'detritus self-shading',         default=0.0_rk, scale_factor=1.0_rk/1000.0_rk )
      call self%get_parameter( self%extdom , 'Exdom',       'm**2/molC', 'dom self-shading',              default=0.0_rk, scale_factor=1.0_rk/1000.0_rk )
   call self%get_parameter( self%BioC(6) , 'rNH4',        'mmolN/m**3', 'NH4 half saturation',             default=0.20_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%BioC(7) , 'rNO3',        'mmolN/m**3', 'NO3 half saturation',             default=0.50_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%BioC(8) , 'psi',         'm**3/mmolN', 'NH4 inhibition',                  default=3.0_rk,   scale_factor=1.0_rk/(redf(1)*redf(6)) )
   call self%get_parameter( self%BioC(9) , 'mPl',         '1/day',      'Pl mortality rate',               default=0.04_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(10), 'mPs',         '1/day',      'Ps mortality rate',               default=0.08_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(11), 'GrZlP',       '1/day',      'Grazing rate Zl on Phyto',        default=0.80_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(12), 'GrZsP',       '1/day',      'Grazing rate Zs on Phyto',        default=1.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(13), 'GrZlZ',       '1/day',      'Grazing rate Zl on Zs',           default=0.50_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(14), 'Rg',          'mmolN/m**3', 'Zs, Zl half saturation',          default=0.50_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%BioC(15), 'mZl',         '1/day',      'Zl mortality rate',               default=0.10_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(16), 'mZs',         '1/day',      'Zs mortality rate',               default=0.20_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(17), 'excZl',       '1/day',      'Zl excretion rate',               default=0.06_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(18), 'excZs',       '1/day',      'Zs excretion rate',               default=0.08_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(19), 'gammaZlp',    '1',          'Zl assim. eff. on plankton',      default=0.75_rk)
   call self%get_parameter( self%BioC(20), 'gammaZsp',    '1',          'Zs assim. eff. on plankton',      default=0.75_rk)
   call self%get_parameter( self%BioC(21), 'gammaZd',     '1',          'Zl & Zs assim. eff. on det',      default=0.75_rk)
   call self%get_parameter( self%BioC(22), 'reminD',      '1/day',      'Detritus remin. rate',            default=0.003_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(23), 'sinkDet',     'm/day',      'Detritus sinking rate',           default=5.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(24), 'Wa',          'm/day',      '???',                             default=1.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(25),  'rPO4',       'mmolP/m**3', 'PO4 half saturation',             default=0.05_rk,  scale_factor=redf(2)*redf(6))
   call self%get_parameter( self%BioC(26),  'rSi',        'mmolSi/m**3','SiO2 half saturation',            default=0.50_rk,  scale_factor=redf(3)*redf(6))
   call self%get_parameter( self%BioC(27),  'regenSi',    '1/day',      'Si regeneration rate',            default=0.015_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(28),  'muBG',       '1/day',      'max growth rate for BG',          default=1.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(29),  'TctrlBG',    '1/degC',     'BG T control beta',               default=1.00_rk)
   call self%get_parameter( self%BioC(30),  'TrefBG',     'degC',       'BG reference temperature',        default=0.00_rk)
   call self%get_parameter( self%BioC(31),  'GrBG',       '1/day',      'BG max grazing rate',             default=0.30_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(32),  'mBG',        '1/day',      'BG mortality rate',               default=0.08_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(33),  'upliftBG',   'm/day',      'BG uplifting rate',               default=0.10_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(34),  'crBotStr',   'N/m**2',     'critic. bot. stress for resusp.', default=0.007_rk)
   call self%get_parameter( self%BioC(35),  'resuspRt',   '1/day',      'resuspension rate',               default=25.0_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(36),  'sedimRt',    'm/day',      'sedimentation rate',              default=3.5_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(37),  'burialRt',   '1/day',      'burial rate',                     default=1E-5_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(38),  'reminSED',   '1/day',      'sediment remineralization rate',  default=0.001_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(39),  'TctrlDenit', '1/degC',     'temp. control denitrification',   default=0.15_rk)
   call self%get_parameter( self%BioC(40),  'RelSEDp1',   'units??',    'P sedim. rel. p1',                default=0.15_rk)
   call self%get_parameter( self%BioC(41),  'RelSEDp2',   'units??',    'P sedim. rel. p2',                default=0.10_rk)
   call self%get_parameter( self%BioC(42),  'reminSEDsi', '1/day',      'sed. remineralization rate Si',   default=0.0002_rk,scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(43),  'sinkOPAL',   'm/day',      'OPAL sinking rate',               default=5.0_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(44),  'sinkBG',     'm/day',      'BG sinking rate',                 default=-1.0_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(45),  'sinkDia',    'm/y',      'Diatom sinking rate',             default=0.0_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%sinkIdet,  'sinkIdet',   'm/day',      'Ice algae sinking rate',          default=0.0_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%T_lim_phy, 'T_lim_phy',  '-',          'Coefficient of temperature dependency relationship',     default=0.0582_rk)
   !  growth fractions
   call self%get_parameter( self%prefZsPs,  'prefZsPs',   '-',          'Grazing preference Zs on Ps',     default=0.70_rk)
   call self%get_parameter( self%prefZsPl,  'prefZsPl',   '-',          'Grazing preference Zs on Pl',     default=0.25_rk)
   call self%get_parameter( self%prefZsD,   'prefZsD',    '-',          'Grazing preference Zs on Det.',   default=0.00_rk)
   call self%get_parameter( self%prefZsBG,  'prefZsBG',   '-',          'Grazing preference Zs on BG',     default=0.30_rk)
   call self%get_parameter( self%prefZlPs,  'prefZlPs',   '-',          'Grazing preference Zl on Ps',     default=0.10_rk)
   call self%get_parameter( self%prefZlPl,  'prefZlPl',   '-',          'Grazing preference Zl on Pl',     default=0.85_rk)
   call self%get_parameter( self%prefZlZs,  'prefZlZs',   '-',          'Grazing preference Zl on Zs',     default=0.15_rk)
   call self%get_parameter( self%prefZlD,   'prefZlD',    '-',          'Grazing preference Zl on Det.',   default=0.00_rk)
   call self%get_parameter( self%prefZlBG,  'prefZlBG',   '-',          'Grazing preference Zl on BG',     default=0.30_rk)
   ! ute Macrobenthos parameter settings
   call self%get_parameter( self%mMB1,      'mMB1',       '1/day', 'Macrobenthos mortality rate',     default=0.001_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%excMB1,    'excMB1',     '1/day', 'Macrobenthos excretion rate',     default=0.025_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%rMB1,      'rMB1',       'mmolN/m**3',  'MB grazing half saturation',      default=0.50_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%GrMB1Z,    'GrMB1Z',     '1/day', 'Grazingrate MB1 on Zooplankton',  default=0.1_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrMB1P,    'GrMB1P',     '1/day', 'Grazingrate MB1 on Phytoplankton',default=0.1_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrMB1Det,  'GrMB1Det',   '1/day', 'Grazingrate MB1 on Detritus',     default=0.1_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrMB1Sed,  'GrMB1Sed',   '1/day', 'Grazingrate MB1 on Sediment',     default=0.1_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%tempcMB1,  'tempcMB1',   '1/day',      'MB1 temperature control excr.',   default=0._rk)
   call self%get_parameter( self%prefMB1Zs, 'prefMB1Zs',  '-',          'Grazing preference MB1 on Zs',    default=0.20_rk)
   call self%get_parameter( self%prefMB1Zl, 'prefMB1Zl',  '-',          'Grazing preference MB1 on Zl',    default=0.30_rk)
   call self%get_parameter( self%prefMB1Det,'prefMB1Det', '-',          'Grazing preference MB1 on Det',   default=0.10_rk)
   call self%get_parameter( self%prefMB1Sed,'prefMB1Sed', '-',          'Grazing preference MB1 on Sed',   default=0.10_rk)
   call self%get_parameter( self%prefMB1P,  'prefMB1P',   '-',          'Grazing preference MB1 on P',     default=0.20_rk)
   call self%get_parameter( self%prefMB1Dom,'prefMB1Dom', '-',          'Grazing preference MB1 on Dom',   default=0.10_rk)
   !parameter for fish1 functional group
   call self%get_parameter( self%mF1,      'mF1',       '1/day', 'fish1 mortality rate',     default=0.001_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%excF1,    'excF1',     '1/day', 'fish1 excretion rate',     default=0.002_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%rF1,      'rF1',       'mmolN/m**3', 'fish1 grazing half saturation',    default=0.70_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%rF1MB,    'rF1MB',     'mmolN/m**3', 'fish1 grazing on MB half saturation',   default=0.90_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%tempcF1,  'tempcF1',   '1/day',      'fish1 temperature control excr.',   default=0.5_rk)
   call self%get_parameter( self%asefF1,   'asefF1',    '-',      'fish1 assimilation efficiency',           default=0.7_rk)
   call self%get_parameter( self%GrF1Zl,   'GrF1Zl',    '1/day', 'Grazingrate fish1 on Mesozoo',           default=0.01_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrF1Zs,   'GrF1Zs',    '1/day', 'Grazingrate fish1 on Microzoo',          default=0.01_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrF1Det,  'GrF1Det',   '1/day', 'Grazingrate fish1 on Detritus',          default=0.05_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrF1MB1,  'GrF1MB1',   '1/day', 'Grazingrate fish1 on macrobenthos1',     default=0.01_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%prefF1Zl, 'prefF1Zl',  '-',          'Grazing preference fish1 on Zl',    default=0.45_rk)
   call self%get_parameter( self%prefF1Zs, 'prefF1Zs',  '-',          'Grazing preference fish1 on Zs',    default=0.25_rk)
   call self%get_parameter( self%prefF1Det,'prefF1Det', '-',          'Grazing preference fish1 on Det',   default=0.05_rk)
   call self%get_parameter( self%prefF1MB1,'prefF1MB1', '-',          'Grazing preference fish1 on MB1',   default=0.25_rk)
   !parameter for fish2 functional group
   call self%get_parameter( self%mF2,      'mF2',       '1/day', 'fish2 mortality rate',     default=0.001_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%excF2,    'excF2',     '1/day', 'fish2 excretion rate',     default=0.002_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%rF2,      'rF2',       'mmolN/m**3', 'fish2 grazing half saturation',    default=0.70_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%rF2MB,    'rF2MB',     'mmolN/m**3', 'fish2 grazing on MB half saturation',   default=0.90_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%tempcF2,  'tempcF2',   '1/day',      'fish2 temperature control excr.',   default=0.5_rk)
   call self%get_parameter( self%asefF2,   'asefF2',    '-',      'fish2 assimilation efficiency',           default=0.7_rk)
   call self%get_parameter( self%GrF2Zl,   'GrF2Zl',    '1/day', 'Grazingrate fish2 on Mesozoo',           default=0.01_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrF2F1,   'GrF2F1',    '1/day', 'Grazingrate fish2 on Microzoo',          default=0.01_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrF2Det,  'GrF2Det',   '1/day', 'Grazingrate fish2 on Detritus',          default=0.05_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrF2MB1,  'GrF2MB1',   '1/day', 'Grazingrate fish2 on macrobenthos1',     default=0.01_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%prefF2Zl, 'prefF2Zl',  '-',          'Grazing preference fish2 on Zl',    default=0.45_rk)
   call self%get_parameter( self%prefF2F1, 'prefF2F1',  '-',          'Grazing preference fish2 on F1',    default=0.25_rk)
   call self%get_parameter( self%prefF2Det,'prefF2Det', '-',          'Grazing preference fish2 on Det',   default=0.05_rk)
   call self%get_parameter( self%prefF2MB1,'prefF2MB1', '-',          'Grazing preference fish2 on MB1',   default=0.25_rk)

   ! chlorophyll-a constants
   call self%get_parameter( self%MINchl2nPs, 'MINchl2nPs', 'mgChl/mmolN', 'minimum Chl to N ratio Ps', default=0.50_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%MAXchl2nPs, 'MAXchl2nPs', 'mgChl/mmolN', 'maximum Chl to N ratio Ps', default=3.83_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%MINchl2nPl, 'MINchl2nPl', 'mgChl/mmolN', 'minimum Chl to N ratio Pl', default=0.50_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%MAXchl2nPl, 'MAXchl2nPl', 'mgChl/mmolN', 'maximum Chl to N ratio Pl', default=2.94_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%MINchl2nBG, 'MINchl2nBG', 'mgChl/mmolN', 'minimum Chl to N ratio BG', default=0.50_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%MAXchl2nBG, 'MAXchl2nBG', 'mgChl/mmolN', 'maximum Chl to N ratio BG', default=3.83_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%alfaPs,     'alfaPs', 'mmolN m2/(mgChl day W)**-1', 'initial slope P-I curve Ps', default=0.0393_rk, scale_factor=redf(1)*redf(6) )
   call self%get_parameter( self%alfaPl,     'alfaPl', 'mmolN m2/(mgChl day W)**-1', 'initial slope P-I curve Pl', default=0.0531_rk, scale_factor=redf(1)*redf(6) )
   call self%get_parameter( self%alfaBG,     'alfaBG', 'mmolN m2/(mgChl day W)**-1', 'initial slope P-I curve BG', default=0.0393_rk, scale_factor=redf(1)*redf(6) )

   ! community sinking constants
   call self%get_parameter( self%sinkDiaD,  'sinkDiaD',    'm/day', 'Detritus originating from diatom sinking rate',       default=5.0_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%sinkFlaD,  'sinkFlaD',    'm/day', 'Detritus originating from flagellates sinking rate',  default=5.0_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%sinkMicD,  'sinkMicD',    'm/day', 'Detritus originating from microzoo sinking rate',     default=5.0_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%sinkMesD,  'sinkMesD',    'm/day', 'Detritus originating from mesozoo sinking rate',      default=5.0_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%sinkBgD,   'sinkBgD',     'm/day', 'Detritus originating from cyanob. sinking rate',      default=5.0_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%sinkFish1D,'sinkFish1D',  'm/day', 'Detritus originating from fish1 sinking rate',        default=5.0_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%sinkFish2D,'sinkFish2D',  'm/day', 'Detritus originating from fish2 sinking rate',        default=5.0_rk, scale_factor=1.0_rk/sedy0)      
   ! add switches
   call self%get_parameter( self%use_cyanos,     'use_cyanos', '', 'switch cyanobacteria', default=.true.)
   call self%get_parameter( self%use_fish,     'use_fish', '', 'switch fish', default=.false.)
   call self%get_parameter( self%couple_co2,     'couple_co2', '', 'switch coupling to carbonate module', default=.false.)
   call self%get_parameter( self%use_chl,     'use_chl', '', 'switch chlorophyll/c dynamics', default=.true.)
   call self%get_parameter( self%community_export, 'community_export','','switch community influenced export parameterization', default=.false.) 
   ! Register state variables
   call self%register_state_variable( self%id_no3,      'no3',     'mgC/m3',    'nitrate',                   minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=5.0_rk*redf(1)*redf(6)  )
   call self%register_state_variable( self%id_nh4,      'nh4',     'mgC/m3',    'ammonium',                  minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=0.1_rk*redf(1)*redf(6)  )
   call self%register_state_variable( self%id_pho,      'pho',     'mgC/m3',    'phosphate',                 minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=0.3_rk*redf(2)*redf(6)  )
   call self%register_state_variable( self%id_sil,      'sil',     'mgC/m3',    'silicate',                  minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=5.0_rk*redf(3)*redf(6)  )
   call self%register_state_variable( self%id_oxy,      'oxy',     'mmolO2/m3', 'oxygen',                    minimum=0.0_rk,   vertical_movement=0.0_rk,  &
                                      initial_value=85.0_rk  )
   call self%register_state_variable( self%id_fla,      'fla',     'mgC/m3',    'small phytoplankton',       minimum=1.0e-7_rk,     vertical_movement=0.0_rk, &
                                      initial_value=1e-4_rk*redf(1)*redf(6))
   call self%register_state_variable( self%id_dia,      'dia',     'mgC/m3',    'large phytoplankton',       minimum=1.0e-7_rk,     vertical_movement=-self%BioC(45) , &
                                      initial_value=1e-4_rk*redf(1)*redf(6) )
   if (self%use_cyanos) then
     call self%register_state_variable( self%id_bg,       'bg',      'mgC/m3',    'cyanobacteria',             minimum=1.0e-14_rk,     vertical_movement=-self%BioC(44) , &
                                      initial_value=1e-4_rk*redf(1)*redf(6) )
     if (self%use_chl) &
       call self%register_state_variable( self%id_bgchl,    'bgchl',   'mgChl/m3',  'cyanobacteria chl-a',       minimum=1.0e-14_rk/20., vertical_movement=-self%BioC(44) , &
                                      initial_value=1e-4_rk*redf(1)*redf(6)/27.)
   end if
   if (self%use_chl) then
     call self%register_state_variable( self%id_diachl,   'diachl',  'mgChl/m3',  'large phytoplankton chl-a', minimum=1.0e-7_rk/27., vertical_movement=-self%BioC(45) , &
                                      initial_value=1e-4_rk*redf(1)*redf(6)/27.)
     call self%register_state_variable( self%id_flachl,   'flachl',  'mgChl/m3',  'small phytoplankton chl-a', minimum=1.0e-7_rk/20., vertical_movement=0.0_rk, &
                                      initial_value=1e-4_rk*redf(1)*redf(6)/20.)
   end if
   call self%register_state_variable( self%id_microzoo, 'microzoo','mgC/m3',    'microzooplankton',          minimum=1.0e-7_rk,     vertical_movement=0.0_rk, &
                                      initial_value=1e-6_rk*redf(1)*redf(6) )
   call self%register_state_variable( self%id_mesozoo,  'mesozoo', 'mgC/m3',    'mesozooplankton',           minimum=1.0e-7_rk,     vertical_movement=0.0_rk, &
                                      initial_value=1e-6_rk*redf(1)*redf(6) )
   call self%register_state_variable( self%id_det1,      'det1',   'mgC/m3',    'fast sinking detritus',     minimum=0.0_rk,        vertical_movement=-self%sinkIdet , &
                                      initial_value=1.0_rk*redf(1)*redf(6)  )
   call self%register_state_variable( self%id_det2,      'det2',   'mgC/m3',    'slow sinking detritus',     minimum=0.0_rk,        vertical_movement=-self%BioC(23) , &
                                      initial_value=1.0_rk*redf(1)*redf(6)  )
   call self%register_state_variable( self%id_opa,      'opa',     'mgC/m3',    'opal',                      minimum=0.0_rk,        vertical_movement=-self%BioC(43) , &
                                      initial_value=2.0_rk*redf(3)*redf(6) )
   call self%register_state_variable( self%id_dom,      'dom',     'mgC/m3',    'labile dissolved om',       minimum=0.0_rk , &
                                      initial_value=3.0_rk*redf(1)*redf(6)   )
   call self%register_state_variable( self%id_sed1,     'sed1',    'mgC/m2',    'sediment detritus',         minimum=0.0_rk , &
                                      initial_value=20.0_rk*redf(1)*redf(6)*redf(18) )
   call self%register_state_variable( self%id_sed2,     'sed2',    'mgC/m2',    'sediment opal',             minimum=0.0_rk , &
                                      initial_value=20.0_rk*redf(3)*redf(6)*redf(20) )
   call self%register_state_variable( self%id_sed3,     'sed3',    'mgC/m2',    'sediment adsorbed pho.',    minimum=0.0_rk , &
                                      initial_value=2.0_rk*redf(2)*redf(6)*redf(19) )
     if (self%use_fish) then
    !  ute register macrobenthos and Fish state variable
   call self%register_state_variable( self%id_mb1,     'mb1',    'mgC/m2',    'macrobenthos',    minimum=0.0_rk,    &      
                                      initial_value=1e-6_rk*redf(1)*redf(6) )
   call self%register_state_variable( self%id_fish1,     'fis1',    'mgC/m3',    'fish functional group 1',    minimum=0.0_rk,    &      
                                      initial_value=1e-6_rk*redf(1)*redf(6) )
   call self%register_state_variable( self%id_fish2,     'fis2',    'mgC/m3',    'fish functional group 2',    minimum=0.0_rk,    &      
                                      initial_value=1e-7_rk*redf(1)*redf(6) )
    end if
   ! community sinking speed state variable
   if (self%community_export) then
      call self%register_state_variable( self%id_dsnk,      'dsnk',    'mgC/m3', 'detritus group 1 sinking speed advector', minimum=0.0_rk, vertical_movement=self%BioC(23) , &
                                      initial_value=2.0_rk*redf(1)*redf(6)*self%BioC(23))
      call self%register_state_variable( self%id_dsnk2,      'dsnk2',    'mgC/m3', 'detritus group 2 sinking speed advector', minimum=0.0_rk, vertical_movement=self%BioC(23) , &
                                      initial_value=2.0_rk*redf(1)*redf(6)*self%BioC(23))
   end if

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_denit,'denit','mmolN/m**3/s', &
         'denitrification rate', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_primprod,'primprod','mgC/m**3/s', &
         'primary production rate', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_secprod,'secprod','mgC/m**3/s', &
         'secondary production rate', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_parmean_diag,'parmean','W/m**2', &
         'daily-mean photosynthetically active radiation', output=output_time_step_averaged)
   if (self%community_export) then
      call self%register_diagnostic_variable(self%id_snkspd,'snkspd','m/d', & 
         'daily-mean detritus group 1 sinking speed', output=output_time_step_averaged) 
      call self%register_diagnostic_variable(self%id_snkspd2,'snkspd2','m/d', & 
         'daily-mean detritus group 2 sinking speed', output=output_time_step_averaged) 
   end if
! CAGLAR: Turned off to avoid excess disk space usage
! calls outputs - simulated Carbon to chlorophyll-a ratio
!   if (self%use_chl) then
!     call self%register_diagnostic_variable(self%id_c2chl_fla,'c2chl_fla','mgC/mgCHL', &
!         'daily-mean C to CHL ratio for flagellates', output=output_time_step_averaged)
!     call self%register_diagnostic_variable(self%id_c2chl_dia,'c2chl_dia','mgC/mgCHL', &
!         'daily-mean C to CHL ratio for diatoms', output=output_time_step_averaged)
!         if (self%use_cyanos) then
!             call self%register_diagnostic_variable(self%id_c2chl_bg,'c2chl_bg','mgC/mgCHL', &
!             'daily-mean C to CHL ratio for cyanobacteria', output=output_time_step_averaged)
!         end if
!   end if
   call self%register_diagnostic_variable(self%id_tbsout,'botstrss','fill_later', &
         'total bottom stress', output=output_time_step_averaged)
    if (self%use_fish) then
 !ute register diagnostics for fish calculation - integration
   call self%register_diagnostic_variable(self%id_meso_integr,'Zlint','mgC/m**2', &
         'depthintegrated mesozooplankton')
   call self%register_diagnostic_variable(self%id_micro_integr,'Zsint','mgC/m**2', &
         'depthintegrated microzooplankton')
   call self%register_diagnostic_variable(self%id_det1_integr,'Det1int','mgC/m**2', &
         'depthintegrated detritus group 1')
   call self%register_diagnostic_variable(self%id_det2_integr,'Det2int','mgC/m**2', &
         'depthintegrated detritus group 2')
   call self%register_diagnostic_variable(self%id_fish_integr,'Fishint','mgC/m**2', &
         'depthintegrated fish1')
   call self%register_diagnostic_variable(self%id_dF1onZl,'F1Zl','mgC/m**2/s', &
         'Fish feeding on Mesozoo')
   call self%register_diagnostic_variable(self%id_dF1onZs,'F1Zs','mgC/m**2/s', &
         'Fish feeding on Microzoo')
   call self%register_diagnostic_variable(self%id_dF1onDet1,'F1De1','mgC/m**2/s', &
         'Fish feeding on Detritus group 1')
   call self%register_diagnostic_variable(self%id_dF1onDet2,'F1De2','mgC/m**2/s', &
         'Fish feeding on Detritus group 2')
!  fish 2
   call self%register_diagnostic_variable(self%id_fish2_integr,'F2int','mgC/m**2', &
         'depthintegrated fish2')
   call self%register_diagnostic_variable(self%id_dF2onZl,'F2Zl','mgC/m**2/s', &
         'Fish2 feeding on Mesozoo')
   call self%register_diagnostic_variable(self%id_dF2onF1,'F2F1','mgC/m**2/s', &
         'Fish2 feeding on Fish1')
   call self%register_diagnostic_variable(self%id_dF2onDet1,'F2De1','mgC/m**2/s', &
         'Fish2 feeding on Detritus group 1')
   call self%register_diagnostic_variable(self%id_dF2onDet2,'F2De2','mgC/m**2/s', &
         'Fish2 feeding on Detritus group 2')
    endif

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_tbs,standard_variables%bottom_stress)
   call self%register_dependency(self%id_sfpar,standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   ! use temporal mean of light for the last 24 hours
   call self%register_dependency(self%id_parmean,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk,missing_value=0.0_rk))
   call self%register_dependency(self%id_meansfpar,temporal_mean(self%id_sfpar,period=86400._rk,resolution=3600._rk))

   if (self%couple_co2) then
     call self%register_state_dependency(self%id_dic, 'dic_target','mmol m-3','dic budget')
     call self%register_state_dependency(self%id_alk, 'alk_target','mmol m-3','alkalinity budget')
   end if

    if (self%use_fish) then
   ! ute dependencies used for fish production 
   call self%register_dependency(self%id_mesoint,vertical_integral(self%id_mesozoo))
   call self%register_dependency(self%id_microint,vertical_integral(self%id_microzoo))
   call self%register_dependency(self%id_det1int,vertical_integral(self%id_det1))
   call self%register_dependency(self%id_det2int,vertical_integral(self%id_det2))
   call self%register_dependency(self%id_fishint,vertical_integral(self%id_fish1))
   call self%register_dependency(self%id_f1zl,'ECO_F1Zl')
   call self%register_dependency(self%id_f1zs,'ECO_F1Zs')
   call self%register_dependency(self%id_f1det1,'ECO_F1De1')
   call self%register_dependency(self%id_f1det2,'ECO_F1De2')
   call self%register_dependency(self%id_fish2int,vertical_integral(self%id_fish2))
   call self%register_dependency(self%id_f2zl,'ECO_F2Zl')
   call self%register_dependency(self%id_f2f1,'ECO_F2F1')
   call self%register_dependency(self%id_f2det1,'ECO_F2De1')
   call self%register_dependency(self%id_f2det2,'ECO_F2De2')
   end if

   return

end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of ECOSMO model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_nersc_ecosmo_ice),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   real(rk) :: no3,nh4,pho,sil,t_sil,oxy,fla,bg,dia
   real(rk) :: flachl,diachl,bgchl,chl2c_fla,chl2c_dia,chl2c_bg
   real(rk) :: microzoo,mesozoo,opa,det1,det2,dom
   real(rk) :: temp,salt,par
   real(rk) :: frem, fremDOM, blight
   real(rk) :: Ts,Tl,Tbg
   real(rk) :: Prod,Ps_prod,Pl_prod,Bg_prod
   real(rk) :: Fs,Fl,ZlonPs,ZlonPl,ZsonD1,ZsonD2,ZlonD1,ZlonD2,ZlonBg,ZsonBg,ZsonPs,ZsonPl,ZlonZs
   real(rk) :: up_no3,up_nh4,up_n,up_pho,up_sil
   real(rk) :: bioom1,bioom2,bioom3,bioom4,bioom5,bioom6,bioom7,bioom8,Onitr
   real(rk) :: rhs,dxxdet1,dxxdet2,rhs_ID
   real(rk) :: Zl_prod, Zs_prod
   real(rk) :: mean_par, mean_surface_par, Bg_fix
   real(rk) :: fla_loss=1.0_rk
   real(rk) :: dia_loss=1.0_rk
   real(rk) :: bg_loss=1.0_rk
   real(rk) :: mic_loss=1.0_rk
   real(rk) :: mes_loss=1.0_rk
   real(rk) :: tbs
   real(rk) :: rhs_oxy,rhs_amm,rhs_nit
! local variables for fish (Ute)   
   real(rk) :: mesoi,microi,deti1,deti2,fishi,F1onZl,F1onZs,F1onDet1,F1onDet2, Fish_prod,ttemp
   real(rk) :: fish1,fish2
   real(rk) :: fishi2,F2onZl,F2onF1,F2onDet1,F2onDet2, Fish2_prod
! local variables for community sinking
   real(rk) :: dsnk,dsnk2,mean_snkspd,mean_snkspd2  
 
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp,temp)
   _GET_(self%id_salt,salt)
   _GET_(self%id_par,par)
   _GET_(self%id_no3,no3)
   _GET_(self%id_nh4,nh4)
   _GET_(self%id_pho,pho)
   _GET_(self%id_sil,sil)
   _GET_(self%id_dia,dia)
   _GET_(self%id_fla,fla)
   if (self%use_cyanos) then
     _GET_(self%id_bg,bg)
   else
     bg=0.0_rk
   end if
   if (self%use_chl) then
     _GET_(self%id_diachl,diachl)
     _GET_(self%id_flachl,flachl)
     if (self%use_cyanos) then
       _GET_(self%id_bgchl,bgchl)
     end if
   end if
   if (self%community_export) then
      _GET_(self%id_dsnk,dsnk)
   end if
   _GET_(self%id_microzoo,microzoo)
   _GET_(self%id_mesozoo,mesozoo)
   _GET_(self%id_det1,det1)
   _GET_(self%id_det2,det2)
   _GET_(self%id_dom,dom)
   _GET_(self%id_opa,opa)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_parmean,mean_par)
   _GET_HORIZONTAL_(self%id_meansfpar,mean_surface_par)
   _GET_HORIZONTAL_(self%id_tbs,tbs)
   if(self%use_fish) then
    ! get horizontal variables for fish estimates
   _GET_HORIZONTAL_(self%id_mesoint,mesoi)
   _GET_HORIZONTAL_(self%id_microint,microi)
   _GET_HORIZONTAL_(self%id_det1int,deti1)
   _GET_HORIZONTAL_(self%id_det2int,deti2)
   _GET_HORIZONTAL_(self%id_fishint,fishi)
   _GET_HORIZONTAL_(self%id_fish2int,fishi2)
   _GET_HORIZONTAL_(self%id_f1zl,F1onZl)
   _GET_HORIZONTAL_(self%id_f1zs,F1onZs)
   _GET_HORIZONTAL_(self%id_f1det1,F1onDet1)
   _GET_HORIZONTAL_(self%id_f1det2,F1onDet2)
   _GET_HORIZONTAL_(self%id_f2zl,F2onZl)
   _GET_HORIZONTAL_(self%id_f2f1,F2onF1)
   _GET_HORIZONTAL_(self%id_f2det1,F2onDet1)
   _GET_HORIZONTAL_(self%id_f2det2,F2onDet2)
   _GET_(self%id_fish2,fish2)
   _GET_(self%id_fish1,fish1)
   end if

   ! CAGLAR
   ! checks - whether the biomass of plankton is below a predefined threshold,
   !          where below the threshold, loss terms are removed from the RHS of
   !          the equations. The idea is to keep plankton safe from extinction.
   ! loss terms are multiplied by the constants below, which can only be set
   ! by the model to 0 or 1.

   fla_loss = max(sign(-1.0_rk,fla-0.01_rk),0.0_rk)       ! flagellates
   dia_loss = max(sign(-1.0_rk,dia-0.01_rk),0.0_rk)       ! diatoms
   bg_loss  = max(sign(-1.0_rk,bg-0.01_rk),0.0_rk)        ! cyanobacteria
   mic_loss = max(sign(-1.0_rk,microzoo-0.001_rk),0.0_rk) !microzooplankton
   mes_loss = max(sign(-1.0_rk,mesozoo-0.001_rk),0.0_rk) ! mesozooplankton
!c-----------------t dependy for fish & MB  ------------------------
        ttemp=temp/((temp+273.15)*273.15)

   ! remineralisation rate
   frem = self%BioC(22) * (1._rk+20._rk*(temp**2/(13._rk**2+temp**2)))
   fremDOM = 10._rk * frem

   ! nutrient limitation factors
   ! k denotes half-saturation values
   up_nh4 = nh4/(self%BioC(6)+nh4)
   up_no3 = no3/(self%BioC(7)+no3)*exp(-self%BioC(8)*nh4)
   up_n = up_nh4+up_no3
   up_pho = pho/(self%BioC(25)+pho)
   t_sil = max(sil-80._rk,0.0_rk)
   up_sil = t_sil/(self%BioC(26)+t_sil)

   ! light limitation
   blight = max(tanh(self%BioC(3)*par),0.0_rk)

   ! temperature dependence
   Ts = exp(self%T_lim_phy * temp)
   Tl = exp(self%T_lim_phy * temp)
   if ((salt<=10.0) .and. (mean_surface_par > self%bg_growth_minimum_daily_rad)) then
     Tbg = 1.0_rk/(1.0_rk + exp(self%BioC(29)*(self%BioC(30)-temp)))
   else
     Tbg = 0.0_rk
   end if

   ! production and nutrient uptake
   Ps_prod = Ts * min(blight, up_n, up_pho)
   Pl_prod = Tl * min(blight, up_n, up_pho, up_sil)
   Prod = self%BioC(1)*Pl_prod*dia + & ! diatoms production
          self%BioC(2)*Ps_prod*fla     ! flagellates production

   if (self%use_cyanos) then
     Bg_prod = Tbg * min(blight, up_n, up_pho)
     if (mean_par > self%nfixation_minimum_daily_par) then
       Bg_fix = Tbg * min(blight, up_pho) - Bg_prod
     end if
     Prod = Prod + self%BioC(28)*Bg_prod*bg ! cyanobacteria production
    else
    Bg_prod =  0.0_rk
    Bg_fix  = 0.0_rk
   end if

   if (self%use_chl) then   
     ! chlorophyll-a to C change
     chl2c_fla = self%MAXchl2nPs * max(0.1,Ps_prod) * self%BioC(2) * sedy0 * fla / &
               (self%alfaPs * par * flachl)
     chl2c_dia = self%MAXchl2nPl * max(0.1,Pl_prod) * self%BioC(1) * sedy0 * dia / &
               (self%alfaPl * par * diachl)

            chl2c_fla = max(self%MINchl2nPs,chl2c_fla)
            chl2c_fla = min(self%MAXchl2nPs,chl2c_fla)
            chl2c_dia = max(self%MINchl2nPl,chl2c_dia)
            chl2c_dia = min(self%MAXchl2nPl,chl2c_dia)
     if (self%use_cyanos) then
       chl2c_bg = self%MAXchl2nBG * max(0.1,Bg_prod) * self%BioC(28) * sedy0 * bg / &
               (self%alfaBG * par * bgchl)
       chl2c_bg  = max(self%MINchl2nBG,chl2c_bg)
       chl2c_bg  = min(self%MAXchl2nBG,chl2c_bg)
     end if
   end if

   ! grazing
   Fs = self%prefZsPs*fla + self%prefZsPl*dia + self%prefZsD*det1 + self%prefZsD*det2 !+ self%prefZsBG*bg
   Fl = self%prefZlPs*fla + self%prefZlPl*dia + self%prefZlZs*microzoo + &
         self%prefZlD*det1 + self%prefZlD*det2 !+ self%prefZlBG*bg
   if (self%use_cyanos) then
    Fs = Fs + self%prefZsBG*bg
    Fl = Fl + self%prefZlBg*bg
   end if

   ZsonPs = fla_loss * self%BioC(12) * self%prefZsPs * fla/(self%BioC(14)  + Fs)
   ZsonPl = dia_loss * self%BioC(12) * self%prefZsPl * dia/(self%BioC(14)  + Fs)
   ZsonD1 =            self%BioC(12) * self%prefZsD  * det1/(self%BioC(14)  + Fs)
   ZsonD2 =            self%BioC(12) * self%prefZsD  * det2/(self%BioC(14) + Fs) ! Add for sympagic system coupling
   

   ZlonPs = fla_loss * self%BioC(11) * self%prefZlPs * fla/(self%BioC(14)  + Fl)
   ZlonPl = dia_loss * self%BioC(11) * self%prefZlPl * dia/(self%BioC(14)  + Fl)
   ZlonD1 =            self%BioC(11) * self%prefZlD  * det1/(self%BioC(14)  + Fl)
   ZlonD2 =            self%BioC(11) * self%prefZlD  * det2/(self%BioC(14) + Fl) ! Add for sympagic system coupling
   ZlonZs = mic_loss * self%BioC(13) * self%prefZlZs * microzoo/(self%BioC(14) + Fl)
   if (self%use_cyanos) then
     ZsonBg = bg_loss  * self%BioC(31) * self%prefZsBG * bg/(self%BioC(14) + Fs)
     ZlonBg = bg_loss  * self%BioC(31) * self%prefZlBG * bg/(self%BioC(14) + Fl)
   else
     ZsonBg=0.0_rk
     ZlonBg=0.0_rk
   end if

   ! nitrification
   Onitr = 0.01_rk * redf(7) !according to Neumann  (Onitr in mlO2/l see also Stigebrand and Wulff)
   bioom1 = 0.0_rk
   bioom2 = 0.0_rk
   bioom3 = 0.0_rk
   bioom4 = 0.0_rk
   bioom5 = 0.0_rk
   bioom6 = 0.0_rk
   bioom7 = 0.0_rk
   bioom8 = 0.0_rk
   if (oxy > 0) then
     bioom1 = 0.1_rk/secs_pr_day * exp(temp*0.11_rk) * oxy/(Onitr+oxy)
     bioom2 = bioom1
     bioom6 = 1.0_rk
   else
     if (no3>0) then
       bioom5 = 5.0_rk
       bioom8 = 1.0_rk
     else
       bioom7 = 1.0_rk
     end if
   end if

! reaction rates

   _SET_ODE_(self%id_fla, (self%BioC(2)*Ps_prod - self%BioC(10)*fla_loss)*fla - ZsonPs*microzoo - ZlonPs*mesozoo)
   _SET_ODE_(self%id_dia, (self%BioC(1)*Pl_prod - self%BioC(9)*dia_loss)*dia - ZsonPl*microzoo - ZlonPl*mesozoo)
   if (self%use_cyanos) then
     _SET_ODE_(self%id_bg,  (self%BioC(28)*(Bg_prod + Bg_fix) - self%BioC(32)*bg_loss)*bg - ZsonBg*microzoo - ZlonBg*mesozoo)
   end if

  ! for chlorophyll-a
   if (self%use_chl) then
     rhs = self%BioC(2)*Ps_prod*chl2c_fla*fla - ( (self%BioC(10)*fla_loss*fla + ZsonPs*microzoo + ZlonPs*mesozoo)*flachl/fla )
     _SET_ODE_(self%id_flachl,rhs)
     rhs = self%BioC(1)*Pl_prod*chl2c_dia*dia - ( (self%BioC(9)*dia*dia_loss + ZsonPl*microzoo + ZlonPl*mesozoo)*diachl/dia )
     _SET_ODE_(self%id_diachl,rhs)
     if (self%use_cyanos) then
       rhs = self%BioC(28)*(Bg_prod + Bg_fix)*chl2c_bg*bg - ((self%BioC(32)*bg*bg_loss + ZsonBg*microzoo + ZlonBg*mesozoo)*bgchl/bg )
       _SET_ODE_(self%id_bgchl,rhs)
     end if  
   end if

   ! microzooplankton

   Zs_prod = self%BioC(20)*(ZsonPs + ZsonPl + ZsonBg) + self%BioC(21)*ZsonD1 + self%BioC(21)*ZsonD2 ! modif sympagic system coupling
   rhs = (Zs_prod - (self%BioC(16) + self%BioC(18) + self%zpr)*mic_loss) * microzoo &
         - ZlonZs * mesozoo
    if (self%use_fish) then
     rhs=rhs - F1onZs*microzoo/microi*fishi  !fish1 predation
!        - F2onZs*microzoo/microi*fishi2  !fish2 predation
    end if
  
    _SET_ODE_(self%id_microzoo, rhs)

   ! mesozooplankton
   Zl_prod = self%BioC(19)*(ZlonPs + ZlonPl + ZlonBg + ZlonZs) + self%BioC(21)*ZlonD1 + self%BioC(21)*ZlonD2 !modif sympagic system coupling
   rhs = (Zl_prod - (self%BioC(15) + self%BioC(17) + self%zpr)*mes_loss) * mesozoo
    if(self%use_fish)then
    rhs=rhs- F1onZl*mesozoo/mesoi*fishi &
          - F2onZl*mesozoo/mesoi*fishi2
    end if
   _SET_ODE_(self%id_mesozoo, rhs)

    if(self%use_fish)then
   !fish1
   Fish_prod=0.
   Fish_prod=self%asefF1*(F1onZs*microzoo/microi  &
             +F1onZl*mesozoo/mesoi  &
             +F1onDet1*det1/deti1 &
 	     +F1onDet2*det2/deti2)*fishi   !6,7,8,9
   rhs=Fish_prod-(self%mF1+self%excF1*exp((self%tempcF1/(8.6173324*10.**(-5.)))*ttemp))*fish1  & 
         - F2onF1*fish1/fishi*fishi2  !fish1 predation
   _SET_ODE_(self%id_fish1, rhs)
   !fish2
   Fish2_prod=0.
   Fish2_prod=self%asefF2*(F2onF1*fish1/fishi  &
             +F2onZl*mesozoo/mesoi  &
             +F2onDet1*det1/deti1 &
	     +F2onDet2*det2/deti2)*fishi2   !6,7,8,9
!     print*,'test fish',Fish2_prod,rhs,fishi2
   rhs=Fish2_prod-(self%mF2+self%excF2*exp((self%tempcF2/(8.6173324*10.**(-5.)))*ttemp))*fish2   
   _SET_ODE_(self%id_fish2, rhs)
!     print*,'test fish',Fish_prod,rhs,fish2
    end if

   ! detritus
   dxxdet1 = (  ((1.0_rk-self%BioC(21))*ZsonD1) * microzoo &
              + ((1.0_rk-self%BioC(21))*ZlonD1) * mesozoo  )

   dxxdet2 = (  ((1.0_rk-self%BioC(20))*(ZsonPs + ZsonPl + ZsonBg) &
              + (1.0_rk-self%BioC(21))*ZsonD2)  * microzoo &
              + ((1.0_rk-self%BioC(19))*(ZlonPs + ZlonPl + ZlonBg + ZlonZs) &
              + (1.0_rk-self%BioC(21))*ZlonD2) * mesozoo &
              + self%BioC(16) * microzoo * mic_loss &
              + self%BioC(15) * mesozoo * mes_loss &
              + self%BioC(9)  * dia * dia_loss &
              + self%BioC(10) * fla * fla_loss)

   if (self%use_cyanos) then
              dxxdet2 = dxxdet2 + (self%BioC(32) * bg * bg_loss )
   end if
   if (self%use_fish)then
   dxxdet1=dxxdet1 + (1.-self%asefF1)*(F1onDet1*det1/deti1)*fishi  &      !Fish 1 detritus
              + (1.-self%asefF2)*(F2onDet1*det1/deti1)*fishi2      !Fish 1 detritus

   dxxdet2=dxxdet2 + (1.-self%asefF1)*(F1onZs*microzoo/microi+F1onZl*mesozoo/mesoi+F1onDet2*det2/deti2)*fishi  &      !Fish 1 detritus
              + self%mF1*fish1 & 
              + (1.-self%asefF2)*(F2onF1*fish1/fishi+F2onZl*mesozoo/mesoi+F2onDet2*det2/deti2)*fishi2  &      !Fish 1 detritus
              + self%mF2*fish2
   end if

   rhs = (1.0_rk-self%frr) * dxxdet1 &
         - ZsonD1 * microzoo &
         - ZlonD1 * mesozoo &
         - frem * det1

    if (self%use_fish)then 
    rhs=rhs - F1onDet1*det1/deti1*fishi &
         - F2onDet1*det1/deti1*fishi2
    end if

   _SET_ODE_(self%id_det1, rhs)

   rhs = (1.0_rk-self%frr) * dxxdet2 &
         - ZsonD2 * microzoo &
         - ZlonD2 * mesozoo &
         - frem * det2

    if (self%use_fish)then 
    rhs=rhs - F1onDet2*det2/deti2*fishi &
         - F2onDet2*det2/deti2*fishi2
    end if

   _SET_ODE_(self%id_det2, rhs)

   if (self%community_export) then
   ! community dependent sinking rate
   ! 
   ! dsnk/det is the calculated (and applied) detritus sinking speed
   ! assumptions: unassimilated food uses detritus sinking rates from the prey,
   ! not the predator. 
     dxxdet1 = (  ((1.0_rk-self%BioC(21))* ZsonD1 * dsnk/det1) * microzoo &
              +   ((1.0_rk-self%BioC(21))* ZlonD1 * dsnk/det1) * mesozoo)
 
     dxxdet2 = (  ((1.0_rk-self%BioC(20))*(ZsonPs * self%sinkFlaD + ZsonPl * self%sinkDiaD + ZsonBg * self%sinkBgD) &
              + (1.0_rk-self%BioC(21)) * ZsonD2 * dsnk2/det2) * microzoo &
              + ((1.0_rk-self%BioC(19))*(ZlonPs * self%sinkFlaD + ZlonPl * self%sinkDiaD + ZlonBg * self%sinkBgD + ZlonZs * self%sinkMicD) &
              + (1.0_rk-self%BioC(21)) * ZlonD2 * dsnk2/det2) * mesozoo &
              + self%BioC(16) * microzoo * mic_loss * self%sinkMicD &
              + self%BioC(15) * mesozoo * mes_loss * self%sinkMesD &
              + self%BioC(10) * fla * fla_loss * self%sinkFlaD &
              + self%BioC(9)  * dia * dia_loss * self%sinkDiaD )
       if (self%use_cyanos) then
              dxxdet2 = dxxdet2 + (self%BioC(32) * bg * bg_loss * self%sinkBgD)
       end if 
       if (self%use_fish)then
          dxxdet1 = dxxdet1 &
                   + (1.-self%asefF1)*(F1onDet1 * dsnk/det1 * det1/deti1)*fishi &
                   + (1.-self%asefF2)*(F2onDet1 * dsnk/det1 * det1/deti1)*fishi2

          dxxdet2 = dxxdet2 &
                   + (1.-self%asefF1)*(F1onZs * self%sinkMicD * microzoo/microi + F1onZl * self%sinkMesD * mesozoo/mesoi + F1onDet2 * dsnk2/det2 * det2/deti2)*fishi &
                   + self%mF1 * fish1 * self%sinkFish1D &
                   + (1.-self%asefF2)*(F2onF1 * self%sinkFish1D * fish1/fishi + F2onZl * self%sinkDiaD * mesozoo/mesoi + F2onDet2 * dsnk2/det2 * det2/deti2)*fishi2 &
                   + self%mF2 * fish2 * self%sinkFish2D
       end if

     rhs = (1.0_rk-self%frr) * dxxdet1 &
         + ( - ZsonD1 * microzoo - ZlonD1 * mesozoo ) * dsnk/det1 &
             - frem * dsnk

        if (self%use_fish) then
           rhs = rhs + ( - F1onDet1*det1/deti1*fishi - F2onDet1*det1/deti1*fishi2 ) * dsnk/det1
        end if

     _SET_ODE_(self%id_dsnk, rhs)

     rhs = (1.0_rk-self%frr) * dxxdet2 &
         + ( - ZsonD2 * microzoo - ZlonD2 * mesozoo ) * dsnk2/det2 &
             - frem * dsnk2

        if (self%use_fish) then
           rhs = rhs + ( - F1onDet2*det2/deti2*fishi - F2onDet2*det2/deti2*fishi2 ) * dsnk2/det2
        end if

     _SET_ODE_(self%id_dsnk2, rhs)

!    dxxdet = (  ((1.0_rk-self%BioC(20))*(ZsonPs + ZsonPl + ZsonBg) &
!              + (1.0_rk-self%BioC(21))*ZsonD) * microzoo * self%sinkMicD &
!              + ((1.0_rk-self%BioC(19))*(ZlonPs + ZlonPl + ZlonBg + ZlonZs) &
!              + (1.0_rk-self%BioC(21))*ZlonD) * mesozoo * self%sinkMesD &
!              + self%BioC(16) * microzoo * mic_loss * self%sinkMicD &
!              + self%BioC(15) * mesozoo * mes_loss * self%sinkMesD &
!              + self%BioC(10) * fla * fla_loss * self%sinkDiaD &
!              + self%BioC(9)  * dia * dia_loss * self%sinkFlaD ) 
!
!     rhs = (1.0_rk-self%frr) * dxxdet &
!         + ( - ZsonD * microzoo - ZlonD * mesozoo ) * dsnk/det &
!             - frem * dsnk
!
!     _SET_ODE_(self%id_dsnk, rhs)


   end if 

   ! labile dissolved organic matter
   _SET_ODE_(self%id_dom, self%frr*(dxxdet1 + dxxdet2) - fremdom * dom)

   ! nitrate
   rhs_nit = -(up_no3+0.5d-10)/(up_n+1.0d-10)*Prod &
         + bioom1 * nh4 &
         - bioom3 * no3 &
         - frem * (det1 + det2) * bioom5 &
         - fremDOM * dom * bioOM5
   _SET_ODE_(self%id_no3, rhs_nit)

   ! ammonium
   rhs_amm = -(up_nh4+0.5d-10)/(up_n+1.0d-10)*Prod &
         + self%BioC(18) * microzoo * mic_loss &
         + self%BioC(17) * mesozoo * mes_loss &
         + frem * (det1 + det2) &
         + fremDOM * dom - bioom1 * nh4
    if(self%use_fish)then
    rhs_amm = rhs_amm + self%excF1*exp((self%tempcF1/(8.6173324*10.**(-5.)))*ttemp)*fish1 &
         + self%excF2*exp((self%tempcF2/(8.6173324*10.**(-5.)))*ttemp)*fish2
    end if

   _SET_ODE_(self%id_nh4, rhs_amm)

   ! phosphate

   rhs = -Prod -self%BioC(28)*bg*Bg_fix &
         + self%BioC(18) * microzoo * mic_loss &
         + self%BioC(17) * mesozoo * mes_loss &
         + frem*(det1+det2) + fremDOM*dom
    if(self%use_fish)then
    rhs = rhs + self%excF1*exp((self%tempcF1/(8.6173324*10.**(-5.)))*ttemp)*fish1 &
         + self%excF2*exp((self%tempcF2/(8.6173324*10.**(-5.)))*ttemp)*fish2
    end if
   _SET_ODE_(self%id_pho, rhs)


   ! silicate
   _SET_ODE_(self%id_sil, -self%BioC(1)*Pl_prod*dia + self%BioC(27)*opa)

   ! opal
   _SET_ODE_(self%id_opa, self%BioC(9)*dia*dia_loss + ZsonPl*microzoo + ZlonPl*mesozoo - self%BioC(27)*opa)

   ! oxygen
  
    if(self%use_fish)then
     rhs_oxy = ((6.625*up_nh4 + 8.125*up_no3+1.d-10)/(up_n+1.d-10)*Prod &
         -bioom6*6.625*(self%BioC(18)*microzoo*mic_loss &
         +self%BioC(17)*mesozoo*mes_loss &
         + self%excF1*exp((self%tempcF1/(8.6173324*10.**(-5.)))*ttemp)*fish1 &
         + self%excF2*exp((self%tempcF2/(8.6173324*10.**(-5.)))*ttemp)*fish2) &
         -frem*(det1+det2)*(bioom6+bioom7)*6.625 &
         -(bioom6+bioom7)*6.625*fremDOM*dom &
         -2.0_rk*bioom1*nh4)*redf(11)*redf(16)
    else
    rhs_oxy = ((6.625*up_nh4 + 8.125*up_no3+1.d-10)/(up_n+1.d-10)*Prod &
         -bioom6*6.625*(self%BioC(18)*microzoo*mic_loss &
         +self%BioC(17)*mesozoo*mes_loss) &
         -frem*(det1+det2)*(bioom6+bioom7)*6.625 &
         -(bioom6+bioom7)*6.625*fremDOM*dom &
         -2.0_rk*bioom1*nh4)*redf(11)*redf(16)
    end if

   _SET_ODE_(self%id_oxy, rhs_oxy)

   ! Carbonate dynamics
   if (self%couple_co2) then
     rhs =  redf(16) *( BioC(18)*microzoo*mic_loss &
            + BioC(17)*mesozoo*mes_loss &
            + frem*(det1+det2) + fremDOM*dom - Prod) 
     _SET_ODE_(self%id_dic, rhs)

     rhs = redf(16)*(rhs_amm-rhs_nit)*redf(11) - 0.5_rk*rhs_oxy*(1._rk-bioom6)
     _SET_ODE_(self%id_alk, rhs)
   end if

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit,(frem*(det1+det2)*bioom5+fremDOM*dom*bioom5)*redf(11)*redf(16))
   _SET_DIAGNOSTIC_(self%id_primprod, Prod + self%BioC(28)*bg*Bg_fix )
   _SET_DIAGNOSTIC_(self%id_secprod, Zl_prod*mesozoo + Zs_prod*microzoo)
   _SET_DIAGNOSTIC_(self%id_parmean_diag, mean_par)
   _SET_DIAGNOSTIC_(self%id_snkspd,dsnk/det1*86400.) ! multiplication by 86400.
!                                                     should not be there, but during development, 
!                                                     it is handy to have the diagnostics
!                                                     in m/days. Will switch to m/sec when development complete.
   _SET_DIAGNOSTIC_(self%id_snkspd2,dsnk2/det2*86400.) 
!   if (self%use_chl) then
!     _SET_DIAGNOSTIC_(self%id_c2chl_fla, 1.0_rk/chl2c_fla)
!     _SET_DIAGNOSTIC_(self%id_c2chl_dia, 1.0_rk/chl2c_dia)
!     if (self%use_cyanos) then
!       _SET_DIAGNOSTIC_(self%id_c2chl_bg, 1.0_rk/chl2c_bg)
!     end if
!   end if
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for the ecosmo model
!
! !INTERFACE:

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_nersc_ecosmo_ice),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   real(rk) :: o2flux, T, tr, S, o2sat, oxy
   real(rk) :: no3flux, phoflux
   real(rk) :: pho,par,bg,blight,tbg,up_pho,prod
!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,T)
   _GET_(self%id_salt,S)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_par,par)
   _GET_(self%id_pho,pho)
   if (self%use_cyanos) then
     _GET_(self%id_bg,bg)
   else
     bg=0.0_rk
   end if

   ! Oxygen saturation micromol/liter__(Benson and Krause, 1984)
   tr = 1.0_rk/(T + 273.15_rk)
   o2sat= exp(- 135.90205_rk              &
       + (1.575701d05 ) * tr               &
       - (6.642308d07 ) * tr**2            &
       + (1.243800d10) * tr**3            &
       - (8.621949d11) * tr**4            &
       - S*(0.017674_rk-10.754_rk*tr+2140.7_rk*tr**2)  )

   o2flux = 5._rk/secs_pr_day * (o2sat - oxy)

   _SET_SURFACE_EXCHANGE_(self%id_oxy,o2flux)

   _SET_SURFACE_EXCHANGE_(self%id_no3,self%surface_deposition_no3)
   _SET_SURFACE_EXCHANGE_(self%id_nh4,self%surface_deposition_nh4)
   _SET_SURFACE_EXCHANGE_(self%id_pho,self%surface_deposition_pho)
   _SET_SURFACE_EXCHANGE_(self%id_sil,self%surface_deposition_sil)

#if 0
   if (self%use_cyanos) then
     ! calculate cyanobacteria surface production
     if (S <= 10.0) then
       tbg = 1.0_rk/(1.0_rk + exp(self%BioC(29)*(self%BioC(30)-T)))
     else
       tbg = 0.0_rk
     end if

     blight=max(tanh(self%BioC(3)*par),0.)
     up_pho = pho/(self%BioC(25)+pho)
     prod = self%BioC(28) * bg * Tbg * min(blight, up_pho) ! cyanobacteria production

   !_SET_ODE_(self%id_bg,  prod)
   !_SET_ODE_(self%id_pho, -prod)
   !_SET_SURFACE_ODE_(id_oxy, ) ! not included in the modular ECOSMO version
   !_SET_SURFACE_ODE_(id_dic, -Prod)
   end if
#endif

   ! Leave spatial loops over the horizontal domain (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Bottom fluxes for the ecosmo model
!
! !INTERFACE:

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_nersc_ecosmo_ice),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk) :: temp, tbs, oxy, no3, det1, det2, opa, sed1, sed2, sed3
   real(rk) :: pho, Rds, Rsd, Rsa, Rsdenit, Rsa_p, yt1, yt2
   real(rk) :: rhs, flux, alk_flux
   real(rk) :: bioom1, bioom2, bioom3, bioom4, bioom5, bioom6, bioom7, bioom8
! ute mb new variables
   real(rk) :: mb1
   real(rk) :: FM,MBonMicro,MBonMeso,MBonDet1,MBonDet2,MBonSed,MBonFla,MBonDia,MBonDom
   real(rk) :: dia, fla, microzoo,mesozoo,dom 
   real(rk) :: ttemp, rcMB,rrcMB,MBflux,MBflux1,MBflux2 
! ute add fish local variables integrated values
   real(rk) :: mesoi,microi,deti1,deti2,fishi,fishi1,fishi2,fishi3
   real(rk) :: Ff1,F1onZl,F1onZs,F1onDet1,F1onDet2,F1onMB1
   real(rk) :: Ff2,F2onZl,F2onF1,F2onDet1,F2onDet2,F2onMB1
! add community sinking local variables
   real(rk) :: dsnk, dsnk2
!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,temp)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_det1,det1)
   _GET_(self%id_det2,det2)
   _GET_(self%id_opa,opa)
   _GET_(self%id_no3,no3)
   if (self%community_export) then
     _GET_(self%id_dsnk,dsnk)
     _GET_(self%id_dsnk2,dsnk2)
   end if
   _GET_HORIZONTAL_(self%id_sed1,sed1)
   _GET_HORIZONTAL_(self%id_sed2,sed2)
   _GET_HORIZONTAL_(self%id_sed3,sed3)
   _GET_HORIZONTAL_(self%id_tbs,tbs)
   if (self%use_fish)then
    !  ute get vars for MB
   _GET_(self%id_microzoo,microzoo)
   _GET_(self%id_mesozoo,mesozoo)
   _GET_(self%id_dia,dia)
   _GET_(self%id_fla,fla)
   _GET_(self%id_dom,dom)
   _GET_HORIZONTAL_(self%id_mb1,mb1)
  !ute get dependencies for fish growth
   _GET_HORIZONTAL_(self%id_mesoint,mesoi)
   _GET_HORIZONTAL_(self%id_microint,microi)
   _GET_HORIZONTAL_(self%id_det1int,deti1)
   _GET_HORIZONTAL_(self%id_det2int,deti2)
   _GET_HORIZONTAL_(self%id_fishint,fishi)
   _GET_HORIZONTAL_(self%id_fish2int,fishi2)
   end if

   bioom1 = 0.0_rk
   bioom2 = 0.0_rk
   bioom3 = 0.0_rk
   bioom4 = 0.0_rk
   bioom5 = 0.0_rk
   bioom6 = 0.0_rk
   bioom7 = 0.0_rk
   bioom8 = 0.0_rk
   if (oxy > 0) then
     bioom1 = 0.1/secs_pr_day * exp(temp*0.11_rk) * oxy/((0.1_rk*redf(7))+oxy)
     bioom2 = bioom1
     bioom6 = 1.0_rk
   else
     if (no3>0) then
       bioom5 = 5.0_rk
       bioom8 = 1.0_rk
     else
       bioom7 = 1.0_rk
     end if
   end if

    if (self%use_fish)then
    !ute grazing fish
    Ff1=self%prefF1Zl*mesoi+self%prefF1Zs*microi+self%prefF1Det*deti1+self%prefF1Det*deti2+self%prefF1MB1*mb1
!    Ff2=self%prefF2Zl*mesoi+self%prefF2Zs*microi+self%prefF2Det*deti+self%prefF2MB1*mb1
    Ff2=self%prefF2Zl*mesoi+self%prefF2F1*fishi+self%prefF2Det*deti1+self%prefF2Det*deti2+self%prefF2MB1*mb1
!
     F1onZl=0.
     F1onZs=0.
     F1onDet1=0.
     F1onDet2=0.
     F1onMB1=0.
    if(oxy .gt. 0.0)then
       if(mesoi.ge.0.) F1onZl=self%prefF1Zl*self%GrF1Zl*mesoi/(self%rF1+Ff1)
       if(microi.ge.0.) F1onZs=self%prefF1Zs*self%GrF1Zs*microi/(self%rF1+Ff1)
       if(deti1.ge.0.) F1onDet1 =self%prefF1Det*self%GrF1Det*deti1/(self%rF1+Ff1)
       if(deti2.ge.0.) F1onDet2 =self%prefF1Det*self%GrF1Det*deti2/(self%rF1+Ff1)
       if(mb1.ge.0.)F1onMB1=self%prefF1MB1*self%GrF1MB1*mb1/(self%rF1MB+Ff1)
       if(mesoi.ge.0.) F2onZl=self%prefF2Zl*self%GrF2Zl*mesoi/(self%rF2+Ff2)
       if(fishi.ge.0.) F2onF1=self%prefF2F1*self%GrF2F1*fishi/(self%rF2+Ff2)
       if(deti1.ge.0.) F2onDet1 =self%prefF2Det*self%GrF2Det*deti1/(self%rF2+Ff2)
       if(deti2.ge.0.) F2onDet2 =self%prefF2Det*self%GrF2Det*deti2/(self%rF2+Ff2)
       if(mb1.ge.0.)F2onMB1=self%prefF2MB1*self%GrF2MB1*mb1/(self%rF2MB+Ff2)
     endif
!c-----------------t dependy for fish & MB  -------------------------
        ttemp=temp/((temp+273.15)*273.15)
!--------------------------------------------------------------
! Ute macrobenthos feeding preferences
     FM=self%prefMB1Zs*microzoo+self%prefMB1Zl*mesozoo   &
        +self%prefMB1Det*det1+self%prefMB1Det*det2+self%prefMB1Sed*sed1        &
        +self%prefMB1P*fla+self%prefMB1P*dia+self%prefMB1Dom*dom
        MBonMicro=0.
        MBonMeso =0.
        MBonDet1 =0.
        MBonDet2 =0.
        MBonSed  =0.
        MBonFla  =0.
        MBonDia  =0.
        MBonDom  =0.

       if(oxy .gt. 0.0)then
       if(microzoo.ge.0.) MBonMicro=self%prefMB1Zs*self%GrMB1Z*microzoo/(self%rMB1+FM)
       if(mesozoo.ge.0.) MBonMeso=self%prefMB1Zl*self%GrMB1Z*mesozoo/(self%rMB1+FM)
       if(det1.ge.0.) MBonDet1 =self%prefMB1Det*self%GrMB1Det*det1/(self%rMB1+FM)
       if(det2.ge.0.) MBonDet2 =self%prefMB1Det*self%GrMB1Det*det2/(self%rMB1+FM)
       if(sed1.ge.0.) MBonSed =self%prefMB1Sed*self%GrMB1Sed*sed1/(self%rMB1+FM)
       if(fla.ge.0.) MBonFla =self%prefMB1P*self%GrMB1P*fla/(self%rMB1+FM)
       if(dia.ge.0.) MBonDia =self%prefMB1P*self%GrMB1P*dia/(self%rMB1+FM)
       if(dom.ge.0.) MBonDom =self%prefMB1Dom*self%GrMB1Det*dom/(self%rMB1+FM)
       endif
      end if
!----citical bottom shear stress
        if (tbs.ge.self%BioC(34)) then
          Rsd=min(self%BioC(35), self%BioC(35) * (tbs**3/(0.1+tbs**3)))
          Rds=0.0_rk
        else if (tbs.lt.self%BioC(34)) then
          Rsd=0.0_rk
          Rds=self%BioC(36)
        end if

!---------------------------------------------------------------
!----denitrification parameter in dependence of available oxygen
        if (oxy .gt. 0.0) then
          Rsa=self%BioC(38)*exp(self%BioC(39)*temp)*1.0_rk
          Rsdenit=0.0_rk
        else if (oxy .le. 0.0) then
          Rsdenit=self%BioC(38)*exp(self%BioC(39)*temp)*2.0_rk
          Rsa=0.0_rk
        end if

        !--- sediment 1 total sediment biomass and nitrogen pool
        rhs = Rds*(det1+det2) - Rsd*sed1 - 2.0_rk*Rsa*sed1 - Rsdenit*sed1 &
              - self%BioC(37)*sed1
        _SET_BOTTOM_ODE_(self%id_sed1, rhs)

        ! community sinking variable exchange
        if (self%community_export) then
           _SET_BOTTOM_EXCHANGE_(self%id_dsnk, Rsd*sed1*dsnk/det1 - Rds*det1*dsnk/det1)   
           _SET_BOTTOM_EXCHANGE_(self%id_dsnk2, Rsd*sed1*dsnk2/det2 - Rds*det2*dsnk2/det2)   
        end if
        ! oxygen
        flux = -(BioOM6*6.625_rk*2.0_rk*Rsa*sed1 &
                 +BioOM7*6.625_rk*Rsdenit*sed1 &
                 +2.0_rk*BioOM1*Rsa*sed1) &
                *REDF(11)*REDF(16)
        _SET_BOTTOM_EXCHANGE_(self%id_oxy, flux)

        ! nitrate
        _SET_BOTTOM_EXCHANGE_(self%id_no3, -BioOM5*Rsdenit*sed1)

        ! detritus
        _SET_BOTTOM_EXCHANGE_(self%id_det1, Rsd*sed1 - Rds*det1)
        _SET_BOTTOM_EXCHANGE_(self%id_det2, Rsd*sed1 - Rds*det2)

        ! ammonium
        _SET_BOTTOM_EXCHANGE_(self%id_nh4, (Rsdenit+Rsa)*sed1)

        if (self%couple_co2) then
          _SET_BOTTOM_EXCHANGE_(self%id_dic, redf(16)*(Rsdenit+2*Rsa)*sed1)
          alk_flux = redf(16)*redf(11)*((Rsdenit+Rsa+bioom5*Rsdenit)*sed1) - 0.5_rk*flux*(1._rk-bioom6)
          _SET_BOTTOM_EXCHANGE_(self%id_alk, alk_flux)
        end if
        !--try out for phosphate ute 2.6.2010
        Rsa_p=self%BioC(38)*exp(self%BioC(39)*temp)*2.0_rk

        if (oxy.gt.0.0) then
          yt2=oxy/375.0_rk   !normieren des wertes wie in Neumann et al 2002
          yt1=yt2**2.0_rk/(self%BioC(41)**2.0_rk+yt2**2.0_rk)

          _SET_BOTTOM_EXCHANGE_(self%id_pho,Rsa_p*(1.0_rk-self%BioC(40)*yt1)*sed3)

          !--sed 3 phosphate pool sediment+remineralization-P release
          _SET_BOTTOM_ODE_(self%id_sed3, 2.0_rk*Rsa*sed1 - Rsa_p*(1.0_rk-self%BioC(40)*yt1)*sed3)

        else if (oxy.le.0.0) then
          _SET_BOTTOM_EXCHANGE_(self%id_pho, Rsa_p*sed3)
          _SET_BOTTOM_ODE_(self%id_sed3, Rsdenit*sed1 - Rsa_p*sed3)
        end if

        ! sediment opal(Si)
        _SET_BOTTOM_ODE_(self%id_sed2, Rds*opa - Rsd*sed2 - self%BioC(42)*sed2 - ( self%BioC(37)*1000.*(sed2**3/(sed2**3 + 1E+12)) )*sed2)
        _SET_BOTTOM_EXCHANGE_(self%id_opa, Rsd*sed2 - Rds*opa)
        _SET_BOTTOM_EXCHANGE_(self%id_sil, self%BioC(42)*sed2)

           _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tbsout, tbs)

     if(self%use_fish)then
       !Ute Macrobenthos & benthic fish
         _SET_BOTTOM_EXCHANGE_(self%id_fish1,+self%asefF1*F1onMB1*fishi)
         _SET_BOTTOM_EXCHANGE_(self%id_fish2,+self%asefF2*F2onMB1*fishi2)
        rrcMB=(MBonSed+MBonDet1+MBonDet2+MBonDom+MBonMicro+MBonMeso+MBonDia+MBonFla) 
        rcMB=(+self%BioC(21)*rrcMB  &
             -self%mMB1 &
           -self%excMB1*exp((self%tempcMB1/(8.6173324*10.**(-4.)))*ttemp))*mb1 &
             -F1onMB1*fishi  &
             -F2onMB1*fishi2  !
         _SET_BOTTOM_ODE_(self%id_mb1,rcMB)     
        MBflux=-MBonSed*mb1+(1.-self%frr)*((1.-self%BioC(21))*rrcMB+self%mMB1)*mb1
         _SET_BOTTOM_ODE_(self%id_sed1,MBflux)     
         _SET_BOTTOM_EXCHANGE_(self%id_microzoo,-MBonMicro*mb1)     
         _SET_BOTTOM_EXCHANGE_(self%id_mesozoo,-MBonMeso*mb1)     
         _SET_BOTTOM_EXCHANGE_(self%id_dia,-MBonDia*mb1)     
         _SET_BOTTOM_EXCHANGE_(self%id_fla,-MBonFla*mb1)     
        MBflux1=-MBonDet1*mb1 
	MBflux2=-MBonDet2*mb1  &
         +(1.-self%frr)*((1.-self%asefF1)*F1onMB1*fishi &
         +(1.-self%asefF2)*F2onMB1*fishi2)
          _SET_BOTTOM_EXCHANGE_(self%id_det2,MBflux2)
	  _SET_BOTTOM_EXCHANGE_(self%id_det1,MBflux1)

        if (self%community_export) then
           MBflux1=-MBonDet1 * mb1 * dsnk/det1 &
                  +(1.-self%frr)*((1.-self%asefF1)*F1onMB1*fishi*dsnk/det1+(1.-self%asefF2)*F2onMB1*fishi2*dsnk/det1) 
           _SET_BOTTOM_EXCHANGE_(self%id_dsnk,MBflux1)

           MBflux=-MBonDet2 * mb1 * dsnk2/det2 &
                  +(1.-self%frr)*((1.-self%asefF1)*F1onMB1*fishi*dsnk2/det2+(1.-self%asefF2)*F2onMB1*fishi2* dsnk2/det2) 
           _SET_BOTTOM_EXCHANGE_(self%id_dsnk2,MBflux2)
        end if     

         MBflux=-MBonDom*mb1+(self%frr)*(((1.-self%BioC(21))*rrcMB+self%mMB1)*mb1 &
        +(1.-self%asefF1)*F1onMB1*fishi  &  !fish    
        +(1.-self%asefF2)*F2onMB1*fishi2)    !fish    
         _SET_BOTTOM_EXCHANGE_(self%id_dom,MBflux)     
        MBflux=+self%excMB1*exp((self%tempcMB1/(8.6173324*10.**(-4.)))*ttemp)*mb1
         _SET_BOTTOM_EXCHANGE_(self%id_nh4,MBflux)     
         _SET_BOTTOM_EXCHANGE_(self%id_pho,MBflux) 
        MBflux =-(BioOM6*6.625_rk*self%excMB1* &
              exp((self%tempcMB1/(8.6173324*10.**(-4.)))*ttemp)*mb1) &
                *REDF(11)*REDF(16)
         _SET_BOTTOM_EXCHANGE_(self%id_oxy,MBflux) 
        fishi1=fishi    
        fishi3=fishi2    
 !ute diagnostics for fish
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_meso_integr,mesoi)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_micro_integr,microi)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_det1_integr,deti1)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_det2_integr,deti2)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fish_integr,fishi1)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fish2_integr,fishi3)
! diagnostics for dependency
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dF1onZl,F1onZl)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dF1onZs,F1onZs)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dF1onDet1,F1onDet1)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dF1onDet2,F1onDet2)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dF2onZl,F2onZl)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dF2onF1,F2onF1)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dF2onDet1,F2onDet1)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dF2onDet2,F2onDet2)
    end if

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
   class (type_nersc_ecosmo_ice), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_

   real(rk)                     :: dom,det1,det2,diachl,flachl,bgchl
   real(rk)                     :: dia,fla,bg
   real(rk)                     :: my_extinction

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.

   _GET_(self%id_det1, det1)
   _GET_(self%id_det2, det2)
   _GET_(self%id_dom, dom)
   _GET_(self%id_dia, dia)
   _GET_(self%id_fla, fla)
     if (self%use_cyanos) then
   _GET_(self%id_bg, bg)
     else
       bg=0.0_rk
     end if

   my_extinction = self%BioC(4)

   if (self%use_chl) then
     _GET_(self%id_diachl, diachl)
     _GET_(self%id_flachl, flachl)
     if (self%use_cyanos) then
       _GET_(self%id_bgchl, bgchl)
     else
       bgchl=0.0_rk
     end if
     my_extinction = my_extinction + self%BioC(5)*(diachl+flachl+bgchl)+ self%extdet*(det1+det2) + self%extdom*dom
   else
     diachl=0.0_rk
     flachl=0.0_rk
     my_extinction = my_extinction + self%BioC(5)*(dia+fla+bg) + self%extdet*(det1+det2) + self%extdom*dom
   end if

   _SET_EXTINCTION_( my_extinction )

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction

! ----- COMMUNITY DEPENDENT VARIABLE SINKING ------------------------------
   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_nersc_ecosmo_ice),intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

      real(rk) :: det1, det2, dsnk, dsnk2, meanspd, meanspd2

      _LOOP_BEGIN_
         if (self%community_export) then
           _GET_(self%id_det1, det1)
           _GET_(self%id_det2, det2)
           _GET_(self%id_dsnk, dsnk)
           _GET_(self%id_dsnk2, dsnk2)

           meanspd = dsnk / det1
           meanspd2 = dsnk2 / det2

           _SET_VERTICAL_MOVEMENT_(self%id_det1,-meanspd)
           _SET_VERTICAL_MOVEMENT_(self%id_det2,-meanspd2)
           _SET_VERTICAL_MOVEMENT_(self%id_dsnk,-meanspd)
           _SET_VERTICAL_MOVEMENT_(self%id_dsnk2,-meanspd2)
         else
           _SET_VERTICAL_MOVEMENT_(self%id_det1,-self%sinkIdet)
           _SET_VERTICAL_MOVEMENT_(self%id_det2,-self%BioC(23))
         endif
      _LOOP_END_
   end subroutine get_vertical_movement
! ------------------------------------------------------------------------- 

   end module fabm_nersc_ecosmo_ice

