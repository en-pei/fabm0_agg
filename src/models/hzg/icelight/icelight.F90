#include "fabm_driver.h"

module hzg_icelight

   use fabm_types
   use fabm_expressions
 
   implicit none

   private

   type,extends(type_base_model),public :: type_hzg_icelight
      ! Identifiers for dependencies [model inputs]
      type (type_horizontal_dependency_id)          :: id_swr0     ! Surface shortwave radiation
      type (type_horizontal_dependency_id)          :: id_iceh     ! Ice thickness
      type (type_horizontal_dependency_id)          :: id_icec     ! Ice concentration
      type (type_horizontal_dependency_id)          :: id_snowh    ! Snow thickness
      type (type_surface_state_variable_id)         :: id_icealgea ! Ice algea concentration
      type (type_horizontal_dependency_id)          :: id_h_bal    ! Ice biological active layer height
      type (type_dependency_id)                     :: id_dz       ! Cell thickness
      type (type_dependency_id)                     :: id_ext      ! Attentuation coefficient for PAR

      ! Identifiers for diagnostic variables [model outputs]
      type (type_diagnostic_variable_id)            :: id_par      ! Photosynthetically active radiation
      type (type_diagnostic_variable_id)            :: id_swr      ! Shortwave radiation
      type (type_horizontal_diagnostic_variable_id) :: id_par0     ! Surface photosynthetically active radiation

      type (type_horizontal_diagnostic_variable_id)    :: id_swr0ice
      type (type_horizontal_diagnostic_variable_id)    :: id_lightice
 
      ! Parameters
      real(rk) :: a,g1,g2
      real(rk) :: gice
      real(rk) :: h_BAL, g_ia
      real(rk) :: C0
      real(rk) :: Catt_ice
      real(rk) :: Catt_sno
      real(rk) :: Catt_alg

      logical  :: couple_ice, couple_bgc

   contains
!     Model procedures
      procedure :: initialize
      procedure :: get_light
   end type type_hzg_icelight

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_hzg_icelight),intent(inout),target :: self
      integer,                  intent(in)           :: configunit
      character(len=265)                             :: icealgea_name,h_bal_name
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%a,          'a',             '-',           'non-visible fraction of shortwave radiation',default=0.58_rk) 
      call self%get_parameter(self%g1,         'g1',            'm',           'e-folding depth of non-visible fraction',    default=0.35_rk)
      call self%get_parameter(self%g2,         'g2',            'm',           'e-folding depth of visible fraction',        default=23.0_rk) 
      call self%get_parameter(self%gice,       'g_ice',         'm',           'e-folding depth of radiation in ice',        default=0.1_rk) 
      call self%get_parameter(self%g_ia,       'g_ia',          'm/(mmolC/m2)','e-folding depth of radiation by ice algea concentration', default=3.0_rk) 
      call self%get_parameter(self%C0,         'C0',            '',            'part of incoming radiation absorbed',        default=0.3_rk)
      call self%get_parameter(self%Catt_ice,   'Catt_ice',      '1/m',         'ice attenuation coefficent',                 default=0.5_rk)
      call self%get_parameter(self%Catt_sno,   'Catt_sno',      '1/m',         'snow attenuation coefficent',                default=0.5_rk)
      call self%get_parameter(self%Catt_alg,   'Catt_alg',      '1/m',         'algae attenuation coefficent',               default=0.5_rk)
      call self%get_parameter(self%couple_ice, 'couple_ice',    'true/false',  'switch use of ice forcing    ',              default=.true.)     
      call self%get_parameter(self%couple_bgc, 'couple_bgc',    'true/false',  'switch use of ice algae model',              default=.true.)     
      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_swr,'swr','W/m^2','shortwave radiation', standard_variable=standard_variables%downwelling_shortwave_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_par,'par','W/m^2','photosynthetically active radiation', standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_par0,'par0','W/m^2','surface photosynthetically active radiation under ice',  standard_variable=standard_variables%surface_downwelling_photosynthetic_radiative_flux,source=source_do_column)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_swr0, standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_ext,  standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_dz,   standard_variables%cell_thickness)
      call self%register_dependency(self%id_iceh, standard_variables%ice_thickness)
      call self%register_dependency(self%id_icec, standard_variables%ice_conc)
      call self%register_dependency(self%id_snowh, standard_variables%snow_thickness)
      call self%register_surface_state_dependency(self%id_icealgea, 'Ialg',  'mmolC/m2', 'ice algea mass in ice') 
      call self%register_dependency(self%id_h_bal,    'Zbot',    'm',        'height of BAL layer in ice')
      call self%register_horizontal_diagnostic_variable(self%id_lightice, 'lightice', '', &
           'light under ice', source=source_do_column)
      call self%register_horizontal_diagnostic_variable(self%id_swr0ice, 'lightocean', '', &
           'total ocean cell light', source=source_do_column)
   end subroutine
   
   subroutine get_light(self,_ARGUMENTS_VERTICAL_)
      class (type_hzg_icelight),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      real(rk) :: swr0,dz,swr,par,z,ext,bioext
      real(rk) :: iceh, icec, snowh
      real(rk) :: Ialg, Zbot
      real(rk) :: swr0ice, par0ice
 
      _GET_HORIZONTAL_(self%id_swr0,swr0)

      if (self%couple_ice) then
      _GET_HORIZONTAL_(self%id_iceh,iceh)
      _GET_HORIZONTAL_(self%id_snowh,snowh)
      _GET_HORIZONTAL_(self%id_icec,icec)
      else
      iceh  = 0.0_rk
      snowh = 0.0_rk
      icec  = 0.0_rk
      end if

      if (self%couple_bgc) then
        _GET_HORIZONTAL_(self%id_icealgea,Ialg)
        _GET_HORIZONTAL_(self%id_h_bal,Zbot)
      else
        Ialg  = 0.0_rk ! no shading by ice algea
        Zbot  = 0.0_rk ! use parameter from namelist instead
      endif

      ! so far, snow thickness not used
      !k_ia scales the algea mass in the biological active layer into an extinction factor [1/m]
      swr0ice = swr0*(1-icec) + icec*swr0*self%C0*exp((-self%Catt_ice * iceh)- (self%Catt_sno*snowh)-(self%g_ia*Ialg*Zbot))
      par0ice = swr0ice*(1-self%a)
      ! so far, snow thickness not used
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_par0,par0ice)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_swr0ice,swr0ice)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_lightice,icec*swr0*self%C0*exp((-self%Catt_ice * iceh)- (self%Catt_sno*snowh)-(self%g_ia*Ialg*Zbot)))
      

      z = 0
      bioext = 0
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_dz,dz)     ! Layer height (m)
         _GET_(self%id_ext,ext)   ! PAR attenuation (m-1)

         ! Set depth to centre of layer
         z = z + dz/2
         bioext = bioext + ext*dz/2

         ! Calculate photosynthetically active radiation (PAR), shortwave radiation, and PAR attenuation.
         par = swr0ice*(1-self%A)*exp(-z/self%g2-bioext)
         swr = par+swr0ice*self%A*exp(-z/self%g1)
         
         ! Move to bottom of layer
         z = z + dz/2
         bioext = bioext + ext*dz/2
         _SET_DIAGNOSTIC_(self%id_swr,swr) ! Shortwave radiation at layer centre
         _SET_DIAGNOSTIC_(self%id_par,par) ! Photosynthetically active radiation at layer centre
      _VERTICAL_LOOP_END_

   end subroutine get_light

end module hzg_icelight
