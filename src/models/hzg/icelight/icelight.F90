#include "fabm_driver.h"

module hzg_icelight

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_hzg_icelight
      ! Identifiers for dependencies [model inputs]
      type (type_horizontal_dependency_id)          :: id_swr0 ! Surface shortwave radiation
      type (type_horizontal_dependency_id)          :: id_iceh ! Ice thickness
      type (type_horizontal_dependency_id)          :: id_icec ! Ice concentration
      type (type_horizontal_dependency_id)          :: id_snowh ! Snow thickness
      type (type_horizontal_dependency_id)          :: id_icealgea ! Ice algea concentration
      type (type_horizontal_dependency_id)          :: id_h_bal ! Ice biological active layer height
      type (type_dependency_id)                     :: id_dz   ! Cell thickness
      type (type_dependency_id)                     :: id_ext  ! Attentuation coefficient for PAR

      ! Identifiers for diagnostic variables [model outputs]
      type (type_diagnostic_variable_id)            :: id_par  ! Photosynthetically active radiation
      type (type_diagnostic_variable_id)            :: id_swr  ! Shortwave radiation
      type (type_horizontal_diagnostic_variable_id):: id_par0 ! Surface photosynthetically active radiation

      ! Parameters
      real(rk) :: a,g1,g2
      real(rk) :: gice
      real(rk) :: h_BAL, g_ia
      logical  :: use_icealgea
      logical  :: use_external_bal
      logical  :: use_snow
      real(rk) :: C0
      real(rk) :: Catt_ice
      real(rk) :: Catt_sno
      real(rk) :: ice_alb
      real(rk) :: sno_alb
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
      call self%get_parameter(self%a, 'a', '-','non-visible fraction of shortwave radiation',default=0.58_rk) 
      call self%get_parameter(self%g1,'g1','m','e-folding depth of non-visible fraction',    default=0.35_rk)
      call self%get_parameter(self%g2,'g2','m','e-folding depth of visible fraction', default=23.0_rk) 
      call self%get_parameter(self%gice,'g_ice','m','e-folding depth of radiation in ice', default=0.1_rk) 
      call self%get_parameter(self%use_snow,'use_snow','true/false','switch use of snow thickness', default=.false.) 
      call self%get_parameter(icealgea_name,'icealgea_name','-','name of icealgea variable', default='')
      self%use_icealgea = icealgea_name .ne. '' 
      call self%get_parameter(self%g_ia,'g_ia','m/(mmolC/m2)','e-folding depth of radiation by ice algea concentration', default=3.0_rk) 
      call self%get_parameter(h_bal_name,'h_bal_name','-','name of BAL variable', default='')
      self%use_external_bal = h_bal_name .ne. '' 
      call self%get_parameter(self%h_bal,'h_bal','m','height of biological active layer in ice', default=0.05_rk) 

      call self%get_parameter(self%C0,      'C0',      '',                       'part of incoming radiation absorbed',            default=0.3_rk)
      call self%get_parameter(self%Catt_ice,'Catt_ice','1/m',                    'ice attenuation coefficent',                     default=0.5_rk)
      call self%get_parameter(self%Catt_sno,'Catt_sno','1/m',                    'snow attenuation coefficent',                    default=0.5_rk)
      call self%get_parameter(self%ice_alb, 'ice_alb', '',                       'albedo in the ice',                              default=0.5_rk)
      call self%get_parameter(self%sno_alb, 'sno_alb', '',                       'albedo in the snow',                             default=0.5_rk)
      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_swr,'swr','W/m^2','shortwave radiation', &
              standard_variable=standard_variables%downwelling_shortwave_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_par,'par','W/m^2','photosynthetically active radiation', &
              standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_par0,'par0','W/m^2','surface photosynthetically active radiation', &
              standard_variable=standard_variables%surface_downwelling_photosynthetic_radiative_flux,source=source_do_column)
      !call self%register_diagnostic_variable(self%id_par0ice,'par0ice','W/m^2','surface photosynthetically active radiation', &
      !        standard_variable=standard_variables%surface_downwelling_photosynthetic_radiative_flux,source=source_do_column)

      ! Register environmental dependencies (temperature, shortwave radiation)
      !call self%register_dependency(self%id_iceh,       'iceh', 'm',    'ice thickness')
      call self%register_dependency(self%id_swr0,standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_ext, standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_dz,  standard_variables%cell_thickness)
      call self%register_dependency(self%id_iceh, standard_variables%ice_thickness)
      call self%register_dependency(self%id_icec, standard_variables%ice_concentration)
      if (self%use_snow) &
        call self%register_dependency(self%id_snowh, standard_variables%snow_thickness)
      if (self%use_icealgea) &
        call self%register_horizontal_dependency(self%id_icealgea, icealgea_name, 'mmolC/m2', 'ice algea mass in ice')
      if (self%use_external_bal) &
        call self%register_horizontal_dependency(self%id_h_bal, h_bal_name, 'm', 'height of BAL layer in ice')
   end subroutine
   
   subroutine get_light(self,_ARGUMENTS_VERTICAL_)
      class (type_hzg_icelight),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      real(rk) :: swr0,dz,swr,par,z,ext,bioext
      real(rk) :: iceh, icec, snowh, par0ice, swr0ice
      real(rk) :: icea, h_bal
      real(rk) :: alb
 
      _GET_HORIZONTAL_(self%id_swr0,swr0)
      _GET_HORIZONTAL_(self%id_iceh,iceh)
      _GET_HORIZONTAL_(self%id_icec,icec)
      if (self%use_icealgea) then
        _GET_HORIZONTAL_(self%id_icealgea,icea)
      else
        icea = 0.0_rk ! no shading by ice algea
      endif
      if (self%use_external_bal) then
        _GET_HORIZONTAL_(self%id_h_bal,h_bal)
      else
        h_bal = self%h_bal ! use parameter from namelist instead
      endif
      if (self%use_snow) then
        _GET_HORIZONTAL_(self%id_snowh,snowh)
      end if
      
      ! so far, snow thickness not used
      !k_ia scales the algea mass in the biological active layer into an extinction factor [1/m]
      swr0ice = swr0*(1-icec) + icec*swr0*exp(-iceh/self%gice-icea/self%g_ia*h_bal)
      par0ice = swr0ice*(1-self%a)
!      if(snowh > 0.0_rk) then 
! 	alb = self%sno_alb
!      else 
! 	alb = self%ice_alb
!      end if 

!      par0ice = (1-icec) * (1.0_rk - alb) * self%C0 * exp((-self%Catt_ice * iceh) - (self%Catt_sno * snowh))

      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_par0,par0ice)
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
       !  _SET_DIAGNOSTIC_(self%id_par0ice,par0ice) ! Photosynthetically active radiation at bottom ice
      _VERTICAL_LOOP_END_

   end subroutine get_light

end module hzg_icelight
