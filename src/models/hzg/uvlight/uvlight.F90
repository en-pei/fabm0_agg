#include "fabm_driver.h"

module hzg_uvlight

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_hzg_uvlight
      ! Identifiers for dependencies [model inputs]
      type (type_horizontal_dependency_id)          :: id_swr0 ! Surface shortwave radiation
      type (type_dependency_id)                     :: id_dz   ! Cell thickness
      type (type_dependency_id)                     :: id_ext  ! Attentuation coefficient for PAR

      ! Identifiers for diagnostic variables [model outputs]
      type (type_diagnostic_variable_id)            :: id_uv  ! UV radiation

      ! Parameters
      real(rk) :: uv_fraction,attenuation_length
   contains
!     Model procedures
      procedure :: initialize
      procedure :: get_light
   end type type_hzg_uvlight

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_hzg_uvlight),intent(inout),target :: self
      integer,                  intent(in)           :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%uv_fraction, 'uv_fraction', '-','UV fraction of shortwave radiation',default=0.04_rk) 
      call self%get_parameter(self%attenuation_length,'attenuation_length','m','e-folding depth of UV fraction', default=1.3_rk)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_uv,'radiation','W/m^2','UV radiation', source=source_do_column)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_swr0,standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_ext, standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_dz,  standard_variables%cell_thickness)
   end subroutine
   
   subroutine get_light(self,_ARGUMENTS_VERTICAL_)
      class (type_hzg_uvlight),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      real(rk) :: swr0,dz,z,ext,bioext
      real(rk) :: uv

      _GET_HORIZONTAL_(self%id_swr0,swr0)

      z = 0
      bioext = 0
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_dz,dz)     ! Layer height (m)
         _GET_(self%id_ext,ext)   ! attenuation from particles use in other
                                  ! submodules (m-1)

         ! Set depth to centre of layer
         z = z + dz/2
         bioext = bioext + ext*dz/2

         ! Calculate UV radiation

         uv = swr0*self%uv_fraction*exp(-z/self%attenuation_length-bioext)

         ! Move to bottom of layer
         z = z + dz/2
         bioext = bioext + ext*dz/2

         _SET_DIAGNOSTIC_(self%id_uv,uv) ! UV radiation at layer centre
      _VERTICAL_LOOP_END_

   end subroutine get_light

end module hzg_uvlight
