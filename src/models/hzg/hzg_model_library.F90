module hzg_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: hzg_model_factory

contains

   subroutine create(self,name,model)

      use fabm_hzg_ecosmo
      use fabm_hzg_omexdia_p
      use fabm_hzg_test
      use fabm_hzg_icealgea
      use hzg_icelight
      ! Add your own model module here...

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('ecosmo');       allocate(type_hzg_ecosmo::model)
         case ('omexdia_p');    allocate(type_hzg_omexdia_p::model)
         case ('test');         allocate(type_hzg_test::model)
         case ('icealgea');     allocate(type_hzg_icealgea::model)
         case ('icelight');     allocate(type_hzg_icelight::model)
         ! Add your new model here
      end select

   end subroutine

end module
