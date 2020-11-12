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

      ! Common models
      use hzg_omexdia_p
      use fabm_hzg_ecosmo

      ! KST specific models
      use fabm_hzg_test
      use fabm_hzg_icealgea
      use hzg_icelight
      use hzg_uvlight
      use fabm_hzg_pops

      ! KSE specific models
      !use hzg_omexdia_p_mpb
      !use hzg_omexdia_cnp
      !use hzg_omexdia_mpb
      use hzg_mpb
      !use hzg_mpb_cnp
      use hzg_jelly
      use hzg_ctenophore_jt
      use hzg_n2pzdq
      use hzg_medmac
      use hzg_maecs
      use hzg_benthic_pool
      use hzg_Ndepoden
      use fabm_hzg_dependencies
      use hzg_agg
      use hzg_lina
      !use hzg_kristineb

      ! Add more models modules here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('ecosmo');       allocate(type_hzg_ecosmo::model)
         case ('omexdia_p');    allocate(type_hzg_omexdia_p::model)
         case ('test');         allocate(type_hzg_test::model)
         case ('icealgea');     allocate(type_hzg_icealgea::model)
         case ('icelight');     allocate(type_hzg_icelight::model)
         case ('uvlight');      allocate(type_hzg_uvlight::model)
         case ('pops');         allocate(type_hzg_pops::model)
         !case ('omexdia_p_mpb'); allocate(type_hzg_omexdia_p_mpb::model)
         !case ('omexdia_cnp'); allocate(type_hzg_omexdia_cnp::model)
         !case ('omexdia_mpb'); allocate(type_hzg_omexdia_mpb::model)
         case ('mpb'); allocate(type_hzg_mpb::model)
         !case ('mpb_cnp'); allocate(type_hzg_mpb_cnp::model)
         case ('jelly'); allocate(type_hzg_jelly::model)
         case ('n2pzdq'); allocate(type_hzg_n2pzdq::model)
         case ('maecs'); allocate(type_hzg_maecs::model)
         case ('medmac'); allocate(type_hzg_medmac::model)
         !case ('kristineb'); allocate(type_hzg_kristineb::model)
         case ('Ndepoden'); allocate(type_hzg_Ndepoden::model)
         case ('benthic_pool'); allocate(type_hzg_benthic_pool::model)
         case ('dependencies'); allocate(type_hzg_dependencies::model)
         case ('agg'); allocate(type_hzg_agg::model)
         case ('ctenophore_jt'); allocate(type_hzg_ctenophore_jt::model)
         case ('lina'); allocate(type_hzg_lina::model)
         ! Add your model type here, and make sure
         ! they appear in the use statement above
      end select

   end subroutine

end module
