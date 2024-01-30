module dyn_core
!========================================================================
! This module contains the routine that computes the PPM flux
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/dyn_core.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, R_GRID
use test_cases, only: calc_winds
use tp_core,    only: fv_tp_1d
use sw_core,    only: time_averaged_cfl
use fv_duogrid, only: ext_scalar_agrid, ext_scalar_cgrid

implicit none

contains
subroutine dy_core(qa, uc, uc_old, bd, gridstruct, time, time_centered, dt, dto2, test_case, hord, lim_fac, dp)
   type(fv_grid_bounds_type), intent(INOUT) :: bd
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   real(R_GRID), intent(in) :: time, dt, dto2
   real(R_GRID), intent(in) :: time_centered
   real(R_GRID), intent(inout) :: lim_fac
   real(R_GRID), intent(inout) :: qa(bd%isd:bd%ied)
   real(R_GRID), intent(inout) :: uc(bd%isd:bd%ied)
   real(R_GRID), intent(inout) :: uc_old(bd%isd:bd%ied)
   integer, intent(IN) :: test_case
   integer, intent(IN) :: hord
   integer, intent(IN) :: dp

   real(R_GRID), pointer :: dx
   real(R_GRID) :: xfx_adv(bd%is:bd%ie+1)
   real(R_GRID) :: crx_adv(bd%is:bd%ie+1)
   real(R_GRID) :: flux(bd%is:bd%ie+1)
   integer :: is, ie, isd, ied, ng
   integer :: i

   is  = bd%is
   ie  = bd%ie
   isd = bd%isd
   ied = bd%ied
   ng  = bd%ng
   dx  => gridstruct%dx

   ! winds
   call calc_winds(uc_old, bd, gridstruct, time         , test_case)
   call calc_winds(uc    , bd, gridstruct, time_centered, test_case)

   ! periodic BC
   call ext_scalar_agrid(qa, bd)
   call ext_scalar_cgrid(uc, bd)
   call ext_scalar_cgrid(uc_old, bd)

   ! compute time averaged cfl
   call time_averaged_cfl(gridstruct, bd, crx_adv, uc_old, uc, dp, dt)

   ! compute adv coeffs
   xfx_adv(is:ie+1)= crx_adv(is:ie+1)

   !call  xppm(flux, qa, crx_adv, hord, is, ie, isd, ied, bd%npx,  lim_fac)
   call fv_tp_1d(qa, crx_adv, bd%npx, hord, flux, xfx_adv, gridstruct, bd, lim_fac)

   ! update scalar field
   qa(is:ie) = qa(is:ie) - (flux(is+1:ie+1)-flux(is:ie))
end subroutine dy_core

end module dyn_core 
