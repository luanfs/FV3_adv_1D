module dyn_core
!========================================================================
! This module contains the routine that computes the PPM flux
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/dyn_core.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, R_GRID
use test_cases, only: calc_winds
use tp_core,    only: xppm
implicit none

contains
subroutine dy_core(qa, uc, uc_old, bd, gridstruct, time, time_centered, dt, dto2, test_case, hord, lim_fac)
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

   real(R_GRID), pointer :: dx
   real(R_GRID) :: xfx_adv(bd%is:bd%ie+1)
   real(R_GRID) :: cx(bd%is:bd%ie+1)
   real(R_GRID) :: flux(bd%is:bd%ie+1)
   integer :: is, ie, isd, ied, ng
   integer :: i

   is  = bd%is
   ie  = bd%ie
   isd = bd%isd
   ied = bd%ied
   ng  = bd%ng
   dx  => gridstruct%dx

   ! periodic BC
   qa(isd:is-1) = qa(ie-ng+1:ie)
   qa(ie+1:ied) = qa(is:is+ng-1)

   ! compute adv coeffs
   cx(is:ie+1) = uc(is:ie+1)*dt/dx
   xfx_adv(is:ie+1)= cx(is:ie+1)

   call  xppm(flux, qa, cx, hord, is, ie, isd, ied, bd%npx,  lim_fac)

   ! update scalar field
   qa(is:ie) = qa(is:ie) - (cx(is+1:ie+1)*flux(is+1:ie+1)-cx(is:ie)*flux(is:ie))
end subroutine dy_core

end module dyn_core 
