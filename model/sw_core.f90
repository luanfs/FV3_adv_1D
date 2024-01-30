module sw_core
!========================================================================
! This module contains the routine that computes the time averaged winds
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/sw_core.F90
!========================================================================
use fv_arrays, only: R_GRID, fv_grid_bounds_type, fv_grid_type
implicit none
contains

subroutine time_averaged_cfl(gridstruct, bd, crx_adv, uc_old, uc, dp, dt)
   !--------------------------------------------------
   ! Compute the time average CFL needed
   ! for the departure point scheme
   !--------------------------------------------------
   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_grid_type), intent(IN), target :: gridstruct 

   real(R_GRID), intent(INOUT), dimension(bd%is:bd%ie+1  ) :: crx_adv
   real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1) :: uc_old
   real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1) :: uc
   real(R_GRID), intent(IN) :: dt

   real(R_GRID), dimension(bd%is-1:bd%ie+2) :: crx_time_centered
   integer, intent(IN):: dp ! departute point method !1 - Euler; 2-RK2

   ! aux
   real(R_GRID) :: a, a1, a2, c1, c2
   integer :: i
   integer :: is,  ie

   real(R_GRID), pointer :: dx

   is  = bd%is
   ie  = bd%ie

   dx => gridstruct%dx
   select case (dp)
     case (1)
        crx_adv(is:ie+1) = uc(is:ie+1)*dt/dx
        !call  compute_cfl(gridstruct, bd, crx_adv, cry_adv, uc, vc, mt, dt, 0, flagstruct%contravariant_wind_init)

     case (2)
        crx_adv(is:ie+1)             = uc_old(is:ie+1)*dt/dx
        crx_time_centered(is-1:ie+2) = uc(is-1:ie+2)*dt/dx
        ! RK2
        ! cfl for dp in x direction
        do i = is, ie+1
           ! Linear interpolation weight
           a = crx_adv(i)*0.5d0
           ! Upwind linear interpolation
           if (a>0.d0) then
              c1 = crx_time_centered(i-1)
              c2 = crx_time_centered(i)
              a1 = a
              a2 = 1.d0-a
           else
              c1 = crx_time_centered(i)
              c2 = crx_time_centered(i+1)
              a1 = 1.d0+a
              a2 = -a
           end if
           crx_adv(i) = a1*c1 + a2*c2
        end do
     case default
        print*, 'ERROR in time_averaged_cfl: invalid departure point,  ', dp
        stop
   end select
end subroutine time_averaged_cfl

end module sw_core 
