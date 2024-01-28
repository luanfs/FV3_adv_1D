module dyn_core
!========================================================================
! This module contains the routine that computes the PPM flux
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/dyn_core.F90
!========================================================================
  use fv_arrays,  only: fv_grid_bounds_type
 implicit none

contains
 subroutine dy_core(time_total, bd)
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(in), optional:: time_total  ! total time (seconds) since start
 end subroutine dy_core

end module dyn_core 
