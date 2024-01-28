module fv_arrays
!========================================================================
! This module contains the routine that computes the PPM flux
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/fv_arrays.F90
!========================================================================

  type fv_grid_bounds_type

     integer :: is,  ie,  js,  je
     integer :: isd, ied, jsd, jed

     integer :: ng = 3 !default

  end type fv_grid_bounds_type
end module fv_arrays
