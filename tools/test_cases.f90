module test_cases
!========================================================================
! Module for ICS
!
! Reference
! https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/tools/test_cases.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, fv_atmos_type, &
                      point_structure, R_GRID, pi
implicit none

contains

!-------------------------------------------------
! compute the ICS
!-------------------------------------------------
subroutine init_case(atm)
   type(fv_atmos_type), intent(INOUT) :: atm

   call init_scalar(atm%qa0, atm%bd, atm%gridstruct, atm%test_case)
   call init_winds (atm%uc0, atm%bd, atm%gridstruct, atm%test_case)
   call calc_winds (atm%uc , atm%bd, atm%gridstruct, atm%time_centered, atm%test_case)
end subroutine init_case

!-------------------------------------------------
! compute initial scalar field
!-------------------------------------------------
subroutine init_scalar(qa, bd, gridstruct, test_case)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN) :: bd
   type(point_structure), pointer, dimension(:) :: agrid
   real(R_GRID), intent(OUT) :: qa(bd%isd:bd%ied)
   real(R_GRID) :: x
   integer, intent(IN) :: test_case
   integer :: is, ie
   integer :: i

   is = bd%is
   ie = bd%ie
 
   agrid => gridstruct%agrid

   if (test_case == 1 .or. test_case==2) then
      do i = is, ie
         x = agrid(i)%x
         qa(i) = dexp(-10*(dcos(pi*x))**2)
      enddo
   else

   endif


end subroutine init_scalar

!-------------------------------------------------
! compute initial winds
!-------------------------------------------------
subroutine init_winds(uc, bd, gridstruct, test_case)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN) :: bd
   type(point_structure), pointer, dimension(:) :: cgrid
   real(R_GRID), intent(OUT) :: uc(bd%isd:bd%ied+1)
   integer, intent(IN) :: test_case

   call calc_winds(uc, bd, gridstruct, 0.d0, test_case)

end subroutine init_winds


!-------------------------------------------------
! compute winds at a given time
!-------------------------------------------------
subroutine calc_winds(uc, bd, gridstruct, time, test_case)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN) :: bd
   type(point_structure), pointer, dimension(:) :: cgrid
   real(R_GRID), intent(OUT) :: uc(bd%isd:bd%ied+1)
   real(R_GRID), intent(IN) :: time
   integer, intent(IN) :: test_case
   integer :: is, ie
   integer :: i

   is = bd%is
   ie = bd%ie
 
   cgrid => gridstruct%cgrid

   do i = is, ie+1
      call compute_wind(uc(i), cgrid(i)%x, time, test_case)
   enddo

end subroutine calc_winds

!-------------------------------------------------
! compute wind at a given time and position
!-------------------------------------------------
subroutine compute_wind(uc, x, t, test_case)
   real(R_GRID), intent(OUT) :: uc
   real(R_GRID), intent(IN)  :: x, t
   integer, intent(IN) :: test_case
   real(R_GRID) :: Tf, u0, u1

   select case (test_case)
      case(1)
         uc = 0.2d0

      case(2)
         Tf = 5
         u0 = 0.2
         u1 = 0.2
         uc = u0*dsin(pi*(x-t/Tf))**2*dcos(pi*t/Tf) + u1

      case default
         print*, 'error in compute_wind: invalid testcase, ', test_case
   end select
end subroutine compute_wind


end module test_cases
