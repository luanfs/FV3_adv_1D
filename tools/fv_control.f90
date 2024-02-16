module fv_control
!========================================================================
! Module for data allocation
!
! Reference
! https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/tools/fv_control.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, fv_atmos_type, &
                      point_structure, R_GRID, erad, pi, pio2
implicit none

contains

!-------------------------------------------------
! define bounds
!-------------------------------------------------
subroutine init_bounds(bd, npx)
   type(fv_grid_bounds_type), intent(INOUT) :: bd
   integer                  , intent(IN)    :: npx

   bd%is  = 1
   bd%ie  = npx
   bd%isd = bd%is - bd%ng
   bd%ied = bd%ie + bd%ng
   bd%npx = npx
end subroutine init_bounds

!-------------------------------------------------
! allocate grid - bd must be filled
!-------------------------------------------------
subroutine init_grid(gridstruct, bd)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN)    :: bd
   type(point_structure), pointer, dimension(:) :: agrid
   type(point_structure), pointer, dimension(:) :: cgrid
   real(R_GRID), pointer :: dx
   real(R_GRID) :: L
   integer :: is, ie, isd, ied
   integer :: i

   is  = bd%is
   isd = bd%isd
   ie  = bd%ie
   ied = bd%ied
   gridstruct%npx = bd%npx

   allocate(gridstruct%agrid(isd:ied))
   allocate(gridstruct%cgrid(isd:ied+1))

   agrid => gridstruct%agrid
   cgrid => gridstruct%cgrid
   dx    => gridstruct%dx

   L = pio2*erad
   dx = L/bd%npx
   
   ! compute c grid local coordinates
   do i = isd, ied+1
      cgrid(i)%x = -L*0.5d0 + (i-1d0)*dx
   enddo

   ! compute a grid local coordinates
   agrid(isd:ied)%x = (cgrid(isd+1:ied+1)%x + cgrid(isd:ied)%x)* 0.5d0
end subroutine init_grid

!-------------------------------------------------
! allocate atmos
!-------------------------------------------------
subroutine init_atmos(atm)
   type(fv_atmos_type), intent(INOUT) :: atm
   integer :: is, ie, isd, ied
   integer :: i
   character (len=60):: n, tc, hord, dp
   is  = atm%bd%is
   isd = atm%bd%isd
   ie  = atm%bd%ie
   ied = atm%bd%ied
   atm%npx = atm%bd%npx

   write(n   ,'(i8)') atm%npx
   write(tc  ,'(i8)') atm%test_case
   write(hord,'(i8)') atm%hord
   write(dp  ,'(i8)') atm%dp
   atm%simulation_name = "tc"//trim(adjustl(tc))//"_N"//trim(adjustl(n))//"_hord"//&
   trim(adjustl(hord))//"_dp"//trim(adjustl(dp))//"_"

   allocate(atm%qa(isd:ied))
   allocate(atm%qa0(isd:ied))

   allocate(atm%ua (isd:ied))
   allocate(atm%uc (isd:ied+1))
   allocate(atm%uc0(isd:ied+1))

   allocate(atm%error_qa(is:ie))
end subroutine init_atmos


!-------------------------------------------------
! init everything
!-------------------------------------------------
subroutine init_model(atm)
   type(fv_atmos_type), intent(inout) :: atm

   call init_bounds(atm%bd, atm%npx)
   call init_grid  (atm%gridstruct, atm%bd)
   call init_atmos (atm)

end subroutine init_model


end module fv_control 
