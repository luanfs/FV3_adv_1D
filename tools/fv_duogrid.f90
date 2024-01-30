module fv_duogrid
!========================================================================
! Module to handle with ghost cell interpolations
!
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type,  &
                      R_GRID
implicit none

contains

!--------------------------------------------------------------------------
! extend agrid field
subroutine ext_scalar_agrid(qa, bd)
   type(fv_grid_bounds_type), intent(IN) :: bd
   real(R_GRID), intent(INOUT) :: qa(bd%isd:bd%ied)
   integer :: is, ie, isd, ied, ng

   is = bd%is
   ie = bd%ie
   isd = bd%isd
   ied = bd%ied
   ng = bd%ng

   ! periodic BC
   qa(isd:is-1) = qa(ie-ng+1:ie)
   qa(ie+1:ied) = qa(is:is+ng-1)
end subroutine ext_scalar_agrid
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! extend cgrid field
subroutine ext_scalar_cgrid(uc, bd)
   type(fv_grid_bounds_type), intent(IN) :: bd
   real(R_GRID), intent(INOUT) :: uc(bd%isd:bd%ied+1)
   integer :: is, ie, isd, ied, ng

   is = bd%is
   ie = bd%ie
   isd = bd%isd
   ied = bd%ied
   ng = bd%ng

   ! periodic BC
   uc(isd:is-1) = uc(ie-ng+1:ie)
   uc(ie+2:ied) = uc(is+1:is+1+ng-1)
end subroutine ext_scalar_cgrid
!--------------------------------------------------------------------------


end module fv_duogrid
