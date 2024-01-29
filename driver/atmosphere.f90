module atmosphere
!========================================================================
!========================================================================
use fv_control, only: init_model
use fv_arrays , only: fv_atmos_type, datadir, R_GRID
use test_cases, only: init_case, calc_winds
use dyn_core  , only: dy_core

implicit none

contains
!--------------------------------------------------------------
! main atmosphere routine
!--------------------------------------------------------------
subroutine atmosphere_main(atm)
   type(fv_atmos_type), intent(inout) :: atm
   integer :: tstep
   logical :: first_step=.false.

   ! init everything
   call atmosphere_init(atm)

   ! loop over time
   do tstep = 1, atm%total_tsteps
      ! compute one time step
      call atmosphere_timestep(atm)

      ! compute diagnostics
      call atmosphere_diag(atm, first_step)

      ! output data
      call atmosphere_output(atm, tstep, atm%total_tsteps)
   enddo

   ! compute error
   call atmosphere_end(atm)
end subroutine atmosphere_main

!--------------------------------------------------------------
! init the model
!--------------------------------------------------------------
subroutine atmosphere_init(atm)
   type(fv_atmos_type), intent(inout) :: atm
   logical :: first_step=.true.

   ! initialize all variables
   call init_model(atm)

   ! time vars
   atm%time          = 0.d0
   atm%time_centered = atm%time + atm%dto2

   ! compute ICs
   call init_case(atm)

   ! cfl
   atm%cfl = maxval(abs(atm%uc0))*atm%dt/atm%gridstruct%dx

   ! update qa and uc
   atm%qa = atm%qa0
   atm%uc_old = atm%uc0

   ! compute initial diagnostics
   call atmosphere_diag(atm, first_step)

   ! plot IC
   call atmosphere_output(atm, 0, atm%total_tsteps)

   first_step=.false.
end subroutine atmosphere_init

!--------------------------------------------------------------
! compute on timestep of atmosphere dynamics
!--------------------------------------------------------------
subroutine atmosphere_timestep(atm)
   type(fv_atmos_type), intent(inout) :: atm
   logical :: first_step=.false.


   ! solves dynamics
   call dy_core(atm%qa, atm%uc, atm%uc_old, atm%bd, atm%gridstruct, atm%time, atm%time_centered,&
                   atm%dt, atm%dto2, atm%test_case, atm%hord, atm%lim_fac)

   ! update times
   atm%time          = atm%time + atm%dt
   atm%time_centered = atm%time + atm%dto2

   ! winds for next timestep
   call calc_winds(atm%uc_old, atm%bd, atm%gridstruct, atm%time         , atm%test_case)
   call calc_winds(atm%uc    , atm%bd, atm%gridstruct, atm%time_centered, atm%test_case)

end subroutine atmosphere_timestep


!--------------------------------------------------------------
! atmosphere output
!--------------------------------------------------------------
subroutine atmosphere_output(atm, step, total_tsteps)
   type(fv_atmos_type), intent(inout) :: atm
   integer, intent(in) :: step, total_tsteps
   integer :: is, ie
   integer :: i, iunit 
   character (len=60):: nplot
   character (len=60):: filename
   is = atm%bd%is
   ie = atm%bd%ie

   
   if(step==0 .or. step==total_tsteps .or. mod(step,10)==0 )then
      write(nplot, '(i8)') atm%nplot
      filename = trim(datadir)//trim(atm%simulation_name)//"t"//trim(adjustl(nplot))//".txt"
      print*, 'saving ', filename
      ! write the data in a text file
      iunit = 19
      open(iunit, file=filename, status='replace')
      write(iunit,*) atm%time
      write(iunit,*) atm%mass_qa_var
      write(iunit,*) atm%cfl
      do i = is, ie
         write(iunit,*) atm%qa(i)
      end do
      close(iunit)
      atm%nplot = atm%nplot + 1
   endif
end subroutine atmosphere_output

!--------------------------------------------------------------
! Compute diagnostics
!--------------------------------------------------------------
subroutine atmosphere_diag(atm, first_step)
   type(fv_atmos_type), intent(inout) :: atm
   logical, intent(in) :: first_step
   integer :: is, ie
   is = atm%bd%is
   ie = atm%bd%ie

   if(first_step) then
      atm%mass_qa0 = sum(abs(atm%qa0(is:ie)))*atm%gridstruct%dx
   else
      atm%mass_qa     = sum(abs(atm%qa(is:ie)))*atm%gridstruct%dx
      atm%mass_qa_var = (atm%mass_qa0-atm%mass_qa)/atm%mass_qa0
   endif
end subroutine atmosphere_diag


!--------------------------------------------------------------
! end the model
!--------------------------------------------------------------
subroutine atmosphere_end(atm)
   type(fv_atmos_type), intent(inout) :: atm
   integer :: is, ie
   is = atm%bd%is
   ie = atm%bd%ie

   ! compute errors of qa
   atm%error_qa(is:ie) = abs(atm%qa(is:ie)-atm%qa0(is:ie))
   atm%linf_error_qa   = maxval(atm%error_qa(is:ie))/maxval(abs(atm%qa0(is:ie)))
   atm%l1_error_qa     = sum   (atm%error_qa(is:ie))/sum   (abs(atm%qa0(is:ie)))
   atm%l2_error_qa     = norm2 (atm%error_qa(is:ie))/norm2 (abs(atm%qa0(is:ie)))

   print*, atm%linf_error_qa, atm%l1_error_qa, atm%l2_error_qa
end subroutine atmosphere_end


end module atmosphere
