program main
!========================================================================
!========================================================================
use atmosphere,  only: atmosphere_main
use fv_arrays,   only: fv_atmos_type

implicit none

type(fv_atmos_type) :: Atm

! Grid size
atm%npx = 48

! Test case
atm%test_case = 2

! Hord (transport scheme)
atm%hord = 8

! Time vars
atm%Tf   = 5.d0
atm%dt   = 0.04d0
atm%dto2 = atm%dt*0.5d0
atm%total_tsteps  = int(atm%Tf/atm%dt)

! Readjust time step
atm%dt  = atm%Tf/atm%total_tsteps
!print*, atm%dt, atm%total_tsteps

call atmosphere_main(Atm)


end program main
