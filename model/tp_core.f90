module tp_core
!========================================================================
! This module contains the routine that computes the PPM flux
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/tp_core.F90
!========================================================================
use fv_arrays, only: fv_grid_type, fv_grid_bounds_type, R_GRID
 implicit none

 private
 public fv_tp_1d

 real, parameter:: ppm_fac = 1.5   ! nonlinear scheme limiter: between 1 and 2
 real, parameter:: r3 = 1./3.
 real, parameter:: near_zero = 1.E-25
 real, parameter:: ppm_limiter = 2.0
 real, parameter:: r12 = 1./12.

! Suresh & Huynh scheme 2.2 (purtabation form)
! The wave-form is more diffusive than scheme 2.1
 real, parameter:: b1 =   0.0375
 real, parameter:: b2 =  -7./30.
 real, parameter:: b3 =  -23./120.
 real, parameter:: b4 =  13./30.
 real, parameter:: b5 = -11./240.
! scheme 2.1: perturbation form
 !real, parameter:: b1 =   1./30.
 !real, parameter:: b2 = -13./60.
 !real, parameter:: b3 = -13./60.
 !real, parameter:: b4 =  0.45
 !real, parameter:: b5 = -0.05

 real, parameter:: t11 = 27./28., t12 = -13./28., t13=3./7.
 real, parameter:: s11 = 11./14., s14 = 4./7.,    s15=3./14.
!----------------------------------------------------
! volume-conserving cubic with 2nd drv=0 at end point:
!----------------------------------------------------
! Non-monotonic
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: p1 =  7./12.     ! 0.58333333
  real, parameter:: p2 = -1./12.
!   q(i+0.5) = p1*(q(i-1)+q(i)) + p2*(q(i-2)+q(i+1))
! integer:: is, ie, js, je, isd, ied, jsd, jed

!
!-----------------------------------------------------------------------
contains
 subroutine fv_tp_1d(q, crx, npx, hord, fx, xfx, gridstruct, bd, lim_fac)
   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(in):: npx
   integer, intent(in)::hord
   real(R_GRID), intent(in)::  crx(bd%is:bd%ie+1)  !
   real(R_GRID), intent(in)::  xfx(bd%is:bd%ie+1)  !
   real(R_GRID), intent(inout):: q(bd%isd:bd%ied)  ! transported scalar
   real(R_GRID), intent(out)::fx(bd%is:bd%ie+1 )    ! Flux in x ( E )

   type(fv_grid_type), intent(IN), target :: gridstruct

   real(R_GRID), intent(in):: lim_fac
! optional Arguments:
! Local:
   real(R_GRID)   fx1(bd%is:bd%ie+1)
   integer i, j

   integer:: is, ie, js, je, isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   isd = bd%isd
   ied = bd%ied

   call  xppm(fx1, q, crx, hord, is, ie, isd, ied, bd%npx,  lim_fac)

   do i=is,ie+1
      fx(i) =  xfx(i) * fx1(i)
   enddo


 end subroutine fv_tp_1d

subroutine xppm(flux, q, c, iord, is, ie, isd, ied, npx, lim_fac)
 integer, INTENT(IN) :: is, ie, isd, ied
 integer, INTENT(IN) :: iord
 integer, INTENT(IN) :: npx
 real(R_GRID)   , INTENT(IN) :: q(isd:ied)
 real(R_GRID)   , INTENT(IN) :: c(is:ie+1) ! Courant   N (like FLUX)
 !real   , intent(IN) :: dxa(isd:ied)
 real(R_GRID)   , intent(IN) :: lim_fac
! !OUTPUT PARAMETERS:
 real(R_GRID)  , INTENT(OUT) :: flux(is:ie+1) !  Flux
! Local
 real(R_GRID), dimension(is-1:ie+1):: bl, br, b0, a4, da1
 real(R_GRID):: q1(isd:ied)
 real(R_GRID), dimension(is:ie+1):: fx0, fx1, xt1
 logical, dimension(is-1:ie+1):: ext5, ext6, smt5, smt6
 logical, dimension(is:ie+1):: hi5, hi6
 real(R_GRID)  al(is-1:ie+2)
 real(R_GRID)  dm(is-2:ie+2)
 real (R_GRID) dq(is-3:ie+2)
 integer:: i, j, ie3, is1, ie1, mord
 real(R_GRID):: x0, x1, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2

 !if ( .not. bounded_domain .and. grid_type<3 ) then
 !   is1 = max(3,is-1);  ie3 = min(npx-2,ie+2)
 !                       ie1 = min(npx-3,ie+1)
 !else
    is1 = is-1;         ie3 = ie+2
                        ie1 = ie+1
 !end if

 mord = abs(iord)


    do i=isd, ied
       q1(i) = q(i)
    enddo

 if ( iord < 7 ) then
! ord = 2: perfectly linear ppm scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6

   do i=is1, ie3
      al(i) = p1*(q1(i-1)+q1(i)) + p2*(q1(i-2)+q1(i+1))
   enddo

   !if ( .not.bounded_domain .and. grid_type<3 ) then
   !  if ( is==1 ) then
   !    al(0) = c1*q1(-2) + c2*q1(-1) + c3*q1(0)
   !    al(1) = 0.5*(((2.*dxa(0,j)+dxa(-1,j))*q1(0)-dxa(0,j)*q1(-1))/(dxa(-1,j)+dxa(0,j)) &
   !          +      ((2.*dxa(1,j)+dxa( 2,j))*q1(1)-dxa(1,j)*q1( 2))/(dxa(1, j)+dxa(2,j)))
   !    al(2) = c3*q1(1) + c2*q1(2) +c1*q1(3)
   !  endif
   !  if ( (ie+1)==npx ) then
   !    al(npx-1) = c1*q1(npx-3) + c2*q1(npx-2) + c3*q1(npx-1)
   !    al(npx) = 0.5*(((2.*dxa(npx-1,j)+dxa(npx-2,j))*q1(npx-1)-dxa(npx-1,j)*q1(npx-2))/(dxa(npx-2,j)+dxa(npx-1,j)) &
   !            +      ((2.*dxa(npx,  j)+dxa(npx+1,j))*q1(npx  )-dxa(npx,  j)*q1(npx+1))/(dxa(npx,  j)+dxa(npx+1,j)))
   !    al(npx+1) = c3*q1(npx) + c2*q1(npx+1) + c1*q1(npx+2)
   !  endif
   !endif

   if ( iord<0 ) then
       do i=is-1, ie+2
          al(i) = max(0., al(i))
       enddo
   endif

   if ( mord==1 ) then  ! perfectly linear scheme
        do i=is-1,ie+1
           bl(i) = al(i)   - q1(i)
           br(i) = al(i+1) - q1(i)
           b0(i) = bl(i) + br(i)
           smt5(i) = abs(lim_fac*b0(i)) < abs(bl(i)-br(i))
        enddo

      do i=is,ie+1
         if ( c(i) > 0. ) then
             fx1(i) = (1.-c(i))*(br(i-1) - c(i)*b0(i-1))
             flux(i) = q1(i-1)
         else
             fx1(i) = (1.+c(i))*(bl(i) + c(i)*b0(i))
             flux(i) = q1(i)
         endif
         if (smt5(i-1).or.smt5(i)) flux(i) = flux(i) + fx1(i)
      enddo

   elseif ( mord==2 ) then  ! perfectly linear scheme

      do i=is,ie+1
         xt = c(i)
         if ( xt > 0. ) then
              qtmp = q1(i-1)
              flux(i) = qtmp + (1.-xt)*(al(i)-qtmp-xt*(al(i-1)+al(i)-(qtmp+qtmp)))
         else
              qtmp = q1(i)
              flux(i) = qtmp + (1.+xt)*(al(i)-qtmp+xt*(al(i)+al(i+1)-(qtmp+qtmp)))
         endif
!        x0 = sign(dim(xt, 0.), 1.)
!        x1 = sign(dim(0., xt), 1.)
!        flux(i,j) = x0*(q1(i-1)+(1.-xt)*(al(i)-qtmp-xt*(al(i-1)+al(i)-(qtmp+qtmp))))     &
!                  + x1*(q1(i)  +(1.+xt)*(al(i)-qtmp+xt*(al(i)+al(i+1)-(qtmp+qtmp))))
      enddo

   elseif ( mord==3 ) then

        do i=is-1,ie+1
           bl(i) = al(i)   - q1(i)
           br(i) = al(i+1) - q1(i)
           b0(i) = bl(i) + br(i)
              x0 = abs(b0(i))
              xt = abs(bl(i)-br(i))
           smt5(i) =    x0 < xt
           smt6(i) = 3.*x0 < xt
        enddo
        do i=is,ie+1
           xt1(i) = c(i)
           if ( xt1(i) > 0. ) then
               if ( smt5(i-1) .or. smt6(i) ) then
                    flux(i) = q1(i-1) + (1.-xt1(i))*(br(i-1) - xt1(i)*b0(i-1))
               else
                    flux(i) = q1(i-1)
               endif
           else
               if ( smt6(i-1) .or. smt5(i) ) then
                    flux(i) = q1(i) + (1.+xt1(i))*(bl(i) + xt1(i)*b0(i))
               else
                    flux(i) = q1(i)
               endif
           endif
        enddo

   elseif ( mord==4 ) then

        do i=is-1,ie+1
           bl(i) = al(i)   - q1(i)
           br(i) = al(i+1) - q1(i)
           b0(i) = bl(i) + br(i)
              x0 = abs(b0(i))
              xt = abs(bl(i)-br(i))
           smt5(i) =    x0 < xt
           smt6(i) = 3.*x0 < xt
        enddo
        do i=is,ie+1
           xt1(i) = c(i)
           hi5(i) = smt5(i-1) .and. smt5(i)   ! more diffusive
           hi6(i) = smt6(i-1) .or.  smt6(i)
           hi5(i) = hi5(i) .or. hi6(i)
        enddo

        do i=is,ie+1
! Low-order only if (ext6(i-1).and.ext6(i)) .AND. ext5(i1).or.ext5(i)()
          if ( xt1(i) > 0. ) then
               fx1(i) = (1.-xt1(i))*(br(i-1) - xt1(i)*b0(i-1))
               flux(i) = q1(i-1)
           else
               fx1(i) = (1.+xt1(i))*(bl(i) + xt1(i)*b0(i))
               flux(i) = q1(i)
           endif
           if ( hi5(i) ) flux(i) = flux(i) + fx1(i)
        enddo

   else

      if ( iord==5 ) then
        do i=is-1,ie+1
           bl(i) = al(i)   - q1(i)
           br(i) = al(i+1) - q1(i)
           b0(i) = bl(i) + br(i)
           smt5(i) = bl(i)*br(i) < 0.
        enddo
      elseif ( iord==-5 ) then
        do i=is-1,ie+1
           bl(i) = al(i)   - q1(i)
           br(i) = al(i+1) - q1(i)
           b0(i) = bl(i) + br(i)
           smt5(i) = bl(i)*br(i) < 0.
           da1(i) = br(i) - bl(i)
           a4(i) = -3.*b0(i)
        enddo
        do i=is-1,ie+1
           if( abs(da1(i)) < -a4(i) ) then
           if( q1(i)+0.25/a4(i)*da1(i)**2+a4(i)*r12 < 0. ) then
             if( .not. smt5(i) ) then
                br(i) = 0.
                bl(i) = 0.
                b0(i) = 0.
             elseif( da1(i) > 0. ) then
                br(i) = -2.*bl(i)
                b0(i) =    -bl(i)
             else
                bl(i) = -2.*br(i)
                b0(i) =    -br(i)
             endif
           endif
           endif
        enddo
      else
        do i=is-1,ie+1
           bl(i) = al(i)   - q1(i)
           br(i) = al(i+1) - q1(i)
           b0(i) = bl(i) + br(i)
           smt5(i) = 3.*abs(b0(i)) < abs(bl(i)-br(i))
        enddo
      endif

      do i=is,ie+1
         if ( c(i) > 0. ) then
              fx1(i) = (1.-c(i))*(br(i-1) - c(i)*b0(i-1))
              flux(i) = q1(i-1)
         else
              fx1(i) = (1.+c(i))*(bl(i) + c(i)*b0(i))
              flux(i) = q1(i)
         endif
         if (smt5(i-1).or.smt5(i)) flux(i) = flux(i) + fx1(i)
      enddo

   endif

 else

! Monotonic constraints:
! ord = 8: PPM with Lin's PPM fast monotone constraint
! ord = 10: PPM with Lin's modification of Huynh 2nd constraint
! ord = 13: positive definite constraint

    do i=is-2,ie+2
          xt = 0.25*(q1(i+1) - q1(i-1))
       dm(i) = sign(min(abs(xt), max(q1(i-1), q1(i), q1(i+1)) - q1(i),  &
                         q1(i) - min(q1(i-1), q1(i), q1(i+1))), xt)
    enddo
    do i=is1,ie1+1
       al(i) = 0.5*(q1(i-1)+q1(i)) + r3*(dm(i-1)-dm(i))
    enddo

    if ( iord==8 ) then
       do i=is1, ie1
          xt = 2.*dm(i)
          bl(i) = -sign(min(abs(xt), abs(al(i  )-q1(i))), xt)
          br(i) =  sign(min(abs(xt), abs(al(i+1)-q1(i))), xt)
       enddo
    elseif ( iord==10 ) then
       do i=is1-2, ie1+1
          dq(i) = 2.*(q1(i+1) - q1(i))
       enddo
       do i=is1, ie1
          bl(i) = al(i  ) - q1(i)
          br(i) = al(i+1) - q1(i)
          if ( abs(dm(i-1))+abs(dm(i))+abs(dm(i+1)) < near_zero ) then
                   bl(i) = 0.
                   br(i) = 0.
          elseif( abs(3.*(bl(i)+br(i))) > abs(bl(i)-br(i)) ) then
                   pmp_2 = dq(i-1)
                   lac_2 = pmp_2 - 0.75*dq(i-2)
                   br(i) = min( max(0., pmp_2, lac_2), max(br(i), min(0., pmp_2, lac_2)) )
                   pmp_1 = -dq(i)
                   lac_1 = pmp_1 + 0.75*dq(i+1)
                   bl(i) = min( max(0., pmp_1, lac_1), max(bl(i), min(0., pmp_1, lac_1)) )
          endif
       enddo
    elseif ( iord==11 ) then
! This is emulation of 2nd van Leer scheme using PPM codes
       do i=is1, ie1
          xt = ppm_fac*dm(i)
          bl(i) = -sign(min(abs(xt), abs(al(i  )-q1(i))), xt)
          br(i) =  sign(min(abs(xt), abs(al(i+1)-q1(i))), xt)
       enddo
    elseif ( iord==7 .or. iord==12 ) then  ! positive definite (Lin & Rood 1996)
       do i=is1, ie1
          bl(i) = al(i)   - q1(i)
          br(i) = al(i+1) - q1(i)
          a4(i) = -3.*(bl(i) + br(i))
           da1(i) = br(i) - bl(i)
          ext5(i) = br(i)*bl(i) > 0.
          ext6(i) = abs(da1(i)) < -a4(i)
       enddo
       do i=is1, ie1
          if( ext6(i) ) then
            if( q1(i)+0.25/a4(i)*da1(i)**2+a4(i)*r12 < 0. ) then
                if( ext5(i) ) then
                   br(i) = 0.
                   bl(i) = 0.
                elseif( da1(i) > 0. ) then
                   br(i) = -2.*bl(i)
                else
                   bl(i) = -2.*br(i)
                endif
            endif
          endif
       enddo
    else
       do i=is1, ie1
          bl(i) = al(i  ) - q1(i)
          br(i) = al(i+1) - q1(i)
       enddo
    endif
! Positive definite constraint:
    if(iord==9 .or. iord==13) call pert_ppm(ie1-is1+1, q1(is1), bl(is1), br(is1), 0)

    !if (.not. bounded_domain .and. grid_type<3) then
    !  if ( is==1 ) then
    !     bl(0) = s14*dm(-1) + s11*(q1(-1)-q1(0))

    !     xt = 0.5*(((2.*dxa(0,j)+dxa(-1,j))*q1(0)-dxa(0,j)*q1(-1))/(dxa(-1,j)+dxa(0,j)) &
    !        +      ((2.*dxa(1,j)+dxa( 2,j))*q1(1)-dxa(1,j)*q1( 2))/(dxa(1, j)+dxa(2,j)))
!   !     if ( iord==8 .or. iord==10 ) then
    !        xt = max(xt, min(q1(-1),q1(0),q1(1),q1(2)))
    !        xt = min(xt, max(q1(-1),q1(0),q1(1),q1(2)))
!        endif
    !     br(0) = xt - q1(0)
    !     bl(1) = xt - q1(1)
    !     xt = s15*q1(1) + s11*q1(2) - s14*dm(2)
    !     br(1) = xt - q1(1)
    !     bl(2) = xt - q1(2)

    !     br(2) = al(3) - q1(2)
    !     call pert_ppm(3, q1(0), bl(0), br(0), 1)
    !  endif
    !  if ( (ie+1)==npx ) then
    !     bl(npx-2) = al(npx-2) - q1(npx-2)
    !
    !    xt = s15*q1(npx-1) + s11*q1(npx-2) + s14*dm(npx-2)
    !    br(npx-2) = xt - q1(npx-2)
    !    bl(npx-1) = xt - q1(npx-1)

    !    xt = 0.5*(((2.*dxa(npx-1,j)+dxa(npx-2,j))*q1(npx-1)-dxa(npx-1,j)*q1(npx-2))/(dxa(npx-2,j)+dxa(npx-1,j)) &
    !       +      ((2.*dxa(npx,  j)+dxa(npx+1,j))*q1(npx  )-dxa(npx,  j)*q1(npx+1))/(dxa(npx,  j)+dxa(npx+1,j)))
!        if ( iord==8 .or. iord==10 ) then
    !       xt = max(xt, min(q1(npx-2),q1(npx-1),q1(npx),q1(npx+1)))
    !!       xt = min(xt, max(q1(npx-2),q1(npx-1),q1(npx),q1(npx+1)))
!        endif
    !     br(npx-1) = xt - q1(npx-1)
    !     bl(npx  ) = xt - q1(npx  )

    !     br(npx) = s11*(q1(npx+1)-q1(npx)) - s14*dm(npx+1)
    !     call pert_ppm(3, q1(npx-2), bl(npx-2), br(npx-2), 1)
    !  endif
    !endif

  endif

  if ( iord==7 ) then
      do i=is-1,ie+1
           b0(i) = bl(i) + br(i)
         smt5(i) = bl(i) * br(i) < 0.
      enddo
      do i=is,ie+1
         if ( c(i) > 0. ) then
              fx1(i) = (1.-c(i))*(br(i-1) - c(i)*b0(i-1))
              flux(i) = q1(i-1)
         else
              fx1(i) = (1.+c(i))*(bl(i) + c(i)*b0(i))
              flux(i) = q1(i)
         endif
         if ( smt5(i-1).or.smt5(i) ) flux(i) = flux(i) + fx1(i)
      enddo
  else
      do i=is,ie+1
         if( c(i)>0. ) then
             flux(i) = q1(i-1) + (1.-c(i))*(br(i-1)-c(i)*(bl(i-1)+br(i-1)))
         else
             flux(i) = q1(i  ) + (1.+c(i))*(bl(i  )+c(i)*(bl(i)+br(i)))
         endif
      enddo
  endif

 end subroutine xppm




 subroutine pert_ppm(im, a0, al, ar, iv)
 integer, intent(in):: im
 integer, intent(in):: iv
 real(R_GRID), intent(in)   :: a0(im)
 real(R_GRID), intent(inout):: al(im), ar(im)
! Local:
 real(R_GRID) a4, da1, da2, a6da, fmin
 integer i

!-----------------------------------
! Optimized PPM in perturbation form:
!-----------------------------------

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
     if ( a0(i) <= 0. ) then
          al(i) = 0.
          ar(i) = 0.
     else
        a4 = -3.*(ar(i) + al(i))
       da1 =      ar(i) - al(i)
      if( abs(da1) < -a4 ) then
         fmin = a0(i) + 0.25/a4*da1**2 + a4*r12
         if( fmin < 0. ) then
             if( ar(i)>0. .and. al(i)>0. ) then
                 ar(i) = 0.
                 al(i) = 0.
             elseif( da1 > 0. ) then
                 ar(i) = -2.*al(i)
             else
                 al(i) = -2.*ar(i)
             endif
         endif
      endif
     endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
       if ( al(i)*ar(i) < 0. ) then
            da1 = al(i) - ar(i)
            da2 = da1**2
            a6da = 3.*(al(i)+ar(i))*da1
! abs(a6da) > da2 --> 3.*abs(al+ar) > abs(al-ar)
            if( a6da < -da2 ) then
                ar(i) = -2.*al(i)
            elseif( a6da > da2 ) then
                al(i) = -2.*ar(i)
            endif
       else
! effect of dm=0 included here
            al(i) = 0.
            ar(i) = 0.
       endif
  enddo
 endif

 end subroutine pert_ppm
end module tp_core 
