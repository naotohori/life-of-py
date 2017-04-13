subroutine drid( nx, nc, x, mask, munuxi)
       
   implicit none
   integer, intent(in)  :: nx, nc
   real*8,  intent(in)  :: x(3,nx)
   integer, intent(in)  :: mask(nx, nc)
   real*8,  intent(out) :: munuxi(3*nc)

   integer :: ia, ic, ix
   integer :: nsum(nc)
   real*8  :: fnsum(nc)
   real*8  :: xc(3), d(3)
   real*8  :: d_inv(nx,nc)
   real*8  :: dd(nx,nc)
   real*8  :: dd_23(nx,nc)
   real*8  :: s
   logical :: flg

   munuxi(:) = 0.0
   nsum(:) = 0
   d_inv(:,:) = 0.0
   ic = 0

   do ix = 1, nx
      xc(:) = x(:,ix)
      flg = .False.

      do ia = 1, nx

         if (mask(ia, ix) < 0) then
            cycle
         endif
         if (.not. flg) then
            flg = .True.
            ic = ic + 1
         endif

         d(:) = x(:,ia) - xc(:)
         d_inv(ia,ic) = 1.0 / sqrt(dot_product(d,d))
         nsum(ic) = nsum(ic) + 1
      enddo
   enddo

   ! Check
   if (ic /= nc) then
      write(*,*) 'Error: ic /= nc in drid.F90'
      stop
   endif

   ! Debug
   !do ic = 1, nc
   !   write(*,*) 'nsum(',ic,')', nsum(ic)
   !enddo

   fnsum(:) = real(nsum(:), kind=8)

   !do ic = 1, nc
   !   do ia = 1, na
   !      munuxi(ic) = munuxi(ic) + d_inv(ia,ic)
   !   enddo
   !   munuxi(ic) = munuxi(ic) / float(nsum(ic))
   !enddo
   munuxi(1:nc) = sum(d_inv,DIM=1) / fnsum(:)

   dd(:,:) = 0.0
   do ic = 1, nc
      do ia = 1, nx
         if (d_inv(ia,ic) > 0.0) then
            dd(ia, ic) = d_inv(ia,ic) - munuxi(ic)
         endif
      enddo
   enddo

   dd_23(:,:) = dd(:,:) ** 2
   munuxi(nc+1:2*nc) = sqrt( sum(dd_23, DIM=1) / fnsum(:) )

   dd_23(:,:) = dd_23(:,:) * dd(:,:)

   !!!! To avoid numerical error
   !!!! (x ** (1./3.) returns NaN for x < 0.0)
   !!!! sum(dd_23,DIM=1) can be negative.
   if (minval( sum(dd_23, DIM=1) ) < 0.0) then
      do ic = 1, nc
         s = sum(dd_23(:,ic))
         if (s < 0.0) then
            munuxi(2*nc+ic) = -( (-s / fnsum(ic) ) ** (1./3.))
         else
            munuxi(2*nc+ic) = (s / fnsum(ic) ) ** (1./3.)
         endif
      enddo
   else
      munuxi(2*nc+1:3*nc) = (sum(dd_23, DIM=1) / fnsum(:) ) ** (1.0/3.0)
   endif

endsubroutine drid
