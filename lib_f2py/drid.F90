subroutine drid( nx, nc, na, x, centroids, atoms, munuxi)
       
   implicit none
   integer, intent(in)  :: nx, nc, na
   real*8,  intent(in)  :: x(3,nx)
   integer, intent(in)  :: centroids(nc)
   integer, intent(in)  :: atoms(na)
   real*8,  intent(out) :: munuxi(3*nc)

   integer :: ia, ic, a, c
   integer :: nsum(nc)
   real*8  :: xc(3), d(3)
   real*8  :: d_inv(na,nc)
   real*8  :: dd(na,nc)
   real*8  :: dd_23(na,nc)

   munuxi(:) = 0.0
   nsum(:) = 0

   do ic = 1, nc
      c = centroids(ic)
      xc(:) = x(:,c)

      do ia = 1, na
         a = atoms(ia)
            
         if (a == c) then
            d_inv(ia,ic) = 0.0
         else
            d(:) = x(:,a) - xc(:)
            d_inv(ia,ic) = 1.0 / sqrt(dot_product(d,d))
            nsum(ic) = nsum(ic) + 1
         endif
      enddo
   enddo

   !do ic = 1, nc
   !   do ia = 1, na
   !      munuxi(ic) = munuxi(ic) + d_inv(ia,ic)
   !   enddo
   !   munuxi(ic) = munuxi(ic) / float(nsum(ic))
   !enddo
   munuxi(1:nc) = sum(d_inv,DIM=1) / nsum(:)

   do ic = 1, nc
      dd(:,ic) = d_inv(:,ic) - munuxi(ic)
   enddo

   dd_23(:,:) = dd(:,:) * dd(:,:)
   munuxi(nc+1:2*nc) = sum(dd_23, DIM=1) / nsum(:)
   munuxi(nc+1:2*nc) = sqrt(munuxi(nc+1:2*nc))

   dd_23(:,:) = dd_23(:,:) * dd(:,:)
   munuxi(2*nc+1:3*nc) = sum(dd_23, DIM=1) / nsum(:)
   munuxi(2*nc+1:3*nc) = munuxi(2*nc+1:3*nc) ** (1.0/3.0)

endsubroutine drid
