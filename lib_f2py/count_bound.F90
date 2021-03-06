subroutine count_bound( xyz, nmp, id0_P, nP, id0_M, nM, r_bound, nbound)
       
   implicit none
   integer, intent(in) :: nmp
   integer, intent(in) :: nP
   integer, intent(in) :: nM
   real*8,  intent(in) :: xyz(3, nmp)
   integer, intent(in) :: id0_P(nP)
   integer, intent(in) :: id0_M(nM)
   real*8,  intent(in) :: r_bound
   integer, intent(out) :: nbound(nP)
   
   integer :: iP, iM
   real*8  :: d2
   real*8  :: xyz_P(3), xyz_M(3), d(3)
   real*8  :: r2_bound

   r2_bound = r_bound ** 2
   nbound(:) = 0

   do iP = 1, nP
      xyz_P = xyz(:, id0_P(iP)+1)

      do iM = 1, nM
         xyz_M = xyz(:, id0_M(iM)+1)

         d = xyz_M - xyz_P
         d2 = dot_product(d,d)
         if (d2 < r2_bound) then
            nbound(iP) = nbound(iP) + 1
         endif
      enddo
   enddo

endsubroutine count_bound
