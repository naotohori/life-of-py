subroutine distance2_hist( xyz, nmp, id0_P, nP, id0_M, nM, r2_hist, &
                          dist2_M, idx_M)
       
   implicit none
   integer, intent(in) :: nmp
   integer, intent(in) :: nP
   integer, intent(in) :: nM
   real*8,  intent(in) :: xyz(3, nmp)
   integer, intent(in) :: id0_P(nP)
   integer, intent(in) :: id0_M(nM)
   real*8,  intent(in) :: r2_hist
   real*8,  intent(out) :: dist2_M(nP*nM)
   integer, intent(out) :: idx_M
   
   integer :: iP, iM
   real*8  :: d2
   real*8  :: xyz_P(3), xyz_M(3), d(3)

   idx_M = 0

   do iP = 1, nP
      xyz_P = xyz(:, id0_P(iP)+1)

      do iM = 1, nM
         if (id0_M(iM) == id0_P(iP)) then
            cycle
         endif
         xyz_M = xyz(:, id0_M(iM)+1)

         d = xyz_M - xyz_P
         d2 = dot_product(d,d)
         if (d2 < r2_hist) then
            idx_M = idx_M + 1
            dist2_M(idx_M) = d2
         endif
      enddo
   enddo

endsubroutine distance2_hist
