subroutine distance2_hist( xyz, nmp, id0_P, nP, id0_Mg, nMg, id0_K, nK, r2_hist, &
                          dist2_Mg, idx_Mg, dist2_K, idx_K)
       
   implicit none
   integer, intent(in) :: nmp
   integer, intent(in) :: nP
   integer, intent(in) :: nMg
   integer, intent(in) :: nK
   real*8,  intent(in) :: xyz(3, nmp)
   integer, intent(in) :: id0_P(nP)
   integer, intent(in) :: id0_Mg(nMg)
   integer, intent(in) :: id0_K(nK)
   real*8,  intent(in) :: r2_hist
   real*8,  intent(out) :: dist2_Mg(nP*nMg)
   real*8,  intent(out) :: dist2_K(nP*nK)
   integer, intent(out) :: idx_Mg
   integer, intent(out) :: idx_K
   
   integer :: iP, iK, iMg
   real*8  :: d2
   real*8  :: xyz_P(3), xyz_K(3), xyz_Mg(3), d(3)

   idx_K = 0
   idx_Mg = 0

   do iP = 1, nP
      xyz_P = xyz(:, id0_P(iP)+1)

      do iK = 1, nK
         xyz_K = xyz(:, id0_K(iK)+1)

         d = xyz_K - xyz_P
         d2 = dot_product(d,d)
         if (d2 < r2_hist) then
            idx_K = idx_K + 1
            dist2_K(idx_K) = d2
         endif
      enddo
      
      do iMg = 1, nMg
         xyz_Mg = xyz(:, id0_Mg(iMg)+1)

         d = xyz_Mg - xyz_P
         d2 = dot_product(d,d)
         if (d2 < r2_hist) then
            idx_Mg = idx_Mg + 1
            dist2_Mg(idx_Mg) = d2
         endif
      enddo
   enddo

endsubroutine distance2_hist
