subroutine distance2_hist_nt( xyz, nmp, id0_P, nP, id0_M, nM, r2_hist, dist2_M, idx_M)
       
   implicit none
   integer, intent(in) :: nmp
   integer, intent(in) :: nP
   integer, intent(in) :: nM
   real*8,  intent(in) :: xyz(3, nmp)
   integer, intent(in) :: id0_P(nP)
   integer, intent(in) :: id0_M(nM)
   real*8,  intent(in) :: r2_hist
   real*8,  intent(out) :: dist2_M(nM,nP)
   integer, intent(out) :: idx_M(nP)
   
   integer :: iP, iM
   real*8  :: d2
   real*8  :: xyz_P(3), xyz_M(3), d(3)

   do iP = 1, nP
      xyz_P = xyz(:, id0_P(iP)+1)

      idx_M(iP) = 0
      do iM = 1, nM
         xyz_M = xyz(:, id0_M(iM)+1)

         d = xyz_M - xyz_P
         d2 = dot_product(d,d)
         if (d2 < r2_hist) then
            idx_M(iP) = idx_M(iP) + 1
            dist2_M(idx_M(iP),iP) = d2
         endif
      enddo
   enddo

endsubroutine distance2_hist_nt

! dist2_m,idx_m = distance2_hist_nt(xyz,id0_p,id0_m,r2_hist,[nmp,np,nm])
! 
! Wrapper for ``distance2_hist_nt``.
! 
! Parameters
! ----------
! xyz : input rank-2 array('d') with bounds (3,nmp)
! id0_p : input rank-1 array('i') with bounds (np)
! id0_m : input rank-1 array('i') with bounds (nm)
! r2_hist : input float
! 
! Other Parameters
! ----------------
! nmp : input int, optional
!     Default: shape(xyz,1)
! np : input int, optional
!     Default: len(id0_p)
! nm : input int, optional
!     Default: len(id0_m)
! 
! Returns
! -------
! dist2_m : rank-2 array('d') with bounds (nm,np)
! idx_m : rank-1 array('i') with bounds (np)
