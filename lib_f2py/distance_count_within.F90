subroutine distance_count_within( xyz, nmp, id0_P, nP, id0_M, nM, r, num)
       
   implicit none
   integer, intent(in) :: nmp
   integer, intent(in) :: nP
   integer, intent(in) :: nM
   real*8,  intent(in) :: xyz(3, nmp)
   integer, intent(in) :: id0_P(nP)
   integer, intent(in) :: id0_M(nM)
   real*8,  intent(in) :: r
   integer, intent(out) :: num(nP)
   
   integer :: iP, iM
   real*8  :: d2, r2
   real*8  :: xyz_P(3), xyz_M(3), d(3)

   r2 = r * r
   num(:) = 0

   do iP = 1, nP
      xyz_P = xyz(:, id0_P(iP)+1)

      do iM = 1, nM
         xyz_M = xyz(:, id0_M(iM)+1)

         d = xyz_M - xyz_P
         d2 = dot_product(d,d)

         if (d2 < r2) then
            num(iP) = num(iP) + 1
         endif
      enddo
   enddo

endsubroutine distance_count_within

! This module 'py_distance_count_within' is auto-generated with f2py (version:2).
! Functions:
!   num =
!  distance_count_within(xyz,id0_p,id0_m,r,nmp=shape(xyz,1),np=len(id0_p),nm=len(id0_m))

