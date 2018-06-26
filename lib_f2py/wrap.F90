subroutine wrap( xyz, nmp, boxsize, boxmax, boxmin )
   implicit none
   integer, intent(in) :: nmp
   real*8,  intent(inout) :: xyz(3, nmp)
   real*8,  intent(in) :: boxsize, boxmax, boxmin

   integer :: imp, k
   real*8 :: x

   do imp = 1, nmp

      do k = 1,3

         x = xyz(k,imp)
         if (x > boxmax) then
            xyz(k,imp) = x - boxsize * (int((x - boxmax)/boxsize) + 1)
         else if (x < boxmin) then
            xyz(k,imp) = x + boxsize * (int((boxmin - x)/boxsize) + 1)
         endif

      enddo
   enddo

endsubroutine wrap
