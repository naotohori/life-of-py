subroutine find_max_min_pbc( xyz, nmp, maxd, boxsize, max_xyz, min_xyz, dxyz)
       
   implicit none
   integer, intent(in) :: nmp
   real*8,  intent(in) :: xyz(3, nmp)
   real*8,  intent(in) :: maxd
   real*8,  intent(in) :: boxsize
   real*8,  intent(out) :: max_xyz(3)
   real*8,  intent(out) :: min_xyz(3)
   real*8,  intent(out) :: dxyz(3)
   
   integer :: imp, k
   real*8 :: delta, x
   real*8 :: pre(3), add(3)

   max_xyz(:) = xyz(:,1)
   min_xyz(:) = xyz(:,1)

   pre(:) = xyz(:,1)
   add(:) = 0.0

   do imp = 1, nmp

      do k = 1,3
         
         x = xyz(k,imp) + add(k)

         delta = x - pre(k)
         if (delta > maxd) then
            x      =   x    - boxsize
            add(k) = add(k) - boxsize
         else if (delta < - maxd) then
            x      =   x    + boxsize
            add(k) = add(k) + boxsize
         endif

         pre(k) = x

         if (max_xyz(k) < x) then
            max_xyz(k) = x
         endif
         if (x < min_xyz(k)) then
            min_xyz(k) = x
         endif

      enddo
   enddo

   dxyz(:) = max_xyz(:) - min_xyz(:)

endsubroutine find_max_min_pbc
