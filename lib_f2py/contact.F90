subroutine contact( xyz, nmp, dist, cmap)
       
   implicit none
   integer, intent(in) :: nmp
   real*8,  intent(in) :: xyz(3, nmp)
   real*8,  intent(in) :: dist
   integer, intent(out) :: cmap(nmp,nmp)
   
   integer :: imp, jmp
   real*8 :: dist2, d2
   real*8 :: xyz_i(3), d(3)

   cmap(:,:) = 0
   dist2 = dist*dist

   do imp = 1, nmp-1

      xyz_i = xyz(:,imp)

      do jmp = imp+1, nmp

         d(:) = xyz_i(:) - xyz(:,jmp) 
         d2 = dot_product(d,d)

         if (d2 <= dist2) then
            cmap(imp,jmp) = 1
         endif
      enddo
   enddo

endsubroutine contact
