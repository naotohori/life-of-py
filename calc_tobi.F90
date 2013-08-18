subroutine calc_tobi( num_type,   &
                      param1,     &
                      param2,     &
                      cut1,       &
                      cut2,       &
                      num_atomI,  &
                      num_atomJ,  &
                      atom2typeI, &
                      atom2typeJ, &
                      xyzI,       &
                      xyzJ,       &
                      ene )
                      
   implicit none

   integer, intent(in)  :: num_type
   real*8,  intent(in)  :: param1(num_type, num_type)
   real*8,  intent(in)  :: param2(num_type, num_type)
   real*8,  intent(in)  :: cut1
   real*8,  intent(in)  :: cut2
   integer, intent(in)  :: num_atomI
   integer, intent(in)  :: num_atomJ
   integer, intent(in)  :: atom2typeI(num_atomI)
   integer, intent(in)  :: atom2typeJ(num_atomJ)
   real*8,  intent(in)  :: xyzI(3, num_atomI)
   real*8,  intent(in)  :: xyzJ(3, num_atomJ)
   real*8,  intent(out) :: ene

   integer :: i,j
   real*8  :: cut1_sq, cut2_sq
   real*8  :: d(3), dist2

   cut1_sq = cut1 ** 2
   cut2_sq = cut2 ** 2

   ene = 0.0e0_8

   do i = 1, num_atomI
      do j = 1, num_atomJ

         d = xyzI(:,i) - xyzJ(:,j)
         dist2 = dot_product(d,d)

         if (dist2 > cut2_sq) then
            cycle
         else if (dist2 <= cut1_sq) then
            ene = ene + param1(atom2typeI(i)+1, atom2typeJ(j)+1)
         else  ! (dist2 <= cut2_sq)
            ene = ene + param2(atom2typeI(i)+1, atom2typeJ(j)+1)
         endif
      enddo
   enddo

   return

endsubroutine calc_tobi
