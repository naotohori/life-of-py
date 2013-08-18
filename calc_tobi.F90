subroutine calc_tobi( num_type,   &
                      param1,     &
                      param2,     &
                      cut1,       &
                      cut2,       &
                      num_chain,  &
                      max_atom,   &
                      num_atom,   &
                      atom2type,  &
                      xyz,        &
                      ene )
                      
   implicit none

   integer, intent(in)  :: num_type
   real*8,  intent(in)  :: param1(num_type, num_type)
   real*8,  intent(in)  :: param2(num_type, num_type)
   real*8,  intent(in)  :: cut1
   real*8,  intent(in)  :: cut2
   integer, intent(in)  :: num_chain
   integer, intent(in)  :: max_atom
   integer, intent(in)  :: num_atom(num_chain)
   integer, intent(in)  :: atom2type(max_atom, num_chain)
   real*8,  intent(in)  :: xyz(3, max_atom, num_chain)
   real*8,  intent(out) :: ene(num_chain,num_chain)

   integer :: i,j, m, n
   real*8  :: e
   real*8  :: cut1_sq, cut2_sq
   real*8  :: d(3), dist2

   cut1_sq = cut1 ** 2
   cut2_sq = cut2 ** 2

   ene(:,:) = 0.0e0

   do m = 1, num_chain
      do n = m+1, num_chain

         e = 0.0e0

         !! pairwise between (atom i of chain m) and (atom j of chain n)
         do i = 1, num_atom(m)
            do j = 1, num_atom(n)
      
               d = xyz(:,i,m) - xyz(:,j,n)
               dist2 = dot_product(d,d)
      
               if (dist2 > cut2_sq) then
                  cycle
               else if (dist2 <= cut1_sq) then
                  e = e + param1(atom2type(i,m)+1, atom2type(j,n)+1)
               else  ! (dist2 <= cut2_sq)
                  e = e + param2(atom2type(i,m)+1, atom2type(j,n)+1)
               endif
            enddo
         enddo

         ene(m,n) = e
         ene(n,m) = e

      enddo
   enddo

   return

endsubroutine calc_tobi
