subroutine ddrid( nc, x1, x2, dd)
       
   implicit none
   integer, intent(in) :: nc
   real*8,  intent(in) :: x1(3*nc)
   real*8,  intent(in) :: x2(3*nc)
   real*8,  intent(out):: dd

   real*8 :: d(3*nc)

   d = x1 - x2
   dd = sqrt ( dot_product(d,d) / (3 * real(nc)) )
      
endsubroutine ddrid
