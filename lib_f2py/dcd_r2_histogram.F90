subroutine dcd_r2_histogram( xyz, nmp, bin_edges, nbin, hist)
       
   implicit none
   integer, intent(in) :: nmp
   integer, intent(in) :: nbin
   real*8,  intent(in) :: xyz(3, nmp)
   real*8,  intent(in) :: bin_edges(nbin+1)
   integer, intent(out) :: hist(nbin)

   integer :: i,j, k, ibin
   integer :: n, n_pre
   real*8 :: xyz_i(3), xyz_j(3), d(3)
   real*8 :: d2(nmp*(nmp-1)/2)

   k = 0 
   do i = 1, nmp
      xyz_i = xyz(:,i)

      do j = i+1, nmp
         k = k + 1
         xyz_j = xyz(:,j)

         d = xyz_j - xyz_i
         d2(k) = dot_product(d,d)
      enddo
   enddo

   n_pre = count( d2 <= bin_edges(1) )
   do ibin = 1, nbin
      n = count( d2 < bin_edges(ibin+1) )
      hist(ibin) = n - n_pre
      n_pre = n
   enddo
      
endsubroutine dcd_r2_histogram
