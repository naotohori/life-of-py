program dcd_traj_shape

   implicit none
   integer, parameter :: PREC=8
   real(PREC), parameter :: PI = 3.1415926535897932384626433832795028841971 

   integer, parameter :: FDCD = 10
   integer, parameter :: FOUT = 20

   ! Threshold of distance of neighboring beads when unwrapping
   real(PREC), parameter :: MAXD = 50.0  

   integer :: nmp_dcd   ! Read from DCD
   integer :: nmp_solute
   real(PREC) :: box_size

   integer :: imp, iframe, nframe
   integer :: istatus
   integer :: idummy
   integer :: iarg
   integer :: iargc

   real(4), allocatable :: xyz_dcd(:,:)
   real(PREC), allocatable :: xyz_solute(:,:)

   real(PREC) :: Rg, D, S
   character(256) :: cfile_dcd, cfile_out, cinp

   iarg = iargc()
   if (iarg /= 4) then
      write(*,*) 'Usage: PROGRAM [DCD file] [box size] [imp solute] [output]'
      stop
   endif

   call getarg(1, cfile_dcd)
   call getarg(2, cinp)
   read(cinp, *) box_size

   call getarg(3, cinp)
   read(cinp, *) nmp_solute

   call getarg(4, cfile_out)

   open(FOUT, file=cfile_out, status='UNKNOWN', action='WRITE', iostat=istatus)
   if (istatus > 0) then
      write(*,*) 'Error in opening output file'
      stop
   endif

   open(FDCD, file=cfile_dcd, status='OLD', action="READ", iostat=istatus, &
              form='UNFORMATTED', access='STREAM') 
   if (istatus > 0) then
      write(*,*) 'Error in opening DCD file'
      stop
   endif

   ! Test with 12_000.dcd
   !nframe = 10000
   !nmp_dcd = 1010
   call dcd_count_frame(FDCD, nframe, nmp_dcd, istatus)
   if (istatus > 0) then
      write(*,   *) 'Warning: the DCD file has an abnormal end'
      write(FOUT,*) '## Warning: the DCD file has an abnormal end'
   endif 

   !write(FOUT,*) '## Number of frames: ',nframe
   !write(FOUT,*) '## The frame number shown below starts from 0.'
   !write(FOUT,*) '## Number of MP in DCD: ',nmp_dcd

   call dcd_skip_header(FDCD)

   allocate(xyz_dcd(3,nmp_dcd))
   allocate(xyz_solute(3,nmp_solute))
 
   !!! iframe counting starts from 0
   do iframe = 0, nframe-1
 
      read (FDCD) idummy
      read (FDCD) (xyz_dcd(1,imp),imp=1,nmp_dcd)
      read (FDCD) idummy
      read (FDCD) idummy
      read (FDCD) (xyz_dcd(2,imp),imp=1,nmp_dcd)
      read (FDCD) idummy
      read (FDCD) idummy
      read (FDCD) (xyz_dcd(3,imp),imp=1,nmp_dcd)
      read (FDCD) idummy

      call unwrap_PBC(box_size, nmp_dcd, xyz_dcd, nmp_solute, xyz_solute)

      call calc_shape(nmp_solute, xyz_solute, Rg, D, S)

      write(FOUT, '(f6.2, 1x, f6.3, 1x, f6.3)') Rg, D, S
      flush(FOUT)

   enddo

   deallocate(xyz_dcd)
   deallocate(xyz_solute)

   stop

contains
   subroutine unwrap_PBC(BOXSIZE, nmp_dcd, d_in, nmp_solute, d_out)

      real(PREC), intent(in) :: BOXSIZE
      integer,    intent(in) :: nmp_dcd
      real(4),    intent(in) :: d_in(3,nmp_dcd)
      integer,    intent(in) :: nmp_solute
      real(PREC), intent(out) :: d_out(3,nmp_solute)

      real(PREC) :: add(3)
      integer :: imp, i
      real(PREC) :: delta, x

      add(1:3) = 0.0e0_PREC
      d_out(1,1) = real(d_in(1,1), kind=PREC)
      d_out(2,1) = real(d_in(2,1), kind=PREC)
      d_out(3,1) = real(d_in(3,1), kind=PREC)

      do imp = 2, nmp_solute
         do i = 1, 3
            x = real(d_in(i,imp), kind=PREC) + add(i)
            delta = x - d_out(i,imp-1)
            if (delta > MAXD) then
               x      = x      - BOXSIZE
               add(i) = add(i) - BOXSIZE
            else if (delta < -MAXD) then
               x      = x      + BOXSIZE
               add(i) = add(i) + BOXSIZE
            endif
            
            d_out(i,imp) = x
         enddo
      enddo

   endsubroutine unwrap_PBC

   subroutine eigen_value(A, w)
      ! Wikipedia "Eigenvalue algorithm"
      ! Given a real symmetric 3x3 matrix A, compute the eigenvalues
      ! Note that acos and cos operate on angles in radians
      real(PREC), intent(in) :: A(3,3)
      real(PREC), intent(out) :: w(3)
      
      real(PREC) :: p1, p2
      real(PREC) :: p, q, detB, r, phi
      real(PREC) :: B(3,3), qI(3,3)

      p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2

      if (p1 == 0.0) then
         ! A is diagonal.
         w(1) = A(1,1)
         w(2) = A(2,2)
         w(3) = A(3,3)
      else
         q = (A(1,1) + A(2,2) + A(3,3)) / 3.0   ! trace(A) is the sum of all diagonal values
         p2 = (A(1,1) - q)**2 + (A(2,2) - q)**2 + (A(3,3) - q)**2 + 2 * p1
         p = sqrt(p2 / 6.0)

         qI(:,:) = 0.0e0_PREC
         qI(1,1) = q
         qI(2,2) = q
         qI(3,3) = q

         !B = (1 / p) * (A - q * I)    ! I is the identity matrix
         B = (1 / p) * (A - qI)
         detB =   B(1,1)*B(2,2)*B(3,3) &
                - B(1,1)*B(2,3)*B(3,2) &
                - B(1,2)*B(2,1)*B(3,3) &
                + B(1,2)*B(2,3)*B(3,1) &
                + B(1,3)*B(2,1)*B(3,2) &
                - B(1,3)*B(2,2)*B(3,1)
         r = 0.5 * detB
      
         ! In exact arithmetic for a symmetric matrix  -1 <= r <= 1
         ! but computation error can leave it slightly outside this range.
         if (r <= -1.0) then
            phi = PI / 3.0
         else if (r >= 1.0) then
            phi = 0.0
         else
            phi = acos(r) / 3.0
         endif
      
         ! the eigenvalues satisfy w(3) <= w(2) <= w(1)
         w(1) = q + 2 * p * cos(phi)
         w(3) = q + 2 * p * cos(phi + (2*PI/3.0))
         w(2) = 3 * q - w(1) - w(3)     ! since trace(A) = eig1 + eig2 + eig3
      endif
   endsubroutine eigen_value

   subroutine calc_shape(N, xyzs, Rg, D, S)
      integer, intent(in) :: N
      real(PREC), intent(in) :: xyzs(3, N)
      real(PREC), intent(out) :: Rg, D, S

      real(PREC) :: w(3)
      real(PREC) :: trT, w_avg
      real(PREC) :: T(3,3)
      integer :: i, j, a, b

      T(1:3, 1:3) = 0.0e0_PREC

      do i = 1, N
         do j = 1, N
            do a = 1, 3
               do b = 1, 3
                  T(a,b) = T(a,b) + (xyzs(a,i) - xyzs(a,j)) * (xyzs(b,i) - xyzs(b,j))
               enddo
            enddo
         enddo
      enddo

      T(:,:) = T(:,:) / (2.0 * N * N)

      call eigen_value(T, w)

      trT = sum(w)
      w_avg = trT / 3.0

      Rg = sqrt(trT)
      D =  1.5 * ( (w(1)-w_avg)**2 + (w(2)-w_avg)**2 + (w(3)-w_avg)**2 ) / (trT**2)
      S = 27.0 * ( (w(1)-w_avg) * (w(2)-w_avg) * (w(3)-w_avg) ) / (trT**3)

   endsubroutine calc_shape

   subroutine dcd_skip_header(f)
      integer, intent(in) :: f
      integer :: i, nblock_size, idummy

      read (f) nblock_size
      do i = 1, nblock_size, 4
        read(f) idummy
      enddo
      read (f) nblock_size

      read (f) nblock_size
      do i = 1, nblock_size, 4
        read(f) idummy
      enddo
      read (f) nblock_size

      read (f) nblock_size
      do i = 1, nblock_size, 4
        read(f) idummy
      enddo
      read (f) nblock_size
   endsubroutine dcd_skip_header

   subroutine dcd_count_frame(f, n, nmp, istatus)
      integer, intent(in) :: f
      integer, intent(out) :: n, nmp, istatus
      real(4) :: xdummy
      logical :: flg_first = .True.

      call dcd_skip_header(f)

      n = 0
      istatus = 0
      L1: do 
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1
         if (flg_first) then
            nmp = idummy / 4
            flg_first = .False.
         endif
         Lx: do imp = 1, nmp
            read (fdcd, iostat=istatus) xdummy
            if (istatus /= 0) exit L1
         enddo Lx
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1

         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1
         Ly: do imp = 1, nmp
            read (fdcd, iostat=istatus) xdummy
            if (istatus /= 0) exit L1
         enddo Ly
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1

         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1
         Lz: do imp = 1, nmp
            read (fdcd, iostat=istatus) xdummy
            if (istatus /= 0) exit L1
         enddo Lz
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1

         n = n + 1
      enddo L1

      rewind(f)

   endsubroutine dcd_count_frame

end program dcd_traj_shape
