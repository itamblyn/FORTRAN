program integrate


implicit none

! integers
integer j, nxmax, nx, ny, nz, i, nl, is,k
integer stride, natoms, natoms_max, ispin, dan_unit
integer ntype, nend, rmdr, cnt, istat, nspin, nstart
integer, allocatable :: type(:)
integer ntype_max, sio, iswitch, origin_atom
parameter(ntype_max=100) ! Maximum number of types of atoms
integer cnt_find

! doubles
double precision, allocatable :: tmp(:)
double precision, allocatable :: density(:,:), xred(:,:)
double precision :: sum, vol, afac, rprim(3,3), dv, v(3), v2(3), a0
double precision :: tol, origin(3), dens1, dens2, delta2, delta
double precision :: dorigin(3)
parameter(a0=0.529177d0,tol=1.d-5)

! characters      

character(100) title1, title2, title3
character(2), allocatable :: label1(:), label(:)

! logicals

logical integrate_file, vasp_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



open(unit=10,file='CHGCAR',status='old',form='formatted')



! standard ouptut
      sio = 6

! format (# of columns for density in input file)
      stride=5
      allocate(tmp(stride))

! number of types of atoms
      write(sio,'(a)') "Number of atom types?"
      read(*,*) ntype
      allocate(type(ntype),label1(ntype))
      do i = 1, ntype
         write(*,'(a,2x,i5,2x,a)') "Label for atom #", i, "type?"
         read(*,*) label1(i)
      enddo
      write(sio,'(a)') "Number of spin states?"
      read(*,*) nspin
      write(sio,'(a)') "Atom number to be at the origin?"
      read(*,*) origin_atom

      read(10,*) title1
      read(10,*) afac
      do i = 1, 3
         read(10,*) (rprim(j,i),j=1,3)
         write(*,*) (rprim(j,i),j=1,3)
         rprim(:,i) = rprim(:,i)*afac/a0
      enddo

! cell volume
      vol = rprim(1,1)*(rprim(2,2)*rprim(3,3)-rprim(3,2)*rprim(2,3))+ &
          rprim(2,1)*(rprim(3,2)*rprim(1,3)-rprim(1,2)*rprim(3,3))+ &
          rprim(3,1)*(rprim(1,2)*rprim(2,3)-rprim(2,2)*rprim(1,3))
      write(sio,'(a,5x,f15.8)') "Volume of unit cell:", vol
      read(10,*) (type(j),j=1,ntype)
      write(*,*) (type(j),j=1,ntype)
      natoms = 0
      do i = 1, ntype
         natoms = natoms + type(i)
      enddo
      allocate(label(natoms))
      k = 0
      do i = 1, ntype
         do j = 1, type(i)
            k = k + 1
            label(k) = label1(i)
         enddo
      enddo

      write(sio,'(a,5x,i5)') "Number of atoms:", natoms
      read(10,*) title2
      allocate(xred(3,natoms))
      if( origin_atom == 0 ) then
         origin(:) = 0.0
      endif
      do i = 1, natoms
         read(10,*) (xred(j,i),j=1,3)
         if( i == origin_atom ) then
            origin(:) = xred(:,i)
            write(*,*) "Origin of the integration:", origin(:)
         endif
      enddo
      read(10,*)
      read(10,*) nx, ny, nz
      do i = 1, 3
         v(i) = rprim(i,1)*origin(1)+rprim(i,2)*origin(2)+rprim(i,3)*origin(3)
      enddo
      origin(:) = v(:)
      write(*,*) "Origin of the integration (Cartesian):", origin(:)
      
      write(sio,'(a,5x,i8)') "Number of real space points:", nx*ny*nz
      allocate(density(nspin,nx*ny*nz),STAT=istat)
      density(:,:) = 0.0
      if( istat /= 0 ) then
         write(sio,'(a)') "Failed to allocate memory for the density."
         write(sio,'(a,5x,f12.8)') "Memory demands (MB):", float(nspin*8*nx*ny*nz)/1.d6
      endif
      dv = vol/float(nx*ny*nz)
      write(sio,'(a,5x,f8.4)') "Integration weight:", dv


      nl = nx*ny*nz/stride
      write(*,*) "Number of input lines:", nl
      sum = 0.0
      cnt = 0
      is = 1
      do i = 1, nl
         read(10,*) (tmp(j),j=1,stride)
         do j = 1, stride
            cnt = cnt + 1
            density(is,cnt) = tmp(j)
            sum = sum + tmp(j)
            if( tmp(j) .lt. 0.0 ) then
               write(sio,'(a)') "Negative density?"
            endif
         enddo
      enddo
      sum = sum/float(nx*ny*nz)
      if( cnt < nx*ny*nz ) then
         read(10,*) (tmp(j-cnt),j=cnt+1,nx*ny*nz)
         do j = cnt+1, nx*ny*nz
            density(is,j) = tmp(j-cnt)
            sum = sum + tmp(j-cnt)
         enddo
      endif
      write(sio,*) "Total electron charge:", sum

! Spin-polarized calculations
      if( nspin == 2 ) then
         is = 2
         vasp_new = .false.
         if( vasp_new ) then
            do i = 1, natoms
               read(10,*) title3
               read(10,*)
            enddo
         endif
         rmdr = natoms
         nend = stride
         if( natoms < nend ) nend = natoms
         do while ( rmdr .gt. 0 )
            read(10,*) (tmp(i),i=1,nend)
            rmdr = rmdr - stride
            if( rmdr .lt. stride ) then
               nend = rmdr
            endif
         enddo
         read(10,*) nx, ny, nz
         nl = nx*ny*nz/stride
         write(*,*) "Number of input lines:", nl
         sum = 0.0
         cnt = 0
         do i = 1, nl
            read(10,*) (tmp(j),j=1,stride)
            do j = 1, stride
               cnt = cnt + 1
               density(is,cnt) = tmp(j)
               sum = sum + tmp(j)
            enddo
         enddo
         if( cnt < nx*ny*nz ) then
            read(10,*) (tmp(j-cnt),j=cnt+1,nx*ny*nz)
            do j = cnt+1, nx*ny*nz
               density(is,j) = tmp(j-cnt)
               sum = sum + tmp(j-cnt)
            enddo
         endif
         sum = sum/float(nx*ny*nz)
         write(sio,'(a,5x,f12.8)') "Total electron moment:", sum
      endif
      close(10)

! produce file for moment real space integration
      integrate_file = .true.
      if( integrate_file ) then
         write(*,*) "Producing a file to be integrated:"
         open(unit=50,file='integrate_moment',status='unknown',form='formatted')
         write(50,*) vol
         write(50,*) (rprim(j,1),j=1,3)
         write(50,*) (rprim(j,2),j=1,3)
         write(50,*) (rprim(j,3),j=1,3)
         cnt = 0
         sum = 0.0

! find a good origin first by shifting the requested atomic site
! onto the 000 mesh point
         delta = 100.0
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  cnt = cnt + 1
                  v(:) = (i-1)*rprim(:,1)/nx+(j-1)*rprim(:,2)/ny+&
                  (k-1)*rprim(:,3)/nz
                  v(:) = v(:) - origin(:)
                  delta2 = sqrt( v(1)*v(1)+v(2)*v(2)+v(3)*v(3) )
                  if( delta2 < delta ) then
                     delta = delta2
                     cnt_find = cnt
                     dorigin(:) = v(:)
                     dens1 = density(1,cnt)/vol
                     if( nspin == 2 ) then
                        dens2 = 0.5*(density(1,cnt)+density(2,cnt))/vol
                     else
                        dens2 = density(1,cnt)/vol
                     endif
                  endif
               enddo
            enddo
         enddo

! write the origin density as the first point, we fake this point
! by putting in "0" and then double counting the point with the
! correct value
         write(50,'((3f12.6,g16.6))') 0.0, 0.0, 0.0, dens1
         sum = sum + dens1
         cnt = 0
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  cnt = cnt + 1
                  if( cnt /= cnt_find ) then
                     v(:) = (i-1)*rprim(:,1)/nx+(j-1)*rprim(:,2)/ny+&
                     (k-1)*rprim(:,3)/nz
                     v(:) = v(:) - origin(:) -dorigin(:)
                     sum = sum + density(1,cnt)/vol
                     write(50,'(3f12.6,g16.6)') v(:), density(1,cnt)/vol
                  endif
               enddo
            enddo
         enddo
         sum = sum*vol/nx/ny/nz
         write(*,*) "Total number of electrons:", sum
         cnt = 0
         sum = 0.0
! write the origin density as the first point
         write(50,'((3f12.6,g16.6))') 0.0, 0.0, 0.0, dens2
         sum = sum + dens2
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  cnt = cnt + 1
                  if( cnt /= cnt_find ) then
                     v(:) = (i-1)*rprim(:,1)/nx+(j-1)*rprim(:,2)/ny+&
                     (k-1)*rprim(:,3)/nz
                     v(:) = v(:) - origin(:) - dorigin(:)
                     write(50,'(3f12.6,g16.6)') v(:), &
                     (density(1,cnt)+density(2,cnt))/vol/2.0
                     sum = sum + 0.5*(density(1,cnt)+density(2,cnt))/vol
                  endif
               enddo
            enddo
         enddo
         sum = sum*vol/nx/ny/nz
         write(*,*) "Total number of up electrons:", sum
         close(50)
      endif


      end
