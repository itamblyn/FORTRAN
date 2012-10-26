! Crystal: A simple program to generate fcc coordinates and write them
! in a coordinate file.  This file can be read by the companion MD programs 
! (and also by the BallRoom plotting program for visualization) and
! provides a starting point for a simulation.
! Crystal is meant to be an interactive program: prompts for the input
! parameters are provided, as well as some suggestions for the two 
! materials simulated by the example MD programs (Lennard-Jonesium and
! Aluminum).
!
! Furio Ercolessi, SISSA, Trieste, May 1995, revised May 1997

program Crystal
implicit none
character*80 :: FileName
double precision :: alat,displac
integer :: nx,ny,nz
call Read_Parameters(FileName,alat,nx,ny,nz,displac)
call Generate_Crystal(FileName(1:len_trim(FileName)), &
                      alat,nx,ny,nz,displac)
end program Crystal


subroutine Read_Parameters(FileName,alat,nx,ny,nz,displac)
!
!  Obtain all parameters from user
!
implicit none
integer, parameter :: crtout = 6
character*(*) :: FileName
double precision :: alat,displac
integer :: nx,ny,nz
logical :: AlreadyThere
double precision, parameter :: cutoffLJ = 2.5d0
double precision, parameter :: cutoffAl = 5.55805441821810d0

write( crtout,'(/,a,/)') 'Crystal: generate fcc coordinates'
write( crtout,'(a,/,a,/)') &
   'Reminder: for LJ truncated at 2.5, a_eq = 1.5496 at P=0', &
   '          for the Al glue system , a_eq = 4.032 A at P=0'
write( crtout,'(a)',advance='NO') 'Lattice spacing a = '
read(*,*,end=900) alat
write( crtout,'(/,a,i3,a,/,a,i3,a,/)') &
   'Reminder: minimum image requires at least', &
    int(2.d0*cutoffLJ/alat) + 1,' cells along x,y,z for LJ', &
   '                                         ', &
    int(2.d0*cutoffAl/alat) + 1,' cells along x,y,z for Al'
write( crtout,'(a)',advance='NO') 'Number of cells along x = '
read(*,*,end=900) nx
write( crtout,'(a)',advance='NO') 'Number of cells along y = '
read(*,*,end=900) ny
write( crtout,'(a)',advance='NO') 'Number of cells along z = '
read(*,*,end=900) nz
write( crtout,'(a)',advance='NO') 'Maximum random displacement = '
read(*,*,end=900) displac
800 continue
write( crtout,'(a)',advance='NO') 'Name of coordinates file = '
read(*,'(a)',end=900) FileName
!
!  immediately check if FileName already exists
!
inquire(file=FileName(1:len_trim(FileName)),exist=AlreadyThere)
if (AlreadyThere) then
   write( crtout ,'(2a)') FileName(1:len_trim(FileName)), &
                        ' already exists, please choose another name.'
   go to 800
endif
return       ! regular exit
900 continue ! incomplete user input
stop
end subroutine Read_Parameters


subroutine Generate_Crystal(FileName,alat,nx,ny,nz,displac)
!
!  Does the actual work of generating the atomic coordinates
!
implicit none
integer, parameter :: crtout = 6
character*(*) :: FileName
double precision :: alat,displac
double precision, dimension(3) :: BoxSize,rands
integer :: nx,ny,nz
integer :: N,i,j,k,L
double precision :: x,y,z
integer, parameter :: nbase = 4                 ! number of atoms in fcc cell
double precision, dimension(3,4), parameter :: rcell = &   ! coords of atoms
     reshape( (/ 0.0d0 , 0.0d0 , 0.0d0 , &
                 0.5d0 , 0.5d0 , 0.0d0 , &
                 0.0d0 , 0.5d0 , 0.5d0 , &
                 0.5d0 , 0.0d0 , 0.5d0 /) , (/ 3,4 /) )

!
!  Open file for coordinates
!
open(unit=1,file=FileName,status='new',form='formatted', &
     action='write',position='rewind',err=700)
!
!  Define and write number of particles N, and box size along x,y,z
!  on the first line of the file
!
N = 4*nx*ny*nz
BoxSize(1) = nx*alat
BoxSize(2) = ny*alat
BoxSize(3) = nz*alat
!
!  The '%' signals to ignore this line to the BallRoom program producing
!  an image from the simulation
!
write(1,'(A1,L2,I7,3E23.15)',err=800) '%',.FALSE.,N,BoxSize
!
!  Generate atomic coordinates and write them
!
do k=0,nz-1
   do j=0,ny-1
      do i=0,nx-1
         do L=1,nbase
            call random_number(rands)       ! get 3 random numbers in (0,1)
            x = alat * ( i + rcell(1,L) ) + 2.d0*displac*( rands(1) - 0.5d0 )
            y = alat * ( j + rcell(2,L) ) + 2.d0*displac*( rands(2) - 0.5d0 )
            z = alat * ( k + rcell(3,L) ) + 2.d0*displac*( rands(3) - 0.5d0 )
            write(1,'(1X,3E23.15)',err=800) x,y,z
         enddo
      enddo
   enddo
enddo
endfile(unit=1) ! this is because we can be overwriting a previous, longer file
close(unit=1)
print*,N,' coordinates written on ',FileName
return
700 continue
   print*,'Generate_Crystal: FATAL: cannot open ',FileName
   stop
800 continue
   print*,'Generate_Crystal: FATAL: write error in ',FileName
   close(unit=1)
   stop
end subroutine Generate_Crystal
