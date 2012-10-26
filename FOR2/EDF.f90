!       program EDF.f90
!***********************************************************
!
!       
!***********************************************************
!

implicit none

integer i, j, k, rho_x, rho_y, rho_z, nl
integer natom, counter, stride
real(8) acell, sum
character(128) fin

double precision, allocatable :: tmp(:)
double precision, allocatable :: coor(:,:)
!double precision, allocatable :: density(:,:,:)
double precision, allocatable :: density(:)

fin='CHGCAR'

stride = 5      ! this is how many columns there are
allocate(tmp(stride))

open(1,file=fin,status='old',ERR=100)

read(1,*)       ! skips title 
read(1,*) acell ! acell in Angstrom
do i=1,3
  read(1,*)     ! skips lattice vectors
end do
read(1,*)  natom
write(*,*) natom
allocate(coor(natom,3))
read(1,*) 

do i=1,natom
  read(1,*) (coor(i,1), j=1,3)
end do

read(1,*)

read(1,*) rho_x, rho_y, rho_z
!allocate(density(rho_x, rho_y, rho_z))
allocate(density(rho_x*rho_y*rho_z))

nl = rho_x*rho_y*rho_z/stride

sum = 0.0
counter = 0

do i = 1, nl
  read(1,*) (tmp(j),j=1,stride)
    do j = 1, stride
      counter = counter + 1
      density(counter) = tmp(j)
      sum = sum + tmp(j)
      if( tmp(j) .lt. 0.0 ) then
        write(*,'(a)') "Negative density?"
      endif
    enddo
enddo
sum = sum/float(rho_x*rho_y*rho_z)
if( counter < rho_x*rho_y*rho_z ) then
  read(1,*) (tmp(j-counter),j=counter+1,rho_x*rho_y*rho_z)
  do j = counter+1, rho_x*rho_y*rho_z
    density(counter) = tmp(j-counter)
    sum = sum + tmp(j-counter)
  enddo
endif

write(*,*) "Total electron charge:", sum


close(1)

100 continue

END
