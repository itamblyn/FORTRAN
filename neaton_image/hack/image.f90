PROGRAM image

!       program image.f90
!***********************************************************
!
!       
!***********************************************************
!

implicit none

integer natom, nsteps, i, j, k, l, counter, maxline
parameter(maxline=10000000)      ! this just sets a max number of lines that can be read in
real(4), allocatable :: trajectory(:,:)
character(128) fin
character(4) dummyc 

integer n_carbon, n_hydrogen, n_nitrogen, first_gold
real(4) carbon_z, surface_z, d, W, QP, sigma

parameter(n_carbon = 6)
parameter(n_hydrogen = 8)
parameter(n_nitrogen = 2)

print *, 'Warning, this program has a CRAZY number of magic numbers!!!'

print *, 'n_carbon', n_carbon, 'n_hydrogen', n_hydrogen, 'n_nitrogen', n_nitrogen

write(6,*) 'Name of data file'
read(5,*) fin

write(6,*) 'Index of first SURFACE gold atom (e.g. 2 in the case of adatom)'
read(5,*) first_gold

open(1,file=fin,status='old',ERR=100)

read(1,*) natom         ! the first line should be the number of atoms

counter = 0

rewind(1)               ! start at the beginning of the file

do i=1,maxline 
  read(1,*,END=100)
  counter = counter + 1
end do

100  continue

if (counter == maxline)  print *, "Error, did not read to end of file" 

nsteps = counter/(natom + 2)

!print *, nsteps

rewind(1)               ! this puts the read head back at top of file

allocate(trajectory(nsteps*natom,3))

do i=1,nsteps

    read(1,*)   !       natom
    read(1,*)   !       timestep or cell, depending on xyz convention
    do j=1,natom
        read(1,*) dummyc, (trajectory((i-1)*natom + j,k),k=1,3)
    end do

end do

! trajectory read in


do i=1,nsteps
  carbon_z = 0.0
  do j=1,n_carbon
    carbon_z = carbon_z + trajectory((i-1)*natom + j,3)
  enddo
  
  surface_z = 0.0
  do j=(n_carbon + n_hydrogen + n_nitrogen + first_gold), (n_carbon + n_hydrogen + n_nitrogen + ((first_gold - 1) + 16))  ! 16 atom surface
    surface_z = surface_z + trajectory((i-1)*natom + j,3)
  end do
  carbon_z = carbon_z/n_carbon
  surface_z = 12.0949700575000 ! 16 atom surface

  print *, 'surface_plane', surface_z

  d = ((carbon_z - surface_z) - 1.0)/0.529177
!  print *, d
  W = -13.605692/(2*d)
  QP = -2.75
  sigma = QP - W 

  print *, "c6_z - surface_z: ", carbon_z - surface_z, "sigma: ", sigma
  print *, "QP: =" , QP
  print *, "W: =" , W
 

end do

END PROGRAM image

SUBROUTINE pbc_round(input_value, i)

implicit none

real(8), intent(IN) :: input_value
real(8), intent(OUT) :: i

i = int(input_value)

if (abs(input_value-i) >= 0.5) then
  if (input_value > 0) i = i + 1.0
  if (input_value < 0) i = i - 1.0
end if

END SUBROUTINE pbc_round
