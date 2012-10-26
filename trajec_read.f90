!       program traj_read.f90
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

!write(6,*) 'Name of data file'
!read(5,*) fin

fin='TRAJEC.xyz'

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

print *, nsteps

rewind(1)               ! this puts the read head back at top of file

allocate(trajectory(nsteps*natom,3))

do i=1,nsteps

    read(1,*)   !       natom
    read(1,*)   !       timestep or cell, depending on xyz convention
    do j=1,natom
        read(1,*) dummyc, (trajectory((i-1)*natom + j,k),k=1,3)
    end do

end do

print *, (trajectory(nsteps*natom, k), k=1,3)

END
