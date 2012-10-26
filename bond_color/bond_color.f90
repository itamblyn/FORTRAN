!       program bond_color.f90
!***********************************************************
! this program reads in a cbn file, and produces an xyz file
! that is colored based on whether or not an atom is molecular
! or atomic
!***********************************************************
!
implicit none

integer n,i,j,k, bonded, nsteps, maxline, counter
parameter (maxline=100000000)
parameter (n=128)
real(4) vector(3), cell(3), dummyf
character(128) fin
character(4) dummyc, element, state(2)

state(1) = "A"
state(2) = "M"

write(6,*) 'Name of data [.cbn] file'
read(5,*) fin

open(1,file=fin,status='old',ERR=100)

do i=1,maxline 
  read(1,*,END=100)
  counter = counter + 1
end do

100  continue

if (counter == maxline)  print *, "Error, did not read to end of file" 

nsteps = (counter - 10)/(n)

print *, nsteps

rewind(1)               ! this puts the read head back at top of file

!       skips header

read(1,*)
   
do i=1,3
  read(1,*) dummyc, dummyc, dummyc, cell(i)
end do

print *, cell

!       skips rest of header

do i=1,6
  read(1,*)
end do

open(2,file="bonding.xyz")

do i=1,nsteps

  write(2,*) n
  write(2,*) i

  do j=1,n

    read(1,*) element, (vector(k),k=1,3), bonded

    if ( bonded == -1 ) then
      write(2,*) state(1), vector(1)*0.529177, vector(2)*0.529177, vector(3)*0.529177 

    else  
      write(2,*) state(2), vector(1)*0.529177, vector(2)*0.529177, vector(3)*0.529177 

    end if

  end do

end do
        

close(1)
close(2)

END
