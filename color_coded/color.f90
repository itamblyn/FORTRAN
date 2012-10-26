!       program color.f90
!***********************************************************
!
!       
!***********************************************************
!

implicit none

integer i, j 
integer n, coordination
parameter (n=256)
real(8) snapshot(3)
character(5) filecounter
character(128) fin
character(2) dummyc, element(n)

element = 'X' ! sets all to X

element(1) = 'H'
element(2) = 'O'
element(3) = 'N'
element(4) = 'C'
element(5) = 'S'

write(6,*) 'Name of data file'
read(5,*) fin

open(1,file=fin,status='old',ERR=100)
open(2,file="color.xyz")

do i=1,1000

!  write(unit=filecounter, fmt='(I5)') (i - 1) + 10000
!  open(2,file="color" // filecounter // ".xyz")

  write(2,*) n
  write(2,*) i

  do j=1,n
    read(1,*) coordination, dummyc, snapshot(1), snapshot(2), snapshot(3)
    write(2,*) element(coordination), snapshot*0.529177
  end do

end do

100    continue

close(1)
close(2)

END
