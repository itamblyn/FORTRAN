!       program vdos.f90
!***********************************************************
!       Calculate normalized VDOS from normalized VACF
!***********************************************************

implicit none

integer i, j, Nsteps
real(8) pi, f, t, X, dt, dummyf
real(8) sum,vacf(:),cosXvacf(:,:)
character(128) fin
allocatable vacf,cosXvacf

! Define pi
pi = 2.0*asin(1.0)

! Retrieve user input
write(*,*) "Output from VACF code"
read(*,*) fin
write(*,*) "Enter timestep in a.u."
read(*,*) dt

! Open necessary files
open(1,file=fin,status='old')
open(2,file="vdos.dat",status='unknown')

! Get the number of steps
Nsteps = 0
100 read(1,*,END=200) 
Nsteps = Nsteps + 1
goto 100
200 continue
REWIND (1) 

! Allocate vacf array and read from file
allocate( vacf(Nsteps) )
do i=1,Nsteps
!  read(1,*) t, vacf(i)
  read(1,*) t, dummyf, vacf(i) 
end do
close(1)

! Convert timestep to fs
dt = dt * 2.418800d-5 * 1000.0
dt = dt / 1000.0

! Allocate cosXvacf array
allocate(cosXvacf(1001,Nsteps))

! Calculate the cosXvacf array
do i=0,1000
  f = i * 0.1
  do j = 1,Nsteps
    cosXvacf(i,j) = dcos(2.0*pi*f*(j-1)*dt)*vacf(j)
  enddo
enddo

! Calculate VDOS and write to file
do i=0,1000
  sum = cosXvacf(i,1) + cosXvacf(i,Nsteps/2)
  X = 4.0
  do j=2,Nsteps/2-1
    sum = sum + X*cosXvacf(i,j)
    X = 6.0 - X
  enddo
  write(2,*) i*0.1, 4.0*sum*dt!/3.0
enddo

close(2)

END
