PROGRAM psd 
!
!       program psd.f90
!***********************************************************
!
!       
!***********************************************************
!

implicit none

include 'fftw3.f'

integer natom, nsteps, i, j, k, l, counter, maxline, dummyi
integer(8) plan
parameter(maxline=10000000)      ! this just sets a max number of lines that can be read in
real(4), allocatable :: trajectory(:,:)
complex, allocatable :: distance_array(:)
complex, allocatable :: fft_distance_array(:)
real(4), allocatable :: freq_array(:)
real(4), allocatable :: psd_array(:)
real(4) rounded, dx, dy, dz, distance, dt, df, time
real(4) lattice_constant, size_of_box_x, size_of_box_y, size_of_box_z 
real(4) power
character(128) fin
character(4) dummyc 

call sfftw_init_threads
call sfftw_plan_with_nthreads(8)

write(6,*) 'Name of data file'
read(5,*) fin

write(6,*) 'Timestep [A.U.]?'
read(5,*) dt

dt = (2.4188e-17)*dt

open(1,file=fin,status='old',ERR=100)
read(1,*) natom         ! the first line should be the number of atoms

write(6,*) 'Lattice constant contained in file? [0=no, 1=yes] '
read(5,*) dummyi

if (dummyi == 0) then
  write(6,*) 'Lattice constant [A]?'
  read(5,*) lattice_constant
  size_of_box_x = lattice_constant
  size_of_box_y = lattice_constant
  size_of_box_z = lattice_constant

else
  read(1,*) size_of_box_x, size_of_box_y, size_of_box_z

end if

rewind(1)               ! start at the beginning of the file

counter = 0

do i=1,maxline 
  read(1,*,END=100)
  counter = counter + 1
end do

100  continue

if (counter == maxline)  print *, "Error, did not read to end of file" 

nsteps = counter/(natom + 2)

write(6,*) 'Detected ', nsteps, ' steps'

df = 1.0/(nsteps*dt)

print *, "Timestep: ", dt, " fs"
print *, "Freqstep: ", df, " 1/fs"
print *, "Detected: ", natom,  " atoms"
print *, "Detected: ", nsteps, " steps"

rewind(1)               ! this puts the read head back at top of file

allocate(trajectory(nsteps*natom,3))
allocate(distance_array(nsteps))
allocate(freq_array(nsteps))
allocate(fft_distance_array(nsteps))
allocate(psd_array(nsteps))

psd_array(:) = 0.0
freq_array(:) = 0.0

do i=1,nsteps

    read(1,*)   !       natom
    read(1,*)   !       timestep or cell, depending on xyz convention
    do j=1,natom
        read(1,*) dummyc, (trajectory((i-1)*natom + j,k),k=1,3)
    end do

    freq_array(i) = (i - 1)*df*(1.0/29979245800.0)

end do

!
! at this point, the trajectory file has been read in
!

! $omp parallel default(private) shared(trajectory, psd_array) 
! $omp do

do j=1,natom

  do k=j+1, natom

    do i=1,nsteps

      dx = trajectory( ((i-1)*natom + j), 1) - trajectory( ((i-1)*natom + k), 1)   
      dy = trajectory( ((i-1)*natom + j), 2) - trajectory( ((i-1)*natom + k), 2)   
      dz = trajectory( ((i-1)*natom + j), 3) - trajectory( ((i-1)*natom + k), 3)

      CALL pbc_round(dx/size_of_box_x, rounded)
      dx = dx - size_of_box_x*rounded

      CALL pbc_round(dy/size_of_box_y, rounded)
      dy = dy - size_of_box_y*rounded

      CALL pbc_round(dz/size_of_box_z, rounded)
      dz = dz - size_of_box_z*rounded

      distance = (dx**2 + dy**2 + dz**2)**(0.5)

      distance_array(i) = cmplx(distance, 0.0)

    end do

    call sfftw_plan_dft_1d ( plan, nsteps, distance_array, fft_distance_array, FFTW_FORWARD, FFTW_ESTIMATE )

    call sfftw_execute ( plan )

    call sfftw_destroy_plan ( plan )
    
    do l=1,nsteps
! $omp critical
      psd_array(l) = psd_array(l) +  (1.0/(2.0*nsteps + 1.0))*(conjg(fft_distance_array(l))*fft_distance_array(l))
! $omp end critical
    end do
! $omp flush (psd_array)

  end do

end do

! $omp end do
! $omp end parallel

power = 0

do i=1,nsteps
  power = power + psd_array(i)
end do

psd_array = psd_array/power

open(2,file="psd.dat")

do i=2,nsteps/2
  write(2,*) freq_array(i), psd_array(i)
end do

close(2)

deallocate(trajectory)
deallocate(distance_array)
deallocate(freq_array)
deallocate(fft_distance_array)
deallocate(psd_array)

END PROGRAM psd 

SUBROUTINE pbc_round(input_value, i)

implicit none

real, intent(IN) :: input_value
real, intent(OUT) :: i

i = int(input_value)

if (abs(input_value-i) >= 0.5) then
  if (input_value > 0) i = i + 1
  if (input_value < 0) i = i - 1
end if

END SUBROUTINE pbc_round
