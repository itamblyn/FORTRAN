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

integer i, l, ndata, maxline, counter
integer(8) plan
parameter(maxline=10000000)      ! this just sets a max number of lines that can be read in
complex, allocatable :: vacf_array(:)
complex, allocatable :: fft_distance_array(:)
real(8), allocatable :: freq_array(:)
real(8), allocatable :: psd_array(:)
real(8) dt, df
real(8) power, dummyf, vacf
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

counter = 0

do i=1,maxline 
  read(1,*,END=100)
  counter = counter + 1
end do

100  continue

if (counter == maxline)  print *, "Error, did not read to end of file" 

ndata = counter

df = 1.0/(ndata*dt)

print *, "Timestep: ", dt, " fs"
print *, "Freqstep: ", df, " 1/fs"

rewind(1)               ! this puts the read head back at top of file

allocate(vacf_array(ndata))
allocate(freq_array(ndata))
allocate(fft_distance_array(ndata))
allocate(psd_array(ndata))

psd_array(:) = 0.0
freq_array(:) = 0.0

do i=1,ndata

    read(1,*) dummyc, vacf
    vacf_array(i) = cmplx(vacf,0.0)
    freq_array(i) = (i - 1)*df*(1.0/29979245800.0)

end do

call sfftw_plan_dft_1d ( plan, ndata, vacf_array, fft_distance_array, FFTW_FORWARD, FFTW_ESTIMATE )
call sfftw_execute ( plan )
call sfftw_destroy_plan ( plan )
    
do l=1,ndata
    psd_array(l) = psd_array(l) +  (1.0/(2.0*ndata + 1.0))*(conjg(fft_distance_array(l))*fft_distance_array(l))
end do

power = 0

do i=1,ndata
  power = power + psd_array(i)
end do

psd_array = psd_array/power

open(2,file="psd.dat")

do i=2,ndata/2
  write(2,*) freq_array(i), psd_array(i)
end do

close(2)

deallocate(vacf_array)
deallocate(freq_array)
deallocate(fft_distance_array)
deallocate(psd_array)

END PROGRAM psd 
