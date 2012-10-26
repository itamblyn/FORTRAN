PROGRAM eig_smear 
!
!       program eig_smear.f90
!***********************************************************
!
!       
!***********************************************************
!

implicit none

integer i, j, k, nkpt, neig, nedos 
integer maxline, counter
parameter(maxline=10000000)      ! this just sets a max number of lines that can be read in
real(8) E, emin, emax, dE, weight, eigenvalue, sigma
real(8) Ha_sum
real(8), allocatable :: Ha_dos(:,:)
real(8), allocatable :: wtk(:)
real(8), allocatable :: Ha_eig_array(:,:)
real, parameter :: pi = 3.14159 ! cannot be changed
character(128) fin, fin2
character(6) dummya, dummyb, dummyc, dummyd, dummye 

nedos = 11001 
sigma = 0.2/27.211383

write(6,*) 'Name of eigenvalue data file '
read(5,*) fin

open(1,file=fin,status='old',ERR=100)

read(1,*) dummya, dummyb, nkpt, dummyc, neig, dummyd
read(1,*)

allocate(Ha_dos(3,nedos))
allocate(wtk(nkpt))
allocate(Ha_eig_array(nkpt,neig))

counter = 2

do i=1,maxline 
  read(1,*,END=100)
  counter = counter + 1
end do

100  continue

if (counter == maxline)  print *, "Error, did not read to end of file" 

wtk(:) = 1.00

rewind(1)               ! this puts the read head back at top of file
read(1,*)               ! skips the header
read(1,*)               ! skips the header

Ha_dos(:,:) = 0.0
Ha_eig_array(:,:) = 0.0

do j=1,nkpt
    do k=1,neig
        read(1,*) Ha_eig_array(j,k)
    end do
end do

emin = minval(Ha_eig_array) - 2.0/27.211383  ! so we can see +/- 2 eV
emax = maxval(Ha_eig_array) + 2.0/27.211383

dE = (emax - emin)/nedos

do i=1,nedos

    E = (i - 1)*dE + emin  ! should I change this to middle..?
                           ! no, this is the left side of the box for
                           ! integration

    do j=1,nkpt

        weight = 1.0*wtk(j)  ! 1 e per state

        do k=1,neig

            eigenvalue = Ha_eig_array(j,k)
            Ha_dos(2,i) = Ha_dos(2,i) + weight*(1.0/2.0)*(erf(sqrt(2.0)*(E + dE -eigenvalue)/(2*sigma)) - erf(sqrt(2.0)*(E - eigenvalue)/(2*sigma)))

        end do 

    end do

end do

open(3,file="eV_dos.dat")

Ha_sum = 0.0

do i=1,nedos
  E = (i-1)*dE + emin + dE/2
  Ha_sum = Ha_sum + Ha_dos(2,i)
  write(3,*) E, Ha_dos(2,i)*27.211383, Ha_sum*27.211383
end do

close(3)

deallocate(wtk)
deallocate(Ha_eig_array)
deallocate(Ha_dos)

END PROGRAM eig_smear 
