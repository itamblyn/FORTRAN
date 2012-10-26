PROGRAM vcf
!
!      program vcf.f90
!***********************************************************
!
!       
!***********************************************************
!

implicit none

integer natom, max_natom, nsteps, max_nsteps, neighbours, max_neighbours
integer i, j, k, l, n, neighbour_dex, bin, numbins, counter
parameter (max_natom=256)
parameter (max_nsteps=10000)
parameter (max_neighbours=127)
parameter (numbins=500)
real(8) velocity(max_natom,3)
real(8) trajectory(max_natom*max_nsteps,3)
real(8) dot_histogram(max_neighbours)
real(8) distance_histogram(numbins + 1)
integer dot_counter(max_neighbours)
integer distance_counter(numbins + 1)
real(8) velocity_a(3), velocity_b(3), snapshot(max_natom,3)
real(8) speed_a, speed_b, dx, dy, dz, rounded, dt, minbin, maxbin, binsize,x 
real(8) p(3), o(3), displacement(3)
integer nn(max_natom*max_nsteps,max_neighbours)
real(8) cell(3), dot, distance, dummyf
character(128) fin, fout_nn, fout_dis
character(4) dummyc 


minbin = 0.0
maxbin = 20.0
binsize = (maxbin - minbin)/numbins

write(6,*) 'Name of data file [.cnn]'
read(5,*) fin
write(6,*) 'Timestep [A.U.]?'
read(5,*) dt 

fout_nn="nn.dat"
fout_dis="dist.dat"

open(1,file=fin,status='old',ERR=200)
open(2,file=fout_nn)
open(3,file=fout_dis)

counter = 0

do i=1,max_nsteps*max_natom

  read(1,*,END=100)
  counter = counter + 1

end do

100    continue

if (counter == max_nsteps*max_natom) print *, "Did not reach end of file, increase max_nsteps"

rewind(1)

read(1,*)
   
do i=1,3
  read(1,*) dummyc, dummyc, dummyc, cell(i)
end do

print *, cell

read(1,*) dummyc, dummyc, dummyc, natom
read(1,*) dummyc, dummyc, dummyc, neighbours

nsteps = (counter - 10)/natom

print *, "Detected ", natom, " atoms"
print *, "Detected ", nsteps, " steps"
print *, "Detected ", neighbours, " neighbours"

dot_histogram(:) = 0.0
distance_histogram(:) = 0.0
dot_counter(:) = 0
distance_counter(:) = 0

!       skips rest of header

do i=1,4
  read(1,*)
end do

print *, "Entering main loop"

do i=1,nsteps

  do j=1,natom

    read(1,*,END=200) dummyc, (trajectory(((i-1)*natom + j),k),k=1,3), (nn(((i-1)*natom + j),k),k=1,neighbours)

    do k=1,neighbours
      nn(((i-1)*natom + j),k) = nn(((i-1)*natom + j),k) + 1
    end do

  end do

end do

print *, "File read complete"

do i=2,nsteps

  do j=1,natom
    do k=1,3
      velocity(j,k) = (trajectory((i-1)*natom + j,k) - trajectory((i-2)*natom + j, k))/dt
    end do
  end do

!  velocities were calculated from unwrapped array

  do j=1,natom
    
!  gets the coordinates of atom j, stores them in vector p
    
    do k=1,3
      p(k) = trajectory((i-1)*natom + j,k)
    end do
    
    do n=1,natom
    
      do k=1,3
!       subtracts of coordinates of atom j from all atoms, including j
!       wraps all atoms into box centered around j
        o(k) = trajectory((i-1)*natom + n,k) - p(k)
        snapshot(n,k) = o(k) - (int(o(k)/cell(k)+nsteps+0.5)-nsteps)*cell(k)
      end do
    
    end do 
    
    do neighbour_dex=1,neighbours

      l = nn((i-1)*natom + j,neighbour_dex) 

      distance = 0.0

      do k=1,3
        displacement(k) = snapshot(l, k)
        distance = distance + snapshot(l,k)**2
      end do

      distance = sqrt(distance) 
      displacement = displacement/distance

      do k=1,3
        velocity_a(k) = velocity(j,k)
        velocity_b(k) = velocity(l,k)
      end do

      velocity_a = velocity_a*displacement
      velocity_b = velocity_b*displacement

!     at this point both velocity vectors have been projected onto the unit
!     displacement vector

      speed_a = 0.0
      speed_b = 0.0

      do k=1,3
        speed_a = speed_a + velocity_a(k)**2
        speed_b = speed_b + velocity_b(k)**2
      end do

      speed_a = sqrt(speed_a)
      speed_b = sqrt(speed_b)

      dot = 0.0

      do k=1,3

        dot = dot + (velocity_a(k)/speed_a)*(velocity_b(k)/speed_b)

      end do


      bin = int((distance - minbin)/binsize)

      if ( speed_a > 0 .and. speed_b > 0) then

        dot_histogram(neighbour_dex) = dot_histogram(neighbour_dex) + dot
        dot_counter(neighbour_dex)   = dot_counter(neighbour_dex) + 1

        if (distance < maxbin) then
 
          distance_histogram(bin) = distance_histogram(bin) + dot
          distance_counter(bin)   = distance_counter(bin) + 1

        end if

      end if

    end do

  end do

end do

do i=1,neighbours

  if (dot_counter(i) > 0) dot_histogram(i) = dot_histogram(i)/dot_counter(i)
  write(2,*) i, dot_histogram(i)

end do

write(3,*) "# distance[Angst], probability (unnormalized for r)"

do i=1,numbins

  x = (i*binsize + 0.5*binsize)*0.529177

  if (distance_counter(i) > 0) distance_histogram(i) = distance_histogram(i)/distance_counter(i)
  write(3,*) x, distance_histogram(i)

end do

200    continue

close(1)
close(2)
close(3)

END PROGRAM vcf

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

