!       program local.f90
!***********************************************************
!
!       
!***********************************************************
!
        implicit none

        integer n, neighbours, i, j, k, l
        parameter (n=128)
        parameter (neighbours=5)
        real(8) snapshot(n,3), origin(3), cell(3), dummyf, unwrapped(3), wrapped(3)
        integer nn(n,neighbours)
        character(128) fin
        character(4) element(neighbours + 1), dummyc 

        element(2) ="A"
        element(3) ="B"
        element(4) ="C"
        element(5) ="D"
        element(6) ="E"

        write(6,*) 'Name of data file'
        read(5,*) fin

        open(1,file=fin,status='old',ERR=100)
        open(2,file="TRAJEC.xyz")

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

        do i=1,1000

          do j=1,128
            read(1,*) element(1), (snapshot(j,k),k=1,3) , (nn(j,k),k=1,5)

!           converts from python to fortran array index 

            nn(j,1) = nn(j,1) + 1
            nn(j,2) = nn(j,2) + 1
            nn(j,3) = nn(j,3) + 1
            nn(j,4) = nn(j,4) + 1
            nn(j,5) = nn(j,5) + 1

!            print *, nn(j,1), nn(j,2), nn(j,3), nn(j,4), nn(j,5)

!           converts to reduced coordinates

            snapshot(j,1) = snapshot(j,1)/cell(1)
            snapshot(j,2) = snapshot(j,2)/cell(2)
            snapshot(j,3) = snapshot(j,3)/cell(3)

!           wraps particles into box

            snapshot(j,1) = snapshot(j,1) - anint(snapshot(j,1))
            snapshot(j,2) = snapshot(j,2) - anint(snapshot(j,2))
            snapshot(j,3) = snapshot(j,3) - anint(snapshot(j,3))

          end do

!        coordinates of origin

          origin(1) = snapshot(1,1)
          origin(2) = snapshot(1,2)
          origin(3) = snapshot(1,3)

          write(2,*) neighbours + 1
          write(2,*) i
          write(2,*) element(1), origin(1) - origin(1), origin(2) - origin(2), origin(3) - origin(3)
       
          do k=1,neighbours

            snapshot(nn(1,k),1) = snapshot(nn(1,k),1) - origin(1)
            snapshot(nn(1,k),2) = snapshot(nn(1,k),2) - origin(2)
            snapshot(nn(1,k),3) = snapshot(nn(1,k),3) - origin(3)

            wrapped(1) = snapshot(nn(1,k),1) - anint(snapshot(nn(1,k),1))
            wrapped(2) = snapshot(nn(1,k),2) - anint(snapshot(nn(1,k),2))
            wrapped(3) = snapshot(nn(1,k),3) - anint(snapshot(nn(1,k),3))

            write(2,*) element(k+1), wrapped(1)*cell(1)*0.529177, wrapped(2)*cell(2)*0.529177, wrapped(3)*cell(3)*0.529177

          end do
        
        end do

100    continue

        close(1)
        close(2)
        END
