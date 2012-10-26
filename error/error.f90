!       program error.f
!***********************************************************
!
!       generate padding for fitting with error bars
!       
!***********************************************************
        implicit none

        integer n, counter, i, j, k
        parameter (n=100000)
        real(8) data(n,4)
        real(8) xvalue, yvalue, xerror, yerror, xstep, ystep, xnew, ynew
        integer numxstep, numystep
!        parameter(numxstep=200,numystep=200)
        character(128) fin

!        write(6,*) 'Name of data file'
        read(5,*) fin

!        write(6,*) 'Number of x steps?'
        read(5,*) numxstep
!        write(6,*) 'Number of y steps?'
        read(5,*) numystep


        open(1,file=fin,status='old',ERR=999)

        open(2,file="error.dat")

        read(1,*) ! skips the header

        counter = 0

        do i=1,n
          read(1,*,END=100) data(i,1), data(i,2), data(i,3), data(i,4)
          counter = counter + 1
        end do

 100    continue

        close(1)

!        do i=1,counter
!          write(2,*) data(i,1), data(i,2), data(i,3), data(i,4)
!        end do 

!        close(2)

        write(2,*) ' I V'

        do i=1, counter
          xvalue = data(i,1)
          yvalue = data(i,2)
          xerror = data(i,3)
          yerror = data(i,4)

          xstep = 2*(xerror/numxstep)
          ystep = 2*(yerror/numystep)

          do j=0,numxstep
          
            do k=0,numystep
          
              xnew = (xvalue + j*xstep) - xerror
              ynew = (yvalue + k*ystep) - yerror
           
              write(2,*) xnew, ynew
          
            end do

          end do

        end do

        close(2)

 999    continue

        END
