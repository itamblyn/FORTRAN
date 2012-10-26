!       program msd.f
!***********************************************************
!       Calculate the mean squared displacement from
!       an xyz file
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer n, m, q, i, j, k
        parameter (n=1000000,m=1024)
        integer natoms(n), timesteps(n)
        real(8) Rx0(m), Ry0(m), Rz0(m)
        real(8) Rx(m), Ry(m), Rz(m)
        real(8) dRx, dRy, dRz, sum_msd, msd(n)
        real(8) D(n)
        character(128) fin, dummy

        write(6,*) 'Name of unwrapped xyz file'
        read(5,*) fin

        open(1,file=fin,status='old',ERR=90)

        read(1,*) natoms(1)
        read(1,*) timesteps(1)
        do j=1,natoms(1)
          read(1,*) dummy, Rx0(j), Ry0(j), Rz0(j)
        enddo        

        open(2,file='msd.dat')

        do i=2,n

          read(1,*,END=100) natoms(i)
          read(1,*) timesteps(i)

          sum_msd = 0.0

          do j=1,natoms(i)

            read(1,*) dummy, Rx(j), Ry(j), Rz(j)

            dRx = Rx(j) - Rx0(j)
            dRy = Ry(j) - Ry0(j)
            dRz = Rz(j) - Rz0(j)
            sum_msd = sum_msd + (dRx**2 + dRy**2 + dRz**2)

          enddo
          
          msd(i) = sum_msd / natoms(i)

          ! Use timestep of 32 au in ps
!          D(i) = msd(i) / (6.0 * timesteps(i)*0.00077404298)

!          D(1) = 0.0

          write(2,*) timesteps(i), msd(i) !, D(i)
!          write(2,*) i - 1, msd(i) !, D(i)

        enddo

 100    continue

        close(1)

        close(2)

 90     continue

        END
