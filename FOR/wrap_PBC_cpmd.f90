!       program wrap_PBC.f
!***********************************************************
!       Wrap trajectories from an xyz file
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer n, i, j
        parameter (n=1000000)
        integer natom, timestep
        real*8 x, y, z
        real*8 ax, ay, az
        character*2 typat
        character*128 fin

        ! Get xyz fname & lattice constants from the user
        write(*,*) 'Name of xyz file'
        read(*,*) fin
        write(*,*) 'Lattice constants (angst): ax, ay, az'
        read(*,*) ax, ay, az

        ! Open input and output xyz files
        open(1,file=fin,status='old',ERR=90)
        open(2,file='wrapped.xyz')

        ! Loop over the file
        do i=1,n

          ! Read and write natom and timestep
          read(1,*,END=100) natom
          read(1,*) 
          write(2,*) natom
          write(2,*) i 

          do j=1,natom

            ! Read in coordinates from xyz file
            read(1,*,END=100) typat, x, y, z

            ! Checks for coords > lattice const
 10         if (x.gt.ax) then
              x = x - ax
              goto 10
            endif
 11         if (x.lt.0.0) then
              x = x + ax
              goto 11
            endif
 20         if (y.gt.ay) then
              y = y - ay
              goto 20
            endif
 21         if (y.lt.0.0) then
              y = y + ay
              goto 21
            endif
 30         if (z.gt.az) then
              z = z - az
              goto 30
            endif
 31         if (z.lt.0.0) then
              z = z + az
              goto 31
            endif

            write(2,*) typat, x, y, z

          enddo

        enddo

 100    continue

        close(1)

 90     continue

        END
