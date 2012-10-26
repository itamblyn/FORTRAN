!       program xyz_util.f
!***********************************************************
!       Perform operations on an xyz file (TRAJEC.xyz)
!       1) replace negative coords with PBC displacements
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer n, m, i, j, k, count
        parameter (n=1000000,m=10000)
        integer natom, timestep
        real(4) r(m,3), threshold
        real(4) dx, dy, dz, dr
        character(2) typat
        character(128) fin

        ! Get the xyz filename from the user
        write(*,*) 'Name of xyz file'
        read(*,*) fin
        write(*,*) 'What dr threshold would you like to use?'
        read(*,*) threshold

        open(1,file=fin,status='old',ERR=90)

        count = 0

        ! Loop over the rest of the file
        do i=1,n

          ! Read in the value of natom and timestep
          read(1,*,END=100) natom
          read(1,*) timestep

          do j=1,natom
            ! Read in coordinates from xyz file
            read(1,*,END=100) typat, r(j,1), r(j,2), r(j,3)
          enddo

          do j=1,natom
            do k=1,natom
              if (j.lt.k) then
                dx = (r(j,1)-r(k,1))
                dy = (r(j,2)-r(k,2))
                dz = (r(j,3)-r(k,3))
                dr = ( dx**2 + dy**2 + dz**2 )**0.5
                if (dr.lt.threshold) then
                  write(*,*) 'WARNING: dr = ',dr,'j=',j,'k=',k
                  count = count + 1
                endif
              endif
            enddo
          enddo

        enddo

 100    continue

        write(*,*) 'Number of configurations checked:',i-1
        write(*,*) "Number of dr's found <",threshold,'was:',count

        close(1)

 90     continue

        END
