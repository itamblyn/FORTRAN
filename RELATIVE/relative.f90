!       program center_of_mass.f90
!***********************************************************
!       Subtract off the COM from an XYZ file
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer n, m, i, j, k, rindex
        parameter (n=1000000)
        parameter (m=1024)
        integer natom, timestep
        real(8) x, y, z, refx, refy, refz, delta
        real(8) sx, sy, sz
        real(8) snapshot(m,3)
        real(8) ax, ay, az
        character(2) typat(m)
        character(128) fin

        delta = 0.01 ! so nothing sits on the edge

        ! Get xyz fname & lattice constants from the user
        write(*,*) 'Name of xyz file'
        read(*,*) fin

        ! index of relative atom
        write(*,*) 'What is the index of the reference atom'
        read(*,*) rindex

        ! Open input and output xyz files
        open(1,file=fin,status='old',ERR=90)
        open(2,file='com.xyz')
        open(3,file='com.dat')

        ! Loop over the file
        do i=1,n

          ! Read and write natom and timestep
          read(1,*,END=100) natom
          read(1,*) 
          write(2,*) natom
          write(2,*) i 

          sx = 0
          sy = 0
          sz = 0

          do j=1,natom

            ! Read in coordinates from xyz file

            read(1,*,END=100) typat(j), (snapshot(j,k), k=1,3) 

            sx = sx + snapshot(j,1)
            sy = sy + snapshot(j,2)
            sz = sz + snapshot(j,3)

          end do

          sx = sx/natom
          sy = sy/natom
          sz = sz/natom

          write(3,*) (sx**2 + sy**2 + sz**2)**0.5
        
          refx = snapshot(rindex,1)
          refy = snapshot(rindex,2)
          refz = snapshot(rindex,3)

          do j=1,natom

            snapshot(j,1) = snapshot(j,1) - refx + delta
            snapshot(j,2) = snapshot(j,2) - refy + delta
            snapshot(j,3) = snapshot(j,3) - refz + delta

            write(2,*) typat(j), (snapshot(j,k),k=1,3) 

          enddo

        enddo

 100    continue

        close(1)
        close(2)
        close(3)

 90     continue

        END
