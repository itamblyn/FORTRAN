!       program bonding_analysis.f
!***********************************************************
!       Translation of Isaac's bonding_analysis.py
!       into fortran for speed
!       However, this assumes constant lattice constants
!***********************************************************
        implicit none
        
        ! n = max # of allowed timesteps
        integer n, i, j, k
        parameter (n=10000000)
        integer natom, tstep, Nsteps
        integer s, p, o, oo, molecule_number
        integer closest_to_p, closest_to_o, pbc_round(3)
        real*8 r(3), alat(3), pbv(4)
        real*4 POS_array(n,4)
        real*8 dx, dy, dz, dr, sanity
        real*8 input_value(3), distance, minimum_distance
        character*2 typat
        character*128 fxyz

        ! Get the xyz filename and alat from the user
        write(*,*) 'Name of xyz file'
        read(*,*) fxyz
        write(*,*) 'Lattice constants a,b,c (bohr)'
        read(*,*) alat(1), alat(2), alat(3)

        ! Open the input and output files
        open(1,file=fxyz,status='old',ERR=90)
        open(2,file='TRAJEC.fbn')
        open(3,file='diss.dat')

        Nsteps = 0

        do i=1,n

          read(1,*,END=100) natom
          read(1,*) tstep
          Nsteps = Nsteps + 1

          do j=1,natom

            read(1,*) typat, r(1), r(2), r(3), dummyvx, dummyvy, dummyvz,(FORCE_array(j+natom*(i-1),k),k=1,3)
            POS_array(j+natom*(i-1),1) = r(1)/0.5291772
            POS_array(j+natom*(i-1),2) = r(2)/0.5291772
            POS_array(j+natom*(i-1),3) = r(3)/0.5291772
            POS_array(j+natom*(i-1),4) = -10
          
          enddo

        enddo

 100    continue

        close(1)
        
        sanity = ( (alat(1)/2)**2 + (alat(2)/2)**2 + (alat(3)/2)**2 )**0.5

        s = 0

        do while (s.lt.Nsteps)

         molecule_number = 0

          ! counts over particles
          p = 0

          do while (p.lt.natom)

            if (POS_array(s*natom+p+1,4).eq.-10) then

              ! counts over other particles
              o = 0
              minimum_distance = sanity
              closest_to_p = p

              do while (o.lt.natom)

                ! Calculate differences
                dx = POS_array(s*natom+p+1,1) - POS_array(s*natom+o+1,1)
                dy = POS_array(s*natom+p+1,2) - POS_array(s*natom+o+1,2)
                dz = POS_array(s*natom+p+1,3) - POS_array(s*natom+o+1,3)

                ! Use the minimum image convention
                input_value(1) = dx/alat(1)
                input_value(2) = dy/alat(2)
                input_value(3) = dz/alat(3)
                pbc_round(1) = int(input_value(1))
                pbc_round(2) = int(input_value(2))
                pbc_round(3) = int(input_value(3))
                if (abs(input_value(1)-pbc_round(1)).ge.0.5) then
                  if (input_value(1).gt.0) pbc_round(1) = pbc_round(1) + 1
                  if (input_value(1).lt.0) pbc_round(1) = pbc_round(1) - 1
                endif
                if (abs(input_value(2)-pbc_round(2)).ge.0.5) then
                  if (input_value(2).gt.0) pbc_round(2) = pbc_round(2) + 1
                  if (input_value(2).lt.0) pbc_round(2) = pbc_round(2) - 1
                endif
                if (abs(input_value(3)-pbc_round(3)).ge.0.5) then
                  if (input_value(3).gt.0) pbc_round(3) = pbc_round(3) + 1
                  if (input_value(3).lt.0) pbc_round(3) = pbc_round(3) - 1
                endif
                
                dx = dx - alat(1)*pbc_round(1)
                dy = dy - alat(2)*pbc_round(2)
                dz = dz - alat(3)*pbc_round(3)

                distance = ( dx**2 + dy**2 + dz**2 )**0.5

                if (distance.gt.sanity) then
                  write(*,*) "Warning, problem with pbc"
                  goto 90
                endif

                if (distance.lt.minimum_distance) then
                  if (distance.ne.0.0) then
                    minimum_distance = distance
                    closest_to_p = o
                  endif
                endif

                o = o + 1

              enddo

              minimum_distance = sanity

              o = closest_to_p
              oo = 0

              if (POS_array(s*natom+o+1,4).eq.-10) then

                do while (oo.lt.natom)

                  ! Calculate differences
                  dx = POS_array(s*natom+o+1,1) - POS_array(s*natom+oo+1,1)
                  dy = POS_array(s*natom+o+1,2) - POS_array(s*natom+oo+1,2)
                  dz = POS_array(s*natom+o+1,3) - POS_array(s*natom+oo+1,3)

                  ! Use the minimum image convention
                  input_value(1) = dx/alat(1)
                  input_value(2) = dy/alat(2)
                  input_value(3) = dz/alat(3)
                  pbc_round(1) = int(input_value(1))
                  pbc_round(2) = int(input_value(2))
                  pbc_round(3) = int(input_value(3))
                  if (abs(input_value(1)-pbc_round(1)).ge.0.5) then
                    if (input_value(1).gt.0) pbc_round(1) = pbc_round(1) + 1
                    if (input_value(1).lt.0) pbc_round(1) = pbc_round(1) - 1
                  endif
                  if (abs(input_value(2)-pbc_round(2)).ge.0.5) then
                    if (input_value(2).gt.0) pbc_round(2) = pbc_round(2) + 1
                    if (input_value(2).lt.0) pbc_round(2) = pbc_round(2) - 1
                  endif
                  if (abs(input_value(3)-pbc_round(3)).ge.0.5) then
                    if (input_value(3).gt.0) pbc_round(3) = pbc_round(3) + 1
                    if (input_value(3).lt.0) pbc_round(3) = pbc_round(3) - 1
                  endif
              
                  dx = dx - alat(1)*pbc_round(1)
                  dy = dy - alat(2)*pbc_round(2)
                  dz = dz - alat(3)*pbc_round(3)

                  distance = ( dx**2 + dy**2 + dz**2 )**0.5

                  if (distance.gt.sanity) then
                    write(*,*) "Warning, problem with pbc"
                    goto 90
                  endif

                  if (distance.lt.minimum_distance) then
                    if (distance.ne.0.0) then
                      minimum_distance = distance
                      closest_to_o = oo
                    endif
                  endif

                  oo = oo + 1

                enddo

                if (closest_to_p.eq.o.and.closest_to_o.eq.p) then

                  molecule_number = molecule_number + 1

                  POS_array(s*natom+p+1,4) = o
                  POS_array(s*natom+o+1,4) = p

                else
                  POS_array(s*natom+p+1,4) = -1

                endif

              else
                POS_array(s*natom+p+1,4) = -1

              endif

            endif

            p = p + 1

          enddo
          
          write(3,*) 2.0*molecule_number/natom

          s = s + 1

        enddo

        close(3)

        ! Write the cbn file
        write(2,*) '# cbn file created with cbn_from_xyz.x'
        write(2,*) '# a =',alat(1)
        write(2,*) '# b =',alat(2)
        write(2,*) '# c =',alat(3)
        write(2,*) '# number_of_particles =',natom
        write(2,*) '# number_of_neighbors =',1
        write(2,*) '#'
        write(2,*) '#'
        write(2,*) '#'
        write(2,*) '# units = bohr'

        do i=1,Nsteps
          do j=1,natom
            write(2,*) typat,POS_array((i-1)*natom+j,1), &
                             POS_array((i-1)*natom+j,2), &
                             POS_array((i-1)*natom+j,3), &
                             int( POS_array((i-1)*natom+j,4) )
          enddo
        enddo

        close(2)

 90     continue

        END
