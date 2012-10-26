program SDF 
!
!***********************************************************
!       Calcuates the 'Spatial Distribution Function'
!       Requires a cbn file.
!***********************************************************
!
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        integer maxBinsR, maxBinsT, counter
        integer nip, nio
        parameter (maxSteps=100000,maxAtoms=1024)
        parameter (maxBinsR=1000,maxBinsT=180)
        integer natom, nsteps, B(maxAtoms), NbinsR, NbinsT, Rbin, Tbin
        integer bonding
        real(8) Rstep, Tstep, all_theta, ve1, ve2, ve, Pi
        real(8) allHist(maxBinsR,maxBinsT,2), nearHist(maxBinsR,maxBinsT,2)
        real(8) Rmax, r, theta, X(maxAtoms), Y(maxAtoms), Z(maxAtoms)
        real(8) alat(3), CX(maxAtoms), CY(maxAtoms), CZ(maxAtoms)
        real(8) allSum, nearSum
        real(8) dx, dy, dz, px, py, pz, ox, oy, oz, oox, ooy, ooz
        real(8) dop, ntp, ndp, nto, ndo
        real(8) do_oo, doop, dooo, ndrp, ndro
        character(8) typat, d1, d2, d3, d4, d5
        character(128) fin

        write(*,*) 'Name of cbn file:'
        read(*,*) fin
        write(*,*) 'Number of R bins:'
        read(*,*) NbinsR
        write(*,*) 'Number of theta bins:'
        read(*,*) NbinsT

        ! Open cbn file and read in header info
        open(1,file=fin,status='old',ERR=90)

        read(1,*) 
        read(1,*) d1, d2, d3, alat(1)
        read(1,*) d1, d2, d3, alat(2)
        read(1,*) d1, d2, d3, alat(3)
        read(1,*) d1, d2, d3, natom
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 

        counter = 0

        do i=1,maxSteps*maxAtoms
          read(1,*,END=80)
          counter = counter + 1
        end do

 80     continue

        nsteps = counter/natom

        print *, nsteps

        rewind(1)

        do i=1,10
          read(1,*)   ! jumps over header
        end do

        ! Calculate variable step sizes
        Rmax = 0.5*((alat(1) + alat(2) + alat(3))/3.0)
        Rstep = Rmax/real(NbinsR)
        Pi = 3.14159265
        Tstep = Pi/real(NbinsT)

        ! Open all histogram output files and write headers
        open(5,file='all_ions.hist')
        open(6,file='all_molecules.hist')
        open(7,file='all_ions_and_molecules.hist')
        open(8,file='nearest_ions.hist')
        open(9,file='nearest_molecules.hist')
        open(10,file='nearest_ions_and_molecules.hist')
        write(5,*)  '# r theta occ Rstep=',Rstep,'Tstep=',Tstep*180.0/Pi
        write(6,*)  '# r theta occ Rstep=',Rstep,'Tstep=',Tstep*180.0/Pi
        write(7,*)  '# r theta occ Rstep=',Rstep,'Tstep=',Tstep*180.0/Pi
        write(8,*)  '# r theta occ Rstep=',Rstep,'Tstep=',Tstep*180.0/Pi
        write(9,*)  '# r theta occ Rstep=',Rstep,'Tstep=',Tstep*180.0/Pi
        write(10,*) '# r theta occ Rstep=',Rstep,'Tstep=',Tstep*180.0/Pi

        ! BEGIN ANALYSIS
        do i=1,maxSteps

          ! Read in current snapshot data from cbn file
          do j=1,natom
            read(1,*,END=100) typat, X(j), Y(j), Z(j), B(j)
            B(j) = B(j) + 1 ! switches to fortran array numbers
          enddo

          do j=1,natom

            ! If particle is bonded
            if (B(j).ne.0) then 

              px = X(j)
              py = Y(j)
              pz = Z(j)

              do k=1,natom
                dx = X(k) - px
                dy = Y(k) - py
                dz = Z(k) - pz
                X(k) = dx - (int(dx/alat(1)+maxSteps+0.5)-maxSteps)*alat(1)  
                Y(k) = dy - (int(dy/alat(2)+maxSteps+0.5)-maxSteps)*alat(2)
                Z(k) = dz - (int(dz/alat(3)+maxSteps+0.5)-maxSteps)*alat(3) 
              enddo

              ! Coords of atom that j is bonded to
              ox = X( B(j) )
              oy = Y( B(j) )
              oz = Z( B(j) )

              ! Do some center shifting of coords
              do k=1,natom
                CX(k) = X(k)  - ox/2.0
                CY(k) = Y(k)  - oy/2.0
                CZ(k) = Z(k)  - oz/2.0
                CX(k) = CX(k) - (int(CX(k)/alat(1)+maxSteps+0.5)-maxSteps)*alat(1)  
                CY(k) = CY(k) - (int(CY(k)/alat(2)+maxSteps+0.5)-maxSteps)*alat(2)  
                CZ(k) = CZ(k) - (int(CZ(k)/alat(3)+maxSteps+0.5)-maxSteps)*alat(3)  
              enddo

              px = CX(j)
              py = CY(j)
              pz = CZ(j)
              ox = CX( B(j) )
              oy = CY( B(j) )
              oz = CZ( B(j) )

              ! Acronyms of variables from Isaac's original code
              ! initialize variables such that they will throw error
              ! if left unaltered
              dop = (px**2 + py**2 + pz**2)**0.5
              ntp = -1
              ndp = 1000000
              nip = maxSteps*natom
              nto = -1
              ndo = 1000000
              nio = maxSteps*natom

              ! Meat of the code
              do k=1,natom

                if (k.ne.j.and.k.ne.B(j)) then

                  oox = CX(k)
                  ooy = CY(k)
                  ooz = CZ(k)

                  do_oo = (oox**2 + ooy**2 + ooz**2)**0.5
                  all_theta = acos( (px*oox + py*ooy + pz*ooz)/(dop*do_oo) )

                  doop = ( (oox-px)**2 + (ooy-py)**2 + (ooz-pz)**2 )**0.5
                  dooo = ( (oox-ox)**2 + (ooy-oy)**2 + (ooz-oz)**2 )**0.5

!                 this looks for the closest particle to either of the
!                 two reference particles

                  if (min(doop,dooo).lt.ndp) then
                    nip  = k
                    ndp  = min(doop,dooo)
                    ntp  = all_theta
                    ndrp = do_oo
                  endif

!                 this looks for the closest partile to the origin

                  if (do_oo.lt.ndo) then
                    nio  = k
                    ndo  = do_oo
                    nto  = all_theta
                    ndro = do_oo
                  endif

                  R = do_oo
                  theta = all_theta
                  bonding = min(B(k) - 1,0) + 2 ! zero means ion, one means bonded

                  Rbin = int(R/Rstep + 1) ! is this correct???
                  Tbin = int(theta/Tstep + 1)  

                  if (Rbin.lt.NbinsR.and.theta.ne.0.and.theta.ne.Pi) then
                    allHist(Rbin,Tbin,bonding) = allHist(Rbin,Tbin,bonding) + 1.0
                  endif

                endif
              enddo

!             if the variable ends in o, it's distance to origin
!             if the variable ends in p, it's distance to nearest
!             particle

              R = ndro
              theta = nto
              bonding = min(B(nio) - 1,0) + 2

              Rbin = int(R/Rstep + 1)  ! is this correct????? 
              Tbin = int(theta/Tstep + 1) 

              if (Rbin.lt.NbinsR.and.theta.ne.0.and.theta.ne.Pi) then
                nearHist(Rbin,Tbin,bonding) = nearHist(Rbin,Tbin,bonding) + 1.0
              endif

            endif

          enddo
       
        enddo

 100    continue

        do Rbin=1,NbinsR
          do Tbin=1,NbinsT
            ve1 = ( (Rbin*Rstep)**3 - ((Rbin-1)*Rstep)**3)
            ve2 = cos((Tbin-1)*Tstep) - cos(Tbin*Tstep)
            ve  = (2.0*Pi/3.0) * ve1 * ve2
            allHist(Rbin,Tbin,1)  = allHist(Rbin,Tbin,1)*(1.0/ve)
            allHist(Rbin,Tbin,2)  = allHist(Rbin,Tbin,2)*(1.0/ve)
            nearHist(Rbin,Tbin,1) = nearHist(Rbin,Tbin,1)*(1.0/ve)
            nearHist(Rbin,Tbin,2) = nearHist(Rbin,Tbin,2)*(1.0/ve)
          enddo
        enddo

        ! Calculate sums of histograms for normalization
        allSum = 0.0
        nearSum = 0.0
        do i=1,NbinsR
          do j=1,NbinsT
            allSum  = allSum  + allHist(i,j,1)
            allSum  = allSum  + allHist(i,j,2)
            nearSum = nearSum + nearHist(i,j,1)
            nearSum = nearSum + nearHist(i,j,2)
          enddo
        enddo

        ! Normalize the histograms
        do i=1,NbinsR
          do j=1,NbinsT
            allHist(i,j,1)  = allHist(i,j,1)  / allSum
            allHist(i,j,2)  = allHist(i,j,2)  / allSum
            nearHist(i,j,1) = nearHist(i,j,1) / nearSum
            nearHist(i,j,2) = nearHist(i,j,2) / nearSum
          enddo
        enddo

        ! Write calculations to output files
        do i=1,NbinsR
          do j=1,NbinsT
            if (NbinsT.eq.1) then
              write(5,*)  (i-0.5)*Rstep, (j-1)*Tstep, allHist(i,j,1)
              write(6,*)  (i-0.5)*Rstep, (j-1)*Tstep, allHist(i,j,2)
              write(7,*)  (i-0.5)*Rstep, (j-1)*Tstep, allHist(i,j,1)+allHist(i,j,2)
              write(8,*)  (i-0.5)*Rstep, (j-1)*Tstep, nearHist(i,j,1)
              write(9,*)  (i-0.5)*Rstep, (j-1)*Tstep, nearHist(i,j,2)
              write(10,*) (i-0.5)*Rstep, (j-1)*Tstep, nearHist(i,j,1)+nearHist(i,j,2)
            endif
            if (NbinsT.ne.1) then
              write(5,*)  (i-0.5)*Rstep, (j-0.5)*Tstep, allHist(i,j,1)
              write(6,*)  (i-0.5)*Rstep, (j-0.5)*Tstep, allHist(i,j,2)
              write(7,*)  (i-0.5)*Rstep, (j-0.5)*Tstep, allHist(i,j,1)+allHist(i,j,2)
              write(8,*)  (i-0.5)*Rstep, (j-0.5)*Tstep, nearHist(i,j,1)
              write(9,*)  (i-0.5)*Rstep, (j-0.5)*Tstep, nearHist(i,j,2)
              write(10,*) (i-0.5)*Rstep, (j-0.5)*Tstep, nearHist(i,j,1)+nearHist(i,j,2)
            endif
          enddo
        enddo

        close(1)
!        close(2)
!        close(3)
!        close(4)
        close(5)
        close(6)
        close(7)
        close(8)
        close(9)
        close(10)

 90     continue

        END
