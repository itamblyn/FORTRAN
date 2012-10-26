!       program stepsize_check.f
!*************************************************
!       Check the stepsize of atoms from an 
!       unwrapped xyz file and make sure they
!       are less than a/2
!*************************************************
        implicit none

        integer N, M
        integer i, j, k
        parameter (N=100000,M=10000)
        integer natom, cur_timestep
        character(2) cur_typat
        real(8) Rx(M),Ry(M),Rz(M)
        real(8) nextRx(M),nextRy(M),nextRz(M)
        real(8) dx,dy,dz,dx_max,dy_max,dz_max
        real(8) ax,ay,az
        character(128) fin
        logical safe

        safe = .true.
        dx_max = 0
        dy_max = 0
        dz_max = 0

        write(*,*) 'Name of unwrapped xyz file'
        read(*,*) fin
        write(*,*) 'Lattice constants a_x,a_y,a_z'
        read(*,*) ax,ay,az

        open(1,file=fin,status='old',ERR=90)

        read(1,*) natom
        read(1,*) cur_timestep
        do j=1,natom
          read(1,*) cur_typat, Rx(j), Ry(j), Rz(j)
        enddo

        do i=1,N

          read(1,*,END=100) natom
          read(1,*) cur_timestep

          do j=1,natom

            read(1,*) cur_typat, nextRx(j), nextRy(j), nextRz(j)

            dx = abs(nextRx(j) - Rx(j))
            dy = abs(nextRy(j) - Ry(j))
            dz = abs(nextRz(j) - Rz(j))

            if (dx.gt.dx_max) then
                dx_max = dx
            endif
            if (dy.gt.dy_max) then
                dy_max = dy
            endif
            if (dz.gt.dz_max) then
                dz_max = dz
            endif

            if (dx.ge.ax/2) then
              print *, 'WARNING: dx =',dx,'>= a_x/2'
              safe = .false.
            endif
            if (dy.ge.ay/2) then
              print *, 'WARNING: dy =',dy,'>= a_y/2'
              safe = .false.
            endif
            if (dz.ge.az/2) then
              print *, 'WARNING: dz =',dz,'>= a_z/2'
              safe = .false.
            endif

          enddo

        enddo

 100    continue

        if (safe) then
          print *,'No warnings were issued in this test'
          print *,'All stepsizes are therefore < a/2'
          print *,'dx_max=',dx_max,'< a_x/2 =',ax/2
          print *,'dy_max=',dy_max,'< a_y/2 =',ay/2
          print *,'dz_max=',dz_max,'< a_z/2 =',az/2
        endif

       close(1)

 90    continue
        
        END
