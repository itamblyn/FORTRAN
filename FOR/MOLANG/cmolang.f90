program MOLANG 

! Note that this program does not like it when a molecules center is completely in line
! with anothers molecular axis. For any type of MD this is almost guarunteed not to
! happen. For the case of solids, a tiny bit of thermal noise should fix the
! issue

use vectors

implicit none

        
integer maxSteps, maxAtoms, i, j, k
integer maxBinsR, maxBinsT, counter
integer nip, nio
parameter (maxSteps=100000,maxAtoms=1024)
parameter (maxBinsR=1000,maxBinsT=180)
integer natom, nsteps, B(maxAtoms), NbinsT, Rbin, bearing_bin, theta_bin, phi_bin, Nbinsphi
integer bonding
real(8) Rstep, Tstep, all_theta, ve1, ve2, ve, Pi, phi_step
real(8) r, bearing, X(maxAtoms), Y(maxAtoms), Z(maxAtoms)
real(8) alat(3), CX(maxAtoms), CY(maxAtoms), CZ(maxAtoms)
real(8) allSum, nearSum
real(8) dx,dy,dz,px,py,pz,ox,oy,oz,oox,ooy,ooz,ooox,oooy,oooz
real(8) dop, ntp, ndp, nto, ndo
real(8) do_oo, doo_ooo, do_oocom, doop, dooo, ndrp, ndro
real(8) o_phi,o_theta 
real(8) avec(3),bvec(3),cvec(3),xvec(3),xprime(3),bprime(3)
real(4) angleHist(maxBinsT,2*maxBinsT,maxBinsT)
character(8) typat, d1, d2, d3, d4, d5
character(128) fin

write(*,*) 'Name of cbn file:'
read(*,*) fin

NbinsT = 9     
Nbinsphi = 18 

! Open cbn file and read in header info
open(1,file=fin,status='old',ERR=90)
read(1,*) 
read(1,*) d1, d2, d3, alat(1)
read(1,*) d1, d2, d3, alat(2)
read(1,*) d1, d2, d3, alat(3)
read(1,*) d1, d2, d3, natom
do i=1,5
  read(1,*)
end do

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
Pi = 3.14159265
Tstep = Pi/real(NbinsT)
phi_step = Pi/real(Nbinsphi)    ! this variable is for the ooo particle

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


          bonding = min(B(k) - 1,0) + 2

          if (bonding.ne.1) then  ! check to see if the nearest particle is bonded

            ooox = CX(B(k))
            oooy = CY(B(k))
            oooz = CZ(B(k))

            doo_ooo = ((ooox - oox)**2 + (oooy-ooy)**2 + (oooz-ooz)**2)**0.5

            if (doo_ooo.lt.2) then   ! this is a hackish way of ignoring
                                     ! molecules that go across a cell boundary
                                     ! to make life easier
         
              avec(1) = (oox+ooox)*0.5
              avec(2) = (ooy+oooy)*0.5  ! this vector points to the COM of the
              avec(3) = (ooz+oooz)*0.5  ! molecule. I reuse the vector later on
                                        ! in the code for something completely different (gets written over)
              do_oocom = length(avec)
  
            endif

          endif

      ! this looks for the closest MOLECULAR CENTER to the origin

          if (do_oocom.lt.ndo) then
            nio  = k
            ndo  = do_oocom
!            nto  = all_theta
!            ndro = do_oo
          endif

        endif

      enddo

!      R = ndro
!      bearing = nto
      bonding = min(B(nio) - 1,0) + 2

      ! now we look at what the molecule is doing (if it is a molecule)...

      if (bonding.ne.1) then  ! check to see if the nearest particle is bonded

        oox = CX(nio)
        ooy = CY(nio)
        ooz = CZ(nio)

        ooox = CX(B(nio))
        oooy = CY(B(nio))
        oooz = CZ(B(nio))

        avec(1) = ox - (px+ox)*0.5
        avec(2) = oy - (py+oy)*0.5
        avec(3) = oz - (pz+oz)*0.5

        ! px + ox, py + oy, pz + oz should return 0, 0, 0 <- it does :)
        
        bvec(1) = (oox+ooox)*0.5 - (px+ox)*0.5
        bvec(2) = (ooy+oooy)*0.5 - (py+oy)*0.5
        bvec(3) = (ooz+oooz)*0.5 - (pz+oz)*0.5

        bearing = acos(dot_product(bvec,avec)/(length(bvec)*length(avec)))  ! I think this deal with the length?
        bearing_bin = int(bearing/Tstep + 1)   ! this is the bin for the CENTER of the ooMolecule

        cvec = cross(avec,bvec)

        xvec(1) = ooox - (oox+ooox)*0.5
        xvec(2) = oooy - (ooy+oooy)*0.5
        xvec(3) = oooz - (ooz+oooz)*0.5

!        print *, "avec: ", avec
!        print *, "bvec: ", bvec
!        print *, "xvec: ", xvec

        o_phi = acos(dot_product(xvec,cvec)/(length(xvec)*length(cvec)))*(180.0/Pi)
        o_phi = 90.0 - o_phi ! this number gives the in/out of planeness...

!        print *, "xvec: ", xvec, "cvec: ", cvec

        xprime = xvec - dot_product(xvec,(cvec/length(cvec)))*(cvec/length(cvec))

!        print *, "xprime: ", xprime

        bprime = bvec - dot_product(bvec,(avec/length(avec)))*(avec/length(avec))

!        print *, "bprime: ", bprime

        o_theta = atan2(dot_product(xprime,(bprime/length(bprime))),dot_product(xprime,(avec/length(avec))))*(180.0/Pi)

        if (o_theta < 0) then
          o_theta = o_theta + 360.0
        endif

        print *, 'pitch: ', NInt(o_phi),'yaw: ', NInt(o_theta), 'bearing: ', NInt(bearing*(180.0/Pi))


        theta_bin = int((o_theta)/(phi_step*(180.0/Pi)) + 1)
        phi_bin =   int((o_phi   +  90)/(phi_step*(180.0/Pi)) + 1)

!        print *, theta_bin, phi_bin

        angleHist(phi_bin, theta_bin, bearing_bin) = angleHist(phi_bin, theta_bin, bearing_bin) + 1

      endif

    endif

  enddo
    
enddo

100    continue

open(3,file='angle.000.hist')
open(4,file='angle.030.hist')
open(5,file='angle.045.hist')
open(6,file='angle.060.hist')
open(7,file='angle.090.hist')
open(8,file='angle.180.hist')

do i=1,Nbinsphi
  do j=1,Nbinsphi*2
    write(3,*) (i-0.5)*(phi_step*(180.0/Pi))-90, (j-0.5)*(phi_step*(180.0/Pi)), angleHist(i,j,1)
    write(4,*) (i-0.5)*(phi_step*(180.0/Pi))-90, (j-0.5)*(phi_step*(180.0/Pi)), angleHist(i,j,Nbinsphi/6)
    write(5,*) (i-0.5)*(phi_step*(180.0/Pi))-90, (j-0.5)*(phi_step*(180.0/Pi)), angleHist(i,j,Nbinsphi/4)
    write(6,*) (i-0.5)*(phi_step*(180.0/Pi))-90, (j-0.5)*(phi_step*(180.0/Pi)), angleHist(i,j,Nbinsphi/3)
    write(7,*) (i-0.5)*(phi_step*(180.0/Pi))-90, (j-0.5)*(phi_step*(180.0/Pi)), angleHist(i,j,Nbinsphi/2)
    write(8,*) (i-0.5)*(phi_step*(180.0/Pi))-90, (j-0.5)*(phi_step*(180.0/Pi)), angleHist(i,j,Nbinsphi)
  enddo
enddo

close(3)
close(4)
close(5)
close(6)
close(7)
close(8)

allSum = sum(angleHist)

open(3,file='com.dat')

do k=1,NbinsT 
  nearSum = 0.0
  do i=1,Nbinsphi
    do j=1,Nbinsphi*2
      nearSum = nearSum + angleHist(i,j,k)
    enddo
  enddo

  write(3,*) (k - 0.5)*(Tstep*(180.0/Pi)), nearSum/allSum

enddo

close(3)

90    continue

end program MOLANG 
