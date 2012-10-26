      program RDF
      implicit none
      integer MaxAts,MaxBin
      real*8 Pi
      parameter(MaxAts=50000)
      parameter(MaxBin=1000)
      parameter(Pi=3.14159265359D0)
      integer i,j,NAtoms,ind,n1,n2
      integer tot,icell
      real*8 g(MaxBin),n(MaxBin),gt(MaxBin),nt(MaxBin)
      real*8 x,y,z,frame,omega,c1,c2
      real*8 r,rmax,rbin,shell,rint
      real*8 a,b,c
      real*8 cx1,cy1,cz1
      real*8 cx2,cy2,cz2
      real*8 Carts(3)
      real*8 Coords(3,MaxAts)
      character*2 At,AtName(MaxAts),a1,a2
      character*128 Line,f1,f2
 
      write(*,*) "Name of input file"
      read(*,*) f1
      write(*,*) "Name of output file"
      read(*,*) f2
      write(*,*) "Atom type 1"
      read(*,*) a1
      write(*,*) "Atom type 2"
      read(*,*) a2
      write(*,*) "Max r"
      read(*,*) rmax
      write(*,*) "Bin size"
      read(*,*) rbin
      write(*,*) "Cell size in xyz file (0=no, 1=yes)"
      read(*,*) icell
      if(icell.eq.0) then
        write(*,*) "Enter cell size a,b,c (angstrom)"
        read(*,*) a,b,c
      endif
       
      open(70,file=f1,status='old')
      open(71,file=f2,status='unknown')
      close(71,status='delete')
      open(71,file=f2,status='new')

      do i=1,MaxBin
        g(i)=0.0D0
        gt(i)=0.0D0
        nt(i)=0.0D0
      enddo
      if(int(rmax/rbin)+1.gt.MaxBin) then
        write(*,*) "Max r too large"
        stop
      endif
      frame=0.0D0
      omega=a*b*c

 1000 continue
      read(70,*,end=2000) NAtoms
      frame=frame+1.0D0
      if(icell.eq.0) then
        read(70,*) 
      else
        read(70,*) a,b,c
      endif
      omega=a*b*c
      if(NAtoms.gt.MaxAts) then
         write(*,*) "Number of atoms too large ",NAtoms,MaxAts
         stop
      endif
      do i=1,NAtoms
        read(70,9000) Line
        call gline(Line,At,Carts)
        AtName(i)=At
        Coords(1,i)=Carts(1)
        Coords(2,i)=Carts(2)
        Coords(3,i)=Carts(3)
      enddo

      tot=0
      n1=0
      do i=1,NAtoms
        if(AtName(i).eq.a1) then
          n1=n1+1
          n2=0
          do j=1,NAtoms
            if(AtName(j).eq.a2.and.j.ne.i) then
              n2=n2+1
              x=Coords(1,i)-Coords(1,j)
              y=Coords(2,i)-Coords(2,j)
              z=Coords(3,i)-Coords(3,j)
              x=x-a*anint(x/a)
              y=y-b*anint(y/b)
              z=z-c*anint(z/c)
              r=sqrt(x*x+y*y+z*z)
              if(r.le.rmax+rbin) then
                ind = int(r/rbin)
                g(ind)=g(ind)+1.0D0
                tot=tot+1
              endif
            endif
          enddo
        endif
      enddo

      rint=0.0D0
      do i=1,int(rmax/rbin)+1
        r = float(i)*rbin
        cx1=0.0D0
        cx2=0.0D0
        cy1=0.0D0
        cy2=0.0D0
        cz1=0.0D0
        cz2=0.0D0

        if(r.gt.a*0.5D0)      cx1=Pi/12.0D0*(a**3-12.0D0*a*r**2+
     >                            16.0D0*r**3)
        if(r+rbin.gt.a*0.5D0) cx2=Pi/12.0D0*(a**3-12.0D0*a*(r+rbin)**2+
     >                           16.0D0*(r+rbin)**3)

        if(r.gt.b*0.5D0)      cy1=Pi/12.0D0*(b**3-12.0D0*b*r**2+
     >                            16.0D0*r**3)
        if(r+rbin.gt.b*0.5D0) cy2=Pi/12.0D0*(b**3-12.0D0*b*(r+rbin)**2+
     >                            16.0D0*(r+rbin)**3)

        if(r.gt.c*0.5D0)      cz1=Pi/12.0D0*(c**3-12.0D0*c*r**2+
     >                            16.0D0*r**3)
        if(r+rbin.gt.c*0.5D0) cz2=Pi/12.0D0*(c**3-12.0D0*c*(r+rbin)**2+
     >                            16.0D0*(r+rbin)**3)

        shell = 4.0D0*Pi/3.0D0*((r+rbin)**3 - r**3) - 
     >          (cx2-cx1)-(cy2-cy1)-(cz2-cz1)

        g(i)=(omega/float(n1*n2))*(g(i)/shell)

        rint=rint+4.0D0*Pi*rbin*g(i)*(r+rbin/2.0D0)**2
        n(i)=rint*float(n2)/omega
        if(tot.gt.0) then
          gt(i)=gt(i)+g(i)
          nt(i)=nt(i)+n(i)
        endif
        g(i)=0.0D0
        n(i)=0.0D0
      enddo

      goto 1000
 2000 continue
 
      do i=1,int(rmax/rbin)
        r = float(i)*rbin+rbin*0.5D0
        write(71,9020) r,gt(i)/frame,nt(i)/frame
      enddo

      close(70)
      close(71)
      stop
 9000 format(A128)
 9010 format(A2,1x,f10.5,1x,f10.5,1x,f10.5)
 9020 format(f10.5,2x,f10.5,2x,f10.5)
      end
