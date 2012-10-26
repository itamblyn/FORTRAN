      program VACF
      implicit none
      integer NAtoms,i,j,k,is,imax
      integer MaxWin,MaxAts,MaxPts
      real au2ps
      parameter (MaxWin=5000)
      parameter (MaxAts=1000)
      parameter (MaxPts=500000)
      parameter (au2ps=0.000024189D0)
      integer spacing
      integer np(MaxPts)
      integer windex(MaxWin),iw
      integer i1,i2,i3,i4,i5
      real*8 dt,norm,x,t,sum,maxtime,MaxT
      real*8 cdata(MaxPts)
      real*8 V(3),Carts(3)
      real*8 tau(3,5,MaxAts)
      real*8 X0(3,MaxAts,MaxWin)
      real*8 Xi(3,MaxAts,MaxWin)
      real*8 C(MaxWin)
      character*128 f1,f2,Line
      character*2 At,AtName(MaxAts)
      character*3 AType
      logical init

      call inp(dt,maxtime,spacing,AType,f1,f2)
      open(70,file=f1,status='old')
      open(71,file=f2,status='unknown')
      close(71,status='delete')
      open(71,file=f2,status='new')
      open(72,file='DiffCoef.data',status='unknown')
      close(72,status='delete')
      open(72,file='DiffCoef.data',status='new')


      do i=1,MaxPts
        np(i)=0
        cdata(i)=0.0D0
      enddo
      do i=1,MaxWin
        windex(i)=0
      enddo
      dt=dt*au2ps
      MaxT=int(min(maxtime/dt,MaxPts))
 
      is=0
      iw=0
      imax=0
      i5=0
 1000 continue
      init=.FALSE.
      read(70,*,end=2000) NAtoms
      if(NAtoms.gt.MaxAts) then
        write(*,*) "Number of atoms too large"
        stop
      endif
      is=is+1
      if(mod(is,1000).eq.0) write(*,9020)
      if(mod(is-3,spacing).eq.0.and.iw.lt.MaxWin) then
        init=.TRUE.
        iw=iw+1
      endif

      i5=mod(i5,5)+1
      i1=mod(i5+0,5)+1
      i2=mod(i5+1,5)+1
      i3=mod(i5+2,5)+1
      i4=mod(i5+3,5)+1

      read(70,*)
      do i=1,NAtoms
        read(70,9000) Line
        call gline(Line,At,Carts)
        AtName(i)=At
        tau(1,i5,i)=Carts(1)
        tau(2,i5,i)=Carts(2)
        tau(3,i5,i)=Carts(3)
 
        if(is.le.2) then
          V(1)=0.0D0
          V(2)=0.0D0          
          V(3)=0.0D0
        elseif(is.eq.3) then
          V(1)=(tau(1,i4,i)-tau(1,i3,i))/dt
          V(2)=(tau(2,i4,i)-tau(2,i3,i))/dt
          V(3)=(tau(3,i4,i)-tau(3,i3,i))/dt
        elseif(is.eq.4) then
          V(1)=0.5D0*(tau(1,i4,i)-tau(1,i2,i))/dt
          V(2)=0.5D0*(tau(2,i4,i)-tau(2,i2,i))/dt
          V(3)=0.5D0*(tau(3,i4,i)-tau(3,i2,i))/dt
        else
          V(1)=(-tau(1,i5,i)+8.0D0*tau(1,i4,i)-
     >           8.0D0*tau(1,i2,i)+tau(1,i1,i))/12.0D0/dt
          V(2)=(-tau(2,i5,i)+8.0D0*tau(2,i4,i)-
     >           8.0D0*tau(2,i2,i)+tau(2,i1,i))/12.0D0/dt
          V(3)=(-tau(3,i5,i)+8.0D0*tau(3,i4,i)-
     >           8.0D0*tau(3,i2,i)+tau(3,i1,i))/12.0D0/dt
        endif
c
c initialize a new window at t_0
c
        if(init) then
          X0(1,i,iw)=V(1)
          X0(2,i,iw)=V(2)
          X0(3,i,iw)=V(3)
        endif
c
c update current velocity in all active windows
c
        do j=1,iw
          Xi(1,i,j)=V(1)
          Xi(2,i,j)=V(2)
          Xi(3,i,j)=V(3)
        enddo
      enddo
c
c compute current average for each active window
c      
      do j=1,iw
        if(windex(j).lt.MaxT) then
          windex(j)=windex(j)+1
          C(j)=0.0D0
          if(AType.eq.'ALL') then 
            do i=1,NAtoms
              C(j)=C(j)+X0(1,i,j)*Xi(1,i,j)+
     >                  X0(2,i,j)*Xi(2,i,j)+
     >                  X0(3,i,j)*Xi(3,i,j)
            enddo
            C(j)=C(j)/float(NAtoms)
          else
            k=0
            do i=1,NAtoms
              if(AtName(i).eq.AType) then
                k=k+1
                C(j)=C(j)+X0(1,i,j)*Xi(1,i,j)+
     >                    X0(2,i,j)*Xi(2,i,j)+
     >                    X0(3,i,j)*Xi(3,i,j)
              endif
            enddo
            if(k.eq.0) then
              write(*,*) "No atoms found of type: ",AType
              stop
            endif
            C(j)=C(j)/float(k)
          endif
          i=windex(j)
          imax=max(imax,i)
          cdata(i)=cdata(i)+C(j)
          np(i)=np(i)+1
        endif
      enddo

      goto 1000
 2000 continue
c
c Average the data
c
      do i=1,imax
        if(np(i).gt.0) cdata(i)=cdata(i)/float(np(i))
      enddo
c
c Dump the results (raw and normalized) to the output file
c
      if(cdata(1).ne.0.0D0) then
        norm=1.0D0/cdata(1)
      else
        norm=1.0D0
      endif
      do i=1,imax
        if(np(i).gt.0) then
          x=cdata(i)
          t=float(i)*dt-dt
          write(71,9010) t,x,x*norm
        endif
      enddo
c
c Compute the diffusion coeficient:
c
c     D = 1/3 * int <v(o).v(t)> dt
c
      sum=0.5D0*(cdata(1)+cdata(imax))*dt
      do i=2,imax-1
        sum=sum+cdata(i)*dt
        write(72,*) i*dt,sum/3.0D0*1.0D-4
      enddo
      write(*,*)
      write(*,9030) sum/3.0D0*1.0D-4
c
c Some final stats
c
      write(*,9040) float(imax)*dt
      write(*,9050) iw
      close(70)
      close(71)
      close(72)
      stop
 9000 format(A128)
 9010 format(F10.5,2X,F15.8,2x,F15.8)
 9020 format(".",$)
 9030 format("Diffusion coef (cm^2/sec) = ",E15.5)
 9040 format("Maximum window length used (ps): ",f10.5)
 9050 format("Number of windows used: ",i8)
      end

      subroutine inp(dt,maxtime,spacing,AType,f1,f2)
      implicit none
      integer iout,input
      parameter (iout=6)
      parameter (input=5)
      integer spacing
      real*8 dt,maxtime
      character*3 AType,string
      character*128 f1,f2
      AType='ALL'
      write(iout,9000)
      read(input,*) f1
      write(iout,9010)
      read(input,*) f2
      write(iout,9020)
      read(input,*) dt
      write(iout,9030)
      read(input,*) spacing
      write(iout,9035)
      read(input,*) maxtime
      write(iout,9040)
      read(input,9050) string
      if(string.ne."") AType=string
      return
 9000 format("Enter name of XYZ file: ",$)
 9010 format("Enter name of output file: ",$)
 9020 format("Enter time between iterations (au): ",$)
 9030 format("Window spacing (iterations): ",$)
 9035 format("Window length (ps): ",$)
 9040 format("Atom type to track [ALL]: ",$)
 9050 format(a3)
      end

