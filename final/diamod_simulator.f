c-------------------------------------------------------------------------------
c     Size distribution Simulator for NEOMOD
c 
c     Based on the parameters listed in the parameter file (diamod.par), 
c     the code generates diameters, albedos and absolute magnitudes for model NEOs.  
c
c     The output is:       H (mag), diameter (km), albedo 
c     
c     The currently available diameter range is 0.001<D<100 km, but note that NEOMOD     
c     was calibrated for 15<H<28 mag, and H<15 and H>28 are extrapolations
c
c     Compile: ifort -o diamod_simulator diamod_simulator.f or gfortran -o diamod_simulator diamod_simulator.f  
c     Run: diamod_simulator < diamod.par > output_file.txt 
c     See diamod.par for input parameters
c------------------------------------------------------------------------------- 

      implicit NONE
      
      character*100 filename,modelname
      integer iseed,j
      integer nd,nseg,npoint,nneo,counter
      real*8 nneo_float,diam1,diam2,maxmod,ran3,diam,wmod
      real*8 mind,maxd,dref,nref,dpoint(100),mind2,maxd2,gamma(100)
      real*8 nlogp(100),y2(100),nlog0,yp1,ypn,nreal,ndiam1,ndiam2
      real*8 dd,ldleft,ldright,ld1,ld2,ldiam,ldref
c      real*8 nh,minh,maxh
      real*8 facd,sig1,sig2,logd,logd1,logd2,logd3,pV,phi,hmag

c     Read input parameters from neomod.par
      read(*,'(a)')filename ! NEO model file (only for bookkeeping)
      read(*,*)iseed        ! seed for the random number generator 
      read(*,*)nneo_float   ! number of NEOs to be generated, <=0 for real number
      read(*,*)diam1,diam2  ! desired diameter range
      if(diam1.ge.diam2) then
         write(*,*)'diam2 in diamod.par must be greater than diam1'
         write(*,*)'Exiting...'
         stop
      end if

c     Read the header of the model file
      open(1,file=filename,status='old')
      read(1,'(a)')modelname
      read(1,*)            !nh,minh,maxh ! absolute magnitude binning: # of H bins, Hmin,Hmax 
      read(1,*)dref,nref        ! reference magnitude, reference pop 
      ldref=log10(dref)
      read(1,*)nseg        ! number of spline segments for absolute magnitude distribution  
      npoint=nseg+1        ! number of boundaries between segments
      read(1,*)(dpoint(j),j=1,npoint)  ! segment boundaries
      mind2=10.d0**dpoint(1)   
      maxd2=10.d0**dpoint(npoint)
      if(diam1.lt.mind2) then
         write(*,*)
     &'The requested diameter range is outside of model domain'
         write(*,*)'diam1 in diamod.par must be >= than ',mind2
         write(*,*)'Exiting...'
         stop
      end if
      if(diam2.gt.maxd2) then
         write(*,*)
     &'The requested diameter range is outside of model domain'
         write(*,*)'diam2 in diamod.par must be <= than ',maxd2
         write(*,*)'Exiting...'
         stop
      end if
      read(1,*)(gamma(j),j=1,nseg) ! segment slopes
      close(1)

c     construct the size distribution
      nlog0=log10(nref)
      nlogp(6)=nlog0+(dpoint(6)-ldref)*gamma(5) ! dref falls into 5nd segment
      nlogp(5)=nlog0+(dpoint(5)-ldref)*gamma(5)
      nlogp(4)=nlogp(5)+(dpoint(4)-dpoint(5))*gamma(4)
      nlogp(3)=nlogp(4)+(dpoint(3)-dpoint(4))*gamma(3)
      nlogp(2)=nlogp(3)+(dpoint(2)-dpoint(3))*gamma(2)
      nlogp(1)=nlogp(2)+(dpoint(1)-dpoint(2))*gamma(1)
      do j=6,npoint
         nlogp(j)=nlogp(j-1)+(dpoint(j)-dpoint(j-1))*gamma(j-1)
      end do
      yp1=1.d30
      ypn=1.d30
      call spline(dpoint,nlogp,npoint,yp1,ypn,y2)

c     Rescale distribution to the desired number of model NEOs (nneo)
      ld1=log10(diam1)
      ld2=log10(diam2)      
      call splint(dpoint,nlogp,y2,npoint,ld1,ndiam1)
      call splint(dpoint,nlogp,y2,npoint,ld2,ndiam2)
      nreal=10.d0**ndiam1-10.d0**ndiam2
      if(nneo_float.ge.0.d0) then
         nneo=nneo_float
      else
         nneo=nreal
      end if

c     determine maxmod (differential!)
c      dd=0.0001d0
c      ldleft=ld1
c      ldright=log10(diam1+dd)
      dd=0.01
      ldleft=ld1
      ldright=ld1+dd
      call splint(dpoint,nlogp,y2,npoint,ldleft,ndiam1)
      call splint(dpoint,nlogp,y2,npoint,ldright,ndiam2)
      maxmod=10.d0**ndiam1-10.d0**ndiam2
      
c     Generate NEOs
      counter=0
 300  continue
c      diam=diam1+(diam2-diam1)*ran3(iseed)
c      ldiam=log10(diam)
c      ldleft=log10(diam-0.5d0*dd)   ! this may become negative if dd is too large
c      ldright=log10(diam+0.5d0*dd)
      ldiam=ld1+(ld2-ld1)*ran3(iseed)
      ldleft=ldiam-0.5d0*dd
      ldright=ldiam+0.5d0*dd
      call splint(dpoint,nlogp,y2,npoint,ldleft,ndiam1)
      call splint(dpoint,nlogp,y2,npoint,ldright,ndiam2)
      wmod=10.d0**ndiam1-10.d0**ndiam2
      if(ran3(iseed)*maxmod.lt.wmod) then   ! rejection method

c     generate albedo & magnitude

c     albedo distribution from Wright+16
c      facd=0.253d0
c      sig1=0.030d0
c      sig2=0.168d0

c     debiased albedo distribution from Model319, unchanging with D
c      facd=0.234d0
c      sig1=0.0287d0
c      sig2=0.170d0

c     debiased albedo distribution from Model330 (1-3 km)
c      facd=0.306d0
c      sig1=0.0251d0
c      sig2=0.162d0

c     size-dependent albedo model
      logd=ldiam
      logd1=log10(sqrt(0.1d0*0.3d0))
      logd2=log10(sqrt(0.3d0*1.0d0)) 
      logd3=log10(sqrt(1.0d0*3.0d0))              
      if(logd.lt.logd1) then    ! Model332
         facd=0.183d0
         sig1=0.0566d0
         sig2=0.182d0   
      end if
      if(logd.ge.logd1.and.logd.lt.logd2) then ! Models 332->331
         facd=0.183d0+(logd-logd1)*(0.212d0-0.183d0)/(logd2-logd1)
         sig1=0.0566d0+(logd-logd1)*(0.0367d0-0.0566d0)/(logd2-logd1)
         sig2=0.182d0+(logd-logd1)*(0.182d0-0.182d0)/(logd2-logd1)  
      end if
      if(logd.ge.logd2.and.logd.lt.logd3) then ! Models 331->330
         facd=0.212d0+(logd-logd2)*(0.306d0-0.212d0)/(logd3-logd2)
         sig1=0.0367d0+(logd-logd2)*(0.0251d0-0.0367d0)/(logd3-logd2)
         sig2=0.182d0+(logd-logd2)*(0.162d0-0.182d0)/(logd3-logd2)  
      end if
      if(logd.ge.logd3) then    ! Model330
         facd=0.306d0
         sig1=0.0251d0
         sig2=0.162d0              
      end if

c     rejection method to select albedo
 400  continue
      pV=ran3(iseed)  ! random albedo between 0 and 1
      phi=facd*(pV*
     &     exp(-pV**2/(2.d0*sig1*sig1))/(sig1*sig1))
      phi=phi+(1.d0-facd)*(pV*
     &     exp(-pV**2/(2.d0*sig2*sig2))/(sig2*sig2))
      if(ran3(iseed)*10.d0.gt.phi) then 
         goto 400  ! reject this albedo value
      end if
      hmag=-5.d0*log10(10.d0**ldiam*sqrt(pV)/1329.d0)
      write(*,8848)hmag,10.d0**ldiam,pV
 8848 format(f7.3,f8.4,f7.4)
      counter=counter+1
      end if
      if(counter.ge.nneo) stop  ! cycle until we generate the desired number of NEOs 
      goto 300

      end

c---------------------------------------------------------------------------

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt.0.99d30) then
         y2(1)=0.d0
         u(1)=0.d0
      else
         y2(1)=-0.5d0
         u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.d0)/p
         u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt.0.99d30) then
         qn=0.d0
         un=0.d0
      else
         qn=0.5d0
         un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      enddo 
      return
      END

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
 1    if(khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
      h=xa(khi)-xa(klo)
c      if (h.eq.0.) pause ’bad xa input in splint’
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
      END

      function ran3(idum)
      implicit real*8(a-h,o-z)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.d-9)
      dimension ma(55)
      data iff /0/
      save
c --------------------------------------------------------------------
      if((idum.ge.0).and.(iff.ne.0)) goto 20
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
20    continue
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      return
      end
