C
C      PROGRAM Ngamma.f (WRITTEN BY D. Dutta and P. Mohanmurthy)
C
C      This program calculates the total number (and difference in number above 
c      and below the orbital plane) of Synchrotron photons emitted by 
c      longitudinally polarized electrons over a horizontal angular 
c      range of dTheta and verticle angular range of +/-alpha = +/-gamma*Psi
c      (where gamma is the Lorentz boost, i.e. +/-alpha approx = +/-1 ), when 
c      traversing a 3 pole wriggler magnet with a field stength of Bwg tesla 
c      and a pole length of Lwg. 
C
c       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT real*8(A-H,O-Z)
       external fn01
       external fn02
       external gn01
       external gn02
       real*8 xheeng(48),xhemu(48)

       parameter (xMe=0.510998902) !electron mass in MeV/c^2
       parameter (GeV2MeV=1000.0) 
       parameter (hbarc=197.3269602) ! MeV*fm
       parameter (xmuB=5.788381749E-11) ! Bohr magnetron MeV/T
       parameter (c=299792458) ! m/s
       parameter (pi=3.141592654)
       parameter (qe=1.602176462E-19) ! coulomb
       parameter (n1=10) ! number of times the integration alg is compounded
       parameter (n2=10)

       Common gamma,y

c       write(6,*)'Enter Ebeam (GeV) and current (micro A)'
c       read(5,*)Ebeam,xIe
c       write(6,*)'Wriggler B-field (T) and Pole length (m)' 
c       read(5,*)Bwg,xLwg
       open(unit=10,file="spinlight_gydep4.dat",status="unknown")
       open(unit=11,file="xenon.dat",status="old")
       do i =1,48
         read(11,*)xheeng(i),xhemu(i)
       enddo
       close(11)
       ymin = 0.01 ! min fractional photon freq (W/Wc)
       ymax=0.02 !initialize
c       Ebeam=4.0+(i-1)*1.0 ! e-beam enengy GeV
       Ebeam=11.0 ! e-beam enengy GeV
       xIe=100.0  ! e-beam current micro A
       Bwg=4.6158    ! B-field in T
c       Bwg=1.0+(i-1)*1.0
c       xLwg=0.066  ! pole length in m
c       write(6,*)'Ebeam =',Ebeam, 'GeV'
       gamma = Ebeam*GeV2MeV/xMe  ! Lorentz boost = E/(Me*c^2)
       R_bend = gamma*hbarc*1.0E-15/(2.*xmuB*Bwg) !bending radius in m
       Omega_o = c/R_bend  ! betatron freq.
       Omega_c = 1.5*gamma**3*Omega_o ! central photon frequency
       E_cent=(Omega_c*hbarc*1.0E-15/c)*1000. ! central photon energy in keV
       xlambda_c=2.*pi*c/Omega_c
       do i=1,1000
       ymax = 0.02+(i-1)*0.01 ! max fractional photon freq (W/Wc)
       y_cent= (ymin+ymax)/2.
       E_min = (ymin*Omega_c*hbarc*1.0E-15/c)*1000. ! min photon energy in keV 
       E_max = (ymax*Omega_c*hbarc*1.0E-15/c)*1000. ! max photon energy in keV
       E_cent =(y_cent*Omega_c*hbarc*1.0E-15/c) ! photon energy in MeV 
       ij=2
       ik=1
       ift=1
       do ii=1,48
        if(E_cent.lt.xheeng(ii).and.ift.eq.1)then
         ij=ii
         ik=ii-1
         ift=ift+1
        endif
       enddo
       xrayabs=xhemu(ik)+
     > (((E_cent-xheeng(ik))/(xheeng(ij)-xheeng(ik)))*
     >                          (xhemu(ij)-xhemu(ik)))
c       absconst=0.023*(2./4.)*(2.*xlambda_c*1.0E+10/y_cent)**2.78
c       if (absconst.lt.0.1) absconst=0.1
c       write(6,*)E_cent,xrayabs
       absconst=xrayabs*(5.9/1000.) 
       Amut=(1.0-exp(-absconst**0.5))
       dTheta = 0.01 ! horizontal angular range
       xLwg = dTheta*R_bend ! pole length in m
c       write(6,*)'Wriggler B-field and pole lengths =',Bwg,xLwg
c       write(6,*)'boost, central freq and bend radius=',gamma,Omega_c,
c     1            R_bend
c       write(6,*)'vertical ang. range, min and max photon energy=',
c     1 dTheta, E_min, E_max,xlambda_c,absconst,y_cent,Amut,E_cent
c       Psi_min = -asin(1./gamma) ! min vertical angle
c       Psi_max = asin(1./gamma) ! max vertical angle
c       alpha_min = gamma*Psi_min ! boosted min vertical angle
c       alpha_max = gamma*Psi_max ! boosted max vertical angle
c       z_min = (ymin/2.)*(1.+alpha_min**2.)**1.5 ! z=(y/2)*(1+alpha)^1.5
c       z_max = (ymax/2.)*(1.+alpha_max**2.)**1.5
c       write(6,*)alpha_min,alpha_max,z_min,z_max
       xNe=xIe*1.0E-06/qe  ! # of electrons
       xi=1.5*gamma**2*hbarc*1.0E-15/xMe/R_bend ! crit parameter
       xHbyHo=gamma*hbarc*1.0E-15/xMe/R_bend
       tau=(8.*sqrt(3.)/15.)*(hbarc*1.0E-15/xMe/c) ! use hbarc/qe^2 =137
     1      *(1./xHbyHo)**3*(1./gamma**2.)*137.0
       sflip=hbarc*(xLwg/c)*(1.+8.*sqrt(3.)/15.)*0.5*1./tau
c       write(6,*)'sflip probability',sflip
       const1=3.*xNe*gamma*dTheta/(4.*pi**2.*137.)
       const2=4.*const1*xi
       call p3pgs(ymin,ymax,n1,fn01,gn01,vint1) ! integrations for N
       call p3pgs(ymin,ymax,n2,fn02,gn02,vint2) ! integraton for Delta N
       xn=const1*vint1
       dn=const2*vint2
       xa=(dn/xn)*sqrt(2.*xn)
       xpow=xn*E_cent*1.6E-19*1.0E+6! power released in W
       dpow1=dn*E_cent*1.6E-19*1.0E+10 ! power of spin light in W
c       write(6,*)'edep',xdpow !Ebeam,gamma,dTheta,gamma*dTheta,const1,vint1
c       write(6,*)'photon energy, #of photons, up/dn diff, assym, 
c     1 analyzing pwr, photons abs'
c       write(6,15)(E_min+E_max)/2000.,xn,dn,dn/xn,xa,xpow,Amut*xn
       write(10,15)y_cent,E_cent,xn,dn,dn/xn,xa,xpow,Amut*xn
       ymin=ymax
       enddo
       close (10)
 15    format(1x,f6.3,1x,f6.3,1x,e15.3,1x,e15.3,1x,e15.3,1x,e15.3,1x,
     1        e15.3,1x,e15.3)
       END

        SUBROUTINE IKV(V,X,VM,BI,DI,BK,DK)
C
C       =======================================================
C       Purpose: Compute modified Bessel functions Iv(x) and
C                Kv(x), and their derivatives
C       Input :  x --- Argument ( x > 0 )
C                v --- Order of Iv(x) and Kv(x)
C                      ( v = n+v0, n = 0,1,2,..., 0 < v0 < 1 )
C       Output:  BI(n) --- In+v0(x)
C                DI(n) --- In+v0'(x)
C                BK(n) --- Kn+v0(x)
C                DK(n) --- Kn+v0'(x)
C                VM --- Highest order computed
C       Routines called:
C            (1) GAMMA for computing the gamma function
C            (2) MSTA1 and MSTA2 to compute the starting
C                point for backward recurrence
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BI(0:*),DI(0:*),BK(0:*),DK(0:*)
        PI=3.141592653589793D0
        X2=X*X
        N=INT(V)
        V0=V-N
        IF (N.EQ.0) N=1
        IF (X.LT.1.0D-100) THEN
           DO 10 K=0,N
              BI(K)=0.0D0
              DI(K)=0.0D0
              BK(K)=-1.0D+300
10            DK(K)=1.0D+300
           IF (V.EQ.0.0) THEN
              BI(0)=1.0D0
              DI(1)=0.5D0
           ENDIF
           VM=V
           RETURN
        ENDIF
        PIV=PI*V0
        VT=4.0D0*V0*V0
        IF (V0.EQ.0.0D0) THEN
           A1=1.0D0
        ELSE
           V0P=1.0D0+V0
           CALL GAMMA(V0P,GAP)
           A1=(0.5D0*X)**V0/GAP
        ENDIF
        K0=14
        IF (X.GE.35.0) K0=10
        IF (X.GE.50.0) K0=8
        IF (X.LE.18.0) THEN
           BI0=1.0D0
           R=1.0D0
           DO 15 K=1,30
              R=0.25D0*R*X2/(K*(K+V0))
              BI0=BI0+R
              IF (DABS(R/BI0).LT.1.0D-15) GO TO 20
15         CONTINUE
20         BI0=BI0*A1
        ELSE
           CA=DEXP(X)/DSQRT(2.0D0*PI*X)
           SUM=1.0D0
           R=1.0D0
           DO 25 K=1,K0
              R=-0.125D0*R*(VT-(2.0D0*K-1.0D0)**2.0)/(K*X)
25            SUM=SUM+R
           BI0=CA*SUM
        ENDIF
        M=MSTA1(X,200)
        IF (M.LT.N) THEN
           N=M
        ELSE
           M=MSTA2(X,N,15)
        ENDIF
        F2=0.0D0
        F1=1.0D-100
        DO 30 K=M,0,-1
           F=2.0D0*(V0+K+1.0D0)/X*F1+F2
           IF (K.LE.N) BI(K)=F
           F2=F1
30         F1=F
        CS=BI0/F
        DO 35 K=0,N
35         BI(K)=CS*BI(K)
        DI(0)=V0/X*BI(0)+BI(1)
        DO 40 K=1,N
40         DI(K)=-(K+V0)/X*BI(K)+BI(K-1)
        IF (X.LE.9.0D0) THEN
           IF (V0.EQ.0.0D0) THEN
              CT=-DLOG(0.5D0*X)-0.5772156649015329D0
              CS=0.0D0
              W0=0.0D0
              R=1.0D0
              DO 45 K=1,50
                 W0=W0+1.0D0/K
                 R=0.25D0*R/(K*K)*X2
                 CS=CS+R*(W0+CT)
                 WA=DABS(CS)
                 IF (DABS((WA-WW)/WA).LT.1.0D-15) GO TO 50
45               WW=WA
50            BK0=CT+CS
           ELSE
              V0N=1.0D0-V0
              CALL GAMMA(V0N,GAN)
              A2=1.0D0/(GAN*(0.5D0*X)**V0)
              A1=(0.5D0*X)**V0/GAP
              SUM=A2-A1
              R1=1.0D0
              R2=1.0D0
              DO 55 K=1,120
                 R1=0.25D0*R1*X2/(K*(K-V0))
                 R2=0.25D0*R2*X2/(K*(K+V0))
                 SUM=SUM+A2*R1-A1*R2
                 WA=DABS(SUM)
                 IF (DABS((WA-WW)/WA).LT.1.0D-15) GO TO 60
55               WW=WA
60            BK0=0.5D0*PI*SUM/DSIN(PIV)
           ENDIF
        ELSE
           CB=DEXP(-X)*DSQRT(0.5D0*PI/X)
           SUM=1.0D0
           R=1.0D0
           DO 65 K=1,K0
              R=0.125D0*R*(VT-(2.0*K-1.0)**2.0)/(K*X)
65            SUM=SUM+R
           BK0=CB*SUM
        ENDIF
        BK1=(1.0D0/X-BI(1)*BK0)/BI(0)
        BK(0)=BK0
        BK(1)=BK1
        DO 70 K=2,N
           BK2=2.0D0*(V0+K-1.0D0)/X*BK1+BK0
           BK(K)=BK2
           BK0=BK1
70         BK1=BK2
        DK(0)=V0/X*BK(0)-BK(1)
        DO 80 K=1,N
80         DK(K)=-(K+V0)/X*BK(K)-BK(K-1)
        VM=N+V0
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function â(x)
C       Input :  x  --- Argument of â(x)
C                       ( x is not equal to 0,-1,-2,úúú)
C       Output:  GA --- â(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward  
C                recurrence such that the magnitude of    
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point   
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END


      SUBROUTINE P3PGS ( A, B, N, FN, GN, VINT )

c*********************************************************************72
C
C  THIS SUBROUTINE USES THE PRODUCT TYPE THREE-POINT GAUSS-
C  LEGENDRE-SIMPSON RULE COMPOUNDED N TIMES TO APPROXIMATE
C  THE INTEGRAL FROM A TO B OF THE FUNCTION FN(X) * GN(X).
C  FN AND GN ARE FUNCTION SUBPROGRAMS WHICH MUST BE SUPPLIED
C  BY THE USER.  THE RESULT IS STORED IN VINT.
C
      DOUBLE PRECISION A, AG, AM(2,3), B, F(2), FN, G(3),
     &  GN, H, VINT, X(2), Y(2), DBLE

      DATA AM(1,1), AM(2,3) / 2 * 1.718245836551854D0 /,
     &  AM(1,2), AM(2,2) / 2 * 1.D0 /, AM(1,3), AM(2,1)
     &  / 2 * -.2182458365518542D0 /

      H = ( B - A ) / DBLE ( FLOAT ( N ) )
      X(1) = A + .1127016653792583D0 * H
      X(2) = A + .8872983346207417D0 * H
      Y(1) = A + H / 2.D0
      Y(2) = A + H
      VINT = 0.D0
      G(3) = GN ( A )
      DO 3 I = 1, N
        AG = FN ( Y(1) )
        G(1) = G(3)
        DO 1 J = 1, 2
          F(J) = FN ( X(J) )
          G(J+1) = GN ( Y(J) )
          X(J) = X(J) + H
1         Y(J) = Y(J) + H
        VINT = VINT + AG * 4.D0 * G(2)
        DO 3 J = 1, 2
          AG = 0.D0
          DO 2 K = 1, 3
2           AG = AG + AM(J,K) * G(K)
3         VINT = VINT + F(J) * AG
      VINT = H * VINT / 9.D0

      RETURN
      END



      function fn01(x)
      implicit none

      integer n

      double precision fn01
      double precision x

      fn01 = x

      return
      end

      function fn02(x)
      implicit none

      integer n

      double precision fn02
      double precision x

      fn02 = x**2.

      return
      end

      function gn01(x)
      implicit real*8(A-H,O-Z)

c      real*8 Psi_min,Psi_max,alpha_min,alpha_max,y,vint3
c      real*8 gamma

c      integer n
      external fn03
      external gn03      

c      double precision gn01
c      double precision x

      parameter (n=20)

      common gamma,y
    
      Psi_min = -asin(1./gamma) ! min vertical angle
      Psi_max = asin(1./gamma) ! max vertical angle
      alpha_min = 2.*gamma*Psi_min ! boosted min vertical angle
      alpha_max = 2.*gamma*Psi_max ! boosted max vertical angle
      alpha_cutoffm=-0.16
      alpha_cutoffp=0.16
      y=x
      call p3pgs(alpha_min,alpha_max,n,fn03,gn03,vint31)
c      call p3pgs(alpha_min,alpha_cutoffm,n,fn03,gn03,vint31)
c      call p3pgs(alpha_cutoffp,alpha_max,n,fn03,gn03,vint32)
      gn01 = vint31 !+vint32

      return
      end

      function gn02(x)
      implicit none

      real*8 Psi_min,Psi_max,alpha_min,alpha_max,y,vint4
      real*8 gamma

      integer n
      external fn04
      external gn04

      double precision gn02
      double precision x
      parameter (n=20)
      common gamma,y
    
      Psi_min = -asin(1./gamma) ! min vertical angle
      Psi_max = asin(1./gamma) ! max vertical angle
      alpha_min = 0. !16 !for the diff the int is from alpha_cutoff to alpha_max
      alpha_max = 2.*gamma*Psi_max ! boosted max vertical angle
      y=x
      call p3pgs(alpha_min,alpha_max,n,fn04,gn04,vint4)
      gn02 = vint4

      return
      end
 
      function fn03(x)
      implicit none


      integer n

      double precision fn03
      double precision x

    
      fn03 = (1+x**2.)**2.

      return
      end

      function fn04(x)
      implicit none


      integer n

      double precision fn04
      double precision x

    
      fn04 = x*(1+x**2.)**1.5

      return
      end

      function gn03(x)
      implicit real*8 (A-H,O-Z)

c      real*8 z,v,k23,k13,gamma,y,vm
c      dimension BI(0:*),DI(0:*),BK(0:*),DK(0:*)
  
c      integer n

      double precision gn03
      double precision x
      common gamma,y
      COMMON BI(0:250),DI(0:250),BK(0:250),DK(0:250)

      xk23=0.
      xk13=0.
      z=(y/2.)*(1+x**2)**1.5 ! z= (omega/2omega_c)*(1+alpha^2)^3/2
c      write(6,*)'gn03 gamma,y,x,z',gamma,y,x,z
      v=2./3.
      CALL IKV(V,z,VM,BI,DI,BK,DK)
      xk23=BK(0)

      v=1./3.
      CALL IKV(V,z,VM,BI,DI,BK,DK)
      xk13=BK(0)
c      write(6,*)'g03 k23 k13',xk23,xk13 
      gn03 = xk23**2. + x**2.*xk13**2/(1+x**2.)

      return
      end

      function gn04(x)
      implicit real*8(A-H,O-Z)

c      real*8 z,v,k23,k13,y,vm,gamma  
c      Dimension BI(0:*),DI(0:*),BK(0:*),DK(0:*)  

c      integer n

c      double precision gn04
c      double precision x
      Common gamma,y
      COMMON BI(0:250),DI(0:250),BK(0:250),DK(0:250)

      z=(y/2.)*(1+x**2)**1.5
      v=2./3.
      CALL IKV(V,z,VM,BI,DI,BK,DK)
      xk23=BK(0)

      v=1./3.
      CALL IKV(V,z,VM,BI,DI,BK,DK)
      xk13=BK(0)
 
      gn04 = xk23*xk13

      return
      end

