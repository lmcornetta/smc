      module mplane2

c***
c  Modified to remove common block CALLIN. Similar to the modification
c  done in module MPLANE1
c***

      use mplane1, only: nap, lap, map, iend, ll, nmax, lmax, mmax,
     .                   ientr, jentr, ifirst, cterm, cnorm, pi,
     .                   twopi, pi3haf, pi5hf2, piquart, nqt, norm
     
      contains

c----------------------------------------------------------------
      subroutine int1nh(s,XA,ALPA,NTYPA,XB,XC,ALPC,NTYPC)
!     subroutine int1nh(s,XA,ALPA,NTYPA,NTYPB,XB,ALPB,XC,ALPC,NTYPC)
c    . NTYPD,XD,ALPD)

      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 PATHOS,PF
c     COMPLEX*16 PFAB,PFCD,S,DAB,PQX(3),T,RN
      COMPLEX*16 S,DAB,PQX,T,RN!,PFAB
      dimension pqx(3)
      COMPLEX*16 RNLM
!fk
      integer :: ntypa, ntypc !, ntypb, ntypb
      real*8 :: xa(3), xb(3), xc(3)!, xd(3)
      real*8 :: alpa, alpc!, alpb, alpd
!     dimension xa(3), xb(3), xc(3) !, xd(3)
C=======================================================================
C
C       HYBRID INTEGRAL
C       (NTYPA,ALPA,XA,YA,ZA)*(1/R(C))*(KX,KY,KZ)
C       R(C)=(XC,YC,ZC)
C       (KX,KY,KZ)=(XB,YB,ZB)
C
C=======================================================================
      COMMON /OVRFLO/PATHOS,PATHO
cccc  COMMON/CALLIN/XA,YA,ZA,ALPA,NTYPA,NTYPB,XB,YB,ZB,ALPB,
cccc * XC,YC,ZC,ALPC,NTYPC,NTYPD,XD,YD,ZD,ALPD
cccc  COMMON /FUNC/NAP(60),LAP(60),MAP(60),IEND(60),CTERM(30),LL(60)
cccc * ,NMAX(60),LMAX(60),MMAX(60),IENTR(60),JENTR(60),CNORM(60)
cccc * ,PI,TWOPI,PI3HAF,PI5HF2,PIQUART,IFIRST
!fk   COMMON/CALLR/L,NM,LM,MM,PQX,PQY,PQZ,MAS
      COMMON/CALLR/L,NM,LM,MM,PQX,MAS
      DIMENSION DAB(8),RN(4),RNLM(4625),NLMP(8)
      COMMON /ERN/RNLM
      EQUIVALENCE (RNLM(1906),RN(1))
      MAS=0
CRVX$ VAX
!fk   PQX=DCMPLX(XA-XC,.5D0*XB/ALPA)
c     PQY=DCMPLX(YA-YC,.5D0*YB/ALPA)
c     PQZ=DCMPLX(ZA-ZC,.5D0*ZB/ALPA)
      PQX(:)=DCMPLX(XA(:)-XC(:),.5D0*XB(:)/ALPA)
CRVX$ END
CRVX$ CRAY
C     PQX=CMPLX(XA-XC,.5D0*XB/ALPA)
C     PQY=CMPLX(YA-YC,.5D0*YB/ALPA)
C     PQZ=CMPLX(ZA-ZC,.5D0*ZB/ALPA)
CRVX$ END
      L=LL(NTYPA)
C=======================================================================
C
C       CHECK, BY SYMMETRY, FOR ZERO-VALUE OF INTEGRAL
C
C=======================================================================
      IF(L.EQ.0) GO TO 1
CRVX$ VAX
      IF((DREAL(PQX(1)).EQ.0.D0).AND.(DIMAG(PQX(1)).EQ.0.D0))MAS=4
      IF((DREAL(PQX(2)).EQ.0.D0).AND.(DIMAG(PQX(2)).EQ.0.D0))MAS=MAS+2
      IF((DREAL(PQX(3)).EQ.0.D0).AND.(DIMAG(PQX(3)).EQ.0.D0))MAS=MAS+1
CRVX$ END
CRVX$ CRAY
C     IF((REAL(PQX).EQ.0.D0).AND.(AIMAG(PQX).EQ.0.D0))MAS=4
C     IF((REAL(PQY).EQ.0.D0).AND.(AIMAG(PQY).EQ.0.D0))MAS=MAS+2
C     IF((REAL(PQZ).EQ.0.D0).AND.(AIMAG(PQZ).EQ.0.D0))MAS=MAS+1
CRVX$ END
      ISYMP=MAS
      IF(IAND(NMAX(NTYPA),1).EQ.0) ISYMP=IAND(ISYMP,3)
      IF(IAND(LMAX(NTYPA),1).EQ.0)  ISYMP=IAND(ISYMP,5)
      IF(IAND(MMAX(NTYPA),1).EQ.0) ISYPM=IAND(ISYMP,6)
      IF(ISYMP.EQ.0) GO TO 1
C=======================================================================
C
C       SYMMETRY REVEALS INTEGRAL IS ZERO
C
C=======================================================================
      S=0.D0
      RETURN
C=======================================================================
C
C       EVALUATE NON-ZERO INTEGRAL
C
C=======================================================================
1     T=ALPA*sum(pqx(:)**2)
CRVX$ VAX
      PF=TWOPI/ALPA*CDEXP(DCMPLX(0.D0,sum(XB(:)*XA(:))))
CRVX$ END
CRVX$ CRAY
C     PF=TWOPI/ALPA*CEXP(CMPLX(0.D0,XB*XA+YB*YA+ZB*ZA))
CRVX$ END
      PATHO=-0.25D0*sum(XB(:)**2)/ALPA
      PATHOS=PATHO-T
      CALL FORMGN2(L,T,ALPA,RN)
      IF(L.EQ.0) GO TO 140
      NM=NMAX(NTYPA)+1
      MM=MMAX(NTYPA)+1
      LM=LMAX(NTYPA)+1

      call expldp2(DAB,NLMP,LP,1,ALPA,NTYPA,XB,ALPC,NTYPC)
!    . XA,ALPA,NTYPA,NTYPB,XB,ALPB,XC,ALPC,NTYPC)
!    . ALPA,NTYPA,XB,ALPC,NTYPC)

      L=L+2
      CALL GETR2
      S=0.D0
      DO 130 I=1,LP
      IND=NLMP(I)
130   S=S+DAB(I)*RNLM(IND)
      S=S*PF
      RETURN
140   S=RN(1)*PF
      RETURN
ccc   END
      end subroutine int1nh
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE FORMGN2(L,T,A,R)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 R,T
cmod  DIMENSION R(1)
      DIMENSION r(*)
      CALL FINT2(L,T,R)
      IF(L.EQ.0) RETURN
      B=-2.D0*A
      S=B
      DO 1 I=1,L
      R(I+1)=R(I+1)*S
1     S=S*B
      RETURN
ccc   END
      end subroutine formgn2
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE FINT2(ITOP,X,FVEC)
      IMPLICIT REAL*8(A-H,O-Z)
cccc  COMPLEX*16   X,Y,TERM,PTLSUM,FVEC    ,X2  ,ARG ,DCERFC2,PATHOS
      COMPLEX*16   X,Y,TERM,PTLSUM,FVEC    ,X2  ,ARG ,PATHOS
C=======================================================================
C
C    SUBROUTINE TO COMPUTE INTEGRALS FROM 0 TO 1 OF
C    (T**2M) * EXP(-X*(T**2)).      STORES VALUES FOR M=0 TO ITOP IN
C      FVEC.
C
C=======================================================================
      COMMON/OVRFLO/PATHOS,PATHO
      DIMENSION FVEC(10)
      DATA SQTPI2/.886226925452758E0/
CRVX$ VAX
      Y = CDEXP(PATHOS)
CRVX$ END
CRVX$ CRAY
C     Y = CEXP(PATHOS)
CRVX$ END
      X2 = 2.D00*X
CRVX$ VAX
      IF (CDABS(X)-10.D00) 10,10,20
CRVX$ END
CRVX$ CRAY
C     IF (CABS(X)-10.D00) 10,10,20
CRVX$ END
C=======================================================================
C
C    SMALL ARGUMENT -- USE ANALYTIC SERIES EXPANSION OF INCOMPLETE
C    GAMMA FUNCTION
C
C=======================================================================
  10  A = ITOP
       A=A+0.5D0
      TERM=1.D00/A
       PTLSUM=TERM
       DO 11 I=2,50
       A=A+1.D00
       TERM=TERM*X/A
       PTLSUM=PTLSUM+TERM
CCRVX$ VAX
C             IF (CDABS(PTLSUM).EQ.0.0D0)GO TO 11
C             IF (CDABS(TERM/PTLSUM)-1.D-12) 12,11,11
CCRVX$ END
CCRVX$ CRAY
CC            IF (CABS(PTLSUM).EQ.0.0D0)GO TO 11
CC            IF (CABS(TERM/PTLSUM)-1.D-12) 12,11,11
CCRVX$ END
CRVX$ VAX
      IF((DREAL(TERM)*DREAL(TERM)+DIMAG(TERM)*DIMAG(TERM)).LE.
     1 (DREAL(PTLSUM)*DREAL(PTLSUM)+DIMAG(PTLSUM)*DIMAG(PTLSUM))
     2 *1.D-24) GO TO 12
CRVX$ END
CRVX$ CRAY
C     IF(( REAL(TERM)* REAL(TERM)+AIMAG(TERM)*AIMAG(TERM)).LE.
C    1 ( REAL(PTLSUM)* REAL(PTLSUM)+AIMAG(PTLSUM)*AIMAG(PTLSUM))
C    2 *1.D-24) GO TO 12
CRVX$ END
  11  CONTINUE
      WRITE(6,999) M,X
  999 FORMAT(27H0CONVERGENCE FAILED IN FINT,I5,2D16.8)
      CALL EXIT (1)
   12 FVEC(ITOP+1) = 0.5D00*PTLSUM*Y
      IF(ITOP.EQ.0) GO TO 15
C=======================================================================
C
C    BACKWARDS RECURSION
C
C=======================================================================
      M = ITOP+1
  14  FVEC(M-1) = (X2*FVEC(M)+Y)/ DFLOAT(2*M-3)
      M = M-1
      IF(M-1) 15,15,14
  15  RETURN
C=======================================================================
C
C    LARGE ARGUMENT - FORWARD RECURSION FROM ERROR FUNCTION
C
C=======================================================================
  20  CONTINUE
CRVX$ VAX
      ARG = CDSQRT(X)
CRVX$ END
CRVX$ CRAY
C     ARG = CSQRT(X)
CRVX$ END
      FVEC(1) = SQTPI2/ARG*( DEXP(PATHO )-             Y*DCERFC2(ARG))
      IF(ITOP.EQ.0) GO TO 17
      I=1
  16  FVEC(I+1) = ( DFLOAT(2*I-1)*FVEC(I)-Y)/X2
      I=I+1
      IF(I-ITOP-1) 16,16,17
  17  RETURN
ccc   END
      end subroutine fint2
c----------------------------------------------------------------

c----------------------------------------------------------------
CRVX$ VAX
      DOUBLE COMPLEX FUNCTION DCERFC2(ZP)
CRVX$ END
CRVX$ CRAY
C     COMPLEX FUNCTION DCERFC2(ZP)
CRVX$ END
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16           ZP,AI,Z,ZA,CERF
      REAL*8 LAMBDA
      LOGICAL B
      AI=(0.D00,1.D00)
      Z=AI*ZP
CRVX$ VAX
      X= DREAL(Z)
      Y=DIMAG(Z)
CRVX$ END
CRVX$ CRAY
C     X= REAL(Z)
C     Y=AIMAG(Z)
CRVX$ END
      IQ=1
      IF (X .LT. 0.D00) IQ=2
      IF (Y .LT. 0.D00) IQ=3
      IF (Y .LT. 0.D00 .AND. X  .GE. 0.D00) IQ=4
      X=DABS(X)
      Y=DABS(Y)
CRVX$ VAX
      ZA=DCMPLX(X,Y)
CRVX$ END
CRVX$ CRAY
C     ZA=CMPLX(X,Y)
CRVX$ END
      IF (Y .LT. 4.29D0 .AND.X .LT. 5.33D0) GO TO 10
      GO TO 15
   10 S=(1.D00-Y/4.29D0)* DSQRT(1.D00-X*X/28.41D0)
      H = 1.6D0*S
      H2=2.D00*H
      K= 6+ 23*S
      LAMBDA=H2**K
      NU=9+ 21*S
      GO TO 20
   15 H=0.D00
      K=0
      NU=8
C=======================================================================
C
C  20 IF(H .GT. 0.D0D0) LAMBDA=H2**K            THIS CARD IS DELETED
C
C=======================================================================
   20 CONTINUE
      B=.FALSE.
      IF(H.EQ. 0.D00 .OR. LAMBDA .EQ. 0.D00) B=.TRUE.
      R1=0.D00
      R2=0.D00
      S1=0.D00
      S2=0.D00
      N=NU+1
   30 N=N-1
      IF(N .LT. 0) GO TO 50
      NP1=N+1
      T1=Y+H+NP1*R1
      T2= X-NP1*R2
      C= .5D0/(T1*T1+T2*T2)
      R1=C*T1
      R2=C*T2
      IF(H.LE.0.D00 .OR.  N .GT. K) GO TO 30
      T1=LAMBDA+S1
      S1=R1* T1-R2*S2
      S2=R2*T1+R1*S2
      LAMBDA=LAMBDA/H2
      GO TO 30
   50 CONTINUE
      IF (Y .NE. 0.D00) GO TO 60
CALPHA
      RE = 0.0D0
      IF (X.LT.20.D0) RE=DEXP(-X*X)/1.12837916709551D0
CAPPN RE=DEXP(-X*X)/1.12837916709551D0

CALPHA
      GO TO 70
   60 RE=S1
      IF (B) RE=R1
   70 AIM=S2
      IF (B) AIM=R2
CRVX$ VAX
      CERF=1.12837916709551D0* DCMPLX(RE,AIM)
CRVX$ END
CRVX$ CRAY
C     CERF=1.12837916709551D0* CMPLX(RE,AIM)
CRVX$ END
      IF (IQ .EQ. 1) GO TO 90
      IF (IQ .NE. 2) GO TO 80
CRVX$ VAX
      CERF=DCONJG(CERF)
CRVX$ END
CRVX$ CRAY
C     CERF=CONJG(CERF)
CRVX$ END
      GO TO 90
CRVX$ VAX
   80 CERF=2.D00*CDEXP(-ZA*ZA)-CERF
      IF (IQ .EQ. 3) GO TO 90
      CERF=DCONJG(CERF)
CRVX$ END
CRVX$ CRAY
C  80 CERF=2.D00*CEXP(-ZA*ZA)-CERF
C     IF (IQ .EQ. 3) GO TO 90
C     CERF=CONJG(CERF)
CRVX$ END
C=======================================================================
C
C  EXP(Z**2) TAKEN OUT OF CERFC 4/22/74 -APPEARS IN FINT
C  90 DCERFC2=CDEXP(Z*Z)*CERF
C
C=======================================================================
  90   DCERFC2 = CERF
      RETURN
ccc   END
      end function dcerfc2
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE GETR2
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ARL,ARM,ARN,RNLM,RDUM,PQX(3),RM,RL,RN
C=======================================================================
C
C     GENERATE R(N+1,L+1,M+1) BY RECURRENCE FOR N+1=1 TO NM,
C     L+1=1 TO LM, M+1=1 TO MM AND (N+1)+(L+1)+(M+1) SMALLER THAN LL+2
C     CALLING PARAMETERS ARE PASSED IN COMMON /CALLR/
C     PQX,PQY,PQZ ARE COORDINATES OF VECTOR PQ
C     MASX IS 1 IF R IS ZERO BY SYMMETRY FOR ODD N, QTHERWISE MASX=0
C     MASY IS 1 IF R IS ZERO BY SYMMETRY FOR ODD L, OTHERWISE MASY=0
C     MASZ IS 1 IF R IS ZERO BY SYMMETRY FOR ODD M, OTHERWISE MASZ=0
C     ZERO VALUES OF R BY SYMMETRY WILL NOT BE COMPUTED OR STORED.
C     REGION RN CONTAINS MODIFIED FN INTEGRALS FROM FORMGN FOR
C     N+1=1 TO LL-1 AT ENTRANCE. RN,RL,RM ARE OVERWRITTEN
C     DURING EXECUTION.
C
C=======================================================================
      COMMON/ERN/RNLM(17,17,16),RDUM
!fk   COMMON /CALLR/LL,NM,LM,MM,PQX,PQY,PQZ,MAS
      COMMON /CALLR/LL,NM,LM,MM,PQX,MAS
      DIMENSION RN(17,2),RL(17,2),RM(17,2)
      DIMENSION FACT(9)
      EQUIVALENCE (RL(1,1),RNLM(2,13,7)),(RM(1,1),RNLM(2,15,7)),
     $(RN(1,1),RNLM(2,11,7))
      DATA FACT/1.D0,1.D0,3.D0,15.D0,105.D0,945.D0,10395.D0,135135.D0,
     * 2027025.D0/
       MASX=IAND(MAS,4)
       MASY=IAND(MAS,2)
       MASZ=IAND(MAS,1)
      IF(MASX.EQ.0) GO TO 2
1      NINC=2
      NCO=0
      GO TO 3
2      NINC=1
3      IF(MASY.EQ.0) GO TO 5
4     LINC=2
      GO TO 10
5     LINC=1
10      CONTINUE
      DO 590 N=1,NM,NINC
      JF1=LL-N
      IF(N-2) 380,400,420
380   CONTINUE
      DO 390 J=1,JF1
390   RL(J,1)=RN(J,1)
      GO TO 440
  400 CONTINUE
      DO 410 J=1,JF1
      RN(J,2)=PQX(1)*RN(J+1,  1)
      RL(J,1)=RN(J,2)
  410 CONTINUE
      GO TO 440
  420 CONTINUE
      IF(MASX.EQ.0) GO TO 421
      NCO=NCO+1
      DO 431 J=1,JF1
431   RL(J,1)=RN(J+NCO,1)*FACT(NCO+1)
      GO TO 440
421    CONTINUE
      FN=DFLOAT(N-2)
      DO 430 J=1,JF1
      ARN=PQX(1)*RN(J+1,2)+FN*RN(J+1,1)
      RN(J+1,1)=RN(J+1,2)
      RN(J,2)=ARN
      RL(J,1)=ARN
  430 CONTINUE
  440 CONTINUE
      LCO=0
      LF=MIN0(LM,JF1)
      LF1=MIN0(JF1+1,LF+MM)
      DO 580 L=1,LF,LINC
      JF2=LF1-L
      IF(L-2) 450,470,490
450   CONTINUE
      DO 460 J=1,JF2
460   RM(J,1)=RL(J,1)
      GO TO 505
  470 CONTINUE
      DO 480 J=1,JF2
      RL(J,2)=PQX(2)*RL(J+1,  1)
      RM(J,1)=RL(J,2)
  480 CONTINUE
      GO TO 505
  490 CONTINUE
      IF(MASY.EQ.0) GO TO 491
      LCO=LCO+1
       DO 492 J=1,JF2
492    RM(J,1)=RL(J+LCO,1)*FACT(LCO+1)
      GO TO 505
  491  CONTINUE
      FL=DFLOAT(L-2)
      DO 500 J=1,JF2
      ARL=PQX(2)*RL(J+1,2)+FL*RL(J+1,1)
      RL(J+1,1)=RL(J+1,2)
      RL(J,2)=ARL
      RM(J,1)=ARL
  500 CONTINUE
  505 CONTINUE
      MF=MIN0(MM,JF2)
      IF(MASZ.NE.0) GO TO 541
      MF1=MF+1
      RNLM(N,L,1)=RM(1,1)
      IF(MF.LT.2) GO TO 580
      JF3=MF1-2
      DO 540 J=1,JF3
      RM(J,2)=PQX(3)*RM(J+1,  1)
  540 CONTINUE
      RNLM(N,L,2)=RM(1,2)
      IF(MF.LT.3) GO TO 580
      DO 570 M=3,MF
      JF3=MF1-M
      FM=DFLOAT(M-2)
      DO 560 J=1,JF3
      ARM=PQX(3)*RM(J+1,2)+FM*RM(J+1,1)
      RM(J+1,1)=RM(J+1,2)
      RM(J,2)=ARM
  560 CONTINUE
      RNLM(N,L,M)=RM(1,2)
  570 CONTINUE
      GO TO 580
541    K=0
      DO 571 M=1,MF,2
      K=K+1
571   RNLM(N,L,M)=RM(K,1)*FACT(K)
  580 CONTINUE
  590 CONTINUE
      RETURN
ccc   END
      end subroutine getr2
c----------------------------------------------------------------

c----------------------------------------------------------------
      subroutine expldp2(D,NLM,K,IT,ALPA,NTYPA,XB,ALPC,NTYPC)
!    . XA,ALPA,NTYPA,NTYPB,XB,ALPB,XC,ALPC,NTYPC)
!    . ALPA,NTYPA,XB,ALPC,NTYPC)

      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 DX,DY,DZ
      COMPLEX*16 D,XAB,YAB,ZAB!,PFAB
!fk   COMPLEX*16 D,XAB(3),PFAB,PFCD
      integer :: ntypa, ntypc!, ntypb!, ntypd
      real*8 :: xb(3)!, xa(3), xc(3)!, xd(3)
      real*8 :: alpa, alpc!, alpb!, alpd
!     dimension xa(3), xb(3), xc(3) !, xd(3)

C=======================================================================
C
C       EVALUATE S AND P FUNCTIONS AND EXPAND GAUSSIANS IN HERMITE
C       POLINOMIALS ABOUT A COMPLEX ORIGIN.
C
C=======================================================================
cccc  COMMON /FUNC/NAP(60),LAP(60),MAP(60),IEND(60),CTERM(30),LL(60)
cccc * ,NMAX(60),LMAX(60),MMAX(60),IENTR(60),JENTR(60),CNORM(60)
cccc * ,PI,TWOPI,PI3HAF,PI5HF2,PIQUART,IFIRST
cccc  COMMON /CALLIN/XA,YA,ZA,ALPA,NTYPA,NTYPB,XB,YB,ZB,ALPB,
cccc * XC,YC,ZC,ALPC,NTYPC,NTYPD,XD,YD,ZD,ALPD

      DIMENSION DX(4),DY(4),DZ(4),IX(4),IY(4),IZ(4)
      DIMENSION D(8),NLM(8)
      MAS=7
      IF(IT.EQ.2) GO TO 4
      A=.5D0/ALPA
      NTY=NTYPA
      IF(XB(1).EQ.0.D0)GO TO 1
CRVX$ VAX
      XAB=DCMPLX(0.D0,A*XB(1))
CRVX$ END
CRVX$ CRAY
C     XAB=CMPLX(0.D0,A*XB)
CRVX$ END
      MAS=IAND(MAS,3)
C=======================================================================
C
C        IAND(,) CORESPONDS TO  .INT. (INTERSECTION = BITWISE AND)
C       IOR(,)  CORESPONDS TO  .UN.  (UNION = BITWISE OR)
C       IAND,IOR ARE FOR VAX-11
C       INT,UN   ARE FOR LRLTRAN
C
C=======================================================================
1     IF(XB(2).EQ.0.D0) GO TO 2
CRVX$ VAX
      YAB=DCMPLX(0.D0,A*XB(2))
CRVX$ END
CRVX$ CRAY
C     YAB=CMPLX(0.D0,A*YB)
CRVX$ END
      MAS=IAND(MAS,5)
2     IF(XB(3).EQ.0.D0) GO TO 8
CRVX$ VAX
      ZAB=DCMPLX(0.D0,A*XB(3))
CRVX$ END
CRVX$ CRAY
C     ZAB=CMPLX(0.D0,A*ZB)
CRVX$ END
      MAS=IAND(MAS,6)
8     IF(IT.EQ.1)GO TO 7
      XAB=-XAB
      YAB=-YAB
      ZAB=-ZAB
      GO TO 7
4     A=.5D0/ALPC
      NTY=NTYPC
cfk   IF(XD.EQ.0.D0) GO TO 5
cfk   IF(XD(1).EQ.0.D0) GO TO 5
CRVX$ VAX
cfk   XAB=DCMPLX(0.D0,-A*XD)
cfk   XAB=DCMPLX(0.D0,-A*XD(1))
CRVX$ END
CRVX$ CRAY
C     XAB=CMPLX(0.D0,-A*XD)
CRVX$ END
      MAS=IAND(MAS,3)
cfk5  IF(YD.EQ.0.D0) GO TO 6
cfk5  IF(XD(2).EQ.0.D0) GO TO 6
CRVX$ VAX
cfk   YAB=DCMPLX(0.D0,-A*YD)
cfk   YAB=DCMPLX(0.D0,-A*XD(2))
CRVX$ END
CRVX$ CRAY
C     YAB=CMPLX(0.D0,-A*YD)
CRVX$ END
      MAS=IAND(MAS,5)
cfk6  IF(ZD.EQ.0.D0) GO TO 7
cfk6  IF(XD(3).EQ.0.D0) GO TO 7
CRVX$ VAX
cfk   ZAB=DCMPLX(0.D0,-A*ZD)
cfk   ZAB=DCMPLX(0.D0,-A*XD(3))
CRVX$ END
CRVX$ CRAY
C     ZAB=CMPLX(0.D0,-A*ZD)
CRVX$ END
      MAS=IAND(MAS,6)
7     CONTINUE
      IF(NTY.GT.4) GO TO 120
C=======================================================================
C
C       CHECK INTEGRALS WITH EXPLICIT FORMS
C
C=======================================================================
      GO TO(10,20,30,40)NTY
10    D(1)=1.D0
      NLM(1)=1
      K=1
      RETURN
C=======================================================================
C
C       TYPE #1 N=L=M=0
C
C=======================================================================
20    D(1)=A
      NLM(1)=2
      K=1
      IF(IAND(MAS,4).NE.0) RETURN
C=======================================================================
C
C       TYPE #2 N=1,L=M=0
C
C=======================================================================
      D(2)=XAB
      NLM(2)=1
      K=2
      RETURN
C=======================================================================
C
C       TYPE #3 L=1,M=N=0
C
C=======================================================================
30    D(1)=A
      NLM(1)=18
      K=1
      IF(IAND(MAS,2).NE.0) RETURN
      D(2)=YAB
      NLM(2)=1
      K=2
C=======================================================================
C
C       TYPE #4 M=1,L=N=0
C
C=======================================================================
      RETURN
40    D(1)=A
      NLM(1)=290
      K=1
      IF(IAND(MAS,1).NE.0) RETURN
      D(2)=ZAB
      NLM(2)=1
      K=2
      RETURN
C=======================================================================
C
C       EVALUATE INTEGRALS WITH NON-EXPLICIT FORMS AND EXPAND X**N IN
C       HERMITE POLYNOMIALS
C
C=======================================================================
120   I=IEND(NTY)
      N=NAP(I)+1
      GO TO (121,122,123,124)N
121   DX(1)=1.D0
      IX(1)=1
      KX=1
      GO TO 130
C=======================================================================
C
C       N=0
C
C=======================================================================
122   DX(1)=A
      IX(1)=2
      IF(IAND(MAS,4).NE.0) GO TO 125
C=======================================================================
C
C       COMPUTE N+1,L+1,M+1 FOR EACH GAUSSIAN
C
C=======================================================================
      DX(2)=XAB
C=======================================================================
C
C       N=1
C
C=======================================================================
      KX=2
      IX(2)=1
      GO TO 130
125   KX=1
      GO TO 130
C=======================================================================
C
C       N=2
C
C=======================================================================
123   DX(1)=A*A
      IX(1)=3
      IF(IAND(MAS,4).NE.0) GO TO 126
      DX(2)=2.D0*XAB*A
      IX(2)=2
      DX(3)=A+XAB*XAB
      IX(3)=1
      KX=3
      GO TO 130
C=======================================================================
C
C       N=3
C
C=======================================================================
126   DX(2)=A
      IX(2)=1
      KX=2
      GO TO 130
124   DX(1)=A*A*A
      IX(1)=4
      IF(IAND(MAS,4).NE.0) GO TO 127
      DX(2)=3.D0*XAB*A*A
      IX(2)=3
      DX(3)=3.D0*(XAB*XAB+A)*A
      IX(3)=2
      DX(4)=XAB*(XAB*XAB+3.D0*A)
      IX(4)=1
      KX=4
      GO TO 130
127   DX(2)=3.D0*A*A
      IX(2)=2
      KX=2
130   N=LAP(I)+1
C=======================================================================
C
C       EXPAND Y**L IN HERMITE POLYNOMIALS
C
C=======================================================================
      GO TO (131,132,133,134)N
131   DY(1)=1.D0
      IY(1)=0
      KY=1
      GO TO 140
132   DY(1)=A
      IY(1)=17
      IF(IAND(MAS,2).NE.0) GO TO 135
      DY(2)=YAB
      IY(2)=0
      KY=2
      GO TO 140
135   KY=1
      GO TO 140
133   DY(1)=A*A
      IY(1)=34
      IF(IAND(MAS,2).NE.0) GO TO 136
      DY(2)=2.D0*YAB*A
      IY(2)=17
      DY(3)=A+YAB*YAB
      IY(3)=0
      KY=3
      GO TO 140
136   DY(2)=A
      IY(2)=0
      KY=2
      GO TO 140
134   DY(1)=A*A*A
      IY(1)=51
      IF(IAND(MAS,2).NE.0) GO TO 137
      DY(2)=3.D0*YAB*A*A
      IY(2)=34
      DY(3)=3.D0*A*(YAB*YAB+A)
      IY(3)=17
      DY(4)=YAB*(YAB*YAB+3.D0*A)
      IY(4)=0
      KY=4
      GO TO 140
137   DY(2)=3.D0*A*A
C=======================================================================
C
C       EXPAND Z**M IN HERMITE POLYNOMIALS
C
C=======================================================================
      IY(2)=17
      KY=2
140   N=MAP(I)+1
      GO TO (141,142,143,144)N
141   DZ(1)=1.D0
      IZ(1)=0
      KZ=1
      GO TO 150
142   DZ(1)=A
      IZ(1)=289
      IF(IAND(MAS,1).NE.0) GO TO 145
      DZ(2)=ZAB
      IZ(2)=0
      KZ=2
      GO TO 150
145   KZ=1
      GO TO 150
143   DZ(1)=A*A
      IZ(1)=578
      IF(IAND(MAS,1).NE.0) GO TO 146
      DZ(2)=2.D0*ZAB*A
      IZ(2)=289
      DZ(3)=A+ZAB*ZAB
      IZ(3)=0
      KZ=3
      GO TO 150
146   DZ(2)=A
      IZ(2)=0
      KZ=2
      GO TO 150
144   DZ(1)=A*A*A
      IZ(1)=867
      IF(IAND(MAS,1).NE.0) GO TO 147
      DZ(2)=3.D0*ZAB*A*A
      IZ(2)=578
      DZ(3)=3.D0*A*(ZAB*ZAB+A)
      IZ(3)=289
      DZ(4)=ZAB*(ZAB*ZAB+3.D0*A)
      IZ(4)=0
      KZ=4
      GO TO 150
147   DZ(2)=3.D0*A*A
      IZ(2)=289
      KZ=2
150   K=0
      DO 151 IXX=1,KX
      DO 151 IYY=1,KY
      DO 151 IZZ=1,KZ
      K=K+1
      D(K)=DX(IXX)*DY(IYY)*DZ(IZZ)
151   NLM(K)=IX(IXX)+IY(IYY)+IZ(IZZ)
      RETURN
ccc   END
      end subroutine expldp2
c----------------------------------------------------------------

      end module mplane2
