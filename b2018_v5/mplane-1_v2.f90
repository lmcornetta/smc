module mplane1

!***
!
!  Modifications introduced by Fábris Kossoski in the begining of 2017.
!
!  In the old code, the routine the performed the two-electron integral
!  (ab|cd) was called for every quartet (a,b,c,d), where a, b and c
!  represent gaussian primitives and d denotes a plane wave.
!  Now, the overlaps (ab) and (cd) are pre-computed in abint and cdint
!  subroutines, while the actual integral is performed in the int2ns 
!  subroutine.
!
!  The part dealing with evaluating the symmetry of the integral was
!  removed, since the direction of the plane waves we employ never coincides
!  with any axis.
!
!  Some variables in int2ns, int1nh and cint were vectorized.
!
!***

!***
!  This is the good old mplane module written in fortran77. To avoid
!  too much coding, the modifications were kept to a minimum. 
!   i) All common  blocks were kept (including those interfacing the onk3dk module); 
!  ii) The block data MPDATA is now obsolete;
! iii) The dimension of array X in subroutine XYZ was increased from 100 to 225 
!      (as in the calling subroutine INT2NS) to avoid violating the declared boundary); 
!  iv) Assumed size array declarations were changed from (1) to (*) to match fortran90 standard
!***

!***
! Modified to include only the subroutines/functions relevant to the
! electronic matrix elements (three gaussians + plane wave), driven by
! INT2NS subroutine.
!***

!***
! Modified to remove the common blocks FUNC, TRANS and DTNORM. The content 
! is now declared in the module head (the common blocks used to be loaded
! in the block data MPDATA. The real variables of FUNC are loaded in
! subroutine SETFN (as in older versions of the code).
!***

!***
! Modified to remove the common block CALLIN. In older versions, this common
! block was used as interface with the main c code (PLANE subroutine in
! ONK3DK module). The variables are now passed (form PLANE to INT2NS) as 
! arguments. These are also passed on to EXPLDP and DD subroutines as 
! arguments, since CALLIN was also used to interface INT2NS and EXPLDP. Note 
! also that useless variables PFAB and PFCD were removed.
!***

!***
! Modified to remove the common block CALLR. In older versions, this common
! block was used as inteface between INT2NS and GETR subroutines. The 
! variables are now passed as arguments.
!***

!***
! Modified to remove the common block CALLEX. In older versions, this common
! block was used as inteface between INT2NS and EXPLDD subroutines. The 
! variables are now passed as arguments.
!***

!***
! Modified to remove the common block OVRFLO. In older versions, this common
! block was used as inteface between INT2NS and FINT subroutines. The 
! variables are now passed as arguments. This is done in a somewhat unelegant
! way since INT2NS calls FORMGN which in turn calls FINT. The variables are
! passed on through FORMGN. This will hopefully make paralelization easier.
!
! Incidentally, the declaration of  had to be changed to either fvec(1)
! or fvec(*) (the latter was employed). Otherwise there would be a compile
!***

!***
! Modified to remove all EQUIVALENCE statments. Since RN(1) becomes FVEC in
! FINT (it is passed from INT2NS to FORMGN and then to FINT), its declaration
! in INT2NS and FORMNG was changed from RN(1) to RN(10) to match the dimension
! of FVEC(10) in INT1 (the former version was based on assumed array sizes).
!
! Also, RNLM(1906:1915)=RN(1:10) is enforced after calling FORMGN from
! INT2NS to mimic the old equivalence statment.
!***

!***
! Modified to remove common block ERN. A few more operations were added to
! subroutine GETR to mimic the olv common block/equivalence functionality
!***

!     BLOCK DATA MPDATA

!***
!     COMMON /FUNC/ NAP(60),LAP(60),MAP(60),IEND(60),CTERM(30),LL(60)
!    1,NMAX(60),LMAX(60),MMAX(60),IENTR(60),JENTR(60),CNORM(60)
!    2,PI,TWOPI,PI3HAF,PI5HF2,PIQUART,IFIRST
!
integer, save :: NAP(60), LAP(60), MAP(60), IEND(60), LL(60),  NMAX(60), LMAX(60), MMAX(60), IENTR(60), JENTR(60), IFIRST
real*8,  save :: CTERM(30), CNORM(60), PI, TWOPI, PI3HAF, PI5HF2, PIQUART
!***

!***
!     COMMON /TRANS/NQT(20)
!
integer, save :: nqt(20)
!***

!***
!     COMMON /DTNORM/NORM(60)
!
integer, save :: norm(60)
!***

!=======================================================================
! EQUIVALENCE (CNORM,NORM)
!=======================================================================

!=======================================================================
! NTYPE SPECIFIES THE FOLLOWING (NLM)=
! 1=(000),2=(100),3=(010),4=(001),5=(110),6=(101),7=(011),
! 8=(200),9=(020),10=(002),
! 11=(111),12=(300),13=(120),14=(102),15=(030),16=(210)
! 17=(012),18=(003),19=(201),20=(021)
!
!  THE VALUE OF PI ON THE FOLLOWING DATA CARD MAY HAVE TO BE CHANGED
!    IF THE COMPILER WILL NOT ACCEPT EXTRA SIGNIFICANT FIGURES
!=======================================================================
DATA IFIRST/0/,PI/3.141592653589793238462D0/
!DATA PI/3.141592653589793238462D0/
DATA NAP /0, 1, 0, 0, 1, 1, 0, 2, 0, 0, 1, 3, 1, 1, 0, 2, 0, 0, &
      2, 0, 0, 2, 0, 3, 1, 1, 3, 1, 1, 4, 2, 0, 0, 0, 2, 4, 24*0/
DATA LAP /0, 0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 0, 3, 1, 1, 0, &
      0, 2, 3, 1, 1, 0, 2, 0, 1, 1, 3, 0, 0, 0, 2, 4, 2, 0, 24*0/
DATA MAP /0, 0, 0, 1, 0, 1, 1, 0, 0, 2, 1, 0, 0, 2, 0, 0, 2, 3, &
      1, 1, 1, 1, 3, 1, 1, 3, 0, 2, 0, 0, 2, 4, 2, 0, 0, 0, 24*0/
DATA CTERM / 1.D0,-1.D0,3.D0,-1.D0,1.D0,-6.D0,1.D0,-1.D0,-1.D0, &
      2.D0,-3.D0,-3.D0,18*0.D0/
DATA IEND /1,2,3,4,5,6,7,8,9,10,9,10,11,12,13,14,15,16,17,18,19,20 &
         ,14,14,17,17,20,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 &
         ,22,23,25,26,28,29,32,34,36, 8*0/
DATA IENTR /1,2,3,4,5,6,7,8,9,10,8,8,11,12,13,14,15,16,17,18,19,20 &
         ,12,13,15,16,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 &
          ,21,22,24,25,27,28,30,32,34, 8*0/
!=======================================================================
!     JENTR CONTAINS THE INDEX OF 1ST COEFFICIENT (CTERM) IN THE FORMULA
!     MINUS THE VALUE OF THE CORRESPONDING IENTR,(SAVING A SUBTRACTION).
!=======================================================================
DATA JENTR/0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-7,0,-10,-11,-12,-13,-14, &
       -15,-16,-17,-18,-19,-2,-12,-5,-15,-8,-18,-20,-21,-22,-23,-24,-25, &
        -26,-27,-28,-29,-30,-31,-32,-33,-34,2*-19,2*-22,3*-25,-27,-29 &
        ,8*0/
DATA LL /0, 3*1, 8*2, 16*3, 24*4, 8*0/
DATA NMAX /0,1,0,0,1,1,0,2,0,0,2,2,1,3,1,1,0,2,0,0,2,0,3,1,2, &
          2,2,2,0,2,0,3,1,1,3,1,1,4,2,0,0,0,2,2,2,3,1,3,1,4,0,4, 8*0/
DATA LMAX /0,0,1,0,1,0,1,0,2,0,2,2,1,0,2,0,3,1,1,0,0,2,2,2,3, &
          1,2,2,3,1,1,0,2,0,1,1,3,0,0,0,2,4,2,3,1,2,2,1,3,0,4,4, 8*0/
DATA MMAX /0,0,0,1,0,1,1,0,0,2,0,2,1,0,0,2,0,0,2,3,1,1,2,2,2, &
          2,3,1,1,1,3,1,1,3,0,2,0,0,2,4,2,0,0,1,3,1,3,2,2,4,4,0, 8*0/
DATA NORM/7*1, 3*3, 4,12,1,15,3,3,15,3,3,15,3,3,60,4,60,4,60,4,15, &
         3,15,15,3,15,15,3,15,105,9,105,9,105,9, 6*24, 3*192, 8*0/
!=======================================================================
! NQT IS USED FOR CORRESPONDANCE TO POLYATOM DEFINITION ON (L,M,N)
!=======================================================================
 DATA NQT/1,2,3,4,8,9,10,5,6,7,12,15,18,16,19,13,20,14,17,11/

! Related to ab overlap
real*8, allocatable, dimension(:), save :: ap, ex
real*8, allocatable, dimension(:,:), save :: px, dab_2d
integer, allocatable, dimension(:), save :: lp_2d, isymab_2d, masab_2d
integer, allocatable, dimension(:), save :: nmp_2d, lmp_2d, mmp_2d, llp
integer, allocatable, dimension(:,:), save :: nlmp_2d, map_ab
integer, save :: n_ab
integer, save :: n_prim2

! Related to cd overlap
complex*16, allocatable, dimension(:), save :: xcxd
complex*16, allocatable, dimension(:,:), save :: qx, dcd_2d
real*8, allocatable, dimension(:), save :: patho_2d, alpc_1d
integer, allocatable, dimension(:), save :: lq_2d, isymcd_2d, mascd_2d, ntypc_1d
integer, allocatable, dimension(:,:), save :: nlmq_2d, map_cd
integer, save :: n_cd

contains

!----------------------------------------------------------------
subroutine setfn

implicit none

!***
!  Local variables
!
integer :: i
!***
 
!***
! The common blocks are now declared in the module head. block data mpdata
! is obsolete.
!     EXTERNAL MPDATA
!***
PI5HF2= DSQRT(DSQRT((2.D0/PI)**3))
do I=1,52
  CNORM(I)=PI5HF2*(2.D0**LL(I))/ DSQRT(DFLOAT(NORM(I)))
end do
TWOPI=2.D0*PI
PI3HAF=PI*DSQRT(PI)
PI5HF2=TWOPI*PI3HAF
PIQUART=TWOPI/DSQRT(PI5HF2)
IFIRST=1

end subroutine setfn
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine int2ns2(s,i_ab,i_cd)

implicit none

!=======================================================================
!     SUBROUTINE TO GIVE INTEGRALS OF WHICH INTEGRAND INCLUDES
!     PLANE WAVE(S)
!=======================================================================
complex*16, intent(out) :: s
integer, intent(in) :: i_ab, i_cd

real*8 :: alpc
integer :: ntypc

!***
!  Local variables
!
!integer :: i , m, n, mn, nla, mfi, ifi, nfi, mla, ila
integer :: l, ii, mm !, ic, ik, j, mk, nk, nn, mn0, iz
!integer :: i0, m0, n0

!=======================================================================
! ELECTRON REPULSION INTEGRAL FOR ONE PLANE WAVE AND THREE GAUSSIANS
! [(NTYPA,ALPA,XA,YA,ZA)*(KX,KY,KZ)*(1/R(1,2))*(NTYPC,ALPC,XC,YX,ZC)*(NTYPD,ALPD,XD,YD,ZD)]
! (KX,KY,KZ)=(XB,YB,ZB)
!=======================================================================

real*8  :: patho
complex*16 :: pathos
complex*16 :: rn(10), rnlm(4625)
integer :: nm, lm
complex*16 :: pqx(3)

integer :: jj
integer :: np, nq
!fk   integer :: mas, nnnn
real*8  :: apq, at, aa
complex*16 :: ac, pf, t

!=======================================================================
!     P=(KX,KY,KZ) PLANE WAVE
!=======================================================================
!     NOTE: THIS PROGRAM USES THE ORIGINAL SUBROUTINES (THANKS ARE
!     DUE TO E.R.DAVIDSON AND T.RESCIGNO) WHICH DO NOT TAKE
!     ACCOUNT OF CONTRACTED BASIS SET EXPLICITLY. THEREFORE,
!     THIS PROGRAM CAN BE MADE FASTER BY THAT MODIFICATION.
!=======================================================================

!     IF(KEYPLR.GT.0)THEN
!     CALL CRD2(CPR,nbfns,max,nbfns,MAX,NTAPE)
!     GO TO 143
!     ENDIF

!=======================================================================
! ONE PLANE WAVE + THREE GAUSSIANS (CPR)
!=======================================================================

PQX(:)=PX(:,i_ab)+QX(:,i_cd)

!fk  Since plane waves never coincide with any plane of symmetry
!    the following part can be skipped, which is related to the
!    symmetry of the integral
!
!     MAS=0
!     IF(pqx(1).eq.0d0)MAS=4
!     IF(pqx(2).eq.0d0)MAS=MAS+2
!     IF(pqx(3).eq.0d0)MAS=MAS+1
!     MAS=IAND(MAS,IAND(MASAB,MASC))
!     MAS=IAND(MAS,IAND(MASAB_2d(i_ab),MASCD_2d(i_cd)))
!     NNNN=IOR(IAND(ISYMAB_2d(i_ab),NOT(ISYMCD_2d(i_cd))),IAND(ISYMCD_2d(i_cd),NOT(ISYMAB_2d(i_ab))))
!     IF(IAND(NNNN,MAS).NE.0) GO TO 210

!=======================================================================
!     ABCD INTEGRAL
!=======================================================================

alpc = alpc_1d(i_cd)
APQ = AP(i_ab) + ALPC
aa = ap(i_ab) * alpc
at = aa / apq
T = AT * sum( pqx(:)*pqx(:) )
PF = EX(i_ab) * xcxd(i_cd) / (aa*DSQRT(APQ))

patho = patho_2d(i_cd)
PATHOS = patho - T
ntypc = ntypc_1d(i_cd)
L = LLP(i_ab) + LL(NTYPC)

! Make it up for the removed EQUIVALENCE statment: RNLM(1906), RN(1)
call formgn(L,T,AT,RN,PATHOS,PATHO)
!     rnlm(1906:1915)=rn(1:10)

IF(L.EQ.0) GO TO 200

rnlm(1906:1915)=rn(1:10)
!fk   rnlm(1906:1912)=rn(1:7)

NM=NMP_2d(i_ab)+NMAX(NTYPC)+1
LM=LMP_2d(i_ab)+LMAX(NTYPC)+1
MM=MMP_2d(i_ab)+MMAX(NTYPC)+1
L=L+2

call getr(L,NM,LM,MM,PQX,rnlm)

if(LP_2d(i_ab).GT.LQ_2d(i_cd)) go to 170

  S=0.D0
  do ii=1,LP_2d(i_ab)
    AC=0.D0
    np = nlmp_2d(ii,i_ab)
    do jj=1,LQ_2d(i_cd)
      nq = nlmq_2d(jj,i_cd)
      AC=AC+DCD_2d(jj,i_cd)*RNLM(NP+NQ)
    end do
    S=S+DAB_2d(ii,i_ab)*AC
  end do
  S=S*PF
  return

170  continue
  S=0.D0
  do jj=1,LQ_2d(i_cd)
    AC=0.D0
    nq = nlmq_2d(jj,i_cd)
    do ii=1,LP_2d(i_ab)
      np = nlmp_2d(ii,i_ab)
      AC=AC+DAB_2d(ii,i_ab)*RNLM(NP+NQ)
    end do
    S=S+DCD_2d(jj,i_cd)*AC
  end do
  S=S*PF
  return

200    S=PF*RN(1)
RETURN

!210   S=0.D0
!RETURN

end subroutine int2ns2
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine formgn(L,T,A,R,PATHOS,PATHO)

implicit none

integer, intent(in) :: l
complex*16, intent(in) :: t
complex*16, intent(out) :: r(10)
real*8, intent(in) :: a, patho
complex*16, intent(in):: pathos

!***
!  Local variables
!
integer :: i
real*8  :: b, s
!***

!     CALL FINT(L,T,R)
call fint(L,T,R,PATHOS,PATHO)

if(l.eq.0) return

B=-2.D0*A
S=B
DO 1 I=1,L
R(I+1)=R(I+1)*S
1     S=S*B

end subroutine formgn
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine fint(ITOP,X,FVEC,PATHOS,PATHO)

implicit none

integer, intent(in) :: itop
complex*16, intent(in)  :: x
complex*16, intent(out) :: fvec(10)
complex*16, intent(in)  :: pathos
real*8, intent(in):: patho

!     COMMON/OVRFLO/PATHOS,PATHO

!     COMPLEX*16    X,Y,TERM,PTLSUM,FVEC    ,X2  ,ARG ,DCERFC,PATHOS
!=======================================================================
!    SUBROUTINE TO COMPUTE INTEGRALS FROM 0 TO 1 OF
!    (T**2M) * EXP(-X*(T**2)).STORES VALUES FOR M=0 TO ITOP IN FVEC.
!=======================================================================

real*8  :: sqtpi2
DATA SQTPI2/.886226925452758E0/

!***
!  Local variables
!
integer :: i, m
real*8  :: a
complex*16 :: y, term, ptlsum, x2, arg
!***
!fk
complex*16 :: x2i
!real*8 :: termr, termi, ptlsumr, ptlsumi

!RVX$ VAX
Y = CDEXP(PATHOS)
!RVX$ END
!RVX$ CRAY
!     Y = CEXP(PATHOS)
!RVX$ END
X2 = 2.D00*X
!RVX$ VAX
IF (CDABS(X)-10.D00) 10,10,20
!RVX$ END
!RVX$ CRAY
!     IF (CABS(X)-10.D00) 10,10,20
!RVX$ END
!=======================================================================
!
!    SMALL ARGUMENT -- USE ANALYTIC SERIES EXPANSION OF INCOMPLETE
!    GAMMA FUNCTION
!
!=======================================================================
  10  A = ITOP
 A=A+0.5D0
TERM=1.D00/A
 PTLSUM=TERM
 DO 11 I=2,50
 A=A+1.D00
 TERM=TERM*X/A
 PTLSUM=PTLSUM+TERM
!CRVX$ VAX
! IF (CDABS(PTLSUM).EQ.0.0D0)GO TO 11
! IF (CDABS(TERM/PTLSUM)-1.D-12) 12,11,11
!CRVX$ END
!CRVX$ CRAY
!CIF (CABS(PTLSUM).EQ.0.0D0)GO TO 11
!CIF (CABS(TERM/PTLSUM)-1.D-12) 12,11,11
!CRVX$ END
!RVX$ VAX
!fk
!     IF((DREAL(TERM)*DREAL(TERM)+DIMAG(TERM)*DIMAG(TERM)).LE.
!    1 (DREAL(PTLSUM)*DREAL(PTLSUM)+DIMAG(PTLSUM)*DIMAG(PTLSUM))
!    2 *1.D-24) GO TO 12
!termr = dreal(term)
!termi = dimag(term)
!ptlsumr = dreal(ptlsum)
!ptlsumi = dimag(ptlsum)
!IF((termr**2+termi**2).LE.(ptlsumr**2+ptlsumi**2)*1.D-24) GO TO 12
IF((DREAL(TERM)**2+DIMAG(TERM)**2).LE.(DREAL(PTLSUM)**2+DIMAG(PTLSUM)**2) *1.D-24) GO TO 12
!fk
!RVX$ END
!RVX$ CRAY
!     IF(( REAL(TERM)* REAL(TERM)+AIMAG(TERM)*AIMAG(TERM)).LE.
!    1 ( REAL(PTLSUM)* REAL(PTLSUM)+AIMAG(PTLSUM)*AIMAG(PTLSUM))
!    2 *1.D-24) GO TO 12
!RVX$ END
  11  CONTINUE
WRITE(6,999) M,X
  999 FORMAT(27H0CONVERGENCE FAILED IN FINT,I5,2D16.8)
CALL EXIT (1)
   12 FVEC(ITOP+1) = 0.5D00*PTLSUM*Y
!IF(ITOP.EQ.0) GO TO 15
IF(ITOP.EQ.0) return
!=======================================================================
!    BACKWARDS RECURSION
!=======================================================================
M = ITOP+1
  14  FVEC(M-1) = (X2*FVEC(M)+Y)/ DFLOAT(2*M-3)
M = M-1
IF(M-1) 15,15,14
  15  RETURN
!=======================================================================
!    LARGE ARGUMENT - FORWARD RECURSION FROM ERROR FUNCTION
!=======================================================================
  20  CONTINUE
!RVX$ VAX
ARG = CDSQRT(X)
!RVX$ END
!RVX$ CRAY
!     ARG = CSQRT(X)
!RVX$ END
FVEC(1) = SQTPI2/ARG*( DEXP(PATHO )-  Y*DCERFC(ARG))
!IF(ITOP.EQ.0) GO TO 17
IF(ITOP.EQ.0) return
I=1
x2i = 1d0/x2
! 16  FVEC(I+1) = ( DFLOAT(2*I-1)*FVEC(I)-Y)/X2
  16  FVEC(I+1) = ( DFLOAT(2*I-1)*FVEC(I)-Y)*x2i
I=I+1
IF(I-ITOP-1) 16,16,17
  17  RETURN

end subroutine fint
!----------------------------------------------------------------

!----------------------------------------------------------------
DOUBLE COMPLEX FUNCTION DCERFC(ZP)

implicit none

complex*16 :: zp

!***
!  Local variables
!
logical :: b
integer :: iq, k, nu, n, np1
real*8  :: lambda, x, y, s, h, h2, r1, r2, s1, s2, t1
real*8  :: t2, c, re, aim
complex*16 :: ai, z, za, cerf
!***

AI=(0.D00,1.D00)
Z=AI*ZP
!RVX$ VAX
X= DREAL(Z)
Y=DIMAG(Z)
!RVX$ END
!RVX$ CRAY
!     X= REAL(Z)
!     Y=AIMAG(Z)
!RVX$ END
IQ=1
IF (X .LT. 0.D00) IQ=2
IF (Y .LT. 0.D00) IQ=3
IF (Y .LT. 0.D00 .AND. X  .GE. 0.D00) IQ=4
X=DABS(X)
Y=DABS(Y)
!RVX$ VAX
ZA=DCMPLX(X,Y)
!RVX$ END
!RVX$ CRAY
!     ZA=CMPLX(X,Y)
!RVX$ END
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
!=======================================================================
!
!  20 IF(H .GT. 0.D0D0) LAMBDA=H2**KTHIS CARD IS DELETED
!
!=======================================================================
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
!ALPHA
RE = 0.0D0
IF (X.LT.20.D0) RE=DEXP(-X*X)/1.12837916709551D0
!APPN RE=DEXP(-X*X)/1.12837916709551D0

!ALPHA
GO TO 70
   60 RE=S1
IF (B) RE=R1
   70 AIM=S2
IF (B) AIM=R2
!RVX$ VAX
CERF=1.12837916709551D0* DCMPLX(RE,AIM)
!RVX$ END
!RVX$ CRAY
!     CERF=1.12837916709551D0* CMPLX(RE,AIM)
!RVX$ END
IF (IQ .EQ. 1) GO TO 90
IF (IQ .NE. 2) GO TO 80
!RVX$ VAX
CERF=DCONJG(CERF)
!RVX$ END
!RVX$ CRAY
!     CERF=CONJG(CERF)
!RVX$ END
GO TO 90
!RVX$ VAX
   80 CERF=2.D00*CDEXP(-ZA*ZA)-CERF
IF (IQ .EQ. 3) GO TO 90
CERF=DCONJG(CERF)
!RVX$ END
!RVX$ CRAY
!  80 CERF=2.D00*CEXP(-ZA*ZA)-CERF
!     IF (IQ .EQ. 3) GO TO 90
!     CERF=CONJG(CERF)
!RVX$ END
!=======================================================================
!  EXP(Z**2) TAKEN OUT OF CERFC 4/22/74 -APPEARS IN FINT
!  90 DCERFC=CDEXP(Z*Z)*CERF
!=======================================================================
  90   DCERFC = CERF
RETURN

end function dcerfc
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine getr0(LL,NM,LM,MM,PQX,MAS,rnlm)

implicit none

integer :: ll, nm, lm, mm, mas
complex*16 :: pqx(3)
!     COMMON /CALLR/LL,NM,LM,MM,PQX,PQY,PQZ,MAS

!=======================================================================
!     GENERATE R(N+1,L+1,M+1) BY RECURRENCE FOR N+1=1 TO NM,
!     L+1=1 TO LM, M+1=1 TO MM AND (N+1)+(L+1)+(M+1) SMALLER THAN LL+2
!     CALLING PARAMETERS ARE PASSED IN COMMON /CALLR/
!     PQX,PQY,PQZ ARE COORDINATES OF VECTOR PQ
!     MASX IS 1 IF R IS ZERO BY SYMMETRY FOR ODD N, QTHERWISE MASX=0
!     MASY IS 1 IF R IS ZERO BY SYMMETRY FOR ODD L, OTHERWISE MASY=0
!     MASZ IS 1 IF R IS ZERO BY SYMMETRY FOR ODD M, OTHERWISE MASZ=0
!     ZERO VALUES OF R BY SYMMETRY WILL NOT BE COMPUTED OR STORED.
!     REGION RN CONTAINS MODIFIED FN INTEGRALS FROM FORMGN FOR
!     N+1=1 TO LL-1 AT ENTRANCE. RN,RL,RM ARE OVERWRITTEN
!     DURING EXECUTION.
!=======================================================================

!     COMMON/ERN/RNLM(17,17,16),RDUM
complex*16 :: rnlm(4625)

complex*16, dimension(17,2) :: rn, rl, rm
!     EQUIVALENCE (RL(1,1),RNLM(2,13,7)),(RM(1,1),RNLM(2,15,7)),
!    $(RN(1,1),RNLM(2,11,7))
!
!     The following mapping will make up for the removed EQUIVALENCE statment:
!   
!     RN: (2,11,7) -> 17*17*6 + 17*10 + 2 = 1906 => RNLM(2,11,7)
!     RL: (2,13,7) -> 17*17*6 + 17*12 + 2 = 1940 => RNLM(2,13,7)
!     RM: (2,15,7) -> 17*17*6 + 17*14 + 2 = 1974 => RNLM(2,13,7)
!
!  INT2NS   GETR    GETR
!     RNLM(1906:1939)    RNLM(2,11:12,7)     RN
!     RNLM(1940:1973)    RNLM(2,13:14,7)     RL
!     RNLM(1974:2007)    RNLM(2,15:16,7)     RM
!   

real*8 :: fact(9)
data fact/1.D0,1.D0,3.D0,15.D0,105.D0,945.D0,10395.D0,135135.D0,2027025.D0/

!***
!  Local variables
!
integer :: masx, masy, masz, ninc, nco, linc, n, jf1, j, lco
integer :: lf, lf1, l, jf2, mf, mf1, jf3, k, m, imln
real*8  :: fn, fl, fm
complex*16 :: arl, arm, arn 
!***

!***
! Missing EQUIVALENCE statment:
!
rn(:,1) = rnlm(1906:1922)
rn(:,2) = rnlm(1923:1939)
rl(:,1) = rnlm(1940:1956)
rl(:,2) = rnlm(1957:1973)
rm(:,1) = rnlm(1974:1990)
rm(:,2) = rnlm(1991:2007)
!***

 MASX=IAND(MAS,4)
 MASY=IAND(MAS,2)
 MASZ=IAND(MAS,1)

IF(MASX.EQ.0) GO TO 2
1 NINC=2
NCO=0
GO TO 3
2 NINC=1
3 IF(MASY.EQ.0) GO TO 5
4     LINC=2
GO TO 10
5     LINC=1
10 CONTINUE

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
!***
!     print*,'getr1: n,l,1',n,l,1
!     RNLM(N,L,1)=RM(1,1)
imln = (l-1)*17 + n
rnlm(imln) = rm(1,1)
!***

IF(MF.LT.2) GO TO 580
JF3=MF1-2
DO 540 J=1,JF3
RM(J,2)=PQX(3)*RM(J+1,  1)
  540 CONTINUE

!***
!     print*,'getr2: n,l,2',n,l,2
!     RNLM(N,L,2)=RM(1,2)
imln = 17*17 + (l-1)*17 + n
rnlm(imln) = rm(1,2)
!***

IF(MF.LT.3) GO TO 580
DO 570 M=3,MF
JF3=MF1-M
FM=DFLOAT(M-2)
DO 560 J=1,JF3
ARM=PQX(3)*RM(J+1,2)+FM*RM(J+1,1)
RM(J+1,1)=RM(J+1,2)
RM(J,2)=ARM
  560 CONTINUE

!***
!     print*,'getr3: n,l,m',n,l,m
!     RNLM(N,L,M)=RM(1,2)
imln = (m-1)*17*17 + (l-1)*17 + n
rnlm(imln) = rm(1,2)
!***

  570 CONTINUE
GO TO 580
541    K=0
DO 571 M=1,MF,2
K=K+1

!***
!     if(m.gt.7) print*,'passou_m:', m
!     if((m.eq.2).or.(m.eq.6)) print*,'achou_m:', m
!     m=1,3,7 apenas
!     if(l.gt.3) print*,'passou_l:', l
!     l=1,2,3 apenas
!     if(n.gt.7) print*,'passou_n:', n
!     n=1,2,3,4,5,6,7 apenas
!     print*,'getr4: n,l,m',n,l,m
!571   RNLM(N,L,M)=RM(K,1)*FACT(K)
imln = (m-1)*17*17 + (l-1)*17 + n
571   rnlm(imln) = rm(k,1) * fact(k)   
!***

  580 CONTINUE
  590 CONTINUE

!***
! Missing EQUIVALENCE statment:
!
rnlm(1906:1922) = rn(:,1)
rnlm(1923:1939) = rn(:,2)
rnlm(1940:1956) = rl(:,1)
rnlm(1957:1973) = rl(:,2)
rnlm(1974:1990) = rm(:,1)
rnlm(1991:2007) = rm(:,2)
!***

end subroutine getr0
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine getr(LL,NM,LM,MM,PQX,rnlm)

implicit none

integer, intent(in) :: ll, nm, lm, mm !, mas
complex*16, intent(in) :: pqx(3)
!     COMMON /CALLR/LL,NM,LM,MM,PQX,PQY,PQZ,MAS

!=======================================================================
!     GENERATE R(N+1,L+1,M+1) BY RECURRENCE FOR N+1=1 TO NM,
!     L+1=1 TO LM, M+1=1 TO MM AND (N+1)+(L+1)+(M+1) SMALLER THAN LL+2
!     CALLING PARAMETERS ARE PASSED IN COMMON /CALLR/
!     PQX,PQY,PQZ ARE COORDINATES OF VECTOR PQ
!     MASX IS 1 IF R IS ZERO BY SYMMETRY FOR ODD N, QTHERWISE MASX=0
!     MASY IS 1 IF R IS ZERO BY SYMMETRY FOR ODD L, OTHERWISE MASY=0
!     MASZ IS 1 IF R IS ZERO BY SYMMETRY FOR ODD M, OTHERWISE MASZ=0
!     ZERO VALUES OF R BY SYMMETRY WILL NOT BE COMPUTED OR STORED.
!     REGION RN CONTAINS MODIFIED FN INTEGRALS FROM FORMGN FOR
!     N+1=1 TO LL-1 AT ENTRANCE. RN,RL,RM ARE OVERWRITTEN
!     DURING EXECUTION.
!=======================================================================

!     complex*16 :: rnlm(4625)
complex*16, intent(out) :: rnlm(4625)

complex*16, dimension(17,2) :: rn, rl, rm
!     EQUIVALENCE (RL(1,1),RNLM(2,13,7)),(RM(1,1),RNLM(2,15,7)),
!    $(RN(1,1),RNLM(2,11,7))
!
!     The following mapping will make up for the removed EQUIVALENCE statment:
!     
!   
!     RN: (2,11,7) -> 17*17*6 + 17*10 + 2 = 1906 => RNLM(2,11,7)
!     RL: (2,13,7) -> 17*17*6 + 17*12 + 2 = 1940 => RNLM(2,13,7)
!     RM: (2,15,7) -> 17*17*6 + 17*14 + 2 = 1974 => RNLM(2,13,7)
!
!  INT2NS   GETR    GETR
!     RNLM(1906:1939)    RNLM(2,11:12,7)     RN
!     RNLM(1940:1973)    RNLM(2,13:14,7)     RL
!     RNLM(1974:2007)    RNLM(2,15:16,7)     RM
!   

real*8  :: fact(9)
data fact/1.D0,1.D0,3.D0,15.D0,105.D0,945.D0,10395.D0,135135.D0,2027025.D0/

!***
!  Local variables
!
!integer :: masx, masy, masz, ninc, nco, linc, n, jf1, j, lco
integer :: ninc, nco, linc, n, jf1, j, lco
integer :: lf, lf1, l, jf2, mf, mf1, jf3, k, m, imln
real*8  :: fn, fl, fm
complex*16 :: arl, arm, arn 
!***

!***
! Missing EQUIVALENCE statment:
!
rn(:,1) = rnlm(1906:1922)
rn(:,2) = rnlm(1923:1939)
rl(:,1) = rnlm(1940:1956)
rl(:,2) = rnlm(1957:1973)
rm(:,1) = rnlm(1974:1990)
rm(:,2) = rnlm(1991:2007)
!***

!fk    MASX=IAND(MAS,4)
!MASY=IAND(MAS,2)
!MASZ=IAND(MAS,1)
!     masx=0
!     masy=0
!     masz=0

!fk   IF(MASX.EQ.0) GO TO 2
go to 2
!fk
1 NINC=2
NCO=0
!fk   GO TO 3
2 NINC=1
!fk3IF(MASY.EQ.0) GO TO 5
go to 5
!fk
4     LINC=2
GO TO 10
5     LINC=1
10 CONTINUE

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
!fk   IF(MASX.EQ.0) GO TO 421
go to 421
!fk
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
!fk   IF(MASY.EQ.0) GO TO 491
go to 491
!fk
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
!fk   IF(MASZ.NE.0) GO TO 541
!fk
MF1=MF+1
!***
!     print*,'getr1: n,l,1',n,l,1
!     RNLM(N,L,1)=RM(1,1)
imln = (l-1)*17 + n
rnlm(imln) = rm(1,1)
!***

IF(MF.LT.2) GO TO 580
JF3=MF1-2
DO 540 J=1,JF3
RM(J,2)=PQX(3)*RM(J+1,  1)
  540 CONTINUE

!***
!     print*,'getr2: n,l,2',n,l,2
!     RNLM(N,L,2)=RM(1,2)
imln = 17*17 + (l-1)*17 + n
rnlm(imln) = rm(1,2)
!***

IF(MF.LT.3) GO TO 580
DO 570 M=3,MF
JF3=MF1-M
FM=DFLOAT(M-2)
DO 560 J=1,JF3
ARM=PQX(3)*RM(J+1,2)+FM*RM(J+1,1)
RM(J+1,1)=RM(J+1,2)
RM(J,2)=ARM
  560 CONTINUE

!***
!     print*,'getr3: n,l,m',n,l,m
!     RNLM(N,L,M)=RM(1,2)
imln = (m-1)*17*17 + (l-1)*17 + n
rnlm(imln) = rm(1,2)
!***

  570 CONTINUE
GO TO 580
541    K=0
DO 571 M=1,MF,2
K=K+1

!***
!     if(m.gt.7) print*,'passou_m:', m
!     if((m.eq.2).or.(m.eq.6)) print*,'achou_m:', m
!     m=1,3,7 apenas
!     if(l.gt.3) print*,'passou_l:', l
!     l=1,2,3 apenas
!     if(n.gt.7) print*,'passou_n:', n
!     n=1,2,3,4,5,6,7 apenas
!     print*,'getr4: n,l,m',n,l,m
!571   RNLM(N,L,M)=RM(K,1)*FACT(K)
imln = (m-1)*17*17 + (l-1)*17 + n
571   rnlm(imln) = rm(k,1) * fact(k)   
!***

  580 CONTINUE
  590 CONTINUE

!***
! Missing EQUIVALENCE statment:
!
rnlm(1906:1922) = rn(:,1)
rnlm(1923:1939) = rn(:,2)
rnlm(1940:1956) = rl(:,1)
rnlm(1957:1973) = rl(:,2)
rnlm(1974:1990) = rm(:,1)
rnlm(1991:2007) = rm(:,2)
!***

!     RETURN
!     END
end subroutine getr
!----------------------------------------------------------------

!----------------------------------------------------------------
!subroutine expldp(D,NLM,K,IT,XA,ALPA,NTYPA,NTYPB,XB,ALPB)
subroutine expldp(D,NLM,K,IT,ALPA,NTYPA,XB)

implicit none

complex*16, intent(out) :: d(8)
integer, intent(out) :: nlm(8)
integer :: k, it
!real*8, intent(in)  :: xa(3), alpa, xb(3), alpb
real*8, intent(in)  :: alpa, xb(3)
integer, intent(in) :: ntypa !, ntypb
!     COMMON /CALLIN/XA,YA,ZA,ALPA,NTYPA,NTYPB,XB,YB,ZB,ALPB,
!    * XC,YC,ZC,ALPC,NTYPC,NTYPD,XD,YD,ZD,ALPD

!=======================================================================
! EVALUATE S AND P FUNCTIONS AND EXPAND GAUSSIANS IN HERMITE
! POLINOMIALS ABOUT A COMPLEX ORIGIN.
!=======================================================================

!***
!  Local variables
!
integer, dimension(4) :: ix, iy, iz
complex*16, dimension(4) :: dx, dy, dz

integer :: mas, nty, i, n, kx, ky, kz, ixx, iyy, izz
real*8  :: a
complex*16 :: xab(3)
!***

MAS=7
!fk   IF(IT.EQ.2) GO TO 4
A=.5D0/ALPA
NTY=NTYPA
!fk   IF(XB.EQ.0.D0)GO TO 1
IF(XB(1).EQ.0.D0)GO TO 1
!RVX$ VAX
!fk   XAB=DCMPLX(0.D0,A*XB)
XAB(1)=DCMPLX(0.D0,A*XB(1))
!RVX$ END
!RVX$ CRAY
!     XAB=CMPLX(0.D0,A*XB)
!RVX$ END
MAS=IAND(MAS,3)
!=======================================================================
!  IAND(,) CORESPONDS TO  .INT. (INTERSECTION = BITWISE AND)
! IOR(,)  CORESPONDS TO  .UN.  (UNION = BITWISE OR)
! IAND,IOR ARE FOR VAX-11
! INT,UN   ARE FOR LRLTRAN
!=======================================================================
!fk1  IF(YB.EQ.0.D0) GO TO 2
1     IF(XB(2).EQ.0.D0) GO TO 2
!RVX$ VAX
XAB(2)=DCMPLX(0.D0,A*XB(2))
!RVX$ END
!RVX$ CRAY
!     YAB=CMPLX(0.D0,A*YB)
!RVX$ END
MAS=IAND(MAS,5)
!fk2  IF(ZB.EQ.0.D0) GO TO 8
2     IF(XB(3).EQ.0.D0) GO TO 8
!RVX$ VAX
!fk   ZAB=DCMPLX(0.D0,A*ZB)
XAB(3)=DCMPLX(0.D0,A*XB(3))
!RVX$ END
!RVX$ CRAY
!     ZAB=CMPLX(0.D0,A*ZB)
!RVX$ END
MAS=IAND(MAS,6)
8     IF(IT.EQ.1)GO TO 7
xab = - xab
GO TO 7
!fk4  A=.5D0/ALPC
!fk   NTY=NTYPC
!fk   IF(XD.EQ.0.D0) GO TO 5
!fk   IF(XD(1).EQ.0.D0) GO TO 5
!RVX$ VAX
!fk   XAB=DCMPLX(0.D0,-A*XD)
!fk   XAB(1)=DCMPLX(0.D0,-A*XD(1))
!RVX$ END
!RVX$ CRAY
!     XAB=CMPLX(0.D0,-A*XD)
!RVX$ END
!fk   MAS=IAND(MAS,3)
!fk5  IF(YD.EQ.0.D0) GO TO 6
!fk5  IF(XD(2).EQ.0.D0) GO TO 6
!RVX$ VAX
!fk   YAB=DCMPLX(0.D0,-A*YD)
!fk   XAB(2)=DCMPLX(0.D0,-A*XD(2))
!RVX$ END
!RVX$ CRAY
!     YAB=CMPLX(0.D0,-A*YD)
!RVX$ END
!fk   MAS=IAND(MAS,5)
!fk6  IF(ZD.EQ.0.D0) GO TO 7
!fk6  IF(XD(3).EQ.0.D0) GO TO 7
!RVX$ VAX
!fk   ZAB=DCMPLX(0.D0,-A*ZD)
!fk   XAB(3)=DCMPLX(0.D0,-A*XD(3))
!RVX$ END
!RVX$ CRAY
!     ZAB=CMPLX(0.D0,-A*ZD)
!RVX$ END
!fk   MAS=IAND(MAS,6)
7     CONTINUE
IF(NTY.GT.4) GO TO 120
!=======================================================================
! CHECK INTEGRALS WITH EXPLICIT FORMS
!=======================================================================
GO TO(10,20,30,40)NTY
10    D(1)=1.D0
NLM(1)=1
K=1
RETURN
!=======================================================================
! TYPE #1 N=L=M=0
!=======================================================================
20    D(1)=A
NLM(1)=2
K=1
IF(IAND(MAS,4).NE.0) RETURN
!=======================================================================
! TYPE #2 N=1,L=M=0
!=======================================================================
!fk   D(2)=XAB
D(2)=XAB(1)
NLM(2)=1
K=2
RETURN
!=======================================================================
! TYPE #3 L=1,M=N=0
!=======================================================================
30    D(1)=A
NLM(1)=18
K=1
IF(IAND(MAS,2).NE.0) RETURN
D(2)=XAB(2)
NLM(2)=1
K=2
!=======================================================================
! TYPE #4 M=1,L=N=0
!=======================================================================
RETURN
40    D(1)=A
NLM(1)=290
K=1
IF(IAND(MAS,1).NE.0) RETURN
D(2)=XAB(3)
NLM(2)=1
K=2
RETURN
!=======================================================================
! EVALUATE INTEGRALS WITH NON-EXPLICIT FORMS AND EXPAND X**N IN
! HERMITE POLYNOMIALS
!=======================================================================
120   I=IEND(NTY)
N=NAP(I)+1
GO TO (121,122,123,124)N
121   DX(1)=1.D0
IX(1)=1
KX=1
GO TO 130
!=======================================================================
! N=0
!=======================================================================
122   DX(1)=A
IX(1)=2
IF(IAND(MAS,4).NE.0) GO TO 125
!=======================================================================
! COMPUTE N+1,L+1,M+1 FOR EACH GAUSSIAN
!=======================================================================
!fk   DX(2)=XAB
DX(2)=XAB(1)
!=======================================================================
! N=1
!=======================================================================
KX=2
IX(2)=1
GO TO 130
125   KX=1
GO TO 130
!=======================================================================
! N=2
!=======================================================================
123   DX(1)=A*A
IX(1)=3
IF(IAND(MAS,4).NE.0) GO TO 126
!fk   DX(2)=2.D0*XAB*A
DX(2)=2.D0*XAB(1)*A
IX(2)=2
!fk   DX(3)=A+XAB*XAB
DX(3)=A+XAB(1)*XAB(1)
IX(3)=1
KX=3
GO TO 130
!=======================================================================
! N=3
!=======================================================================
126   DX(2)=A
IX(2)=1
KX=2
GO TO 130
124   DX(1)=A*A*A
IX(1)=4
IF(IAND(MAS,4).NE.0) GO TO 127
!fk   DX(2)=3.D0*XAB*A*A
DX(2)=3.D0*XAB(1)*A*A
IX(2)=3
!fk   DX(3)=3.D0*(XAB*XAB+A)*A
DX(3)=3.D0*(XAB(1)*XAB(1)+A)*A
IX(3)=2
!fk   DX(4)=XAB*(XAB*XAB+3.D0*A)
DX(4)=XAB(1)*(XAB(1)*XAB(1)+3.D0*A)
IX(4)=1
KX=4
GO TO 130
127   DX(2)=3.D0*A*A
IX(2)=2
KX=2
130   N=LAP(I)+1
!=======================================================================
! EXPAND Y**L IN HERMITE POLYNOMIALS
!=======================================================================
GO TO (131,132,133,134)N
131   DY(1)=1.D0
IY(1)=0
KY=1
GO TO 140
132   DY(1)=A
IY(1)=17
IF(IAND(MAS,2).NE.0) GO TO 135
DY(2)=XAB(2)
IY(2)=0
KY=2
GO TO 140
135   KY=1
GO TO 140
133   DY(1)=A*A
IY(1)=34
IF(IAND(MAS,2).NE.0) GO TO 136
DY(2)=2.D0*XAB(2)*A
IY(2)=17
DY(3)=A+XAB(2)*XAB(2)
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
DY(2)=3.D0*XAB(2)*A*A
IY(2)=34
DY(3)=3.D0*A*(XAB(2)*XAB(2)+A)
IY(3)=17
DY(4)=XAB(2)*(XAB(2)*XAB(2)+3.D0*A)
IY(4)=0
KY=4
GO TO 140
137   DY(2)=3.D0*A*A
!=======================================================================
! EXPAND Z**M IN HERMITE POLYNOMIALS
!=======================================================================
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
DZ(2)=XAB(3)
IZ(2)=0
KZ=2
GO TO 150
145   KZ=1
GO TO 150
143   DZ(1)=A*A
IZ(1)=578
IF(IAND(MAS,1).NE.0) GO TO 146
DZ(2)=2.D0*XAB(3)*A
IZ(2)=289
DZ(3)=A+XAB(3)*XAB(3)
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
DZ(2)=3.D0*XAB(3)*A*A
IZ(2)=578
DZ(3)=3.D0*A*(XAB(3)*XAB(3)+A)
IZ(3)=289
DZ(4)=XAB(3)*(XAB(3)*XAB(3)+3.D0*A)
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
end subroutine expldp
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine expldd(D,NLM,K,a,pd1,pd4,n)

implicit none

real*8, intent(out)  :: d(105)
integer, intent(out) :: nlm(8)
integer, intent(out) :: k
real*8, intent(in)   :: a, pd1(3), pd4(3)
integer, intent(in)  :: n
!     COMMON/CALLEX/A,P(6),N

!=======================================================================
! SUBROUTINE TO FIND EXPLICIT FORMULAE FOR EXPANDING PRODUCTS OF
! TWO GAUSSIANS IN HERMITE POLYNOMIALS ABOUT A REAL CENTER.
!=======================================================================

real*8  :: p(6)
integer :: m(15), n1(15), n2(11), n3(9), l(9)
DATA M /1, 2, 3, 4*4, 4*5, 4*6/
DATA N1 /1,17,289,1,2,18,290,17,18,34,306,289,290,306,578/
DATA N2 /1, 17, 289, 0, 1, 17, 289, 0, 1, 17, 289/
DATA N3 /1, 1, 0, 17, 0, 17, 0, 289, 289/
DATA L /2, 3, 0, 1, 0, 3, 0, 1, 2/

!***
!  Local variables
integer :: i
!***

!***
! Load array P(1:6). This is necessar because common block CALLEX
! no longer exists.
p(1:3) = pd1(1:3)
p(4:6) = pd4(1:3)
!***

GO TO (10,20,20,20,20,30,40,40,20,40,30,40,20,40,40,30), N
!=======================================================================
!     THIS SUBROUTINE BRANCHES TO 4 DIFFERENT SETS OF FORMULAS FOR D
!=======================================================================
   10 CONTINUE
!=======================================================================
!     S  S
!=======================================================================
D(1)=1.D0
NLM(1)=0
   15 CONTINUE
K=1
RETURN
   20 CONTINUE
!=======================================================================
!     P S OR S P (ALL P)
!=======================================================================
D(1)=A
NLM(1)=N1(N-1)
K=M(N-1)
IF(P(K).EQ.0.D0) GO TO 15
D(2)=P(K)
NLM(2)=0
K=2
RETURN
   30 CONTINUE
!=======================================================================
!     PX PX, PY PY, OR PZ PZ
!=======================================================================
D(1)=A*A
NLM(1)=N1(N-1)
NLM(2)=0
K=M(N-1)
IF(P(K).EQ.0.D0) GO TO 35
D(2)=P(K)*P(K-3)+DABS(A)
D(3)=A*(P(K)+P(K-3))
NLM(3)=N2(N-5)
K=3
RETURN
   35 CONTINUE
D(2)=DABS(A)
K=2
RETURN
   40 CONTINUE
!=======================================================================
!     PY PX, PZ PX, PX PY, PZ PY, PX PZ, OR PY PZ
!=======================================================================
D(1)=A*A
I=L(N-6)
NLM(1)=N1(N-1)
K=M(N-1)
IF(P(K).EQ.0.D0) GO TO 45
D(2)=A*P(K)
NLM(2)=N2(N-5)
IF(P(I).EQ.0.D0) GO TO 50
D(3)=A*P(I)
D(4)=P(K)*P(I)
NLM(3)=N3(N-6)
NLM(4)=0
K=4
RETURN
   45 CONTINUE
IF(P(I).EQ.0.D0) GO TO 15
D(2)=A*P(I)
NLM(2)=N3(N-6)
   50 CONTINUE
K=2
RETURN
!     END
end subroutine expldd
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine XYZ(X,NAM,NBM,A,PAX,PBX)

!=======================================================================
! EXPAND XA**N,XB**N IN HERMITE POLYNOMIALS ABOUT THIRD CENTER
!=======================================================================

implicit none

real*8, intent(out)  :: x(225)
integer, intent(in) :: nam, nbm
real*8, intent(in)  :: a, pax, pbx

!***
!  Local variables
!
integer :: i, n, k, nnu, nn, nup, j, nb, nnup
!***

!
! Modificado para evitar que X viole o limite 100. No codigo original,
! o programador parece estar consciente de que isso aconteceria. A dimensao
! 225 corresponde ao dimensionamento na rotina de chamada (INT2NS), embora
! haja risco potencial relacionado ao equivalence em INT2NS.

!     print*,'X',X
!     print*,'NAM,NBM',NAM,NBM
!     print*,'A,PAX,PBX',A,PAX,PBX
!***

X(1)=1.D0
 IF(PAX.EQ.0.D0) GO TO 110
120 IF(NAM.LT.1) GO TO 10
X(10)=PAX
X(11)=A
IF(NAM.LT.2) GO TO 10
X(19)=A+PAX*PAX
X(20)=2.D0*A*PAX
X(21)=A*A
IF(NAM.LT.3) GO TO 10
X(28)=X(20)+PAX*X(19)
X(29)=2.D0*X(21)+PAX*X(20)+A*X(19)
X(30)=PAX*X(21)+A*X(20)
X(31)=A*X(21)
IF(NAM.LT.4) GO TO  10
X(37)=X(29)+PAX*X(28)
X(38)=2.D0*X(30)+PAX*X(29)+A*X(28)
X(39)=3.D0*X(31)+PAX*X(30)+A*X(29)
X(40)= PAX*X(31)+A*X(30)
X(41)= A*X(31)
10    IF(NBM.LT.1) RETURN
X(46)=PBX
X(47)=A
IF(NAM.LT.1) GO TO 40
X(55)=A+PBX*PAX
X(56)=A*(PAX+PBX)
X(57)=A*A
IF(NAM.LT.2) GO TO 40
I=19
DO 30 N=2,NAM
!     print*,'aqui1: I+45',i+45
!     print*,'aqui1: I+46',i+46
X(I+45)=X(I+1)+PBX*X(I)
K=I+N
!     print*,'aqui2: K+45',K+45
!     print*,'aqui2: K+46',K+46
X(K+46)=A*X(K)
X(K+45)=PBX*X(K)+A*X(K-1)
NNU=N-1
DO 20 NN=1,NNU
K=I+NN
!     print*,'aqui3: K+45',K+45
X(K+45)=DFLOAT(NN+1)*X(K+1)+PBX*X(K)+A*X(K-1)
20    CONTINUE
I=I+9
 30   CONTINUE
40    IF(NBM.LT.2) RETURN
X(91)=A+PBX*PBX
X(92)=2.D0*A*PBX
X(93)=A*A
IF(NAM.LT.1) GO TO 70
I=55
DO 60 N=1,NAM
!     print*,'aqui4: I+45',I+45
X(I+45)=X(I+1)+PBX*X(I)
K=I+N
!     print*,'aqui5: K+46',K+46
!     print*,'aqui5: K+47',K+47
X(K+47)=A*X(K+1)
X(K+46)=PBX*X(K+1)+A*X(K)
DO 50 NN=1,N
K=I+NN
!     print*,'aqui6: K+45',K+45
X(K+45)=DFLOAT(NN+1)*X(K+1)+PBX*X(K)+A*X(K-1)
50    CONTINUE
I=I+9
60    CONTINUE
70    CONTINUE
IF(NBM.LT.3) RETURN
NUP=NAM+1
J=91
DO 100 NB=3,NBM
I=J
DO 90 N=1,NUP
X(I+45)=X(I+1)+PBX*X(I)
NNUP=N+NB-3
K=I+NNUP
X(K+47)=A*X(K+1)
X(K+46)=PBX*X(K+1)+A*X(K)
DO 80 NN=1,NNUP
K=I+NN
X(K+45)=DFLOAT(NN+1)*X(K+1)+PBX*X(K)+A*X(K-1)
80    CONTINUE
90    I=I+9
100   J=J+45
RETURN
110   IF(PBX.NE.0.D0) GO TO 120
 N=NAM+NBM
 IF(N.LT.1) RETURN
 X(11)=A
 IF(N.LT.2) RETURN
 X(19)=A
 X(21)=A*A
 IF(N.LT.3) RETURN
 X(29)=2.D0*X(21)+A*X(19)
 X(31)=A*X(21)
 IF(N.LT.4) RETURN
 X(37)=X(29)
  X(39)=3.D0*X(31)+A*X(29)
 X(41)=A*X(31)
 IF(N.LT.5) RETURN
 X(47)=2.D0*X(39)+A*X(37)
 X(49)=4.D0*X(41)+A*X(39)
 X(51)=A*X(41)
 IF(N.LT.6) RETURN
 X(55)=X(47)
 X(57)=3.D0*X(49)+A*X(47)
  X(59)=5.D0*X(51)+A*X(49)
 X(61)=A*X(51)
 IF(N.LT.7) RETURN
 X(65)=2.D0*X(57)+A*X(55)
 X(67)=4.D0*X(59)+A*X(57)
 X(69)=6.D0*X(61)+A*X(59)
 X(71)=A*X(61)
   IF(N.LT.8) RETURN
 X(73)=X(65)
 X(75)=3.D0*X(67)+A*X(65)
 X(77)=5.D0*X(69)+A*X(67)
 X(79)=7.D0*X(71)+A*X(69)
 X(81)=A*X(71)
 RETURN
end subroutine xyz
!----------------------------------------------------------------

!----------------------------------------------------------------
!subroutine dd(D,NLM,K,ISIGN,MASK,ISYM,X,Y,Z,xa,alpa,ntypa,ntypb,xb,alpb,XC,ALPC,NTYPC,NTYPD,XD,ALPD)
subroutine dd(D,NLM,K,ISIGN,MASK,ISYM,X,Y,Z,ntypa,ntypb,ntypc,ntypd)

implicit none

real*8, intent(out) :: d(105)
!!!!!!!!!!!!!!!!!!!!!!!!!
!! ATTENTION !!!
!!!!!!!!!!!!!!!!!!!!!!!!!
!fk   integer, intent(out) :: nlm(8)
integer, intent(out) :: nlm(105)
integer, intent(out) :: k
integer, intent(in) :: isign, mask, isym
real*8, intent(in) :: x(225), y(225), z(225)
integer, intent(in) :: ntypa, ntypb, ntypc, ntypd
!real*8, intent(in) :: xc(3), alpc, xd(3), alpd
!real*8, intent(in) :: xa(3), alpa, xb(3), alpb

!     COMMON /CALLIN/XA,YA,ZA,ALPA,NTYPA,NTYPB,XB,YB,ZB,ALPB
!    1  ,XC,YC,ZC,ALPC,NTYPC,NTYPD,XD,YD,ZD,ALPD

!=======================================================================
! MULTIPLY EXPANSIONS OF X,Y, AND Z TO OBTAIN FINAL POLYNOMIALS
!=======================================================================

!=======================================================================
! THIS STATEMENT REMOVED:X,Y,Z-MATRICES PASSED THRU ARGUMENT. 
! TLG  4/19/83.
!     COMMON /ERN/ X(225),Y(225),Z(225),R(3950)
!=======================================================================

!***
!  Local variables
!
integer :: n1, n2, i1, i1f, j1, i2, i2f, nf, lf, mf, nls
integer :: ninc, nfac, nnf, linc, nlinc, lfac, llf, minc
integer :: nlminc, mmf, na, na1, la, la1, ma, ma1, nb, lb
integer :: mb, j2, mfac, n, l, m, ii, nn, nnll, ll1, nnllmm
integer :: mm
!***

IF(ISIGN.EQ.0) GO TO 300
N1=NTYPC
N2=NTYPD
GO TO 310
  300 CONTINUE
N1=NTYPA
N2=NTYPB
  310 CONTINUE
I1=IENTR(N1)
I1F=IEND(N1)
J1=JENTR(N1)
I2=IENTR(N2)
I2F=IEND(N2)
J2=JENTR(N2)
NF=NMAX(N1)+NMAX(N2)+1
LF=LMAX(N1)+LMAX(N2)+1
MF=MMAX(N1)+MMAX(N2)+1
NLS=-1
IF(IAND(MASK,4).NE.4) GO TO 2
1     NINC=2
NFAC=9
IF(IAND(ISYM,4).NE.4) GO TO 4
3     NNF=2
GO TO 5
4 NNF=1
 GO TO 5
2     NINC=1
NNF=1
NFAC=45
5     IF(IAND(MASK,2).NE.2) GO TO 7
6     LINC=2
NLINC=34
LFAC=9
IF(IAND(ISYM,2).NE.2) GO TO 9
8     LLF=2
NLS=NLS+17
GO TO 10
9 LLF=1
GO TO 10
7 LINC=1
NLINC=17
LLF=1
LFAC=45
10     IF(IAND(MASK,1).NE.1) GO TO 12
11    MINC=2
NLMINC=578
 MFAC=9
 IF(IAND(ISYM,1).NE.1) GO TO 14
13    MMF=2
NLS=NLS+289
GO TO 15
14    MMF=1
GO TO 15
12    MINC=1
NLMINC=289
MMF=1
MFAC=45
15    CONTINUE
NA=NAP(I1)*9
NA1=NA
LA=LAP(I1)*9
LA1=LA
MA=MAP(I1)*9
MA1=MA
NB=NAP(I2)*NFAC
LB=LAP(I2)*LFAC
MB=MAP(I2)*MFAC
N=NA+NB
L=LA+LB
M=MA+MB
II=I2
K=0
DO 370 NN=NNF,NF,NINC
NNLL=NN+NLS
!     DO 360 LL=LLF,LF,LINC
DO 360 ll1=LLF,LF,LINC
NNLLMM=NNLL
DO 350 MM=MMF,MF,MINC
IF(II.NE.I2F) GO TO 100
IF(I1.NE.I1F) GO TO 100
!     D(K+1)=X(NN+N)*Y(LL+L)*Z(MM+M)
D(K+1)=X(NN+N)*Y(ll1+L)*Z(MM+M)
GO TO 330
100   D(K+1)=0.D0
  320 CONTINUE
IF(MM.GT.(MAP(I1)+MAP(II)+1)) GO TO 110
!     IF(LL.GT.(LAP(I1)+LAP(II)+1)) GO TO 110
IF(ll1.GT.(LAP(I1)+LAP(II)+1)) GO TO 110
IF(NN.GT.(NAP(I1)+NAP(II)+1)) GO TO 110
!     D(K+1)=D(K+1)+CTERM(I1+J1)*CTERM(II+J2)*X(NN+N)*Y(LL+L)*Z(MM+M)
D(K+1)=D(K+1)+CTERM(I1+J1)*CTERM(II+J2)*X(NN+N)*Y(ll1+L)*Z(MM+M)
110   CONTINUE
IF(II.EQ.I2F) GO TO 325
II=II+1
N=NA+NAP(II)*NFAC
L=LA+LAP(II)*LFAC
M=MA+MAP(II)*MFAC
GO TO 320
  325 CONTINUE
IF(I1.EQ.I1F) GO TO 331
I1=I1+1
II=I2
NA=NAP(I1)*9
LA=LAP(I1)*9
MA=MAP(I1)*9
N=NA+NB
L=LA+LB
M=MA+MB
GO TO 320
331   CONTINUE
I1=IENTR(N1)
II=I2
NA=NA1
LA=LA1
MA=MA1
N=NA+NB
L=LA+LB
M=MA+MB
  330 CONTINUE
IF(D(K+1).EQ.0.D0) GO TO 340
NLM(K+1)=NNLLMM
K=K+1
  340 CONTINUE
NNLLMM=NNLLMM+NLMINC
  350 CONTINUE
NNLL=NNLL+NLINC
  360 CONTINUE
  370 CONTINUE
RETURN
end subroutine dd
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine abint

use legacy, only: nbfns, ntype, nfirst, nlast, eta, none1, n_prim

implicit none

!***
!  Local variables:
!
integer :: nlmp(105)
real*8  :: dab(105)
real*8  :: x(225), y(225), z(225)

real*8  :: a, pax(3), pbx(3)
integer :: masab, isymab, nmp, lmp, mmp, lp, isign
real*8  :: api, xdif(3)
integer :: nbfn, nbmax, max
integer :: NTYPA, NTYPB, NTYPC, NTYPD
real*8  :: XA(3), ALPA, XB(3), ALPB!, XC(3), ALPC, XD(3), ALPD
integer :: m, mn, mfi, mla, nla, nfi, mk, nk
integer :: i_ab, i_mk, i_nk, nm, nm_prim
integer :: n, nn

nbfn=nbfns
nbmax=nbfns
max=nbmax*(nbmax+1)/2

n_prim2 = n_prim**2

allocate( map_ab(max,n_prim2) )
n_ab = 0
do m=1,nbfns
  do n=1,m
    nm = none1(m) + n
    do mk=nfirst(m),nlast(m)
      i_mk = mk - nfirst(m) + 1
     do nk=nfirst(n),nlast(n)
       i_nk = nk - nfirst(n) + 1
       nm_prim = (i_mk-1) * n_prim + i_nk
       n_ab = n_ab + 1
       map_ab(nm,nm_prim) = n_ab
      end do
    end do
  end do
end do

allocate( lp_2d(n_ab) )
allocate( llp(n_ab) )
allocate( ap(n_ab), ex(n_ab) )
allocate( px(3,n_ab) )
allocate( dab_2d(105,n_ab) )
allocate( nlmp_2d(105,n_ab) )
!fk   allocate( isymab_2d(n_ab) )
!fk   allocate( masab_2d(n_ab) )
allocate( nmp_2d(n_ab), lmp_2d(n_ab), mmp_2d(n_ab) )

!=======================================================================
!     AB CHARGE DISTRIBUTION
!=======================================================================
do M=1,nbfns
  NTYPA=NQT(NTYPE(M))
  MFI=NFIRST(M)
  XA(1:3)=ETA(MFI,1:3)
  MLA=NLAST(M)
  do N=1,M
    NTYPB=NQT(NTYPE(N))
    NFI=NFIRST(N)
    XB(1:3)=ETA(NFI,1:3)
    NLA=NLAST(N)
    MN=NONE1(M)+N
    nm = none1(m) + n
    do MK=MFI,MLA
      i_mk = mk - mfi + 1
      ALPA=ETA(MK,4)
      do NK=NFI,NLA
        i_nk = nk - nfi + 1
        ALPB=ETA(NK,4)
        nm_prim = (i_mk-1) * n_prim + i_nk
        i_ab = map_ab(nm,nm_prim)
        AP(i_ab)=ALPA+ALPB
        API=1.D0/AP(i_ab)
        A=-ALPB*API
        MASAB=0
        XDIF(:)=XA(:)-XB(:)
        IF(XDIF(1).EQ.0.D0) MASAB=4
        IF(XDIF(2).EQ.0.D0) MASAB=MASAB+2
        IF(XDIF(3).EQ.0.D0) MASAB=MASAB+1
!fk     masab_2d(i_ab) = masab
        ISYMAB=MASAB
        pax(:) = xdif(:) * a
        pbx(:) = pax(:) + xdif(:)
        px(:,i_ab) = pax(:) + xa(:)
!       LLP=LL(NTYPA)+LL(NTYPB)
        llp(i_ab) = LL(NTYPA)+LL(NTYPB)
!=======================================================================
! CHECK FOR SYMMETRY OF INTEGRAL
!=======================================================================
        NMP=NMAX(NTYPA)+NMAX(NTYPB)
        IF(IAND(NMP,1).EQ.0) ISYMAB=IAND(ISYMAB,3)
        LMP=LMAX(NTYPA)+LMAX(NTYPB)
        IF(IAND(LMP,1).EQ.0) ISYMAB=IAND(ISYMAB,5)
        MMP=MMAX(NTYPA)+MMAX(NTYPB)
        IF(IAND(MMP,1).EQ.0) ISYMAB=IAND(ISYMAB,6)
        nmp_2d(i_ab) = nmp
        lmp_2d(i_ab) = lmp
        mmp_2d(i_ab) = mmp
!fk     isymab_2d(i_ab) = isymab
!       EX=A*ALPA*(XDIF*XDIF+YDIF*YDIF+ZDIF*ZDIF)
!       EX(i_ab)=A*ALPA*sum(xdif(:)**2)
!       EX(i_ab)=exp(A*ALPA*sum(xdif(:)**2))
        EX(i_ab)=PI5HF2*exp(A*ALPA*sum(xdif(:)**2))
        A=0.5D0*API
        if((ntypa.le.4).and.(ntypb.le.4)) then
!fk       N=NTYPA+4*(NTYPB-1)
          nn=NTYPA+4*(NTYPB-1)
          call expldd(DAB,NLMP,LP,A,PAX,PBX,nn)
        else
          CALL XYZ(X,NMAX(NTYPA),NMAX(NTYPB),A,PAX(1),PBX(1))
          CALL XYZ(Y,LMAX(NTYPA),LMAX(NTYPB),A,PAX(2),PBX(2))
          CALL XYZ(Z,MMAX(NTYPA),MMAX(NTYPB),A,PAX(3),PBX(3))
          ISIGN=1
!=======================================================================
! THE FOLLOWING CARD WAS CHANGED SO THAT THE X,Y,Z-MATRICES ARE
! PASSED INTO SUBROUTINE DD VIA THE ARGUMENT.  TLG  4/19/83.
!=======================================================================
!         call dd(DAB,NLMP,LP,ISIGN,MASAB,ISYMAB,X,Y,Z,xc,alpc,ntypc,ntypd,xd,alpd,XA,ALPA,NTYPA,NTYPB,XB,ALPB)
          call dd(DAB,NLMP,LP,ISIGN,MASAB,ISYMAB,X,Y,Z,ntypc,ntypd,ntypa,ntypb)
        end if
        lp_2d(i_ab) = lp
        dab_2d(:,i_ab) = dab(:)
        nlmp_2d(:,i_ab) = nlmp(:)
      end do
    end do
  end do
end do

end subroutine abint
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine cdint(xd)

use legacy, only: nbfns, ntype, nfirst, nlast, eta, n_prim

implicit none

real*8, intent(in)  :: xd(3)

!***
!  Local variables:
!
integer :: nlmq(8)
complex*16 :: dcd(8)
integer :: masc
integer :: isymcd, lq, i, ifi, ila, ik
integer :: nbfn, nbmax 
integer :: NTYPC !, NTYPD
real*8  :: XC(3), ALPC !, ALPD
integer :: i_cd, i_ik
logical, save :: first_call_cd = .true.

nbfn=nbfns
nbmax=nbfns

if( first_call_cd ) then
  allocate( map_cd(nbmax,n_prim) )
  n_cd = 0
  do i=1,nbmax
    do ik=nfirst(i),nlast(i)
      i_ik = ik - nfirst(i) + 1
      n_cd = n_cd + 1
      map_cd(i,i_ik)= n_cd
    end do
  end do
  allocate( lq_2d(n_cd) )
  allocate( qx(3,n_cd) )
  allocate( dcd_2d(8,n_cd) )
  allocate( nlmq_2d(8,n_cd) )
!fk   allocate( mascd_2d(n_cd) )
!fk   allocate( isymcd_2d(n_cd) )
  allocate( xcxd(n_cd) )
  allocate( patho_2d(n_cd) )
  allocate( alpc_1d(n_cd) )
  allocate( ntypc_1d(n_cd) )
  first_call_cd = .false.
end if

!=======================================================================
!   CD CHARGE DISTRIBUTION
!=======================================================================
do I=1,nbmax
  NTYPC=NQT(NTYPE(I))
  IFI=NFIRST(I)
  XC(1:3)=ETA(IFI,1:3)
  ILA=NLAST(I)
  do IK=IFI,ILA
    ALPC=ETA(IK,4)
    i_cd = map_cd(i,ik-ifi+1)
    ntypc_1d(i_cd) = ntypc
    alpc_1d(i_cd) = alpc
    MASC=0
    qx(:,i_cd)=DCMPLX(-XC(:),-.5D0*XD(:)/ALPC)
    IF(XD(1).EQ.0.D0)MASC=4
    IF(XD(2).EQ.0.D0)MASC=MASC+2
    IF(XD(3).EQ.0.D0)MASC=MASC+1
    ISYMCD=MASC
!fk mascd_2d(i_cd) = masc
    IF(IAND(NMAX(NTYPC),1).EQ.0) ISYMCD=IAND(ISYMCD,3)
    IF(IAND(LMAX(NTYPC),1).EQ.0) ISYMCD=IAND(ISYMCD,5)
    IF(IAND(MMAX(NTYPC),1).EQ.0) ISYMCD=IAND(ISYMCD,6)
!fk  isymcd_2d(i_cd) = isymcd
!   call expldp(DCD,NLMQ,LQ,0,XC,ALPC,NTYPC,NTYPD,XD,ALPD)
    call expldp(DCD,NLMQ,LQ,0,ALPC,NTYPC,XD)
    lq_2d(i_cd) = lq
    dcd_2d(:,i_cd) = dcd(:)
    nlmq_2d(:,i_cd) = nlmq(:)
!fk xcxd(i_cd) = sum( xc(:) * xd(:) )
    xcxd(i_cd) = exp( cmplx(0d0,sum( xc(:) * xd(:) )) )
    if(mod(ll(ntypc),2).eq.1) xcxd(i_cd) = - xcxd(i_cd)
    patho_2d(i_cd) = -0.25d0 * sum( xd(:)*xd(:) ) / alpc
  end do
end do

end subroutine cdint
!----------------------------------------------------------------

end module mplane1
!----------------------------------------------------------------
