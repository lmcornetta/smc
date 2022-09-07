module poly_conf 

use param
use precision

implicit none

!***
!  Modifications introduced by Fábris Kossoski in the begining of 2017.
!
!  elmbmod subroutine is now performed in part A. Since the computed
!  matrix eaqqd2 is linear in energy, it can be simply rescaled in part B.
!  This matrix is now printed in WFN.FILE, just after the energy-independent 
!  eaqqd matrix.
!
!  e_corr variable keeps track of the so-called correlation functions,
!  which is used further in the code.
!
!  The square matrixes haqqd and eaqqd were introduced, which
!  makes assignments in apresult  much faster than if dealing with 
!  triangular forms.
!
!  Now the product hci(k)*hci(l) is performed only once and stored in hci2.
!
!  The indexes of vvd and vvdg were exchanged, such that their later use
!  when updating eaqqd becomes much faster. Also, the previous 
!  vvd(g)(...,3) elements were properly joined with vvd(g)(...,2).
!
!  Multiplicative constants were moved from inner to outer loops.
!
!  Orbital eigenvalues (epsilon) are now printed, which is usefull when
!  an energy criteria is used to build the CSFs space.
!
!  Hartree to eV convertion factor modified to the more accurate value
!  of 27.21138386.
!***

!***
!  Fuction RPMOL modified to calculate 2-e integrals on demand
!  as requested by NTWOEL=3
!
!  Modified to use the NROUTE flag explained in module param
!
!  Modified to include solvent effects employing the COSMO
!  PCM approach as implemented in GAMESS. Controlled by the
!  flag SOLVAT.
!
!  Modified to write down the operator 0.5*(PV+VP). This option
!  is related to NROUTE=5.
!
!  Modified to improve the interface with module legacy of the new part B.
!  The arrays are now written to unit n_wvfun in a different order for
!  convenience. 
!
!  Modified to work with large integers as defined by the parameter
!  I_LARGE in module PRECISION.
!
!  Modified to keep track of the corrections to the normalization of the
!  correlation configurations, [NA-I](I), in case the scattering Hamiltonian
!  is saved to be diagonalized. These configurations should be normalized
!  with 1.0 (formally), but they are actually normalized with 1/sqrt(2.0d0).
!  Note that this does not affect the scattering amplitude in view of the
!  arbitrary (variational) coefficients multiplying the configurations in
!  the expansion of the trial wave function. The modifications were
!  implemented in the IDENT and ELMDN0MOD subroutine and the configuration
!  maps are loaded into C_CORR and N_CORR arrays.
!***

!***
! Then following commen blocks of the original code are
! defined here
!
!     COMMON /PH/NG,NS,NGT(NGX),NHT(NSX),NPT(NSX),NG2(NSX),NGT2(NGX,NSX)
!    1,NUMSD(NSBFX,NSX,2),NUMSQ(NSBFX,NSX)
integer, save :: NG, NS
integer, allocatable, save :: NGT(:),NHT(:),NPT(:),NG2(:), NGT2(:,:), NUMSD(:,:,:),NUMSQ(:,:), nums(:,:)
!
integer, save :: ngmax

!     COMMON /PROJCT/NEXC2(NCHLDX)
integer, allocatable, save :: NEXC2(:)

!     COMMON /SPIN/NSP(3,3),SP(3,3,2),NCSPN2(NCHLDX),MST(2,NPSCX),
!    1NAD(NPSCX)
integer, save :: NSP(3,3)
real*8,  save :: SP(3,3,2)
integer, allocatable, save :: NCSPN2(:),MST(:,:), NAD(:)

!     COMMON /CHNLTR/NEXC(NPSCX),NCSPN(NPSCX),NCTR(NPCTX),
!    1NCTRIN(2,NCHLX)
integer, allocatable, save :: nexc(:), ncspn(:)
integer, allocatable, save :: nctr(:), nctrin(:,:)

!     COMMON /HAMIL/ENE,EG(NCHLX),DEG(NCHLX),F(NBFNX,NBFNX),
!    1F2(NBFNX,NBFNX),NFC,XNFC2,NFC3
integer, save :: nfc, nfc3
real*8,  save :: ene, xnfc, xnfc2
real*8, allocatable, save :: deg(:), f(:,:), f2(:,:)
!   EG declared in poly-gaus module

!     COMMON /PRESV/EAQQD(DDX,NDDX),EAQQQ(NDQX,NDQX),
!    1HAQQD(NDDX,NDDX),HAQQQ(NDQX,NDQX)
! positron
real*8,  allocatable, save :: eaqq(:,:)

real*8,  allocatable, save :: eaqqd_tri(:), eaqqq_tri(:)
real*8,  allocatable, save :: eaqqd(:,:), eaqqq(:,:)
real*8,  allocatable, save :: eaqqd2_tri(:), eaqqq2_tri(:)
real*8,  allocatable, save :: eaqqd2(:,:), eaqqq2(:,:)
real*8,  allocatable, save :: haqqd_tri(:), haqqq_tri(:)
real*8,  allocatable, save :: haqqd(:,:), haqqq(:,:)

integer(i_large), save :: ij_ed = 0
integer(i_large), save :: ij_hd = 0
integer(i_large), save :: ij_eq = 0
integer(i_large), save :: ij_hq = 0

!     COMMON /NGTIN/ NGTINV(NBFNX),NUMINV(NOCCX,NMOFX),NGTIN2(NSBFX,NSX)
integer, allocatable, save :: ngtinv(:)
integer, allocatable, save :: numinv(:,:), ngtin2(:,:)

!     COMMON /PV/VVDG(NBFNX,NDDX),VVQG(NBFNX,NDQX),
!    1VVD(NBFNX,NDDX,NSAX,3),VVQ(NBFNX,NDQX,NSAX,3)
real*8, allocatable, save :: vvdg(:,:), vvqg(:,:)
real*8, allocatable, save :: vvd(:,:,:,:), vvq(:,:,:,:), vv(:,:,:)
!***

!***
!  Other important variables defined here
!
integer, save :: NELC, NOCC, NCHL, NCHLD, JSCHW
integer, save :: ND, NDD, NDQ, NSA, NCHLS, NCHLT
integer, save :: NPSC, NPCT
integer, save :: NPSCMAX, NPCTMAX
integer, allocatable :: NSNGLT(:), NTRPLT(:)
integer, allocatable :: NOIN(:,:)
! positron
integer, save :: ndop

real*8, save  :: eg0
real*8, allocatable, save  :: hci(:,:,:), hci2(:,:,:), egn(:)

!***
!  This array will keep track of the correlation configurations and
!  saved to disk to correct their normalization (to obtain the correct
!  pseudo spectrum of the scattering Hamiltonian).
!
integer, save :: n_corr
integer,  allocatable, save :: c_corr(:)
logical,  allocatable, save :: e_corr(:)
!***

contains

!----------------------------------------------------------------
subroutine sctin_new

implicit none

!***
!    This subroutine reads the wave function block only to load
!    the variables and arrays used in poly. The other arrays are
!    deallocated to save memory and then reallocated in sctin0
!***

!***
!     Local variables
!
integer :: NSING,NTRIP

integer :: i, k, kk
integer :: stat
!***

call ini_flag(flag_w1)

read(5,*,iostat=stat) NELC,NOCC,NG,NS,NCHL,JSCHW,NSING,NTRIP
call write_stat('NELC,NOCC,NG,NS,NCHL,JSCHW,NSING,NTRIP',stat)
read(5,*,iostat=stat) NSA,NCHLT,NCHLS
call write_stat('NSA,NCHLT,NCHLS',stat)
NCHLD=NCHL-1

!***
!  Allocate arrays related to configuration space
!
allocate( ngt(ng) )
ngmax = ng

if(nchld .ge. 1) allocate( nexc2(nchld), ncspn2(nchld) )

if(ns .ge. 1) then
  allocate( nsnglt(ns), ntrplt(ns) )
  allocate( nht(ns), npt(ns) )
  allocate( ng2(ns) )
end if
!***

read(5,*,iostat=stat) (NGT(I),I=1,NG)
call write_stat('ngt',stat)

if(ns.ne.0) then

if(NSING.EQ.0.AND.NS.GT.1) then
  read(5,*,iostat=stat) (NSNGLT(I),I=1,NS)
  call write_stat('nsnglt',stat)
else
  NSNGLT(:)=NSING
end if

if(NTRIP.EQ.0.AND.NS.GT.1) then
  read(5,*,iostat=stat) (NTRPLT(I),I=1,NS)
  call write_stat('ntrplt',stat)
else
  NTRPLT(:)=NTRIP
end if

read(5,*,iostat=stat) (NHT(I),I=1,NS)
call write_stat('nht',stat)
read(5,*,iostat=stat) (NPT(I),I=1,NS)
call write_stat('npt',stat)
read(5,*,iostat=stat) (NG2(I),I=1,NS)
call write_stat('ng2',stat)

do i = 1,ns
  ngmax = max( ngmax,ng2(i) )
end do
allocate( ngt2(ngmax,ns) )

do K=1,NS
  if(ng2(k).ne.0) then
    read(5,*,iostat=stat) (NGT2(I,K),I=1,NG2(K))
    call write_stat('ngt2',stat)
  else
    KK=K-1
    NG2(K)=NG2(KK)
    do I=1,NG2(K)
      NGT2(I,K)=NGT2(I,KK)
    end do
  end if
end do

if(nchl.gt.1) then
  read(5,*,iostat=stat) (NEXC2(I),I=1,NCHLD)
  call write_stat('nexc2',stat)
  read(5,*,iostat=stat) (NCSPN2(I),I=1,NCHLD)
  call write_stat('ncspn2',stat)
end if

end if ! ns.ne.0

call end_flag(flag_w2)

!***
!  Not necessary for the moment...
!
deallocate(ngt)
if( allocated(nsnglt) ) deallocate(nsnglt)
if( allocated(ntrplt) ) deallocate(ntrplt)
if( allocated(ng2) )    deallocate(ng2)
if( allocated(ngt2) )   deallocate(ngt2)
if( allocated(ncspn2) ) deallocate(ncspn2)
!***

!***
! To avoid problems in poly these must be allocated. They
! will not be used at all since ns=0.
!
if( ns .eq. 0) then
  allocate( npt(nchl), nexc2(nchl) )
  npt(1)   = 0
  nexc2(1) = 0
end if
!***

!***
!  Inclusion of solvent effects is currently not implemented 
!  with electronic excitation. Avoid problems ahead...
!
if( (solvat.ne.'VACUUM') .and. (nchl.gt.1) ) then
  write(6,"(//,'SOLVATION EFFECTS NOT IMPLEMENTED',' WITH ELECTRONIC EXCITATION')")
  stop'ABORT...'
end if
!***
  
end subroutine sctin_new
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine sctin0

implicit none

!***
!     Local variables
!
integer :: NSING,NTRIP
integer :: i, k, kk
integer :: stat
!***

call ini_flag(flag_w1)

read(5,*,iostat=stat) NELC,NOCC,NG,NS,NCHL,JSCHW,NSING,NTRIP
call write_stat('NELC,NOCC,NG,NS,NCHL,JSCHW,NSING,NTRIP',stat)
read(5,*,iostat=stat) NSA,NCHLT,NCHLS
call write_stat('NSA,NCHLT,NCHLS',stat)
!
!roma 16/09/2002
!=======================================================================
! IF NSING=0 AND NS>1(MORE THAN ONE EXCITATED STATE)THEN READ NSNGLT(NS)
! IF NTRIP=0 AND NS>1(MORE THAN ONE EXCITATED STATE)THEN READ NTRPLT(NS)
!=======================================================================
NCHLD=NCHL-1

!***
!  Allocate arrays related to configuration space
!  NHT, NPT and NEXC2 have been allocated in sctin_new
!
allocate( ngt(ng) )

if(nchld .ge. 1) allocate( ncspn2(nchld) )

if(ns .ge. 1) then
  allocate( nsnglt(ns), ntrplt(ns) )
  allocate( ng2(ns) )
end if

read(5,*,iostat=stat) (NGT(I),I=1,NG)
call write_stat('ngt',stat)

if(ns.ne.0) then

if(NSING.EQ.0.AND.NS.GT.1) then
  read(5,*,iostat=stat) (NSNGLT(I),I=1,NS)
  call write_stat('nsnglt',stat)
else
  NSNGLT(:)=NSING
end if

if(NTRIP.EQ.0.AND.NS.GT.1) then
  read(5,*,iostat=stat) (NTRPLT(I),I=1,NS)
  call write_stat('ntrplt',stat)
else
  NTRPLT(:)=NTRIP
end if

read(5,*,iostat=stat) (NHT(I),I=1,NS)
call write_stat('nht',stat)
read(5,*,iostat=stat) (NPT(I),I=1,NS)
call write_stat('npt',stat)
read(5,*,iostat=stat) (NG2(I),I=1,NS)
call write_stat('ng2',stat)

allocate( ngt2(ngmax,ns) )
ngt2 = 0

do K=1,NS
  if(ng2(k).ne.0) then
    read(5,*,iostat=stat) (NGT2(I,K),I=1,NG2(K))
    call write_stat('ngt2',stat)
  else
    KK=K-1
    NG2(K)=NG2(KK)
    do I=1,NG2(K)
      NGT2(I,K)=NGT2(I,KK)
    end do
  end if
end do

!=======================================================================
! SELECT THE EXCITED CONFIGURATIN FROM THE ABOVE NS CONFIGURATIONS
!=======================================================================
if(nchl.gt.1) then
  read(5,*,iostat=stat) (NEXC2(I),I=1,NCHLD)
  call write_stat('nexc2',stat)
  read(5,*,iostat=stat) (NCSPN2(I),I=1,NCHLD)
  call write_stat('ncspn2',stat)
end if

end if ! ns.ne.0

call end_flag(flag_w2)

!=======================================================================
!
! NCSPN2 MUST BE -1 (FOR SINGLET) OR 1 (FOR TRIPLET).
!
! AUTOMATIC EXTENSION OF BASIS SET FOR OPEN CHANNEL ORBITALS:
! ALL THE PARTICLE STATES DEFINING EXCITED STATES ARE USED
! AS BASIS FUNCTIONS FOR OTHER OPEN CHANNEL ORBITALS.
!
!#    DO 450 III=1,NCHLD
!#    II=NEXC2(III)
!#    M=0
!#    DO 460 KKK=1,NCHLD
!#    K=NPT(NEXC2(KKK))
!#    DO 470 JJ=1,NG2(II)
!#    J=NGT2(JJ,II)
!#470IF(J.EQ.K) GO TO 460
!#    M=M+1
!#    NADD(M)=K
!#460 CONTINUE
!#    IF(M.EQ.0) GO TO 450
!#    IM=NG2(II)
!#    DO 480 L=1,M
!#    IM=IM+1
!#480NGT2(IM,II)=NADD(L)
!#    NG2(II)=NG2(II)+M
!#450 CONTINUE
!
!=======================================================================
! 250 CONTINUE

!***

!***
!=======================================================================
! PREPARATION FOR ORTHOGONALITY CONDITIONS. SEE 'GINTS'.
!=======================================================================
!   Moved to subroutine POLY
!     DO 400 I=1,NOCC
!     400   NOH(I)=I
!     NOOR=NOCC
!   IF(NCHL.EQ.1) GO TO 430
!   DO 431 I=1,NCHLD
!   II=NPT(NEXC2(I))
!   DO 432 M=1,NOOR
!     432IF(II.EQ.NOH(M)) GO TO 431
!   NOOR=NOOR+1
!  NOH(NOOR)=II
!     431 CONTINUE
!     430 CONTINUE
!***

write(6,520)
  520 FORMAT(//'**** SCHWINGER TYPE FUNCTIONAL *********************')
if(NCHLD.NE.0) write(6,305)
  305 FORMAT(//'**** MULTI-CHANNEL SCATTERING *****************')
if(NS.EQ.0) write(6,310)
  310 FORMAT(//'**** STATIC-EXCHANGE APPROXIMATION FOR ELASTIC SCATTERING **********')
if(NCHLD.EQ.0.AND.NS.NE.0) write(6,306)
  306 FORMAT(//'**** ELASTIC SCATTERING WITH POLARIZATION EFFECT ***')
ND=NG
write(6,500)
  500 FORMAT(///'**** INPUT DATA FOR SCATTERING BASIS AND CONFIGURATIONS ****')
write(6,600) NELC,NOCC
  600 FORMAT(//'    NELC=',I2,'  NOCC=',I2)
!mbformat
!mb   30(I2,',') para 20(I4)
!mbformat

!     WRITE(6,801) NG,(NGT(I),I=1,NG)
! 801 FORMAT(//'  NGT ',I3,/,20(1X,I3))
!     IF(NS.EQ.0) RETURN
!     WRITE(6,802) NS,(NHT(I),I=1,NS)
! 802 FORMAT(//'  NHT ',I3,'/ ',20(1X,I3))
!     WRITE(6,803) NS,(NPT(I),I=1,NS)
! 803 FORMAT(//'  NPT ',I3,'/ ',20(1X,I3))
!     DO 300 K=1,NS
! 300WRITE(6,320) K,NG2(K),(NGT2(I,K),I=1,NG2(K))
! 320 FORMAT(//'  NGT2 ',I3,'  NG2=',I3,/,20(1X,I3))
!mb
!     DO 777 K=1,NS
!     IF(NSNGLT(K).NE.0) WRITE(6,760)NHT(K),NPT(K)
!     IF(NTRPLT(K).NE.0) WRITE(6,761)NHT(K),NPT(K)
! 777 CONTINUE

!***
write(6,"(//,10x,'*** MOBSCI INFORMATION ***')")
write(6,"(5x,'NSA   = ',i4)") nsa
write(6,"(5x,'NCHLT = ',i4)") nchlt
write(6,"(5x,'NCHLS = ',i4)") nchls
!***

! 760 FORMAT(//'&&&& SINGLET CONFIGURATION IS NOT USED FOR CLOSED CHANNELS &
!                    &FOR THE PRIMITIVE CONFIGURATION : HOLE =',I3,'PARTICLE =',I3)
! 761 FORMAT(//'&&&& TRIPLET CONFIGURATION IS NOT USED FOR CLOSED CHANNELS &
!                    &FOR THE PRIMITIVE CONFIGURATION : HOLE =',I3,'PARTICLE =',I3)
if(NCHL.LE.1) return
write(6,810) NCHL,(NEXC2(I),I=1,NCHLD)
  810 FORMAT(//'  NUMBER OF PLANE WAVES OF DIFFERENT MAGNITUDE = ',I3//&
      '     EXCITED STATES OF TARGET=',30I3)
write(6,811) (NCSPN2(I),I=1,NCHLD)
  811 FORMAT(//'     SPIN MULTIPLICITY OF CHANNELS=',30I3)
!=======================================================================
!NOTE: SO FAR THE NUMBER OF OPEN CHANNELS IS NOT COUNTED
!CORRECTLY. THEREFORE 'NCHL' DOES NOT MEAN THE NUMBER OF OPEN
!CHANNEL, IF TRIPLET STATE IS OPEN. THIS COUNT WILL BE MADE IN
!SUBROUTINE 'IDENT'.
!=======================================================================

end subroutine sctin0
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine ident_electron

!***
! This subroutine has been modified to load a map of the correlation
! configurations into the C_CORR array in case the Hamiltonian 
! matrix is saved
!***

use poly_gaus, only: nbfns, nmof_tar

implicit none

!***
!  Local variables
integer :: i, j, k, l, ii, jj, kk, ll, i1, na, ma, nstwo
integer :: NDP, JKD, JKQ, nnnx

integer, allocatable :: nnn(:,:)
real*8, allocatable  :: nopen(:,:)
!***

!romaJK=NG
JKD=NG
JKQ=0
!roma
NDP=NG
!=======================================================================
!
! WE HAVE (NS+1) PRIMITIVE COMFIGULATIONS WHICH MAKE EXACT OR PSEUDO
! STATES. HERE WE DEFINE SERIAL NUMBERS FOR "PHYSICALLY" INDEPENDENT
! CONFIGULATIONS,  TAKING ACCOUNT OF SPIN FUNCTIONS.
! // PRIMITIVE EXCITATION (NA-I) CAN HAVE THREE  CONFIGULATIONS;
!    (NA-I)(AB-BA)/RT2,  (NA-I)(AB+BA)/RT2, AND (NA-I)(AA).
!    SO, IN THE DEFINITION  NCTR(M)=K
!    M:  SERIAL NUMBER OF "PHYSICALLY" INDEPENDENT CONFIG., WHERE
!  M=1 STANDS FOR THE GROUND STATE.
!    K:  SERIAL NUMBER OF PRIMITIVE CONFINGULATIONS (NA-I), WHERE
!  K=1 IS FOR THE GROUND STATE.
!
!   IN THE DEFINITIONS OF NAD(N),NEXC(N),NCSPN(N),MST(1,N).....,
!     N INDICATES THE SERIAL NUMBER OF "PHYSICALLY" INDEPENDENT CONF.,
!     WHERE N=0  IS RESERVED FOR THE GROUND STATE.
!
!     /////  SPIN FUNCTIONS FOR OPEN SHELL PART (NA-I)J ARE DEFINED
!IN THIS PROGRAM BY (AB-BA)A/RT2, (AB+BA)A/RT2, AAB, WHERE
!A=ALPHA SPIN,  B=BETA SPIN.
!
!=======================================================================
!     Write (6,*)'Estou na poly '
!     Write (6,*)'NNNX=3*NSX*NGX ',NNNX,NSX,NGX
!     Write (6,*)'ND= ',ND
!     Write (6,*)'NDD= ',NDD
!     Write (6,*)'NDQ= ',NDQ

! fk
write(*,'(/a,i0)') '*** NUMBER OF OCCUPIED ORBITALS = ', nocc
write(*,'(/a,i0)') '*** NUMBER OF VIRTUAL ORBITALS  = ', nmof_tar - nocc
write(*,'(/a,i0)') '*** NUMBER OF EXCITATIONS       = ', ns

!***
!  This used to be done before calling ident (main)
!
allocate( nctrin(2,nchl) )
nd=ng
ndd=ng
ndq=0
npsc = 0
nctrin(1,1) = 1
nctrin(2,1) = 0

if(ns .eq. 0) then
!  This used to be done after calling ident (main)
  if(keyq.eq.0) ndq=1
  npct = npsc + 1
  return
end if
!***

!***
!  Allocate NUMS and NNN
nnnx = 3*ns*ngmax

allocate( numsd(nbfns,ns,2), numsq(nbfns,ns) )
numsd = 0
numsq = 0
allocate( nnn(4,nnnx) )
!***

!IF(NCHLD.EQ.0) GO TO 100
if(nchld.ne.0) then

!***
!  Allocate arrays related to electronic excitation
!
!mal
nstwo=2*ns
!mal
!mal  allocate( nopen(ns,2) )
allocate( nopen(nstwo,2) )

!   Get NPSC and NPCT first...
do i = 1,nchld
  if( ncspn2(i) .eq. 0) stop 'NCSPN2.EQ.0 IN IDENT'
  if( ncspn2(i) .gt. 0 ) then
    npsc = npsc + 1
  else
    npsc = npsc + 2
  end if
end do
npct = npsc + 1

npscmax=nchls+2*nchlt
npctmax=npscmax+1
allocate( nexc(npscmax), ncspn(npscmax) )
allocate( nad(npscmax), mst(2,npscmax) )
allocate( nctr(npctmax) )

!   This used to be done before calling ident (main)
NCTR(1)=1
!***

do II=1,NCHLD
  I=NEXC2(II)
  NOPEN(I,1)=0
  NOPEN(I,2)=0
end do

J=0
do I=1,NCHLD
  I1=I+1
!IF(NCSPN2(I)) 230,999,220
  if(ncspn2(i).lt.0) then
    NOPEN(NEXC2(I),1)=1
    J=J+1
    NCTR(J+1)=I1
    NCTRIN(1,I1)=J+1
    NCTRIN(2,I1)=0
    NAD(J)=1
    NEXC(J)=NEXC2(I)
    NCSPN(J)=1
  else if(ncspn2(i).gt.0) then
    NOPEN(NEXC2(I),2)=1
    J=J+1
    NCTR(J+1)=I1
    NCTRIN(1,I1)=J+1
    NAD(J)=2
    NEXC(J)=NEXC2(I)
    NCSPN(J)=-1
!=======================================================================
!   BE CAREFUL ABOUT THE SIGN OF NCSPN : 1 FOR SINGLET,  -1 FOR TRIPLET.
!=======================================================================
    J=J+1
    NCTR(J+1)=I1
    NCTRIN(2,I1)=J+1
    NAD(J)=3
    NEXC(J)=NEXC2(I)
    NCSPN(J)=0
  else
    goto 999
  end if
end do
NPSC=J

do K=1,NPSC
  MST(1,K)=1
  MST(2,K)=NCSPN(K)
end do

!  100 CONTINUE
end if ! nchld.ne.0

!roma 08/10/2002 - Aqui comecam as modificacoes: rodada apenas para dubletos
!
do II=1,NS
  NA=NHT(II)
  I=NPT(II)
  ld1: do JJ=1,NG2(II)
    J=NGT2(JJ,II)
    if(II.ne.1) then
      ld2: do KK=1,II-1
        MA=NHT(KK)
        IF(MA.NE.NA) cycle ld2
        K=NPT(KK)
        IF(I.EQ.K) GO TO 999
        IF(K.NE.J) cycle ld2
        do LL=1,NG2(KK)
          L=NGT2(LL,KK)
          IF(L.EQ.I) cycle ld1
        end do
      end do ld2
    end if
    if(NSNGLT(II).eq.0) then
      JKD=JKD+1
      NUMSD(JJ,II,1)=JKD
    end if
    if(NTRPLT(II).eq.0) then
! fk This next line must be before numsd(jj,ii,2)!
      if(J.EQ.I) cycle ld1
      JKD=JKD+1
      NUMSD(JJ,II,2)=JKD
    end if
! roma      JK=JK+1
! roma      NUMS(JJ,II,3)=JK
  end do ld1
end do

!roma Modificacao: rodada apenas para quartetos
!
do II=1,NS
  NA=NHT(II)
  I=NPT(II)
  lq1: do JJ=1,NG2(II)
    J=NGT2(JJ,II)
    IF(II.ne.1) then
      lq2: do KK=1,II-1
        MA=NHT(KK)
        IF(MA.NE.NA) cycle lq2
        K=NPT(KK)
        IF(I.EQ.K) GO TO 999
        IF(K.NE.J) cycle lq2
        do LL=1,NG2(KK)
          L=NGT2(LL,KK)
          IF(L.EQ.I) cycle lq1
        end do
      end do lq2
    end if
    IF(NTRPLT(II).NE.0) cycle lq1
    IF(J.EQ.I) cycle lq1
    JKQ=JKQ+1
    NUMSQ(JJ,II)=JKQ
  end do lq1
end do
!roma  440 ND=JK
NDD=JKD
NDQ=JKQ

ND=NDD+NDQ
write(6,400) ND,NDP,NPSC

400 FORMAT(//'*** INDEPENDENT SLATER DETERMINANTS = ',I8//&
             '*** DIMENSION OF DIRECT(OPEN CHANNEL) SPACE = ',I8//&
             '*** NPSC = ',I8)

write(6,401) NDD,NDQ
401 FORMAT(//'*** IND. SLATER DETERMINANTS IN DUBLET SPACE  = ',I8//&
             '*** IND. SLATER DETERMINANTS IN QUARTET SPACE = ',I8)

K=0
do I=1,NS
  NA=NHT(I)
  II=NPT(I)
  do J=1,NG2(I)
    JJ=NGT2(J,I)
    K=K+1
    NNN(1,K)=NUMSD(J,I,1)
    NNN(2,K)=NA
    NNN(3,K)=II
    NNN(4,K)=JJ
    K=K+1
    NNN(1,K)=NUMSD(J,I,2)
    NNN(2,K)=NA
    NNN(3,K)=II
    NNN(4,K)=JJ
    K=K+1
    NNN(1,K)=NUMSQ(J,I)
    NNN(2,K)=NA
    NNN(3,K)=II
    NNN(4,K)=JJ
  end do
end do

! fk
write(*,'(/a,i0)') '*** K = ', k

write(6,402)
402 FORMAT(//'**  NUMS, NA-I, J **')

call printn(NNN,4,NNNX,4,K)

!***
! This used to be done after coming back from ident (main)
!
if(keyq.eq.0) ndq=1
npct = npsc + 1
!***

!***
! Keep track of correlation configurations
!
!fk e_corr is required in elmb routine
!     If( nwhaqq .ne. 0 ) Then
  n_corr = 0
  if( ns .gt. 0 ) then
  allocate( c_corr(ns) )
  allocate( e_corr(ndd) )
  e_corr = .false.
  end if
  do i = 1, ns
    ii = npt(i)
    do j = 1,ng2(i)
      jj = ngt2(j,i)
      if( (numsd(j,i,1).ne.0) .and. (ii.eq.jj) ) then
        n_corr = n_corr + 1
        c_corr(n_corr) = numsd(j,i,1)
        e_corr(numsd(j,i,1)) = .true.
      end if
    end do
  end do
!     End If
!***

!***
! Deallocate those no longer needed
!
!mal
!     print*,'vou dealocar nnn'
!mal
!mal  deallocate( nnn )
!mal  if( allocated(nopen) ) deallocate( nopen )
!mal  if( keyq .eq.0 ) deallocate( numsq )
!***
!mal
!     print*,'vou sair da ident'
!mal

return
  999 write(6,800) II,KK
  800   FORMAT(//'  THE SAME HOLE-PARTICLE STRUCTURE HAS APPEARED(IDENT) / EXECUTION STOP',2I3)
stop

end subroutine ident_electron
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine fock

use poly_gaus, only: nbfns, nsbf, tket, atrt, eg, d

implicit none

!***
!  Local variables
!
integer :: i, j, ij, ia, ib, m, nn2, nnf
integer :: nsbf1, insbf, nbfn
character(12) :: pcm_label
real*8  :: x, y, xf, def
real*8  :: gtot, ges, vnuc
real*8, allocatable :: f3(:), fip(:)
character(len=*), parameter :: orbitals_energies_file = 'orbitals_energies'
integer :: nunit
integer :: stat

!  Computes Fock matrix F(I,J) and auxiliary F2(I,J):
!
!  F(I,J) = TKET(I,J) + ATRT(I,J) + SUM_M [ 2*RPMOL(I,J,M,M) + RPMOL(I,M,J,M) ]
!  For electrons: F2(I,J) = F(I,J) - TKET(I,J)
!  For positrons: F2(I,J) = ATRT(I,J) - SUM_M [ 2*RPMOL(I,J,M,M) ]
!  where  RPMOL(I,J,M,M) = ( I(1) J(1) / M(2) M(2) )

!=======================================================================
! NOCC= # OF OCCUPIED ORBITALS
! NFC,XNFC2,NFC3 ARE DETERMINED HERE.
! NOCC=NELC/2(CLOSED SYSTEM).
!=======================================================================

!***
!  Allocate fock matrices F (usual definition) and F2,
!  and also DEG
 
allocate( fip(nocc) )
allocate( f(nbfns,nbfns) )
allocate( deg(nchl) )
if( nrfock .eq. 0) then
  allocate( f2(nbfns,nbfns) )
else
  nnf = nbfns * (nbfns+1) / 2
  allocate( f2(nnf,1) )
  allocate( f3(nnf) )
end if
nbfn=nbfns

NFC=1
XNFC2=1.0D0/DSQRT(DFLOAT(NELC+1))
XNFC=1.0D0/DFLOAT(NELC+1)
NFC3=1

if (NSBF.ne.NBFN) then
  NSBF1=NSBF+1
  do i=NSBF1,NBFN
    do j=NSBF1,NBFN
      f(j,i) = 0.0d0
    end do
  end do
end if

if(nrfock .eq. 0) then

if( projectile.eq.-1 ) then
  do I=1,NBFN
    INSBF=MIN0(NSBF,I)
    do J=1,INSBF
      X=0.0D0
      Y=0.0D0
      do M=1,NOCC
        X=X+RPMOL(I,J,M,M)
        Y=Y+RPMOL(I,M,J,M)
      end do
      Y=2d0*X-Y
      f(j,i) = y+tket(j,i)+atrt(j,i)
      f(i,j) = f(j,i)
      f2(j,i) = y+atrt(j,i)
      f2(i,j) = f2(j,i)
    end do
  end do
else
  do I=1,NBFN
    INSBF=MIN0(NSBF,I)
    do J=1,INSBF
      X=0.0D0
      Y=0.0D0
      do M=1,NOCC
        X=X+RPMOL(I,J,M,M)
        Y=Y+RPMOL(I,M,J,M)
      end do
      Y=2d0*X-Y
      f(j,i) = y+tket(j,i)+atrt(j,i)
      f(i,j) = f(j,i)
      f2(j,i) = -(atrt(j,i)+2d0*x)
      f2(i,j) = f2(j,i)
    end do
  end do
end if

else

!   Use F2 as work space first
  open(unit=n_fockm, file='FOCK.FILE', form='FORMATTED', iostat=stat)
  call open_stat('FOCK.FILE',stat)
  read(n_fockm,*,end=999) nn2
  if(nn2 .ne. nnf) stop'WRONG DIMENSION IN FOCK.FILE' 
  read(n_fockm,*,end=999) ( f2(i,1), i=1,nn2 )

  if(solvat .eq. 'COSPCM') then
    read(n_fockm,"(a12)",end=999) pcm_label
    if(pcm_label .ne. 'PCM ENERGIES') then
      write(6,"(//,'WRONG LABEL IN FOCK.FILE: ','PCM_LABEL = ',a7)") pcm_label
      stop'ABORT...'
    end if
    read(n_fockm,*,end=999) gtot, ges
    read(n_fockm,*,end=999) vnuc
    write(6,"(//,10x,'***  SOLVATION ENERGY DATA ***')") 
    write(6,"('FREE ENERGY IN SOLVENT =',' <PSI| H(0)+V/2 |PSI> = ',f16.8)") gtot 
    write(6,"('INTERNAL ENERGY IN SOLVENT = ','<PSI| H(0) |PSI> = ', f16.8)") gtot-ges
    write(6,"('ELECTROSTATIC INTERACTION','   = ', f16.8)") ges
    gtot = gtot - vnuc
  end if

  ij = 0
  do i = 1,nbfn
    do j = 1,i
      ij = ij + 1
      f(i,j) = f2(ij,1)
      f(j,i) = f2(ij,1)
    end do
  end do
!   Transform to MO basis
  deallocate( f2 )
  allocate( f2(nbfn,nbfn) )
  ij = 0
  do i = 1,nbfn
    do j = 1,i 
      ij = ij + 1
      do ia = 1,nbfn
        xf = 0d0
        do ib = 1,nbfn
          xf = xf + d(ib,j) * f(ia,ib)
        end do
        f2(ia,j) = xf
      end do 
      xf = 0d0 
      do ia = 1,nbfn
        xf = xf + f2(ia,j) * d(ia,i)
      end do
      f3(ij) = xf
    end do
  end do
!   If NRFOCK<>0, build F2 as F2 = F - TKET
!   to account for solvation effects in COSPCM
  ij = 0
  do i = 1,nbfn
    do j = 1,i
      ij = ij + 1
      f(i,j)  = f3(ij)
      f2(i,j) = f3(ij) - tket(i,j)
      f(j,i)  = f3(ij)
      f2(j,i) = f2(i,j)
    end do
  end do
  deallocate( f3 )
end if
!***

open(newunit=nunit,file=orbitals_energies_file,form="formatted",iostat=stat)
call open_stat(orbitals_energies_file,stat)
do i=1,nbfn
  write(nunit,*) f(i,i)
end do
close(nunit)

if(solvat .eq. 'COSPCM' ) then
  do i = 1,nocc
    eg0=eg0+f(i,i)+tket(i,i)+atrt(i,i)
  end do
  write(6,"(/'INPUT EG0 OVERWRITTEN WITH THE CORRECT',' COS-PCM ELECTRONIC ENERGY',/,'EG0 = ',1pd15.8)") gtot 
  write(6,"(/'CALCULATED VALUE WITHOUT PROPER ','ONE-ELECTRON CORRECTIONS IN ATRT:',/,'EG0 = ',1pd15.8)") eg0
  eg0 = gtot
  deg(1) = eg0
  eg(1)  = eg0
  return
end if


eg0 = 0.0d0
do i=1,nocc
  eg0=eg0+f(i,i)+tket(i,i)+atrt(i,i)
end do

!mod
! open(unit=1,file='fock.lis',form='formatted')
! write(1,"('EG0 = ',1pe15.8)") eg0
! do i=1,nbfn
!   do j=1,nbfn
!     write(1,"(2(i4),2x,1pe15.8)") i,j,f(i,j)
!   end do
! end do
! close(1)
!mod

write(6,500)
write(6,"(5x,'INPUT FILE EG = ', 1pe15.8)") eg(1)
write(6,"(5x,'CALCULATED EG = ', 1pe15.8)") eg0
!mal
write(6,501)
  501 FORMAT(//'** OCCUPIED ENERGY ORBITALS (IN HARTREE) & 
                &IONIZATION ENERGIES (IN EV)**'/)
write(6,"(5x,'oeo(',i2,') = ', 1pe15.8)") (i,f(i,i),i=1,nocc)
do i=1,nocc
  fip(i) = -f(i,i) * toHartree
end do
write(6,"(5x,'ip(',i2,') = ', 1pe15.8)") (i,fip(i),i=1,nocc)
!mal
def = abs( eg0-eg(1) )
if(def .gt. 100d0*discrep) then
  write(6,101) def
  stop'ABORT...'
end if

if(def .gt. discrep) write(6,102) def

deg(1) = eg0
eg(1)  = eg0

return

  101 FORMAT(//'???? BIG DISCREPANCY IS FOUND BETWEEN INPUT EG0 AND CALCULATED EG0',D15.8/&
               '    EXECUTION STOP')
  102 format(//,20x,'***** WARNING ****',/,'MILD DISCREPANCY ',&
                    'BETWEEN INPUT EG0 AND CALCULATED EG0',1pd15.8,//)
!***
  500 FORMAT(//'*** GROUND STATE ENERGY OF THE TARGET (IN FOCK) ***'/)
 write(6,5555)
 5555 FORMAT(//'&&& F-MATRIX (FOCK) &&&')
  call print(F,nbfns,nbfns,nbfns,nbfns)
!     RETURN

 999  write(6,"('END OF FILE READING FOCK.FILE')")
stop 'ABORT...'

end subroutine fock
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine load_spin

implicit none

integer :: nnspin(3,3)
DATA NNSPIN/1,-1,-1,1,-1,1,1,1,-1/
integer :: i, j

!***
!  In the f77 verion, NSPIN is a dummy for NSP of common block SPIN.
!  RT2, RT3, etc. are set in module param
!fk
!   The normalization RT3/RRT2 is never used in practice since the
!   correlation configurations are stored in the NUMS(j,i,1) and
!   SP(i,j,1) by convention. 

do I=1,3
  do J=1,3
    nsp(j,i) = nnspin(j,i)
  end do
end do

!fk   XNELC1=DFLOAT(NELC)+1.0D0
SP(1,1,1)=RRT2
SP(1,2,1)=RRT2
SP(1,3,1)=0.0D0
SP(2,1,1)=-RRT6
SP(2,2,1)=RRT6
SP(2,3,1)=-2.0D0*RRT6
SP(3,1,1)=-RRT3
SP(3,2,1)=RRT3
SP(3,3,1)=RRT3
SP(1,1,2)=0.0D0
SP(1,2,2)=RT3/RRT2
SP(1,3,2)=0.0D0
SP(2,1,2)=0.0D0
SP(2,2,2)=0.0D0
SP(2,3,2)=0.0D0
SP(3,1,2)=0.0D0
SP(3,2,2)=0.0D0
SP(3,3,2)=0.0D0

end subroutine load_spin
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine elmdn0mod_electron

implicit none

!***
! Local variables
!
integer :: i, j, k, l, m, n, ii, jj, kk, ll, im
integer :: jkd1, jkd2, jkq, key1, key2, lmax, jm
integer :: na, nb, nam, nbn, kn, ln, mx, nx, imx, knx
integer(i_large) :: n_haqqd, n_haqqq
real*8  :: RES1, RES12, RES22, RES32
real*8  :: x, y, rx(3,3)
integer :: alloc_stat, stat

!***
!  This subroutine has been debugged by Fabris Kossoski in Semptember
!  2016. The modifications have been kept to a minimum (see !fk comments). 
!  The map of correlation configurations (C_CORR) is saved along with the
!  scattering Hamiltonian.
!***

n_haqqd = int8(ndd) * int8( (ndd+1) ) / int8(2)
allocate( haqqd_tri(n_haqqd), stat=alloc_stat ) 
call allocation_stat('haqqd_tri',alloc_stat)

if(keyq .ne. 0) then
  n_haqqq = ndq*(ndq+1)/2
  allocate( haqqq_tri(n_haqqq), stat=alloc_stat )
  call allocation_stat('haqqq_tri',alloc_stat)
end if
!=======================================================================
!     NOTE: HPART=(PV+VP)/2+(K-(KP+PK)/(N+1)*2/(N+1)
!    =(PV+VP)-P(K0)P+K/(N+1)
!     WHERE K=E-H, K0=E-H0
!     NOTE 2: E-HN IS CONSIDERED IN ELMDN1, WHERE HN IS THE HAMILTONIAN
! OF THE TARGET.    SO, ACTUALLY K0 HERE IS JUST T(N+1).
! SEE ALSO THE NOTE OF ELMDN1.
!=======================================================================
!=======================================================================
!     <0>(I) -H <0>(J)
!=======================================================================
do II=1,NG
  I=NGT(II)
  do JJ=1,II
    J=NGT(JJ)
    ij_hd = ij_hd + 1
    CALL HGG(I,J,RES1)
    haqqd_tri(ij_hd)=-RES1
  end do
end do

!=======================================================================
!     H(N) PART IN P(K0)P IS CONSIDERED IN ELMDN1
!=======================================================================
IF(NS.EQ.0) GO TO 1200
!=======================================================================
!     <0>(I) H <NA-K>(L)
!     PROJECTION OF EXCITED CONFIGULATIONS <NA-K>(L) ON TO <0> SPACE
!     IS ZERO.   SO, ONLY -H/(N+1) HAS NON-ZERO VALUE.
!=======================================================================
 
do JJ=1,NS
  NA=NHT(JJ)
  J=NPT(JJ)
  do KK=1,NG2(JJ)
    K=NGT2(KK,JJ)
    JKD1=NUMSD(KK,JJ,1)
    JKD2=NUMSD(KK,JJ,2)
    IF(KEYQ.ne.0) JKQ=NUMSQ(KK,JJ)
    if((jkd1.ne.0).or.(jkd2.ne.0).or.(jkq.ne.0)) then
      do II=1,NG
        I=NGT(II)
        CALL HGS(-NA,-J,K,I,RES22)
        CALL HGS( NA, J,K,I,RES12)
        IF(JKD1.ne.0) then
          ij_hd = n_diag(ii,jkd1)
          haqqd_tri(ij_hd)=-(RES12+RES22)*RRT2
        end if
        IF(JKD2.ne.0) then
          CALL HGS(-NA,J,-K,I,RES32)
          ij_hd = n_diag(ii,jkd2)
          haqqd_tri(ij_hd)=-(-(RES12-RES22)-2.0D0*RES32)*RRT6
        end if
        if((keyq.ne.0) .and. (jkq.ne.0)) then
          ij_hq = n_diag(ii,jkq)
          haqqq_tri(ij_hq)=-(-(RES12-RES22)+RES32)*RRT3
        end if
      end do
    end if
  end do
end do

!=======================================================================
!     <NA-I>(J) -H <NB-K>(L)
!=======================================================================

!***
! There were a number of useless assignments to KEY1 and KEY2. These are
! now carried outside the DO 120 loop. Note that always KEY1=KEY2=1 in
! practice and the correlation configurations are always normalized with
! 1/sqrt(2.0d0).
key1 = 1
key2 = 1

do II=1,NS
  NA=NHT(II)
  I=NPT(II)
  do JJ=1,NG2(II)
    J=NGT2(JJ,II)
    do KK=1,II
      NB=NHT(KK)
      K=NPT(KK)
      LMAX=NG2(KK)
      IF(KK.EQ.II) LMAX=JJ
      loop_ll: do LL=1,LMAX
        L=NGT2(LL,KK)
        do M=1,3
          NAM=NA*nsp(M,1)
          IM=I*nsp(M,2)
          JM=J*nsp(M,3)
          IF(IM.eq.JM) then
            RX(M,:)=0.0D0
            cycle 
          end if
          do N=1,3
            NBN=NB*nsp(N,1)
            KN=K*nsp(N,2)
            LN=L*nsp(N,3)
            IF(KN.eq.LN) then
              RX(M,N)=0.0D0
              cycle 
            end if
            CALL HSS(NAM,IM,JM,NBN,KN,LN,RX(M,N))
          end do 
        end do 
 
!       Rodada no subspaco dos dubletos
        do MX=1,2
          IMX=NUMSD(JJ,II,MX)
          IF(IMX.EQ.0) cycle 
          do NX=1,2
            KNX=NUMSD(LL,KK,NX)
            IF(KNX.EQ.0) cycle 
            X=0.0D0
            do M=1,3
              Y=0.0D0
              do N=1,3
                Y=Y+SP(MX,N,KEY1)*RX(N,M)
              end do
              X=X+Y*SP(NX,M,KEY2)
            end do
            ij_hd = n_diag(imx,knx)
            haqqd_tri(ij_hd)=-X
          end do 
        end do 
 
!       Rodada no subspaco dos quartetos
        if(keyq .ne. 0) then
          MX=3
          IMX=NUMSQ(JJ,II)
          IF(IMX.EQ.0) cycle loop_ll
          NX=3
          KNX=NUMSQ(LL,KK)
          IF(KNX.EQ.0) cycle loop_ll
          X=0.0D0
          do M=1,3
            Y=0.0D0
            do N=1,3
              Y=Y+SP(MX,N,KEY1)*RX(N,M)
            end do
            X=X+Y*SP(NX,M,KEY2)
          end do
          ij_hq = n_diag(imx,knx)
          haqqq_tri(ij_hq)=-X
        end if
      end do loop_ll
    end do
  end do
end do

 1200 CONTINUE

!fk
! do i=1,n_corr
! ij_hd = n_diag(c_corr(i),c_corr(i))
! haqqd_tri(ij_hd) = haqqd_tri(ij_hd) * rt2
! do j=1,ndd
! ij_hd = n_diag(c_corr(i),j)
! haqqd_tri(ij_hd) = haqqd_tri(ij_hd) * rt2
! end do
! end do
!fk

if(nwhaqq .ne. 0) then
  open(unit=n_hamil, file='HAM.FILE', form='unformatted', iostat=stat)
  call open_stat('HAM.FILE',stat)
  write(n_hamil) ndd
  write(n_hamil) haqqd_tri
  if( (ns .gt. 0) .and. ( n_corr .gt. 0) ) then
    write(n_hamil) n_corr 
    write(n_hamil) ( c_corr(i), i = 1,n_corr )
  end if
  if(keyq .ne. 0) then
    write(n_hamil) ndq
    write(n_hamil) haqqq_tri
  end if
end if

allocate( haqqd(ndd,ndd), stat=alloc_stat ) 
call allocation_stat('haqqd',alloc_stat)
l = 0
do i=1,ndd
  do j = 1,i
    l = l + 1
    haqqd(i,j) = haqqd_tri(l) 
    haqqd(j,i) = haqqd(i,j)
  end do
end do
deallocate( haqqd_tri )

if(keyq .ne. 0) then
  allocate( haqqq(ndq,ndq), stat=alloc_stat ) 
  call allocation_stat('haqqq',alloc_stat)
  l = 0
  do i=1,ndq
    do j = 1,i
      l = l + 1
      haqqq(i,j) = haqqq_tri(l)   
      haqqq(j,i) = haqqq(i,j)
    end do
  end do
  deallocate( haqqq_tri )
end if

end subroutine elmdn0mod_electron
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine hgg(I,J,RES)

      use poly_gaus, only: eg

      implicit none

!=======================================================================
!
!       <0>(I) -H <0>(J)
!      THE COEFFICIENT OF THE ANTISYMMETRIZER FOR N+1 BODY SYSTEM
!      IS 1.0D0/DSQRT((N+1)!)
!
!=======================================================================

      integer, intent(in) :: i, j
      real*8, intent(out) :: res
      
      RES=-F(IABS(I),IABS(J))
      IF(I.EQ.J) RES=RES-EG(1)

      end subroutine hgg
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine hgs(NA,I,J,K,RES)

      implicit none

!=======================================================================
!
!       <NA-I>(J)  -H  <0>(K)
!
!=======================================================================

      integer, intent(in) :: na, i, j, k
      real*8, intent(out) :: res

!***
!  Local variables
!
      integer :: ia, ja, ka, naa, iprt, ix
!***

      RES=0.0D0
      IF(I.EQ.J) RETURN
        IF(K.EQ.J) GO TO 300
        IF(K.EQ.I) GO TO 200
      IA=IABS(I)
      JA=IABS(J)
      KA=IABS(K)
      NAA=IABS(NA)
      IF(I*NA.GT.0) RES=-RPMOL(IA,NAA,JA,KA)
      IF(I*K.GT.0) RES=RES+RPMOL(IA,KA,JA,NAA)
      RETURN
!=======================================================================
!
!      <NA-J>(K) & <0>(K)
!
!=======================================================================
  200   IPRT=-1
        IA=IABS(J)
      IX=J
        GO TO 420
!=======================================================================
!
!       <0>(K) & <NA-I>(K)
!
!=======================================================================
  300   IPRT=1
! 400   IA=IABS(I)
        IA=IABS(I)
      IX=I
  420    NAA=IABS(NA)
      KA=IABS(K)
      RES=-(F(NAA,IA)+RPMOL(NAA,IA,KA,KA))
      IF(IX*K.GT.0) RES=RES+RPMOL(NAA,KA,IA,KA)
! 500   RES=RES*IPRT
        RES=RES*IPRT

      end subroutine hgs
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine hss(NA,II,JJ,NB,K,L,RES)

      use poly_gaus, only: eg

      implicit none

!=======================================================================
!
!      <NA-II>(JJ) -H <NB-K>(L)
!
!=======================================================================

      integer, intent(in) :: na, ii, jj, nb, k, l
      real*8, intent(out) :: res

!***
!  Local variables
!
      integer :: iter, iprt, j, i, kk, ia, nba
      integer :: naa, ka, la, ja
!***

      RES=0.0D0
      IF(II.EQ.JJ.OR.K.EQ.L) RETURN
      ITER=0
      IPRT=1
      J=JJ
      I=II
  111 CONTINUE
        IF(J.EQ.L) GO TO 100
        IF(J.EQ.K) GO TO 200
      IF(ITER.NE.0) GO TO 112
      J=II
      I=JJ
      IPRT=-IPRT
      ITER=1
      GO TO 111
  112 CONTINUE
        IF(NA.EQ.NB) GO TO 50
        RETURN
  200   IPRT=-IPRT
        KK=L
      GO TO 220
  100      KK=K
!=======================================================================
!
!       <NA-I>(J) & <NB-KK>(J)
!
!=======================================================================
  220      IF(I.EQ.KK) GO TO 300
        IF(NA.EQ.NB) GO TO 350
!=======================================================================
!
!      <NA-I>(J) & <NB-KK>(J)
!
!=======================================================================
      IA=IABS(I)
      NBA=IABS(NB)
      NAA=IABS(NA)
      KA=IABS(KK)
      IF(I*NA.GT.0) RES=RPMOL(IA,NAA,NBA,KA)
      IF(I*KK.GT.0) RES=RES-RPMOL(IA,KA,NBA,NAA)
      GO TO 800
!=======================================================================
!
!      <NA-I>(J) & <NA-KK>(J)
!
!=======================================================================
  350      IA=IABS(I)
      KA=IABS(KK)
      NAA=IABS(NA)
      JA=IABS(J)
      RES=F(IA,KA)-RPMOL(IA,KA,NAA,NAA)+RPMOL(IA,KA,JA,JA)
      IF(I*J.GT.0) RES=RES-RPMOL(IA,JA,JA,KA)
      IF(I*NA.GT.0) RES=RES+RPMOL(IA,NAA,NAA,KA)
      GO TO 800
   50   CONTINUE
!=======================================================================
!
!       <NA-I>(J) & <NA-K>(L)
!
!=======================================================================
        IA=IABS(I)
        JA=IABS(J)
        KA=IABS(K)
      LA=IABS(L)
        IF(I*K.GT.0) RES=RES+RPMOL(IA,KA,JA,LA)
        IF(I*L.GT.0) RES=RES-RPMOL(IA,LA,JA,KA)
        GO TO 800
  300   IF(NA.EQ.NB) GO TO 400
!=======================================================================
!
!       <NA-I>(J) & <NB-I>(J)
!
!=======================================================================
        IPRT=-IPRT
!=======================================================================
!
!       (NA,I,J) & (NB,I,J)
!
!=======================================================================
        NAA=IABS(NA)
        NBA=IABS(NB)
        IA=IABS(I)
        JA=IABS(J)
        RES=F(NAA,NBA)+RPMOL(NAA,NBA,IA,IA)+RPMOL(NAA,NBA,JA,JA)
        IF(NA*I.GT.0) RES=RES-RPMOL(NAA,IA,NBA,IA)
        IF(NA*J.GT.0) RES=RES-RPMOL(NAA,JA,NBA,JA)
        GO TO 800
  400   CONTINUE
!=======================================================================
!
!       <NA-I>(J) & <NA-I>(J)
!
!=======================================================================
        NAA=IABS(NA)
        IA=IABS(I)
        JA=IABS(J)
          RES=EG(1)-F(NAA,NAA)+F(IA,IA)-RPMOL(IA,IA,NAA,NAA) &
           +F(JA,JA)-RPMOL(JA,JA,NAA,NAA)+RPMOL(IA,IA,JA,JA)
        IF(I*J.GT.0) RES=RES-RPMOL(IA,JA,IA,JA)
      IF(I*NA.GT.0) RES=RES+RPMOL(IA,NAA,IA,NAA)
      IF(J*NA.GT.0) RES=RES+RPMOL(JA,NAA,JA,NAA)
  800   RES=-RES*IPRT

      end subroutine hss
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine apresult

use poly_gaus, only: nbfns

implicit none

!***
!  Local variables
!
integer :: i, j, k, l, m, n, ii, jj, kk, na, nc, ll, klx, klx1, klx2, imx, mx
real*8  :: APRESS, APREST, APS, APT, APS1,APS2,APRESS1,APRESS2, APT1,APT2,APREST1,APREST2
integer :: alloc_stat
!***

!=======================================================================
!     NSA -> DIMENSION OF THE ACTIVE SPACE
!=======================================================================

allocate( hci2(nsa,nsa,2) )
if( nchls.ne.0 ) then
  do l=1,nsa
    do k=1,nsa
      hci2(k,l,1) = sum( hci(k,1:nchls,1)*hci(l,1:nchls,1) )
    end do
  end do
end if
if( nchlt.ne.0 ) then
  do l=1,nsa
    do k=1,nsa
      hci2(k,l,2) = sum( hci(k,1:nchlt,2)*hci(l,1:nchlt,2) )
    end do
  end do
end if

allocate( eaqqd(ndd,ndd), stat=alloc_stat )
call allocation_stat('eaqqd',alloc_stat)
eaqqd = 0d0

if( nroute .eq. 5 ) return

allocate( ngtinv(nbfns) )
ngtinv=0
do I=1,NG
  II=NGT(I)
  NGTINV(II)=I
end do
if( ns .gt. 0) then 
  allocate( numinv(nocc,nbfns), ngtin2(nbfns,ns) )
  numinv=0
  ngtin2=0
  do I=1,NS
    NA=NHT(I)
    II=NPT(I)
    NUMINV(NA,II)=I
    do J=1,NG2(I)
      JJ=NGT2(J,I)
      NGTIN2(JJ,I)=J
    end do
  end do
end if

!***
! This part was modified to work with triangular matrices
! according to the convention: do i=1,n / do j=1,i
!
!     DO 1 I=1,NDD
!     DO 5 J=NG+1,NDD
!   ij_ed = n_diag(i,j)
!   eaqqd_tri(ij_ed)=eaqqd_tri(ij_ed)-HAQQD_tri(ij_ed)
!   5 CONTINUE
!   1 CONTINUE

!=======================================================================
!     INICIO DA RODADA NO SUB-ESPACO DOS DUBLETOS (ELASTICO E EXCITACOES)
!=======================================================================

do j=ng+1,ndd
  eaqqd(:,j)=eaqqd(:,j)-haqqd(:,j)
end do

!-----------------------------------------------------------------------
! Caracterizacao da configuracao | 0(N) > : estado fundamental
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Caracterizacao da excitacao:  P0 = | 0 >< 0 |
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Caracterizacao das excitacoes:  P = | NA-K >< NB-L |
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Inicio do loop sobre os canais abertos tipo singleto
!-----------------------------------------------------------------------
if(nchls.ne.0) then
  do m = 1,ng
    N=NGT(M)
    do K=1,NSA
      do L=1,NSA
        CALL APSG(L,K,N,APRESS,J)
        if(apress.eq.0d0 .or. j.eq.0) cycle 
        aps = apress * hci2(k,l,1)
        eaqqd(:,m)=eaqqd(:,m)+haqqd(:,j)*APS
      end do 
    end do
  end do
end if

!-----------------------------------------------------------------------
! Inicio do loop sobre os canais abertos tipo tripleto
!-----------------------------------------------------------------------
if(nchlt.ne.0) then
  do m = 1,ng
    N=NGT(M)
    do K=1,NSA
      do L=1,NSA
        CALL APTG(L,K,N,APREST,J)
        if(aprest.eq.0d0 .or. j.eq.0) cycle
        apt = aprest * hci2(k,l,2)
        eaqqd(:,m)=eaqqd(:,m)+haqqd(:,j)*apt
      end do
    end do
  end do
end if

!-----------------------------------------------------------------------
! Caracterizacao das configuracoes | NC-N(M) > : estados excitados
!-----------------------------------------------------------------------
do II=1,NS
  NC=NHT(II)
  N=NPT(II)
  do J=1,NG2(II)
    M=NGT2(J,II)
    do MX=1,2
      IMX=NUMSD(J,II,MX)
      IF(IMX.EQ.0) cycle
!-----------------------------------------------------------------------
! Caracterizacao das excitacoes:  P = | NA-K >< NB-L |
!-----------------------------------------------------------------------
      do KK=1,NSA
        do LL=1,NSA
          if(nchls.ne.0) then
!-----------------------------------------------------------------------
! Inicio do loop sobre os canais abertos tipo singleto
!-----------------------------------------------------------------------
            CALL APSD(II,LL,KK,NC,N,M,MX,APRESS1,APRESS2,KLX1,KLX2,IMX)
            aps = hci2(kk,ll,1)
            if((KLX1.NE.0).and.(apress1.ne.0d0)) then
              APS1=APS *APRESS1
              eaqqd(:,imx)=eaqqd(:,imx)+haqqd(:,klx1)* aps1
            endif
            if((KLX2.NE.0).and.(apress2.ne.0d0)) then
              APS2=APS *APRESS2
              eaqqd(:,imx)=eaqqd(:,imx)+haqqd(:,klx2)* aps2
            endif
          end if
!-----------------------------------------------------------------------
! Inicio do loop sobre os canais abertos tipo tripleto
!-----------------------------------------------------------------------
          if(nchlt.ne.0) then
            CALL APTD(II,LL,KK,NC,N,M,MX,APREST1,APREST2,KLX1,KLX2,IMX)
            apt = hci2(kk,ll,2)
            if((KLX1.NE.0).and.(aprest1.ne.0d0)) then
              APT1=APT *APREST1
              eaqqd(:,imx)=eaqqd(:,imx)+haqqd(:,klx1)* apt1
            endif
            if((KLX2.NE.0).and.(aprest2.ne.0d0)) then
              APT2=APT *APREST2
              eaqqd(:,imx)=eaqqd(:,imx)+haqqd(:,klx2)* apt2
            endif
          end if
        end do
      end do
    end do
  end do
end do

deallocate( haqqd )

!***
!fk
!symmetrize
!!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j)
!     do i = 1,ndd
! do j = 1,i
!EAQQD_2d(I,J)=0.5d0*(EAQQD_2d(I,J)+EAQQD_2d(J,I))
!EAQQD_2d(J,I)=EAQQD_2d(I,J)
! end do
!     end do
!!$OMP END PARALLEL DO
!  
!!  The result of this loop was 1) Keep diagonal elements as before;
!!  2) keep elements within (1:NG) and (NG+1:NDD) blocks as before; 
!!  3) add up contributions from the off-diagonal elemets (those
!!  coupling the two blocks) and divide by 2. In this version, the off
!!  diagonal elements are already summed, so we only divide by 2. 
!     DO 95 I=1,NDD
!     DO 96 J=1,I
!     EAQQD(I,J)=(EAQQD(I,J)+EAQQD(J,I))/2.0D0
!     EAQQD(J,I)=EAQQD(I,J)
!  96 CONTINUE
!  95 CONTINUE

!=======================================================================
! INICIO DA RODADA NO SUB-ESPACO DOS QUARTETOS (EXCITACOES)
!=======================================================================

if(keyq.ne.0 .and. nchlt.ne.0) then

allocate( eaqqq(ndq,ndq), stat=alloc_stat )
call allocation_stat('eaqqq',alloc_stat)
eaqqq = - haqqq
 
do II=1,NS
  NC=NHT(II)
  N=NPT(II)
  do J=1,NG2(II)
    M=NGT2(J,II)
    IMX=NUMSQ(J,II)
    if(IMX.EQ.0) cycle
    do KK=1,NSA
      do LL=1,NSA
        CALL APTQ(II,LL,KK,NC,N,M,APREST,KLX,IMX)
        if(KLX.NE.0) then
          apt = aprest * hci2(kk,ll,2)
          eaqqq(:,imx)=eaqqq(:,imx)+haqqq(:,klx)* apt
        endif
      end do
    end do
  end do
end do

deallocate( haqqq )

end if

! The factor FAC_Q was introduced to account for the correction 
! introduced by the loops DO 195 and DO 196 below
!     EAQQQ(I,IMX)=EAQQQ(I,IMX)+HAQQQ(I,KLX)*APT*XNFC
!     DO 170 I=1,NDQ
!     ij_eq = n_diag(i,imx)
!     ij_hq = n_diag(i,klx)
!!    eaqqq_tri(ij_eq)=eaqqq_tri(ij_eq)+haqqq_tri(ij_hq)*APT*XNFC
!     eaqqq_tri(ij_eq)=eaqqq_tri(ij_eq)+haqqq_tri(ij_hq)*APT
! 170 CONTINUE

!  This was accounted for by the factor fac_q above
!  (i.e., the factor 0.5 in the elements coupling
!  the (1:NG) and (NG+1:NDQ) blocks)
!  DO 195 I=1,NDQ
!  DO 196 J=1,I
!  EAQQQ(I,J)=(EAQQQ(I,J)+EAQQQ(J,I))/2.0D0
!  EAQQQ(J,I)=EAQQQ(I,J)
!    196 CONTINUE
!    195 CONTINUE

!  Deallocate local arrays
!fk Well, not yet, it will be used in elmb
!     deallocate( ngtinv, numinv )
!     if( allocated(ngtin2) ) deallocate( ngtin2 )

end subroutine apresult
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine apsg(LL,KK,N,APRESS,J)

      implicit none

      integer, intent(in) :: ll, kk, n
      integer, intent(out) :: j
      real*8, intent(out)  :: apress

!***
!  Local variables
!
      integer :: k, l
      integer :: na, nb
!***

!=======================================================================
!
!   A_N+1 x P1 x | 0(N) >  =  A_N+1  | NA-K > < NB-L | 0(N) >
!
!=======================================================================
!
      NA=NHT(KK)
      K=NPT(KK)
!
      NB=NHT(LL)
      L=NPT(LL)
!
      IF(N.NE.L) GO TO 200
      IF(NA.EQ.NB) GO TO 260
      IF(NA.NE.NB) GO TO 200
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | 0(N) >  =  A_N+1  | NA-K > < NB-L | 0(N) >
!
!-----------------------------------------------------------------------      
!
  260 APRESS=0.5D0
      J=NGTINV(K)
      GO TO 300
 
  200 APRESS=0.0D0
 
  300 CONTINUE

      end subroutine apsg
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine aptg(LL,KK,N,APREST,J)

      implicit none

      integer, intent(in) :: ll, kk, n
      integer, intent(out) :: j
      real*8, intent(out)  :: aprest

!***
!  Local variables
!
      integer :: k, l
      integer :: na, nb
!***

!=======================================================================
!
!   A_N+1 x (P2+P3) x | 0(N) >  =  A_N+1  | NA-K > < NB-L | 0(N) >
!
!=======================================================================
!      

      NA=NHT(KK)
      K=NPT(KK)
!
      NB=NHT(LL)
      L=NPT(LL)
!
      IF(N.NE.L) GO TO 400
      IF(NA.EQ.NB.AND.K.EQ.L) GO TO 460
      IF(NA.EQ.NB.AND.K.NE.L) GO TO 470
      IF(NA.NE.NB.AND.K.EQ.L) GO TO 400
      IF(NA.NE.NB.AND.K.NE.L) GO TO 400
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | 0(N) >  =  A_N+1  | NA-K > < NA-L | 0(N) >
!
!-----------------------------------------------------------------------      
!
  470 APREST=1.5D0
      J=NGTINV(K)
      GO TO 500
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | 0(N) >  =  A_N+1  | NA-K > < NA-K | 0(N) >
!
!-----------------------------------------------------------------------      
!
  460 APREST=1.5D0
      J=NGTINV(K)
      GO TO 500
!
  400 APREST=0.0D0
!
  500 CONTINUE

      end subroutine aptg
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine apsd(I,LL,KK,NC,N,M,MX,APRESS1,APRESS2,KLX1,KLX2,IMX)

      implicit none

      integer, intent(in) :: i, ll, kk, nc, n, m, mx, imx
      integer, intent(out) :: klx1, klx2
      real*8, intent(out)  :: apress1, apress2

!***
!  Local variables
!
      integer :: j, k, l
      integer :: na, nb, kkinv, jinv
!***

!=======================================================================  
!      
      NA=NHT(KK)
      K=NPT(KK)
!
      NB=NHT(LL)
      L=NPT(LL)
!
      IF(NB.NE.NC) GO TO 777
      IF(NA.EQ.NB.AND.K.EQ.L) GO TO 460
      IF(NA.EQ.NB.AND.K.NE.L) GO TO 470
      IF(NA.NE.NB.AND.K.EQ.L) GO TO 480
      IF(NA.NE.NB.AND.K.NE.L) GO TO 490 
!
!=======================================================================
!     NA.NE.NB   AND   K.NE.L
!=======================================================================
!
!    490 write(6,*)'Passei pelo loop 490'
  490 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 777 
 
      IF(L.EQ.M) GO TO 400
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-L(M) >  =  A_N+1  | NA-K > < NB-L | NB-L(M) >
!
!-----------------------------------------------------------------------      
!
      J=NGTIN2(M,KK)
 
      GO TO (70,80),MX
 
   70 APRESS1=1.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=-0.5D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.5D0*RT3
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
 
   80 APRESS1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-N(L) >  =  A_N+1  | NA-K > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
  400 J=NGTIN2(N,KK)
 
      GO TO (71,81),MX
 
   71 APRESS1=-0.5D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=0.25D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=-0.25D0*RT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
 
   81 APRESS1=0.5D0*RT3
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=-0.75D0*RRT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.75D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 600
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-L(L) >  =  A_N+1  | NA-K > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
      J=NGTIN2(L,KK)
 
      GO TO (72,82),MX
 
   72 APRESS1=0.5D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=-0.25D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.25D0*RT3
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
 
   82 APRESS1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=0.0D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.0D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
!
!=======================================================================
!     NA.NE.NB   AND   K.EQ.L
!=======================================================================
!
!    480 write(6,*)'Passei pelo loop 480'
  480 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 777 
 
      IF(L.EQ.M) GO TO 410
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-L(M) >  =  A_N+1  | NA-L > < NB-L | NB-L(M) >
!
!-----------------------------------------------------------------------      
!
      J=NGTIN2(M,KK)
 
      GO TO (73,83),MX
 
   73 APRESS1=1.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=-0.5D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.5D0*RT3
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
 
   83 APRESS1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-N(L) >  =  A_N+1  | NA-L > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
  410 J=NGTIN2(N,KK)
 
      GO TO (74,84),MX
 
   74 APRESS1=-0.5D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=0.25D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=-0.25D0*RT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
 
   84 APRESS1=0.5D0*RT3
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=-0.75D0*RRT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.75D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 600
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-L(L) >  =  A_N+1  | NA-L > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
      J=NGTIN2(L,KK)
 
      GO TO (75,85),MX
 
   75 APRESS1=0.5D0
      KLX1=NUMSD(J,KK,1)
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
      GO TO 600
 
   85 APRESS1=0.0D0
      KLX1=NUMSD(J,KK,1)
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
      GO TO 600
!
!=======================================================================
!     NA.EQ.NB   AND   K.NE.L
!=======================================================================
!
!    470 if(imx.eq.6) write(6,*)'Passei pelo loop 470'
  470 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 777 
 
      IF(L.EQ.M) GO TO 420
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-L(M) >  =  A_N+1  | NB-K > < NB-L | NB-L(M) >
!
!-----------------------------------------------------------------------      
!
!      if(imx.eq.6) write(6,*)'Passei pelo CASO (III-B-singleto)'
!
      J=NGTIN2(M,KK)
 
      GO TO (76,86),MX
 
   76 APRESS1=1.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=-0.5D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.5D0*RT3
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
 
   86 APRESS1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-N(L) >  =  A_N+1  | NB-K > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
!    420 if(imx.eq.6) write(6,*)'Passei pelo CASO (III-A-singleto)'
  420 continue
 
      J=NGTIN2(N,KK)
 
      GO TO (77,87),MX
 
   77 APRESS1=-0.5D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=0.25D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=-0.25D0*RT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 600
      GO TO 425
 
   87 APRESS1=0.5D0*RT3
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=-0.75D0*RRT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(N.EQ.K) GO TO 415
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.75D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
  415 IF(L.NE.N) GO TO 600
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-L(L) >  =  A_N+1  | NB-K > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
!    425 if(imx.eq.6) write(6,*)'Passei pelo CASO (III-C-singleto)'
  425 continue
 
      J=NGTIN2(L,KK)
 
      GO TO (78,88),MX
 
   78 APRESS1=0.5D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=-0.25D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.25D0*RT3
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
 
   88 APRESS1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=0.0D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.0D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
!
!=======================================================================
!     NA.EQ.NB   AND   K.EQ.L
!=======================================================================
!
!    460 if(imx.eq.6) write(6,*)'Passei pelo loop 460'
  460 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 777 
 
      IF(L.EQ.M) GO TO 430
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-L(M) >  =  A_N+1  | NA-K > < NB-L | NC-N(M) >
!
!   CASO B: NA=NB, K=L, K=N
!
!-----------------------------------------------------------------------      
!
!        if(imx.eq.6) write(6,*)'Passei pelo CASO (IV-B-singleto)'
!
      J=NGTIN2(M,KK)
 
      GO TO (79,89),MX
 
   79 APRESS1=1.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=-0.5D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.5D0*RT3
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
 
   89 APRESS1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 600
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-N(L) >  =  A_N+1  | NA-K > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
!    430 if(imx.eq.6) write(6,*)'Passei pelo CASO (IV-A-singleto)'
  430 continue
 
      J=NGTIN2(N,KK)
 
      GO TO (90,100),MX
 
   90 APRESS1=-0.5D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=0.25D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=-0.25D0*RT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 600
      GO TO 440
 
  100 APRESS1=0.5D0*RT3
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APRESS1=-0.75D0*RRT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APRESS2=0.75D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 600
!
!-----------------------------------------------------------------------
!
!   A_N+1 x P1 x | NB-L(L) >  =  A_N+1  | NB-L > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
!    440 if(imx.eq.6) write(6,*)'Passei pelo CASO (IV-C-singleto)'
  440 continue
 
      J=NGTIN2(L,KK)
 
      GO TO (91,101),MX
 
   91 APRESS1=0.5D0
      KLX1=NUMSD(J,KK,1)
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      GO TO 600
 
  101 APRESS1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      APRESS2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      GO TO 600
!
!=======================================================================
!
  777 APRESS1=0.0D0
      APRESS2=0.0D0
 
  600 CONTINUE
 
      end subroutine apsd
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine aptd(I,LL,KK,NC,N,M,MX,APREST1,APREST2,KLX1,KLX2,IMX)

      implicit none

      integer, intent(in) :: i, ll, kk, nc, n, m, mx, imx
      integer, intent(out) :: klx1, klx2
      real*8, intent(out)  :: aprest1, aprest2

!***
!  Local variables
!
      integer :: j, k, l
      integer :: na, nb, kkinv, jinv
!***

!=======================================================================
!
      NA=NHT(KK)
      K=NPT(KK)
!
      NB=NHT(LL)
      L=NPT(LL)
!
      IF(NC.NE.NB) GO TO 888
      IF(NA.EQ.NB.AND.K.EQ.L) GO TO 560
      IF(NA.EQ.NB.AND.K.NE.L) GO TO 570
      IF(NA.NE.NB.AND.K.EQ.L) GO TO 580
      IF(NA.NE.NB.AND.K.NE.L) GO TO 590 
!
!=======================================================================
!     NA.NE.NB   AND   K.NE.L
!=======================================================================
!
!    590 write(6,*)'Passei pelo loop 590'
  590 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 888 
 
      IF(L.EQ.M) GO TO 500
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(M) >  =  A_N+1  | NA-K > < NB-L | NB-L(M) >
!
!-----------------------------------------------------------------------      
!
      J=NGTIN2(M,KK)
 
      GO TO (170,180),MX
 
  170 APREST1=0.0D0
 
      IF(M.EQ.K) THEN
      APREST1=0.0D0
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.0D0
 
      IF(M.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
 
  180 APREST1=0.0D0
 
      IF(M.EQ.K) THEN
      APREST1=RT3
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=1.5D0*RRT3
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=1.0D0
 
      IF(M.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.5D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
 
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-N(L) >  =  A_N+1  | NA-K > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
  500 J=NGTIN2(N,KK)
 
      GO TO (171,181),MX
 
  171 APREST1=0.0D0
 
      IF(N.EQ.K) THEN
      APREST1=1.5D0
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.75D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.5D0*RT3
 
      IF(N.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.25D0*RT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
 
  181 APREST1=0.0D0
 
      IF(N.EQ.K) THEN
      APREST1=0.5D0*RT3
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.75D0*RRT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.5D0
 
      IF(N.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.25D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 700
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(L) >  =  A_N+1  | NA-K > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
      J=NGTIN2(L,KK)
 
      GO TO (172,182),MX
 
  172 APREST1=0.0D0
 
      IF(L.EQ.K) THEN
      APREST1=1.5D0
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.75D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.5D0*RT3
 
      IF(L.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST1=0.25D0*RT3
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
 
  182 APREST1=0.0D0
 
      IF(L.EQ.K) THEN
      APREST1=0.0D0
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.0D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.0D0
 
      IF(L.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST1=0.0D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
!
!=======================================================================
!     NA.NE.NB   AND   K.EQ.L
!=======================================================================
!
!    580 write(6,*)'Passei pelo loop 580'
  580 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 888 
 
      IF(L.EQ.M) GO TO 510
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(M) >  =  A_N+1  | NA-L > < NB-L | NB-L(M) >
!
!-----------------------------------------------------------------------      
!
      J=NGTIN2(M,KK)
 
      GO TO (173,183),MX
 
  173 APREST1=0.0D0
 
      IF(M.EQ.K) THEN
      APREST1=0.0D0
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.0D0
 
      IF(M.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
 
  183 APREST1=0.0D0
 
      IF(M.EQ.K) THEN
      APREST1=RT3
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=1.5D0*RRT3
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=1.0D0
 
      IF(M.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.5D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-N(L) >  =  A_N+1  | NA-L > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
  510 J=NGTIN2(N,KK)
 
      GO TO (174,184),MX
 
  174 APREST1=0.0D0
 
      IF(N.EQ.K) THEN
      APREST1=1.5D0
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.75D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.5D0*RT3
 
      IF(N.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2) 
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.25D0*RT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
 
  184 APREST1=0.0D0
 
      IF(N.EQ.K) THEN
      APREST1=0.5D0*RT3
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.75D0*RRT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.5D0
 
      IF(N.EQ.K) THEN
      APREST1=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.25D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 700
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(L) >  =  A_N+1  | NA-L > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
      J=NGTIN2(L,KK)
 
      GO TO (175,185),MX
 
  175 APREST1=1.5D0
      KLX1=NUMSD(J,KK,1)
 
      APREST2=0.0D0
      KLX2=NUMSD(J,KK,2) 
 
      GO TO 700
 
  185 APREST1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      APREST2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      GO TO 700
!
!=======================================================================
!     NA.EQ.NB   AND   K.NE.L
!=======================================================================
!
!    570 if(imx.eq.6) write(6,*)'Passei pelo loop 570'
  570 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 888 
 
      IF(L.EQ.M) GO TO 520
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(M) >  =  A_N+1  | NB-K > < NB-L | NB-L(M) >
!
!-----------------------------------------------------------------------      
!
      if(imx.eq.6) write(6,*)'Passei pelo CASO (III-B-tripleto)'
 
      J=NGTIN2(M,KK)
 
      GO TO (176,186),MX
 
  176 APREST1=0.0D0
 
      IF(M.EQ.K) THEN
      APREST1=0.0D0
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.0D0
 
      IF(M.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
 
  186 APREST1=0.0D0
 
      IF(M.EQ.K) THEN
      APREST1=RT3
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=1.5D0*RRT3
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=1.0D0
 
      IF(M.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.5D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-N(L) >  =  A_N+1  | NB-K > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
  520 if(imx.eq.6) write(6,*)'Passei pelo CASO (III-A-tripleto)'
 
      J=NGTIN2(N,KK)
 
      GO TO (177,187),MX
 
  177 APREST1=0.0D0
 
      IF(N.EQ.K) THEN
      APREST1=1.5D0
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.75D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.5D0*RT3
 
      IF(N.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.25D0*RT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 700
      GO TO 525 
 
  187 APREST1=0.0D0
 
      IF(N.EQ.K) THEN
      APREST1=0.5D0*RT3
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.75D0*RRT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.5D0
 
      IF(N.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.25D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 700
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(L) >  =  A_N+1  | NB-K > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
!    525 if(imx.eq.6) write(6,*) 'Passei pelo CASO (III-C-tripleto)'
  525 continue
 
      J=NGTIN2(L,KK)
 
      GO TO (178,188),MX
 
  178 APREST1=0.0D0
 
      IF(L.EQ.K) THEN
      APREST1=1.5D0
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.75D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.5D0*RT3
 
      IF(L.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.25D0*RT3
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
 
  188 APREST1=0.0D0
 
      IF(L.EQ.K) THEN
      APREST1=0.0D0
      ENDIF
 
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.0D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.0D0
 
      IF(L.EQ.K) THEN
      APREST2=0.0D0
      ENDIF
 
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.0D0
      KKINV=NUMINV(NA,L)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
!
!=======================================================================
!     NA.EQ.NB   AND   K.EQ.L
!=======================================================================
!
!    560 if(imx.eq.6) write(6,*)'Passei pelo loop 560'
  560 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 888 
 
      IF(L.EQ.M) GO TO 530
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(M) >  =  A_N+1  | NB-L > < NB-L | NB-L(M) >
!
!-----------------------------------------------------------------------      
!
!        if(imx.eq.6) write(6,*)'Passei pelo CASO (IV-B-tripleto)'
!
      J=NGTIN2(M,KK)
 
      GO TO (179,189),MX
 
  179 APREST1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
 
  189 APREST1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=1.5D0*RRT3
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=1.0D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.5D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      GO TO 700
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-N(L) >  =  A_N+1  | NB-L > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
!    530 if(imx.eq.6) write(6,*)'Passei pelo CASO (IV-A-tripleto)'
  530 continue
 
      J=NGTIN2(N,KK)
 
      GO TO (190,200),MX
 
  190 APREST1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.75D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.5D0*RT3
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.25D0*RT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 700
      GO TO 540
 
  200 APREST1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      IF(KLX1.EQ.0) THEN
      APREST1=0.75D0*RRT3
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX1=NUMSD(JINV,KKINV,1)
      ENDIF
      ENDIF
 
      APREST2=0.5D0
      KLX2=NUMSD(J,KK,2)
 
      IF(KLX2.EQ.0) THEN
      APREST2=0.25D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX2=NUMSD(JINV,KKINV,2)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 700
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(L) >  =  A_N+1  | NB-L > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
!    540 if(imx.eq.6) write(6,*)'Passei pelo CASO (IV-C-tripleto)'
  540 continue
 
      J=NGTIN2(L,KK)
 
      GO TO (191,201),MX
 
  191 APREST1=1.5D0
      KLX1=NUMSD(J,KK,1)
 
      APREST2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      GO TO 700
 
  201 APREST1=0.0D0
      KLX1=NUMSD(J,KK,1)
 
      APREST2=0.0D0
      KLX2=NUMSD(J,KK,2)
 
      GO TO 700
 
!=======================================================================
  
  888 APREST1=0.0D0
      APREST2=0.0D0
 
  700 CONTINUE
 
      end subroutine aptd
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine aptq(I,LL,KK,NC,N,M,APREST,KLX,IMX)

      implicit none

      integer, intent(in) :: i, ll, kk, nc, n, m, imx
      integer, intent(out) :: klx
      real*8, intent(out)  :: aprest

!***
!  Local variables
!
      integer :: j, k, l
      integer :: na, nb, jinv, kkinv
!***
!=======================================================================
!
      NA=NHT(KK)
      K=NPT(KK)
!
      NB=NHT(LL)
      L=NPT(LL)
!
      IF(NC.NE.NB) GO TO 999
      IF(NA.EQ.NB.AND.K.EQ.L) GO TO 660
      IF(NA.EQ.NB.AND.K.NE.L) GO TO 670
      IF(NA.NE.NB.AND.K.EQ.L) GO TO 680
      IF(NA.NE.NB.AND.K.NE.L) GO TO 690 
!
!=======================================================================
!     NA.NE.NB   AND   K.NE.L
!=======================================================================
!
!    690 write(6,*)'Passei pelo loop 690'
  690 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 999 
 
      IF(L.EQ.M) GO TO 600
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(M) >  =  A_N+1  | NA-K > < NB-L | NB-L(M) >
!
!-----------------------------------------------------------------------      
!
      J=NGTIN2(M,KK)
 
      APREST=1.0D0
 
      IF(M.EQ.K) THEN
      APREST=0.0D0
      ENDIF
 
      KLX=NUMSQ(J,KK)
 
      IF(KLX.EQ.0) THEN
      APREST=-1.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX=NUMSQ(JINV,KKINV)
      ENDIF
      ENDIF
 
      GO TO 800
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-N(L) >  =  A_N+1  | NA-K > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
  600 J=NGTIN2(N,KK)
 
      APREST=-1.0D0
 
      IF(N.EQ.K) THEN
      APREST=0.0D0
      ENDIF
 
      KLX=NUMSQ(J,KK)
 
      IF(KLX.EQ.0) THEN
      APREST=1.0D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX=NUMSQ(JINV,KKINV)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 800
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(L) >  =  A_N+1  | NA-K > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
      J=NGTIN2(L,KK)
 
      APREST=0.0D0
      KLX=NUMSQ(J,KK)
      GO TO 800
!
!=======================================================================
!     NA.NE.NB   AND   K.EQ.L
!=======================================================================
!
!    680 write(6,*)'Passei pelo loop 680'
  680 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 999 
 
      IF(L.EQ.M) GO TO 610
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(M) >  =  A_N+1  | NA-L > < NB-L | NB-L(M) >
!
!-----------------------------------------------------------------------      
!
      J=NGTIN2(M,KK)
 
      APREST=1.0D0
 
      IF(M.EQ.K) THEN
      APREST=0.0D0
      ENDIF
 
      KLX=NUMSQ(J,KK)
 
      IF(KLX.EQ.0) THEN
      APREST=-1.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX=NUMSQ(JINV,KKINV)
      ENDIF
      ENDIF
 
      GO TO 800
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-N(L) >  =  A_N+1  | NA-L > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
  610 J=NGTIN2(N,KK)
 
      APREST=-1.0D0
 
      IF(N.EQ.K) THEN
      APREST=0.0D0
      ENDIF
 
      KLX=NUMSQ(J,KK)
 
      IF(KLX.EQ.0) THEN
      APREST=1.0D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX=NUMSQ(JINV,KKINV)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 800
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(L) >  =  A_N+1  | NA-L > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
      J=NGTIN2(L,KK)
 
      APREST=0.0D0
      KLX=NUMSQ(J,KK)
      GO TO 800
!
!=======================================================================
!     NA.EQ.NB   AND   K.NE.L
!=======================================================================
!
!    670 if(imx.eq.6) write(6,*)'Passei pelo loop 670'
  670 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 999 
 
      IF(L.EQ.M) GO TO 620
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(M) >  =  A_N+1  | NB-K > < NB-L | NB-L(M) >
!
!-----------------------------------------------------------------------      
!
      J=NGTIN2(M,KK)
 
      APREST=1.0D0
 
      IF(M.EQ.K) THEN
      APREST=0.0D0
      ENDIF
 
      KLX=NUMSQ(J,KK)
 
      IF(KLX.EQ.0) THEN
      APREST=-1.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX=NUMSQ(JINV,KKINV)
      ENDIF
      ENDIF
 
      GO TO 800
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-N(L) >  =  A_N+1  | NB-K > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
  620 J=NGTIN2(N,KK)
 
      APREST=-1.0D0
 
      IF(N.EQ.K) THEN
      APREST=0.0D0
      ENDIF
 
      KLX=NUMSQ(J,KK)
 
      IF(KLX.EQ.0) THEN
      APREST=1.0D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX=NUMSQ(JINV,KKINV)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 800
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(L) >  =  A_N+1  | NB-K > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
      J=NGTIN2(L,KK)
 
      APREST=0.0D0
      KLX=NUMSQ(J,KK)
 
      GO TO 800
!
!=======================================================================
!     NA.EQ.NB   AND   K.EQ.L
!=======================================================================
!
!    660 if(imx.eq.6) write(6,*)'Passei pelo loop 660'
  660 continue
 
      IF(L.NE.N.AND.L.NE.M) GO TO 999 
 
      IF(L.EQ.M) GO TO 630
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(M) >  =  A_N+1  | NB-L > < NB-L | NB-L(M) >
!
!-----------------------------------------------------------------------      
!
      J=NGTIN2(M,KK)
 
      APREST=1.0D0
      KLX=NUMSQ(J,KK)
 
      IF(KLX.EQ.0) THEN
      APREST=-1.0D0
      KKINV=NUMINV(NA,M)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX=NUMSQ(JINV,KKINV)
      ENDIF
      ENDIF
 
      GO TO 800
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-N(L) >  =  A_N+1  | NB-L > < NB-L | NB-N(L) >
!
!-----------------------------------------------------------------------
!
  630 J=NGTIN2(N,KK)
 
      APREST=-1.0D0
      KLX=NUMSQ(J,KK)
 
      IF(KLX.EQ.0) THEN
      APREST=1.0D0
      KKINV=NUMINV(NA,N)
!mal
      IF(KKINV.NE.0)THEN 
!mal
      JINV=NGTIN2(K,KKINV)
      KLX=NUMSQ(JINV,KKINV)
      ENDIF
      ENDIF
 
      IF(L.NE.N) GO TO 800
!
!-----------------------------------------------------------------------
!
!   A_N+1 x (P2+P3) x | NB-L(L) >  =  A_N+1  | NB-L > < NB-L | NB-L(L) >
!
!-----------------------------------------------------------------------
!
      J=NGTIN2(L,KK)
 
      APREST=0.0D0
      KLX=NUMSQ(J,KK)
 
      GO TO 800
 
!=======================================================================
  
  999 APREST=0.0D0
 
  800 CONTINUE
 
      end subroutine aptq
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine vva00smod_electron

use poly_gaus, only: nbfns

implicit none

!=======================================================================
!     PREPARATION FOR VGV ELEMENTS OF SHORT-RANGE FUNCTIONS.
!=======================================================================
!
!=======================================================================
! Elementos de matriz VVDG(NBFN,NDD) e VVQG(NBFN,NDQ)
!=======================================================================

!***
! Local variables
!
integer :: i, j, m, ii, jj, i1, i2, i3, na1, na2, na3
integer :: ij1, ij2, ij3
real*8  :: r, r1, r2, r3
integer :: nbfn
integer :: alloc_stat, stat
!***

nbfn = nbfns

!=======================================================================
! Rodada no sub-espaco dos dubletos
!=======================================================================

allocate( vvdg(ndd,nbfns), stat=alloc_stat )
call allocation_stat('vvdg',alloc_stat)
vvdg = 0.0d0

!=======================================================================
!     < 0[M] | V | 0(I) >
!=======================================================================
do II=1,NG
  I=NGT(II)
  do M=1,NBFN
    call PVGG(M,I,R)
    VVDG(ii,m)=R
  end do
end do

if(ns.ne.0) then

!=======================================================================
!     < 0[M] | V | NA-K(L) >
!=======================================================================
do II=1,NS
  NA1=NHT(II)
  NA2=-NHT(II)
  NA3=-NHT(II)
  I1=NPT(II)
  I2=-NPT(II)
  I3=NPT(II)
  do JJ=1,NG2(II)
    J=NGT2(JJ,II)
    IJ1=NUMSD(JJ,II,1)
    IJ2=NUMSD(JJ,II,2)
    if((ij1.ne.0).or.(ij2.ne.0)) then
      do M=1,NBFN
        call PVGS(M,NA1,I1,J,R1)
        call PVGS(M,NA2,I2,J,R2)
        call PVGS(M,NA3,I3,-J,R3)
        IF(IJ1.NE.0) VVDG(ij1,m)=(R1+R2)*RRT2
        IF(IJ2.NE.0) VVDG(ij2,m)=(R2-R1-2.0D0*R3)*RRT6
      end do
    end if
  end do
end do

!=======================================================================
! Rodada no sub-espaco dos quartetos
!=======================================================================

if(keyq.ne.0) then
  allocate( vvqg(ndq,nbfns), stat=alloc_stat )
  call allocation_stat('vvdq',alloc_stat)
  vvqg = 0.0d0
  do II=1,NS
    NA1=NHT(II)
    NA2=-NHT(II)
    NA3=-NHT(II)
    I1=NPT(II)
    I2=-NPT(II)
    I3=NPT(II)
    do JJ=1,NG2(II)
      J=NGT2(JJ,II)
      IJ3=NUMSQ(JJ,II)
      do M=1,NBFN
        call PVGS(M,NA1,I1,J,R1)
        call PVGS(M,NA2,I2,J,R2)
        call PVGS(M,NA3,I3,-J,R3)
        IF(IJ3.NE.0) VVQG(ij3,m)=(R2-R1+R3)*RRT3
      end do
    end do
  end do
end if

end if

if(nroute .ne. 0) then
  open(unit=n_vvd_1,file='VVDG.FILE',form='unformatted',iostat=stat)
  call open_stat('VVDG.FILE',stat)
  write(n_vvd_1) nbfns,ndd
  write(n_vvd_1) vvdg
  close(n_vvd_1)
end if

!mod
!     open(unit=1,file='vvdg.lis',form='formatted')
!     do i=1,nbfn
! do j=1,ndd
! ij1 = ij1 + 1
! write(1,"(2(i4),2x,1pe15.8)") i,j,vvdg(i,j)
! end do
!     end do
!     close(1)
!mod
 
end subroutine vva00smod_electron
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine vvamulmod_electron

use poly_gaus, only: nbfns

implicit none

!=======================================================================
! Elementos de matriz VVD(NBFN,NDD,NSA,3) e VVQ(NBFN,NDQ,NSA,3)
!=======================================================================
! MULTI-CHANNEL VERSION OF VGV00S.
!=======================================================================

!***
!  Local variables
!
integer :: i, j, k, l, m, ii, na1, na2, na3, k1, k2, k3
integer :: nb1, nb2, nb3, jj1, jj2, jj3, jl1, jl2, jl3, ll
integer :: i1, i2, i3
real*8  :: r1, r2, r3, r11, r21, r31, r12, r22, r32, r13
real*8  :: r23, r33
integer :: nbfn
integer :: alloc_stat
!***

nbfn = nbfns

allocate ( vvd(ndd,nbfns,nsa,2), stat=alloc_stat )
call allocation_stat('vvd',alloc_stat)
if(keyq .ne. 0) then
  allocate( vvq(ndq,nbfns,nsa,2), stat=alloc_stat )
  call allocation_stat('vvq',alloc_stat)
end if

!=======================================================================
!< NA-K[M] | V | 0(I) >
!=======================================================================

do K=1,NSA
  NA1=NHT(K)
  NA2=-NHT(K)
  NA3=-NHT(K)
  K1=NPT(K)
  K2=-NPT(K)
  K3=NPT(K)
  do M=1,NBFN
    do I=1,NG
      II=NGT(I)
      call PVSG(NA1,K1,M,II,R1)
      call PVSG(NA2,K2,M,II,R2)
      call PVSG(NA3,K3,-M,II,R3)
      VVD(i,m,K,1)=(R1+R2)*RRT2
      VVD(i,m,K,2)=(R1-R2)*RRT2 + r3 * rt2
    end do
  end do
end do

!=======================================================================
!< NA-I[M] | V | NB-J(L) >
!=======================================================================

do I=1,NSA
  NA1=NHT(I)
  NA2=-NHT(I)
  NA3=-NHT(I)
  I1=NPT(I)
  I2=-NPT(I)
  I3=NPT(I)
  do M=1,NBFN
    do J=1,NS
      NB1=NHT(J)
      NB2=-NHT(J)
      NB3=-NHT(J)
      JJ1=NPT(J)
      JJ2=-NPT(J)
      JJ3=NPT(J)
      do L=1,NG2(J)
        JL1=NUMSD(L,J,1)
        JL2=NUMSD(L,J,2)
!=======================================================================
! Rodada no sub-espaco dos dubletos
!=======================================================================
        if((jl1.ne.0).or.(jl2.ne.0)) then
          LL=NGT2(L,J)
          call PVSS(NA1,I1,M,NB1,JJ1,LL,R11)
          call PVSS(NA2,I2,M,NB1,JJ1,LL,R21)
          call PVSS(NA3,I3,-M,NB1,JJ1,LL,R31)
          call PVSS(NA1,I1,M,NB2,JJ2,LL,R12)
          call PVSS(NA2,I2,M,NB2,JJ2,LL,R22)
          call PVSS(NA3,I3,-M,NB2,JJ2,LL,R32)
          call PVSS(NA1,I1,M,NB3,JJ3,-LL,R13)
          call PVSS(NA2,I2,M,NB3,JJ3,-LL,R23)
          call PVSS(NA3,I3,-M,NB3,JJ3,-LL,R33)
          IF(JL1.ne.0) then
            VVD(jl1,m,I,1)=(R11+R12+R21+R22)/2.0D0
            VVD(jl1,m,I,2)=((R11+R12)-(R21+R22))/2.0D0 + (r31+r32)
          end if
          IF(JL2.ne.0) then
            VVD(jl2,m,I,1)=(-R11+R12-2.0D0*R13-R21+R22-2.0D0*R23)*RRT3/2.0D0
            VVD(jl2,m,I,2)=(-R11+R12-2.0D0*R13+R21-R22+2.0D0*R23)*RRT3/2.0D0 + (-R31+R32-2.0D0*R33) * rrt3
          end if
        end if
!=======================================================================
! Rodada no sub-espaco dos quartetos
!=======================================================================
        if(keyq.ne.0) then
          jl3=numsq(l,j)
          IF(JL3.ne.0) then
            VVQ(jl3,m,I,1)=(-R11+R12+R13-R21+R22+R23)*RRT6
            VVQ(jl3,m,I,2)=(-R11+R12+R13+R21-R22-R23)*RRT6
            VVQ(jl3,m,I,3)=(-R31+R32+R33)*RRT3
          end if
        end if
      end do
    end do
  end do
end do

!mal
! if(nroute .ne. 0) then
! open(unit=n_vvd_2,file='VVD.FILE',form='unformatted')
! write(n_vvd_2) nbfns,ndd,nsa
! close(n_vvd_2)
! end if
!mal

end subroutine vvamulmod_electron
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine pvresult

implicit none

!=======================================================================
! THIS SUBROUTINE GIVES THE RESULT OF PROJECTION
! OPERATION ON CONFIGURATIONS
!=======================================================================

!***
!  Local variables
!
integer :: i, j, k, l, m, n, nc, jg, mx, imx, kk, ll, jd, jq
real*8  :: rgs, rgt, rds, rdt, rqt, pvs, pvt
integer :: stat
!-----------------------------------------------------------------------
! Caracterizacao da configuracao < 0(N) | : estado fundamental
!-----------------------------------------------------------------------

! The factor FAC_D was introduced to account for the correction
! introduced by the loops DO 195 and DO 196 below. It is necessary
! because EAQQD already has values assigned in APRESULT, so 
! multiplying by 0.5d0 afterwards would not work properly

do m = 1,ng
  n=ngt(m)
!-----------------------------------------------------------------------
! Caracterizacao das excitacoes:  P0 = | 0 >< 0 |
!-----------------------------------------------------------------------
  eaqqd(:,m)=eaqqd(:,m)+VVDG(:,n)
!-----------------------------------------------------------------------
! Caracterizacao das excitacoes:  P = | NA-K >< NB-L |
!-----------------------------------------------------------------------
  do K=1,NSA
    call FACEPG(K,N,RGS,RGT,JG)
!-----------------------------------------------------------------------
! Inicio do loop sobre os canais abertos tipo singleto
!-----------------------------------------------------------------------
    IF((NCHLS.ne.0) .and. (rgs.ne.0d0)) then
      do L=1,NSA
        pvs = rgs * hci2(k,l,1)
        eaqqd(:,m)=eaqqd(:,m)+VVD(:,jg,l,1)*pvs
      end do
    end if
!-----------------------------------------------------------------------
! Inicio do loop sobre os canais abertos tipo tripleto 1 e 2
!-----------------------------------------------------------------------
    IF((NCHLT.ne.0) .and. (rgt.ne.0d0)) then
      do L=1,NSA
        pvt = rgt * hci2(k,l,2) 
        eaqqd(:,m)=eaqqd(:,m)+VVD(:,jg,L,2)*PVT
      end do
    end if
  end do
end do

!-----------------------------------------------------------------------
! Caracterizacao das configuracoes < NC-N(M) | : estados excitados
!-----------------------------------------------------------------------
do I=1,NS
  NC=NHT(I)
  N=NPT(I)
  do J=1,NG2(I)
    M=NGT2(J,I)
    do MX=1,2
      IMX=NUMSD(J,I,MX)
      IF(IMX.EQ.0) cycle
!-----------------------------------------------------------------------
! Caracterizacao da excitacao:  P0 = | 0 >< 0 |
! A superposicao entre uma configuracao primitiva < NC-N(M) | e o estado
! fundamental do alvo | 0 > eh sempre nula, ou seja, < NC-N(M) | 0 > = 0
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Caracterizacao das excitacoes:  P = | NA-K >< NB-L |
!-----------------------------------------------------------------------
      do KK=1,NSA
        call FACEPD(NC,N,M,MX,KK,RDS,RDT,JD)
!-----------------------------------------------------------------------
! Inicio do loop sobre os canais abertos tipo singleto
!-----------------------------------------------------------------------
        IF((NCHLS.ne.0) .and. (rds.ne.0d0)) then
          do LL=1,NSA
            pvs = rds * hci2(kk,ll,1) 
            eaqqd(:,imx)=eaqqd(:,imx)+VVD(:,jd,ll,1)*pvs
          end do
        end if
!-----------------------------------------------------------------------
! Inicio do loop sobre os canais abertos tipo tripleto 1 e 2
!-----------------------------------------------------------------------
        IF((NCHLT.ne.0) .and. (rdt.ne.0d0)) then
          do LL=1,NSA
            pvt = rdt * hci2(kk,ll,2)
            eaqqd(:,imx)=eaqqd(:,imx)+VVD(:,jd,ll,2)*PVT
          end do
        end if
      end do
    end do
  end do
end do

!!  This was accounted for by the factor FAC_D above
!      DO 95 I=1,NDD
!      DO 96 J=1,I
!      EAQQD(I,J)=(EAQQD(I,J)+EAQQD(J,I))/2.0D0
!      EAQQD(J,I)=EAQQD(I,J)
!   96 CONTINUE
!   95 CONTINUE
!***
!
!WRITE(6,358)
!  358 FORMAT(//' ##### EAQQD MATRIX (PVRESULT) #######')
!CALL PRINT(EAQQD,NDD,NDD,NDD,NDD)
!=======================================================================

!=======================================================================
! INICIO DA RODADA NO SUB-ESPACO DOS QUARTETOS (EXCITACOES)
!=======================================================================
if(keyq.ne.0) then
do I=1,NS
  NC=NHT(I)
  N=NPT(I)
  do J=1,NG2(I)
    M=NGT2(J,I)
    IMX=NUMSQ(J,I)
    IF(IMX.EQ.0) cycle
    do KK=1,NSA
      call FACEPQ(NC,N,M,KK,RQT,JQ)
      IF(rqt.eq.0d0) cycle
      do LL=1,NSA
        pvt = rqt * hci2(kk,ll,2)
        eaqqq(:,imx)=eaqqq(:,imx)+(RT2*VVQ(:,jq,LL,2)-VVQ(:,jq,LL,3))*PVT
      end do
    end do
  end do
end do
end if

!  This was accounted for by the factor fac_q above
!  (i.e., the factor 0.5 in the elements coupling
!  the (1:NG) and (NG+1:NDQ) blocks)
!  DO 195 I=1,NDQ
!  DO 196 J=1,I
!  EAQQQ(I,J)=(EAQQQ(I,J)+EAQQQ(J,I))/2.0D0
!  EAQQQ(J,I)=EAQQQ(I,J)
!    196 CONTINUE
!    195 CONTINUE
!WRITE(6,360)
!  360 FORMAT(//' ##### EAQQQ MATRIX (PVRESULT) #######')
!CALL PRINT(EAQQQ,NDQ,NDQ,NDQ,NDQ)

!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j)
do i = 1,ndd
  do j = 1,i
    eaqqd(i,j)=0.5d0*(eaqqd(i,j)+eaqqd(j,i))
    eaqqd(j,i)=eaqqd(i,j)
  end do
end do
!$OMP END PARALLEL DO

if(npeaqq.ne.0) then
  write(6,778)
  778 FORMAT(//'#### EAQQD MATRIX (ELMDN0)')
  call print(eaqqd,NDD,NDD,NDD,NDD)
  if(keyq.eq.0) return
  write(6,779)
  779 FORMAT(//'#### EAQQD MATRIX (ELMDN0)')
  call print(eaqqq,NDQ,NDQ,NDQ,NDQ)
endif

!   Save 0.5*(PV+VP) in case NROUTE=5
if( nroute .eq. 5 ) then
  open(unit=n_pvvp2,file='PVVP.FILE',form='unformatted',iostat=stat)
  call open_stat('PVVP.FILE',stat)
  write(n_pvvp2) size(eaqqd)
  write(n_pvvp2) eaqqd
  close(n_pvvp2)
end if

end subroutine pvresult
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine facepg(KK,N,RGS,RGT,JG)

      implicit none

!=======================================================================
!
!      OVERLAP BETWEEN TARGET FUNCTIONS(N-BODY) AND GROUND STATE BASIS
!      FUNCTION (DUBLET SUBSPACE OF CONFIGURATIONS)
!
!=======================================================================

      integer, intent(in)  :: kk, n
      integer, intent(out) :: jg
      real*8, intent(out)  :: rgs, rgt

!***
!  Local variables
!
      integer :: k, na
!***

!=======================================================================
!
!      OVERLAP BETWEEN < 0(N) | AND | 0 >
!      THIS IS TRIVIAL; NOT WRITTEN HERE.
!
!      OVERLAP BETWEEN < 0(N) | AND | NA-K >
!
!=======================================================================
!
      RGS=0.0D0
      RGT=0.0D0
!***
!  The original code used to rely on NHT(1)=NPT(1)=0 in a static-
!  exchange run. This would cause problems since the allocation of
!  these arrays is now conditional (NS>0).
!
      if(ns .eq. 0) return
!***
 
      NA=NHT(KK)
      K=NPT(KK)
 
!     IF(N.NE.K) GO TO 10
      if(n.ne.k) return
 
!-----------------------------------------------------------------------
!
!     < 0(K) | NA-K >
!
!-----------------------------------------------------------------------
 
      RGS=-RRT2
      RGT=-RRT2
      JG=NA
 
!     10 CONTINUE

      end subroutine facepg
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine facepd(NC,N,M,MX,KK,RDS,RDT,JD)

      implicit none

!=======================================================================
!
!      OVERLAP BETWEEN TARGET FUNCTIONS(N-BODY) AND BASIS FUNCTIONS
!      AT THE DUBLET SUBSPACE OF CONFIGURATIONS
!
!=======================================================================

      integer, intent(in)  :: nc, n, m, mx, kk
      integer, intent(out) :: jd
      real*8, intent(out)  :: rds, rdt

!***
!  Local variables
!
      integer :: k, na
!***

!=======================================================================
!
      RDS=0.0D0
      RDT=0.0D0
!
      NA=NHT(KK)
!fk   K=NPT(KK)
!
      IF(NA.NE.NC) GO TO 100
      K=NPT(KK)
!
!=======================================================================
!
!      OVERLAP BETWEEN < NC-N(M) | AND | NA-K >
!      NOTE: ALL THE VALUES INDICATING ORBITALS SUCH AS NA,I, AND J
!            ARE POSITIVE IN THIS SUBROUTINE.
!            NEITHER < NA-K | NOR | NC-N(M) > IS A SINGLE DETERMINANT.
!
!=======================================================================
!
      IF(M.EQ.K.AND.N.EQ.K) GO TO 70
      IF(M.EQ.K) GO TO 10
      IF(N.EQ.K) GO TO 40
      GO TO 100
!
!-----------------------------------------------------------------------
!
!     < NA-N(K) | NA-K >
!
!-----------------------------------------------------------------------
!
   10 GO TO (20,30),MX
   20 RDS=-0.5d0
      RDT=-0.5d0
      JD=N
      GO TO 100
   30 RDS=RRT3_x3d2
      RDT=-RRT3_d2
      JD=N
      GO TO 100
!
!-----------------------------------------------------------------------
!
!     < NA-K(M) | NA-K >
!
!-----------------------------------------------------------------------
!
   40 GO TO (50,60),MX
   50 RDS=1d0
      RDT=0.0D0
      JD=M
      GO TO 100
   60 RDS=0.0D0
      RDT=-RRT3
      JD=M
      GO TO 100
!
!-----------------------------------------------------------------------
!
!     < NA-K(K) | NA-K >
!
!-----------------------------------------------------------------------
!
   70 GO TO (80,90),MX
   80 RDS=0.5d0
      RDT=-0.5d0
      JD=K
      GO TO 100
   90 RDS=RRT3_x3d2
      RDT=-RRT3_x3d2
      JD=K
!
!-----------------------------------------------------------------------
!
 100  CONTINUE

      end subroutine facepd
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine facepq(NC,N,M,KK,RQT,JQ)

      implicit none

!=======================================================================
!
!      OVERLAP BETWEEN TARGET FUNCTIONS(N-BODY) AND BASIS FUNCTIONS
!      AT THE QUARTET SUBSPACE OF CONFIGURATIONS
!
!=======================================================================

      integer, intent(in)  :: nc, n, m, kk
      integer, intent(out) :: jq
      real*8, intent(out)  :: rqt

!***
!  Local variables
!
      integer :: k, na
      real*8  :: rdt
!***

      RDT=0.0D0
!
      NA=NHT(KK)
      K=NPT(KK)
!
      IF(NA.NE.NC) GO TO 200
!
!=======================================================================
!
!      OVERLAP BETWEEN < NC-N(M) | AND | NA-K >
!      NOTE: ALL THE VALUES INDICATING ORBITALS SUCH AS NA,I, AND J
!            ARE POSITIVE IN THIS SUBROUTINE.
!            NEITHER < NA-K | NOR | NC-N(M) > IS A SINGLE DETERMINANT.
!
!=======================================================================
!
      IF(M.EQ.K.AND.N.EQ.K) GO TO 130
      IF(M.EQ.K) GO TO 110
      IF(N.EQ.K) GO TO 120
      GO TO 200
!
!-----------------------------------------------------------------------
!
!     < NA-N(K) | NA-K >
!
!-----------------------------------------------------------------------
!
  110 RQT=RRT3
      JQ=N
      GO TO 200
!
!-----------------------------------------------------------------------
!
!     < NA-K(M) | NA-K >
!
!-----------------------------------------------------------------------
!
  120 RQT=-RRT3
      JQ=M
      GO TO 200
!
!-----------------------------------------------------------------------
!
!     < NA-K(K) | NA-K >
!
!-----------------------------------------------------------------------
!
  130 RQT=0.0D0
      JQ=K
!
!-----------------------------------------------------------------------
!
  200 CONTINUE

      end subroutine facepq
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine pvgg(I,J,RES)

      implicit none

!=======================================================================
!
!       <0><I> V <0>(J)
!       <0><I> MEANS THAT I IS NOT ANTISYMMETRIZED WITH <0>.
!       <0>(J) MEANS THAT J IS ANTISYMMETRIZED WITH <0>.
!       I SHOULD NOT BE ONE OF OCCUPIED ORBITALS.
!      THE COEFFICIENT OF ANTISYMMETRIZER = 1.0D0/DSQRT(N!) FOR <0>,
!      THAT OF <0>(J) = 1.0D0/DSQRT((N+1)!).
!
!=======================================================================

      integer, intent(in) :: i, j
      real*8, intent(out) :: res

      RES=F2(IABS(I),IABS(J))

      end subroutine pvgg
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine pvsg(NA,I,J,K,RES)

      implicit none

!=======================================================================
!
!       <NA-I><J> V <0>(K)
!       LOOK AT PVGG
!
!=======================================================================

      integer, intent(in) :: na, i, j, k
      real*8, intent(out) :: res

!***
!  Local variables
!
      integer :: naa, ia, ja, ka
!***

      RES=0.0D0
      IF(J.EQ.0.OR.NA.EQ.0) RETURN
        NAA=IABS(NA)
        IA=IABS(I)
        JA=IABS(J)
        KA=IABS(K)
      IF(I.EQ.K) GO TO 100
        IF(J*K.GT.0) RES=RPMOL(JA,KA,IA,NAA)
        IF(J*NA.GT.0) RES=RES-RPMOL(JA,NAA,IA,KA)
      RETURN
!=======================================================================
!
!      <NA-I><J> V <0>(I)
!
!=======================================================================
  100 RES=-(F2(JA,NAA)+RPMOL(JA,NAA,IA,IA))
      IF(I*J.GT.0) RES=RES+RPMOL(IA,JA,IA,NAA)

      end subroutine pvsg
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine pvgs(I,NA,J,K,RES)

      implicit none

!=======================================================================
!
!       <0><I> V <NA-J>(K)
!       LOOK AT PVGG
!
!=======================================================================

      integer, intent(in) :: i, na, j, k
      real*8, intent(out) :: res

!***
!  Local variables
!
      integer :: naa, ia, ja, ka
!***

      RES=0.0D0
      IF(J.EQ.K) RETURN
      IA=IABS(I)
      NAA=IABS(NA)
      JA=IABS(J)
      KA=IABS(K)
      IF(I*K.GT.0) RES=RPMOL(IA,KA,NAA,JA)
      IF(I*J.GT.0) RES=RES-RPMOL(IA,JA,NAA,KA)

      end subroutine pvgs
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine pvss(NA,I,J,NB,K,L,RES)

      implicit none

!=======================================================================
!
!       <NA-I><J> V <NB-K>(L)
!       SEE PVGG.
!       J CAN BE ONE OF OCCUPIED ORBITALS
!
!=======================================================================

      integer, intent(in) :: na, i, j, nb, k, l
      real*8, intent(out) :: res

!***
!  Local variables
!
      integer :: ia, ja, ka, la, naa, nba
!***

      RES=0.0D0
      IF(NA.EQ.0) RETURN
      IF(J.EQ.0.OR.K.EQ.L) RETURN
      IF(NA.EQ.NB) GO TO 100
      IF(I.EQ.K) GO TO 150
      IF(I.EQ.L) GO TO 160
      RETURN
!=======================================================================
!
!      <NA-I><J> V <NB-I>(L)
!
!=======================================================================
  150 JA=IABS(J)
      LA=IABS(L)
      NBA=IABS(NB)
      NAA=IABS(NA)
      IF(J*L.GT.0) RES=-RPMOL(JA,LA,NBA,NAA)
      IF(J*NA.GT.0) RES=RES+RPMOL(JA,NAA,NBA,LA)
      GO TO 800
!=======================================================================
!
!      <NA-I><J> V <NB-K>(I)
!
!=======================================================================
  160 JA=IABS(J)
      KA=IABS(K)
      NAA=IABS(NA)
      NBA=IABS(NB)
      IF(J*K.GT.0) RES=RPMOL(JA,KA,NBA,NAA)
      IF(J*NA.GT.0) RES=RES-RPMOL(JA,NAA,NBA,KA)
      GO TO 800
!=======================================================================
!
!      <NA-I><J> V <NA-K>(L)
!
!=======================================================================
  100 IA=IABS(I)
      JA=IABS(J)
      KA=IABS(K)
      LA=IABS(L)
      NAA=IABS(NA)
      IF(I.NE.K) GO TO 300
      RES=F2(JA,LA)-RPMOL(JA,LA,NAA,NAA)+RPMOL(JA,LA,IA,IA)
      IF(I*J.GT.0) RES=RES-RPMOL(JA,IA,IA,LA)
      IF(J*NA.GT.0) RES=RES+RPMOL(JA,NAA,LA,NAA)
      GO TO 800
  300    IF(I.NE.L) GO TO 350
      RES=-(F2(JA,KA)-RPMOL(JA,KA,NAA,NAA)+RPMOL(JA,KA,IA,IA))
      IF(J*I.GT.0) RES=RES+RPMOL(JA,IA,IA,KA)
      IF(J*NA.GT.0) RES=RES-RPMOL(JA,NAA,KA,NAA)
      GO TO 800
  350 CONTINUE
      IF(J*L.GT.0) RES=RES+RPMOL(JA,LA,IA,KA)
      IF(J*K.GT.0) RES=RES-RPMOL(JA,KA,IA,LA)
  800 continue

      end subroutine pvss
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine ivo

implicit none

!***
!  Local variables
!
integer :: i, k, n, na, nb, ins, jci, kci
integer :: ich, ichs, icht, nss
real*8  :: res

real*8, allocatable :: del(:)
integer :: stat
!***

!***
if(nsa.ne.0) then
 
  allocate( hci(nsa,nsa,2), egn(nsa) )

  if(nsa.eq.1) then
    HCI(1,1,1)=1.0D0
  else
    open(unit=n_hci_1,file='HCI1.FILE',form='unformatted',iostat=stat)
    call open_stat('HCI1.FILE',stat)
!mal
! rewind(n_hci_1)
!mal
    read(n_hci_1,iostat=stat)NSS
    call write_stat('nss in HCI1.FILE',stat)
    if(nss.ne.nsa) stop 'NSA DIFERE EM HCI1.FILE'
    read(n_hci_1,iostat=stat) egn
    call write_stat('egn in HCI1.FILE',stat)
    read(n_hci_1,iostat=stat) hci(:,:,1)
    call write_stat('hci in HCI1.FILE',stat)
  endif
 
  if(nsa.eq.1) then
    hci(1,1,2)=1.0d0
  else
    open(unit=n_hci_2,file='HCI2.FILE',form='unformatted',iostat=stat)
    call open_stat('HCI2.FILE',stat)
!mal
! rewind(n_hci_2)
!mal
    read(n_hci_2,iostat=stat)NSS
    call write_stat('nss in HCI2.FILE',stat)
    if(nss.ne.nsa) stop 'NSA DIFERE EM HCI2.FILE'
    read(n_hci_2,iostat=stat) egn
    call write_stat('egn in HCI2.FILE',stat)
    read(n_hci_2,iostat=stat) hci(:,:,2)
    call write_stat('hci in HCI2.FILE',stat)
  endif
 
! WRITE(6,777)
! 777   FORMAT(//'#### HCI MATRIX (IVO)')
! call print(HCI,NSA,NSA,NSA,NSA)
 
end if
!***

deg(2:nchl) = 0.0d0

if(nchls.ne.0) then

  N=1
  do JCI=1,NSA
    NA=NHT(JCI)
    I=NPT(JCI)
    do KCI=1,NSA
      NB=NHT(KCI)
      K=NPT(KCI)
      call ALVOC(NA,I,NB,K,RES,N)
      do ICHS=2,NCHLS+1
        INS=ICHS-1
        DEG(ICHS)=DEG(ICHS)+HCI(JCI,INS,1)*HCI(KCI,INS,1)*RES
      end do
    end do
  end do
end if
 
if(nchlt.ne.0) then
 
  N=2
  do JCI=1,NSA
    NA=NHT(JCI)
    I=NPT(JCI)
    do KCI=1,NSA
      NB=NHT(KCI)
      K=NPT(KCI)
      call ALVOC(NA,I,NB,K,RES,N)
      do ICHT=NCHLS+2,NCHL
        INS=ICHT-NCHLS-1
        DEG(ICHT)=DEG(ICHT)+HCI(JCI,INS,2)*HCI(KCI,INS,2)*RES
      end do
    end do
  end do
end if
 

allocate( del(nchl) )

write(6,600)
  600 FORMAT(//'******* CALCULATED TARGET STATES  *******'/,&
     12X,'CHANNEL',1X,'ENERGY(HARTREE)',2X,'DELTA(eV)')
do ich=1,nchl
  del(ich) = (deg(ich)-deg(1)) * toHartree
  write(6,700) ich, deg(ich), del(ich)
!MTV
!mal  EG(ICH)=DEG(ICH)
!MTV
end do
  700 FORMAT(3X,I3,3X,F15.8,2X,F8.3)

deallocate( del )

end subroutine ivo
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine alvoc(NA,I,NB,K,RES,N)

      use poly_gaus, only: eg

      implicit none

!=======================================================================
!
!     IN THIS SUBROUTINE IS CALCULATED THE FOLLOWING MATRIX ELEMENT
!
!
!                    <NA-I> (J) HN <NB-K> (L)
!
!     WHERE:
!
!     <NA-I> <J> ==>> MEANS THAT J IS NOT ANTISYMMETRIZED WITH <NA-I>
!
!     <NB-K> (L) ==>> MEANS THAT L IS ANTISYMMETRIZED WITH <NB-K>
!
!     HN = TARGET HAMILTONIAN => SUM i h(i) + SUM i e j  1/rij
!
!=======================================================================

      integer, intent(in) :: na, i, nb, k, n
      real*8, intent(out) :: res

!***
!  Local variables
!
      integer :: ia, ka, naa, nba
!***

      RES=0.0D0
      IF(NA.EQ.NB) GO TO 350
      GO TO 440
!
  350 IF(I.EQ.K) GO TO 360
!
!=======================================================================
!
!                <NA-I>  HN  <NA-K>
!
!======================================================================
!
      IA=IABS(I)
      KA=IABS(K)
      NAA=IABS(NA)
!
      IF(N.EQ.2) GO TO 10
      RES = F(IA,KA)-RPMOL(IA,KA,NAA,NAA)+2.0D0*RPMOL(IA,NAA,KA,NAA)
      RETURN
!
   10 RES = F(IA,KA)-RPMOL(IA,KA,NAA,NAA)
      RETURN
!
!=======================================================================
!
!                <NA-I>  HN  <NA-I>
!
!=======================================================================
!
  360 NAA=IABS(NA)
      IA=IABS(I)
!
      IF(N.EQ.2) GO TO 20
      RES = EG(1)-F(NAA,NAA)+F(IA,IA)-RPMOL(IA,IA,NAA,NAA) &
          +2.0D0*RPMOL(IA,NAA,IA,NAA)
      RETURN
!
   20 RES = EG(1)-F(NAA,NAA)+F(IA,IA)-RPMOL(IA,IA,NAA,NAA)
      RETURN
!
  440 IF(I.EQ.K) GO TO 370
!
!=======================================================================
!
!                <NA-I>  HN  <NB-K>
!
!=======================================================================
!
      IA=IABS(I)
      KA=IABS(K)
      NBA=IABS(NB)
      NAA=IABS(NA)
!
      IF(N.EQ.2) GO TO 30
      RES = 2.0D0*RPMOL(NAA,IA,NBA,KA)-RPMOL(NAA,NBA,IA,KA)
      RETURN
!
   30 RES =  -RPMOL(NAA,NBA,IA,KA)
      RETURN
!
!=======================================================================
!
!                <NA-I>  HN  <NB-I>
!
!=======================================================================
!
  370 NAA=IABS(NA)
      NBA=IABS(NB)
      IA=IABS(I)
!
      IF(N.EQ.2) GO TO 40
      RES = -F(NBA,NAA)-RPMOL(NBA,NAA,IA,IA)+2.0D0*RPMOL(NBA,IA,NAA,IA)
      RETURN
!
   40 RES = -F(NBA,NAA)-RPMOL(NBA,NAA,IA,IA)

      end subroutine alvoc
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine read_all

use poly_gaus, only: nbfns

implicit none

!***
!  Local variables
!
integer :: nconf, nfun
integer(i_large) :: ndd2
integer :: alloc_stat, stat
!***

open(unit=n_hamil, file='HAM.FILE', form='unformatted',iostat=stat)
call open_stat('HAM1.FILE',stat)
read(n_hamil,iostat=stat) nconf
call write_stat('nconf in HAM.FILE',stat)
if(nconf .ne. ndd) stop'WRONG NDD IN READ_ALL (N_HAMIL)'
ndd2 = int8(ndd) * int8( (ndd + 1 ) ) / int8(2)
allocate( haqqd_tri(ndd2), stat=alloc_stat )
call allocation_stat('haqqd_tri',alloc_stat)
read(n_hamil,iostat=stat) haqqd_tri
call write_stat('haqqd_tri in HAM.FILE',stat)
close(n_hamil)

open(unit=n_vvd_1, file='VVDG.FILE', form='unformatted',iostat=stat)
call open_stat('VVDG.FILE',stat)
read(n_vvd_1,iostat=stat) nfun,nconf
call write_stat('nfun, nconf in VVDG.FILE',stat)
if(nconf .ne. ndd)   stop'WRONG NDD IN READ_ALL   (N_VVD_1)'
if(nfun  .ne. nbfns) stop'WRONG NBFNS IN READ_ALL (N_VVD_1)'
allocate( vvdg(nbfns,ndd), stat=alloc_stat )
call allocation_stat('vvdg',alloc_stat)
read(n_vvd_1,iostat=stat) vvdg
call write_stat('vvdg in VVDG.FILE',stat)
!
close(n_vvd_1)

if(nchl .ge. 2) then
  open(unit=n_vvd_2, file='VVD.FILE', form='unformatted',iostat=stat)
  call open_stat('VVD1.FILE',stat)
  read(n_vvd_2,iostat=stat) nfun,nconf
  call write_stat('nfun, nconf in VVD.FILE',stat)
  if(nconf .ne. ndd)   stop'WRONG NDD IN READ_ALL   (N_VVD_2)'
  if(nfun  .ne. nbfns) stop'WRONG NBFNS IN READ_ALL (N_VVD_2)'
  allocate ( vvd(nbfns,ndd,nsa,2), stat=alloc_stat )
  call allocation_stat('vvd',alloc_stat)
  read(n_vvd_2,iostat=stat) vvd
  call write_stat('vvd in VVD.FILE',stat)
close(n_vvd_2)
  
end if

end subroutine read_all
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine wrdown_wfn!(title)

use poly_gaus, only : eg

implicit none

integer :: stat
    
!real*8, intent(in) :: title(18)

open( unit=n_wvfun, file='WFN.FILE', form='unformatted', iostat=stat)
call open_stat('WFN.FILE',stat)

!***
!  NBLOCK1 / NBLOCK2 / NSNGLT / NTRPLT
write(n_wvfun) projectile

write(n_wvfun) nwhaqq, npeaqq, nprall, nphaqq, keyq, jschw
write(n_wvfun) ng, ngmax, ns, nelc, nocc, nchl, nd, ndd, ndq
write(n_wvfun) nsa, nchls, nchlt, npct, npsc

!  XNFC2 / TITLE
!write(n_wvfun) title

!  BLOCK3 / BLOCK6
write(n_wvfun) eg
write(n_wvfun) sp
! deallocate( eg, deg, f, f2 )
! deallocate( eg )

write(n_wvfun) ngt
if(ns .gt. 0) then
  write(n_wvfun) nsnglt, ntrplt
  write(n_wvfun) nht, npt
  write(n_wvfun) ng2, ngt2

  write(n_wvfun) numsd
! deallocate( numsd )

  if(keyq .ne. 0) then
    write(n_wvfun) numsq
!   deallocate( numsq )
  end if

! deallocate( nht, npt, ng2, ngt2)
end if
!     deallocate( ngt )

!  NBLOCK7 / NBLOC12 / NBLOC13 
write(n_wvfun) nctrin
write(n_wvfun) nsp

if(nchl .ge. 2) then
  write(n_wvfun) nexc, ncspn, nctr
  write(n_wvfun) mst, nad, ncspn2
! deallocate( nexc, ncspn, nctr, nctrin )
! deallocate( mst, nad, ncspn2 )
end if

close(n_wvfun)

end subroutine wrdown_wfn
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine wrdown_den

use poly_gaus, only : wr2

implicit none

!***
!  Local variables
!
integer :: i, j, ij
integer(i_large) :: ndd2
real*8, allocatable :: eaqq(:,:)
integer :: stat, alloc_stat

open( unit=n_apden, file='DEN.FILE', form='unformatted', iostat=stat)
call open_stat('DEN.FILE',stat)

!mal inicio da simetrizacao do denominador
!     do i = 1,ndd
! do j = 1,i-1
!   ij_ed = n_diag(i,j)
!   eaqqd_tri(ij_ed) = 0.5d0 * eaqqd_tri(ij_ed)
! end do
!     end do
!mal fim da simetrizacao do denominador
ndd2 = int8(ndd) * int8( (ndd + 1 ) ) / int8(2)
allocate( eaqqd_tri(ndd2), stat=alloc_stat )
call allocation_stat('eaqqd_tri',alloc_stat)
do i = 1,ndd
  do j = 1,i
    ij_ed = n_diag(i,j)
    eaqqd_tri(ij_ed) = 0.5d0*(eaqqd(i,j)+eaqqd(j,i))
  end do
end do
deallocate( eaqqd )

eaqqd_tri = - 2d0 * pi * eaqqd_tri 
write(n_apden) eaqqd_tri
deallocate( eaqqd_tri )

allocate( eaqqd2_tri(ndd2), stat=alloc_stat )
call allocation_stat('eaqqd2_tri',alloc_stat)
ij = 0
do i=1,ndd
  do j = 1,i
    ij = ij + 1
    eaqqd2_tri(ij) = 0.5d0*(eaqqd2(i,j)+eaqqd2(j,i))
  end do
end do
deallocate( eaqqd2 )

eaqqd2_tri = - 2d0 * pi * eaqqd2_tri
write(n_apden) eaqqd2_tri
deallocate( eaqqd2_tri )

if( keyq.ne.0) then
!mal inicio da simetrizacao do denominador
  do i = 1,ndq
    do j = 1,i-1
      ij_ed = n_diag(i,j)
      eaqqq_tri(ij_ed) = 0.5d0 * eaqqq_tri(ij_ed)
    end do
  end do
!mal fim da simetrizacao do denominador
  eaqqq_tri = - 2d0 * pi * eaqqq_tri 
  write(n_apden) eaqqq_tri
  deallocate( eaqqq_tri )

  ij = 0
  do i=1,ndq
    do j = 1,i
      ij = ij + 1
      eaqqq(i,j) = eaqqq_tri(ij)
      eaqqq(j,i) = eaqqq_tri(ij)
    end do
  end do
  deallocate( eaqqq_tri )
! CALL WR2(eaqqq,ndq,ndq,ndq,ndq,n_apden)
! CALL WR2(eaqqq2,ndq,ndq,ndq,ndq,n_apden)

  eaqqq2_tri = - 2d0 * pi * eaqqq2_tri 
  write(n_apden) eaqqq2_tri
  deallocate( eaqqq2_tri )

end if

!***
close(n_apden)

end subroutine wrdown_den
!----------------------------------------------------------------

!----------------------------------------------------------------
integer(i_large) function n_diag(ii,jj)

implicit none

integer, intent(in) :: ii, jj

!***
!  Local variables
!
integer :: i, j
!***

!   Assumes (i .ge. j), i.e., do i=1,n / do j=1,i
i = max(ii,jj)
j = min(ii,jj)

!     n_diag = i*(i-1)/2 + j
n_diag = int8(i) * int8( (i-1) ) / int8(2) + int8(j)

end function 
!----------------------------------------------------------------

!-----------------------------------------------------------------------
DOUBLE PRECISION FUNCTION RPMOL(I,J,K,L)

use poly_gaus, only:  none1, ntwo, erpt, erpt_fun

implicit none

integer :: i, j, k, l

!***
!  Local variables
!
integer :: ma, mi, ij, kl
integer(i_large) :: ijkl
!***

!=======================================================================
!
!<I(1)J(1) / K(2)L(2)>.
!
!=======================================================================

!***
if(ntwoel .eq. 3) then
  rpmol = erpt_fun(i,j,k,l)
  return
end if
!***

MA=MAX0(I,J)
MI=MIN0(I,J)
IJ=NONE1(MA)+MI
MA=MAX0(K,L)
MI=MIN0(K,L)
KL=NONE1(MA)+MI
MA=MAX0(IJ,KL)
MI=MIN0(IJ,KL)
IJKL=NTWO(MA)+MI
RPMOL=ERPT(IJKL)
!mod
!     write(98,"(4(i3,1x),1pe13.6)") i,j,k,l,rpmol
!     print*,'ntwo:', mi, ma, ntwo(ma)
!mod

end function rpmol
!-----------------------------------------------------------------------

!----------------------------------------------------------------
subroutine elmbmod

implicit none

!***
!  Local variables
!
integer :: j, k, l, m, n, nc, mx, imx, ii, kk, ll, klx1, klx2, klx
real*8  :: aps, apress, aprest, apt, apress1, apress2, aps1, aps2, aprest1, aprest2, apt1, apt2
integer :: alloc_stat
!***

!=======================================================================
!     THIS SUBROUTINE GIVES THE RESULT OF PROJECTION AND
!     ANTISSYMETRIZATION OPERATIONS ON CONFIGURATIONS
!     AND ADDS THE ELECTRON INCIDENT ENERGY TERM ON THE
!     DIAGONAL ELEMENTS
!=======================================================================
!     ENEN=ENE/(NELC+1)
!=======================================================================
! INICIO DA RODADA NO SUB-ESPACO DOS DUBLETOS (ELASTICO E EXCITACOES)
!=======================================================================
allocate( eaqqd2(ndd,ndd), stat=alloc_stat )
call allocation_stat('eaqqd2',alloc_stat)
eaqqd2 = 0d0

!-----------------------------------------------------------------------
!     Caracterizacao da configuracao | 0(N) > : estado fundamental
!-----------------------------------------------------------------------
 
do M=1,NG
  N=NGT(M)
!-----------------------------------------------------------------------
!    Caracterizacao da excitacao: E - P0 = E - | 0 >< 0 |
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!    Caracterizacao das excitacoes: E - P = E - | NA-K >< NB-L |
!-----------------------------------------------------------------------
  do K=1,NSA
    do L=1,NSA
!-----------------------------------------------------------------------
!     Inicio do loop sobre os canais abertos tipo singleto
!-----------------------------------------------------------------------
      If( nchls .ne. 0 ) Then
        call APSG(L,K,N,APRESS,J)
        if((j.ge.1).and.(j.le.ng).and.(apress.ne.0)) then
          aps = apress * hci2(k,l,1)
          eaqqd2(j,M)=eaqqd2(j,M)-APS
        end if
      End If
!-----------------------------------------------------------------------
!     Inicio do loop sobre os canais abertos tipo tripleto
!-----------------------------------------------------------------------
      If( nchlt .ne. 0 ) Then
        call APTG(L,K,N,APREST,J)
        if((j.ge.1).and.(j.le.ng).and.(aprest.ne.0)) then
          apt = aprest * hci2(k,l,2)
          eaqqd2(j,M)=eaqqd2(j,M)-APT
        end if
      End If
   end do
 end do
end do

!-----------------------------------------------------------------------
!     Caracterizacao das configuracoes | NC-N(M) > : estados excitados
!-----------------------------------------------------------------------
!
do II=1,NS
  NC=NHT(II)
  N=NPT(II)
  do J=1,NG2(II)
    M=NGT2(J,II)
    do MX=1,2
      IMX=NUMSD(J,II,MX)
      if(IMX.EQ.0) cycle 
      if(N.EQ.M) then
        eaqqd2(IMX,IMX)=eaqqd2(IMX,IMX)+0.5d0
      else
        eaqqd2(IMX,IMX)=eaqqd2(IMX,IMX)+1d0
      endif
!-----------------------------------------------------------------------
!     Caracterizacao das excitacoes: E - P = E - | NA-K >< NB-L |
!-----------------------------------------------------------------------
      do KK=1,NSA
        do LL=1,NSA
!-----------------------------------------------------------------------
!     Inicio do loop sobre os canais abertos tipo singleto
!-----------------------------------------------------------------------
          If( nchls .ne. 0 ) Then
            call APSD(II,LL,KK,NC,N,M,MX,APRESS1,APRESS2,KLX1,KLX2,IMX)
            aps = hci2(kk,ll,1)
            if((KLX1.NE.0).and.(apress1.ne.0d0)) then
              APS1=APS *APRESS1
              if(e_corr(klx1)) then
                eaqqd2(klx1,IMX)=eaqqd2(klx1,IMX)-APS1*0.5d0
              else 
                eaqqd2(klx1,IMX)=eaqqd2(klx1,IMX)-APS1
              end if
            endif
            if((KLX2.NE.0).and.(apress2.ne.0d0)) then
              APS2=APS *APRESS2
              eaqqd2(klx2,IMX)=eaqqd2(klx2,IMX)-APS2
            endif
          End If
!-----------------------------------------------------------------------
!     Inicio do loop sobre os canais abertos tipo tripleto
!-----------------------------------------------------------------------
          If( nchlt .ne. 0 ) Then
            call APTD(II,LL,KK,NC,N,M,MX,APREST1,APREST2,KLX1,KLX2,IMX)
            apt = hci2(kk,ll,2)
            if((KLX1.NE.0).and.(aprest1.ne.0d0)) then
              APT1=APT *APREST1
              if(e_corr(klx1)) then
                eaqqd2(klx1,IMX)=eaqqd2(klx1,IMX)-APT1*0.5d0
              else
                eaqqd2(klx1,IMX)=eaqqd2(klx1,IMX)-APT1
              end if
            endif
            if((KLX2.NE.0).and.(aprest2.ne.0d0)) then
              APT2=APT *APREST2
              eaqqd2(klx2,IMX)=eaqqd2(klx2,IMX)-APT2
            endif
          End If
        end do
      end do
    end do 
  end do
end do
!      IF(kEYSTOP.EQ.6)CALL ELMBSYM(EAQQD,NDD,EAQQQ,NDQ,KEYQ,NCHLT)

if( keyq .eq. 0 ) return
!=======================================================================
!     INICIO DA RODADA NO SUB-ESPACO DOS QUARTETOS (EXCITACOES)
!=======================================================================
allocate( eaqqq2(ndq,ndq), stat=alloc_stat )
call allocation_stat('eaqqq2',alloc_stat)
 
do II=1,NS
  NC=NHT(II)
  N=NPT(II)
  do J=1,NG2(II)
    M=NGT2(J,II)
    IMX=NUMSQ(J,II)
    if(IMX.eq.0) cycle 
    eaqqq2(IMX,IMX)=eaqqq2(IMX,IMX)+1d0
    do KK=1,NSA
      do LL=1,NSA
        call APTQ(II,LL,KK,NC,N,M,APREST,KLX,IMX)
        if((KLX.NE.0).and.(aprest.ne.0d0)) then
          apt = aprest * hci2(kk,ll,2) 
          eaqqq2(klx,IMX)=eaqqq2(klx,IMX)-APT
        endif
        end do
    end do
  135 CONTINUE
  end do
end do
 
deallocate( ngtinv )
if( allocated(numinv) ) deallocate( numinv )
if( allocated(ngtin2) ) deallocate( ngtin2 )

end subroutine elmbmod
!----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine pvpgg(J,L,RES)
!
! =======================================================================
!
!       ESTA SUBROUTINA CALCULA O SEGUINTE ELEMENTO DE MATRIZ
!
!                       DO  DENOMINADOR
!
!
!       [I] [0] PVP [0] [J] ==>> [I] [0] ( VPN - VPE ) [0] [J]
!
!
!       [0] [I] ==>> SIGNIFICA QUE [I] NAO EH ANTISSIMETRIZADO COM [0]
!
!       [0] [J] ==>> SIGNIFICA QUE [J] NAO EH ANTISSIMETRIZADO COM [0]
!
!
! =======================================================================
!
      implicit none

      integer, intent(in) :: j, l
      real*8, intent(out) :: res

      integer :: ja, la
!=======================================================================
!=======================================================================
!     COMMON /HAMIL/EG(NCHLX),DEG(NCHLX),F(NBFNX,NBFNX),F2(NBFNX,NBFNX)
!     COMMON /MOLONE/TKET(NBFNX,NBFNX),ATRT(NBFNX,NBFNX)
!
      JA=IABS(J)
      LA=IABS(L)
      RES = F2(JA,LA)

      end subroutine pvpgg
!----------------------------------------------------------------
 
!----------------------------------------------------------------
      subroutine qvqss(NA,I,J,NB,K,L,RES)
!
!=======================================================================
!
!       ESTA SUBROUTINA CALCULA O SEGUINTE ELEMENTO DE MATRIZ
!
!
!                    [J] [NA-I] QVQ [NB-K] [L]
!
!       ONDE:
!
!       [J] [NA-I] ==>> SIGNIFICA QUE [J] NAO EH ANTISSIM. COM [NA-I]
!
!       [NB-K] [L] ==>> SIGNIFICA QUE [L] NAO EH ANTISSIM. COM [NB-K]
!
!       V = VPN - VPE
!
!       Q ==>>  CANAL FECHADO ==>> ! [NA-I] > < [NA-I] ! .....
!
!
!=======================================================================
!
      implicit none

      integer, intent(in) :: na, i, j, nb, k, l
      real*8, intent(out) :: res

      integer :: naa, nba, ia, ja, la, ka
!=======================================================================
!=======================================================================
!     COMMON /HAMIL/EG(NCHLX),DEG(NCHLX),F(NBFNX,NBFNX),F2(NBFNX,NBFNX)
!     COMMON /MOLONE/TKET(NBFNX,NBFNX),ATRT(NBFNX,NBFNX)
!
! =====================================
!     CARGA = +1
! =====================================
!
      RES=0.0D0
      IF(NA.EQ.NB) GO TO 100
      IF(I.EQ.K) GO TO 150
      RETURN
!
! =======================================================================
!
!     <  [J] [NA-I] ! V ! [NB-I] [L]  >
!
! =======================================================================
!
  150 JA = IABS(J)
      LA = IABS(L)
      NBA = IABS(NB)
      NAA = IABS(NA)
!
      RES = RPMOL(JA,LA,NBA,NAA)
      RETURN
!
  100 IA=IABS(I)
      JA=IABS(J)
      KA=IABS(K)
      LA=IABS(L)
      NAA=IABS(NA)
!
      IF(I.NE.K) GO TO 300
!
! =======================================================================
!
!     <  [J] [NA-I] ! V ! [NA-I] [L]  >
!
! =======================================================================
!
      RES = F2(JA,LA) + RPMOL(JA,LA,NAA,NAA) - &
            RPMOL(JA,LA,IA,IA)
!
      RETURN
!
! =======================================================================
!
!     <  [J] [NA-I] ! V ! [NA-K] [L]  >
!
! =======================================================================
!
  300 RES = -RPMOL(JA,LA,IA,KA)
 
      end subroutine qvqss
!----------------------------------------------------------------
 
!----------------------------------------------------------------
      subroutine qtqss(NA,I,J,NB,K,L,RES)
!
      use poly_gaus, only: tket
!=======================================================================
!
!       ESTA SUBROUTINA CALCULA O SEGUINTE ELEMENTO DE MATRIZ
!
!
!                    [J] [NA-I] QTQ [NB-K] [L]
!
!       ONDE:
!
!       [J] [NA-I] ==>> SIGNIFICA QUE [J] NAO EH ANTISSIM. COM [NA-I]
!
!       [NB-K] [L] ==>> SIGNIFICA QUE [L] NAO EH ANTISSIM. COM [NB-K]
!
!       T = ENERGIA CINETICA DA PARTICULA INCIDENTE
!
!       Q ==>>  CANAL FECHADO ==>> ! [NA-I] > < [NA-I] ! + .....
!
!=======================================================================
!
      implicit none

      integer, intent(in) :: na, i, j, nb, k, l
      real*8, intent(out) :: res

      integer :: ja, la
!=======================================================================
!=======================================================================
!=======================================================================
!     COMMON /HAMIL/EG(NCHLX),DEG(NCHLX),F(NBFNX,NBFNX),F2(NBFNX,NBFNX)
!     COMMON /MOLONE/TKET(NBFNX,NBFNX),ATRT(NBFNX,NBFNX)
!
      RES=0.0D0
      IF(NA.NE.NB.OR.I.NE.K)RETURN
!
      JA=IABS(J)
      LA=IABS(L)
!
! =======================================================================
!
!     <  [J] [NA-I] ! T ! [NA-I] [L]  >
!
! =======================================================================
!
      RES = TKET(JA,LA)

      end subroutine qtqss
!----------------------------------------------------------------
 
!----------------------------------------------------------------
      subroutine qhqss(NA,I,J,NB,K,L,RES)
!
      use poly_gaus, only: eg
! =======================================================================
!
!       ESTA SUBROUTINA CALCULA O SEGUINTE ELEMENTO DE MATRIZ
!
!
!                    [J] [NA-I] QHQ [NB-K] [L]
!
!       ONDE:
!
!       [J] [NA-I] ==>> SIGNIFICA QUE [J] NAO EH ANTISSIM. COM [NA-I]
!
!       [NB-K] [L] ==>> SIGNIFICA QUE [L] NAO EH ANTISSIM. COM [NB-K]
!
!       H = HAMILTONIANO DO ALVO => SOMA i h(i) + SOMA i e j  1/rij
!
!       Q ==>>  CANAL FECHADO ==>> ! [NA-I] > < [NA-I] ! + .....
!
! =======================================================================
!
      implicit none

      integer, intent(in) :: na, i, j, nb, k, l
      real*8, intent(out) :: res

      integer :: ia, naa, nba, ka
!=======================================================================
! =======================================================================
! =======================================================================
!     COMMON /HAMIL/EG(NCHLX),DEG(NCHLX),F(NBFNX,NBFNX),F2(NBFNX,NBFNX)
!     COMMON /MOLONE/TKET(NBFNX,NBFNX),ATRT(NBFNX,NBFNX)
!
      RES=0.0D0
      IF(J.EQ.L) GO TO 220
      RETURN
!
  220 IF(NA.EQ.NB) GO TO 350
      GO TO 440
!
  350 IF(I.EQ.K) GO TO 360
!
! =======================================================================
!
!     <  [J] [NA-I] ! HN ! [NA-K] [J]  >
!
! =======================================================================
!
      IA=IABS(I)
      KA=IABS(K)
      NAA=IABS(NA)
!
      RES=F(IA,KA)-RPMOL(IA,KA,NAA,NAA)+2*RPMOL(IA,NAA,KA,NAA)
!
      RETURN
!
! =======================================================================
!
!     <  [J] [NA-I] ! HN ! [NA-I] [J]  >
!
! =======================================================================
!
  360 NAA=IABS(NA)
      IA=IABS(I)
!
      RES=EG(1)-F(NAA,NAA)+F(IA,IA)-RPMOL(IA,IA,NAA,NAA) &
          +2*RPMOL(IA,NAA,IA,NAA)
!
      RETURN
!
  440 IF(I.EQ.K) GO TO 370
!
! =======================================================================
!
!     <  [J] [NA-I] ! HN ! [NB-K] [J]  >
!
! =======================================================================
!
      IA=IABS(I)
      KA=IABS(K)
      NBA=IABS(NB)
      NAA=IABS(NA)
!
      RES = 2*RPMOL(NAA,IA,NBA,KA)-RPMOL(NAA,NBA,IA,KA)
!
      RETURN
!
! =======================================================================
!
!     <  [J] [NA-I] ! HN ! [NB-I] [J]  >
!
! =======================================================================
!
  370 NAA=IABS(NA)
      NBA=IABS(NB)
      IA=IABS(I)
 
      RES = -F(NBA,NAA)-RPMOL(NBA,NAA,IA,IA)+2*RPMOL(NBA,IA,NAA,IA)
 
      end subroutine qhqss
!----------------------------------------------------------------
 
!----------------------------------------------------------------
      subroutine qvpsg(NA,I,J,L,RES)
!
! =======================================================================
!
!       ESTA SUBROUTINA CALCULA O SEGUINTE ELEMENTO DE MATRIZ
!
!
!                    [J] [NA-I] QVP [0] [L]
!
!       ONDE:
!
!       [J] [NA-I] ==>> SIGNIFICA QUE [J] NAO EH ANTISSIM. COM [NA-I]
!
!       [0] [L] ==>> SIGNIFICA QUE [L] NAO EH ANTISSIMETRIZADO COM [0]
!
!       V = VPN - VPE
!
!       Q ==>>  CANAL FECHADO ==>> ! [NA-I] > < [NA-I] !
!
!       P ==>>  CANAL ABERTO  ==>> ! [0] > < [0] !
!
! =======================================================================
!
      implicit none

      integer, intent(in) :: na, i, j, l
      real*8, intent(out) :: res

      integer :: ja, la, ia, naa
!=======================================================================
!=======================================================================
!=======================================================================
!     COMMON /HAMIL/EG(NCHLX),DEG(NCHLX),F(NBFNX,NBFNX),F2(NBFNX,NBFNX)
!     COMMON /MOLONE/TKET(NBFNX,NBFNX),ATRT(NBFNX,NBFNX)
!
!     RT2=DSQRT(2.0D0)
      JA=IABS(J)
      LA=IABS(L)
      IA=IABS(I)
      NAA=IABS(NA)
 
      RES = -RT2*RPMOL(JA,LA,IA,NAA)

      end subroutine qvpsg
!----------------------------------------------------------------
 
!----------------------------------------------------------------
      subroutine pvqss(NA,I,J,NB,K,L,RES)
!
!=======================================================================
!
!       ESTA SUBROUTINA CALCULA O SEGUINTE ELEMENTO DE MATRIZ
!
!
!                    [J] [NA-I] PVQ [NB-K] [L]
!
!       ONDE:
!
!       [J] [NA-I] ==>> SIGNIFICA QUE [J] NAO EH ANTISSIM. COM [NA-I]
!
!       [NB-K] [L] ==>> SIGNIFICA QUE [L] NAO EH ANTISSIM. COM [NB-K]
!
!       V = VPN - VPE
!
!       Q ==>>  CANAL FECHADO ==>> ! [NA-I] > < [NA-I] ! .....
!
!
!=======================================================================
!
      implicit none

      integer, intent(in) :: na, i, j, nb, k, l
      real*8, intent(out) :: res

      integer :: ja, la, ia, naa, nba, ka
!=======================================================================
!=======================================================================
!=======================================================================
!     COMMON /HAMIL/EG(NCHLX),DEG(NCHLX),F(NBFNX,NBFNX),F2(NBFNX,NBFNX)
!     COMMON /MOLONE/TKET(NBFNX,NBFNX),ATRT(NBFNX,NBFNX)
!
! =====================================
!     CARGA = +1
! =====================================
!
      RES=0.0D0
      IF(NA.EQ.NB) GO TO 100
      IF(I.EQ.K) GO TO 150
      RETURN
!
! =======================================================================
!
!     <  [J] [NA-I] ! V ! [NB-I] [L]  >
!
! =======================================================================
!
  150 JA = IABS(J)
      LA = IABS(L)
      NBA = IABS(NB)
      NAA = IABS(NA)
!
      RES = RPMOL(JA,LA,NBA,NAA)
      RETURN
!
  100 IA=IABS(I)
      JA=IABS(J)
      KA=IABS(K)
      LA=IABS(L)
      NAA=IABS(NA)
!
      IF(I.NE.K) GO TO 300
!
! =======================================================================
!
!     <  [J] [NA-I] ! V ! [NA-I] [L]  >
!
! =======================================================================
!
      RES = F2(JA,LA) + RPMOL(JA,LA,NAA,NAA) - &
            RPMOL(JA,LA,IA,IA)
!
      RETURN
!
! =======================================================================
!
!     <  [J] [NA-I] ! V ! [NA-K] [L]  >
!
! =======================================================================
 
  300 RES = -RPMOL(JA,LA,IA,KA)

      end subroutine pvqss
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine elmdn0_positron

implicit none

integer :: i, j, ii, jj, k, na, ka, kk
integer :: m, js, nb, lmax, l, ll, jj1, ll1
real*8  :: res
integer :: alloc_stat

!=======================================================================
!=======================================================================
!     COMMON /PH/NG,NS,NGT(NGX),NHT(NSX),NPT(NSX),NG2(NSX),NGT2(NGX,NSX)
!    1     ,NUMS(NSBFX,NSX)
!     COMMON /HAMIL/EG(NCHLX),DEG(NCHLX),F(NBFNX,NBFNX),F2(NBFNX,NBFNX)
!     COMMON /ORTGNL/NOOR,NOH(NSBFX),NOHIN(NSBFX)
!     COMMON /MOLONE/TKET(NBFNX,NBFNX),ATRT(NBFNX,NBFNX)
!
! =======================================================================
!  ELEMENTOS DE MATRIX PARA O DENOMINADOR. PARTE RELATIVA AO
!  HAMILTONIANO E AO POTENCIAL .
!  AQUI SO CONSIDERAREMOS A PARTE ELASTICA NESSA SUBROUTINA
!  OPERADOR DE PROJECAO EH [0] [0]
! =======================================================================
!     CALL READER(NOCC,NBFN,NSBF,NELC)
! ======================================================================
!       +           ^
!      A  =  PVP + QHQ - VG V
!                          P
!             ^   ^
!      ONDE : H = H0 - V
!             ^
!             H0 = E - H(N) - T(N+1)
! =======================================================================
!    NG ==>>>  # DE FUNCOES BASE P/ ORBITAIS DE ESPALHAMENTO
!    NGT ( I=1,NG ) ==>>> FUNCAO BASE P/ A FUNCAO DE ESPALHAMENTO
!

allocate( eaqq(nd,nd), stat=alloc_stat )
call allocation_stat('eaqq',alloc_stat)
eaqq = 0.0d0
 
do II=1,NG
  I=NGT(II)
  do JJ=1,II
    J=NGT(JJ)
!   SUBROUTINE PVPGG = < [I] [O] / PVP / [O] [J] >
    call PVPGG(I,J,RES)
    EAQQ(II,JJ) = RES + EAQQ(II,JJ)
    EAQQ(JJ,II) = EAQQ(II,JJ)
  end do
end do

if(ns.ne.0) then

NCHLD=NCHL-1
if(nchld.gt.0) then

  do K=1,NCHLD
    NA = NHT(K)
    KA = NPT(K)
    do JJ=1,NG2(K)
      M=NGT2(JJ,K)
      JS=NUMS(JJ,K)
      do II=1,NG
        I=NGT(II)
        call QVPSG(NA,KA,M,I,res)
        EAQQ(JS,II)=EAQQ(JS,II)+res
        EAQQ(II,JS)=EAQQ(JS,II)
      end do
    end do
  end do

  do II=1,NCHLD
    NA=NHT(II)
    I=NPT(II)
    do JJ=1,NG2(II)
      J=NGT2(JJ,II)
      do KK=1,II
        NB=NHT(KK)
        K=NPT(KK)
        LMAX=NG2(KK)
        IF(KK.EQ.II) LMAX=JJ
        do LL=1,LMAX
          L=NGT2(LL,KK)
          call QVQSS(NA,I,J,NB,K,L,RES)
          JJ1=NUMS(JJ,II)
          LL1=NUMS(LL,KK)
          EAQQ(JJ1,LL1)=EAQQ(JJ1,LL1)+RES
          EAQQ(LL1,JJ1)=EAQQ(JJ1,LL1)
        end do
      end do
    end do
  end do

end if
 
if(ns.ne.nchld) then

  do II=NCHL,NS
    NA=NHT(II)
    I=NPT(II)
    do JJ=1,NG2(II)
      J=NGT2(JJ,II)
      do KK=NCHL,II
        NB=NHT(KK)
        K=NPT(KK)
        LMAX=NG2(KK)
        IF(KK.EQ.II) LMAX=JJ
        do LL=1,LMAX
          L=NGT2(LL,KK)
          call QVQSS(NA,I,J,NB,K,L,RES)
          JJ1=NUMS(JJ,II)
          LL1=NUMS(LL,KK)
          EAQQ(JJ1,LL1)=EAQQ(JJ1,LL1)-RES
          EAQQ(LL1,JJ1)=EAQQ(JJ1,LL1)
        end do
      end do
    end do
  end do
 
  do II=NCHL,NS
    NA=NHT(II)
    I=NPT(II)
    do JJ=1,NG2(II)
      J=NGT2(JJ,II)
      do KK=NCHL,II
        NB=NHT(KK)
        K=NPT(KK)
        LMAX=NG2(KK)
        IF(KK.EQ.II) LMAX=JJ
        do LL=1,LMAX
          L=NGT2(LL,KK)
          call QTQSS(NA,I,J,NB,K,L,RES)
          JJ1=NUMS(JJ,II)
          LL1=NUMS(LL,KK)
          EAQQ(JJ1,LL1)=EAQQ(JJ1,LL1)-RES
          EAQQ(LL1,JJ1)=EAQQ(JJ1,LL1)
        end do
      end do
    end do
  end do

  do II=NCHL,NS
    NA=NHT(II)
    I=NPT(II)
    do JJ=1,NG2(II)
      J=NGT2(JJ,II)
      do KK=NCHL,II
        NB=NHT(KK)
        K=NPT(KK)
        LMAX=NG2(KK)
        IF(KK.EQ.II) LMAX=JJ
        do LL=1,LMAX
          L=NGT2(LL,KK)
          call QHQSS(NA,I,J,NB,K,L,RES)
          JJ1=NUMS(JJ,II)
          LL1=NUMS(LL,KK)
          EAQQ(JJ1,LL1)=EAQQ(JJ1,LL1)-RES
          EAQQ(LL1,JJ1)=EAQQ(JJ1,LL1)
        end do
      end do
    end do
  end do

end if

end if

!     call PRINT(eaqq,ND,ND,ND,ND)

end subroutine elmdn0_positron
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine vva00s_positron
 
use poly_gaus, only: nbfns

! =======================================================================
!     PREPARACAO PARA ELEMENTOS VGV DE FUNCOES DE CURTO ALCANCE
! =======================================================================
 
implicit none

integer :: imal, jmal, kmal, i, ii, m, j, ij, na, jj
real*8 :: r, res
integer :: nbfn
integer :: alloc_stat

!=======================================================================
!     COMMON /PH/NG,NS,NGT(NGX),NHT(NSX),NPT(NSX),NG2(NSX),NGT2(NGX,NSX)
!    1     ,NUMS(NSBFX,NSX)
!     COMMON /HAMIL/EG(NCHLX),DEG(NCHLX),F(NBFNX,NBFNX),F2(NBFNX,NBFNX)
!     COMMON /ORTGNL/NOOR,NOH(NSBFX),NOHIN(NSBFX)
!     COMMON VV(NBFNX,NDX,NCHLX)
!=======================================================================

! =======================================================================
!  VV(NBFN,ND,NCHL)
! =======================================================================
 
nbfn = nbfns

allocate ( vv(nbfn,nd,nsa), stat=alloc_stat )
call allocation_stat('vv',alloc_stat)
vv = 0.0d0
 
! =======================================================================
! < [M] [O] ! V ! [O] [I] >
! =======================================================================
 
do ii=1,ng
  i=ngt(ii)
  do m=1,nbfn
    call pvpgg(m,i,r)
    vv(m,ii,1) = r
  end do
end do

! =======================================================================
!  < [M] [O] ! V ! [NA-I] [J] >
! =======================================================================

if(ns.ne.0) then
 
do ii=1,ns
  na=nht(ii)
  i=npt(ii)
  do jj=1,ng2(ii)
    j=ngt2(jj,ii)
    ij=nums(jj,ii)
    do m=1,nbfn
      call qvpsg(na,i,j,m,res)
      vv(m,ij,1) = res
    end do
  end do
end do

end if

end subroutine vva00s_positron
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine vvamul_positron

use poly_gaus, only: nbfns

implicit none

integer :: i, ii, m, j, na, ka, kp, ip, jj, nb, jl, ll, k, l
real*8 :: r1
integer :: nbfn

!=======================================================================
!     COMMON /PH/NG,NS,NGT(NGX),NHT(NSX),NPT(NSX),NG2(NSX),NGT2(NGX,NSX)
!    1     ,NUMS(NSBFX,NSX)
!     COMMON /HAMIL/EG(NCHLX),DEG(NCHLX),F(NBFNX,NBFNX),F2(NBFNX,NBFNX)
!     COMMON /ORTGNL/NOOR,NOH(NSBFX),NOHIN(NSBFX)
!     COMMON VV(NBFNX,NDX,NCHLX)
!=======================================================================

! =======================================================================
!  VV(NBFN,ND,NCHL)
!  VERSAO MULTICANAL DE VGV00S
! =======================================================================

nbfn = nbfns

! =======================================================================
!   <NA-K><M> V <0>(I)
! =======================================================================

nchld = nchl-1

do k=1,nchld
  na = nht(k)
  ka = npt(k)
  kp = k + 1
  do m=1,nbfn
    do i=1,ng
      ii=ngt(i)
      call qvpsg(na,ka,m,ii,r1)
      vv(m,i,kp) = r1
    end do
  end do
end do

! =======================================================================
!   <NA-I><M> V <NB-J>(L)
! =======================================================================
do i=1,nchld
  ip = i + 1
  na = nht(i)
  ii = npt(i)
  do m=1,nbfn
    do j=1,ns
      jj=npt(j)
      nb=nht(j)
      do l=1,ng2(j)
        jl=nums(l,j)
        ll=ngt2(l,j)
        call qvqss(na,ii,m,nb,jj,ll,r1)
        vv(m,jl,ip) = r1
      end do
    end do
  end do
end do

end subroutine vvamul_positron
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine ident_positron

use poly_gaus, only: nbfns

!=======================================================================
!  PROGRAMA PARA DETERMINAR AS CONFIGURACOES INDEPENDENTES
!=======================================================================
implicit none

integer :: id, i, j
integer :: nbfn

!=======================================================================
!     COMMON /PH/NG,NS,NGT(NGX),NHT(NSX),NPT(NSX),NG2(NSX),NGT2(NGX,NSX)
!    1     ,NUMS(NSBFX,NSX)
!=======================================================================

nbfn = nbfns

allocate( nums(nbfn,ns) )

id=ng
nchld=nchl-1
if(nchl.ge.2) then
  do i=1,nchld
    do j=1,ng2(i)
      id=id+1
      nums(j,i)=id
      write(6,100) nums(j,i), j, i
    end do
  end do
end if

ndop=id
nd=id
if(nchld.eq.ns) return

do i=nchl,ns
  do j=1,ng2(i)
    id=id+1
    nums(j,i)=id
!   WRITE(6,200)NUMS(J,I),J,I
  end do
end do

nd=id
write(6,300) ndop, nd-ndop, nd

  100 FORMAT(/,13H  OPEN NUMS= ,I5,7H    J= ,I5,7H    I= ,I5)
! 200 FORMAT(/,15H  CLOSED NUMS= ,I5,7H    J= ,I5,7H    I= ,I5)
  300 FORMAT('//,20H Open Space: NDOP=  ,I5,//,25H Closed Space: ND-NDOP=,I5,&
             &//,31H Number of Configurations: ND =    ,I5,//')

end subroutine ident_positron
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine printn(KK,MD,ND,M,N)

implicit none

integer, intent(in) :: md, nd, m, n
integer, intent(in) :: kk(md,nd)

!***
!  Local variables
!
integer :: n20, ni, nn, nl, in, i, j
!***

N20=N/20
NI=1
if(N20.ge.1) then
  NN=N-N20*20
  NL=20
  do IN=1,N20
    WRITE(6,300) ((KK(I,J),J=NI,NL),I=1,M)
    NI=NI+20
    NL=NL+20
  end do
  IF(NN.EQ.0) RETURN
end if

write(6,400)
400 FORMAT(1H0)
do I=1,M
  WRITE(6,600) (KK(I,J),J=NI,N)
end do

600 FORMAT(20(1X,I6))
300 FORMAT(//20(1X,I6))

end subroutine printn
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine print(A,MD,ND,M,N)

implicit none

integer, intent(in) :: md, nd, m, n
real*8, intent(in) :: A(md,md)

!***
!  Local variables
!
integer :: n7, ni, nl, in, nn, i, j
!**

N7=N/7
NI=1
if(N7.ge.1) then
  NN=N-N7*7
  NL=7
  do IN=1,N7
    write(6,300) ((A(I,J),J=NI,NL),I=1,M)
    NI=NI+7
    NL=NL+7
  end do
  IF(NN.EQ.0) RETURN
end if

write(6,400)
400 FORMAT(1H0)
do I=1,M
  write(6,600) (A(I,J),J=NI,N)
end do
600 FORMAT(1H ,7(1X,1P1D10.3))
300 FORMAT(//7(1X,1P1D10.3))

end subroutine print
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module poly_conf
