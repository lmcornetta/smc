      module poly_gaus 

      use param
      use precision

      implicit none

c***
c  Modifications introduced by Fábris Kossoski in the begining of 2017.
c
c  The maximum number of primitives in a single contracted function
c  is assigned to the n_prim variable, and written in ORB.FILE,
c  as it is required in calculation of two-electron integrals 
c  with plane waves.
c***

c***
c Subroutine MINTS was modified to save a little memory in
c the transformation of 2-E integrals. The auxiliary arrays 
c RIN and RIN2 are now triangular with respect to two indexes: 
c RIN(N,N,N) became RIN(N,N2) and RIN2(N,N) became RIN(N2), 
c were N2=N*(N+1)/2
c
c Modified to use a Fock matrix calculated elsewhere (controlled
c by te flag NRFOCK)
c
c Modified to read molecular orbitals in GAMESS' and ALCHEMY's
c format (flag ORBFMT)
c
c  Modified to automatically generate dummy orbitals if
c  spherical d-shells (5 components) are used. Controlled
c  by flag NDSHEL
c
c  Modified to read nuclear attraction integrals in AO basis (ATR)
c  from VNUC.FILE. The code proceeds to calculation if this file
c  is not found
c
c  Modified to automatically generate dummy orbitals if
c  linear dependent MOs are present. Controlled by flag NLIDEP
c
c  Modified to set NCPP=NON (number of atomic centers). The older 
c  version (NCPP=30) was related to the fortran77 implementation 
c
c  Modified to write ncpp in PSE.FILE to improve the interface 
c  with partB (module legacy)
c
c  Modified to work with large integers as defined by the parameter
c  I_LARGE in module PRECISION.
c
c  Modified to keep the real NMOF value to be used in the CIS_TARGET
c  module. To avoid problems in GINTS, NMOF=NMOFF but this actually
c  causes other problems in CIS_TARGET. NMOF_TAR was created to avoid
c  that. This solution os not very elegant but GINTS is even less.
c
c  Modified to remove the COMMON block SPECS used to interfacing with
c  the AMPOLY module (indexes for 2-electron integral computation).
c  The GAMMA common block has also been removed along with the 
c  EQUIVALENCE statment between this common block and the array F(13).
c
c  Modified to use the new 2-e integrals procedures. Now there is only
c  one MINTS subroutines with several ao-to-mo transformation options
c  controlled by the NTWOEL parameter. The ERP array (AO-basis 2-e
c  integrals) is no longer declared locally. 
c***

c***
c  These are necessary elsewhere
c
      integer, save :: NGAUS, NMOF, NSBF, NBFNS, nmof_tar
      integer, save :: NON, NOC
      integer(i_large), save :: RPMOLX
c***

c***
c  Parameter NTMX: number of types of atomic orbitals
c  Parameter NCPP: pseudopotential arrays
c
      integer, parameter :: NTRNX = 10
      integer, parameter :: NTMX  = 20
      integer, save :: ITYPE(NTMX), NR(NTMX,3)
c***
!fk
!     integer, save :: MAXTYP, NLIST
      integer, save :: maxtyp, ntrn

c***
c  PP block
c     COMMON /VAI/ IPST,LGUIMA,LPSEUD(30),BLKC1(30),BLKC2(60),BLKC3(60),
c    #             BLKI1(270),BLKI2(540)
c
c  Sendo N a dimensao de LPSEUD (numero de centros), as dimensoes
c  devem seguir: BLKC1(N), BLKC2(2*N), BLKC3(2*N), BLKI1(9*N),
c                BLKI2(18*N) 
      integer, save :: NCPP
      integer, save :: IPST, LGUIMA
      integer, allocatable, save :: LPSEUD(:)
      real*8,  allocatable, save :: BLKC1(:),BLKC2(:),
     .         BLKC3(:), BLKI1(:),BLKI2(:)

c    COMMON /RESUME/NRESGAU,IRESUM(NGMX,2)
      integer, save :: NRESGAU
      integer, save, allocatable :: IRESUM(:,:)

c    COMMON /ERGNUC/  ENERGY
      real*8, save :: ENERGY

c    COMMON /OVRPG/OVPM(NBFNX,NBFNX)
      real*8, allocatable, save :: OVPM(:,:)

c     COMMON /NUMB/NONE1(NBFNX),NONE2(NSBFX),NONE3(NBFNX),NTWO(NTWOX)
      integer, allocatable, save :: none1(:), none2(:), 
     .  none3(:)
      integer(i_large), allocatable, save :: ntwo(:)

c     COMMON /ORTGNL/NOOR,NOH(NSBFX),NOHIN(NSBFX)
      integer, save :: NOOR
      integer, allocatable, save :: NOH(:),NOHIN(:)

c     COMMON /HAMIL/ENE,EG(NCHLX),DEG(NCHLX),F(NBFNX,NBFNX),
c    1F2(NBFNX,NBFNX),NFC,XNFC2,NFC3
      real*8, allocatable, save :: eg(:)

c     COMMON /MOLONE/TKET(NBFNX,NBFNX),ATRT(NBFNX,NBFNX)
      real*8, allocatable, save :: tket(:,:), atrt(:,:)

c     COMMON /MOLTWO/ERPT(RPMOLX)
      real*8, allocatable, save :: erp(:), erpt(:)

c     COMMON /SYM/MTRANS(NBFNX,NTRNX),LIST(NBFNX,NTRNX),ITEMP(NBFNX),
c    1NOS(NBFNX)
      integer, allocatable :: mtrans(:,:), list(:,:), itemp(:), nos(:)
c***

c***
c  The following arrays are defined here
c
      integer, save :: JCON(10), ILBL(18), ILAB(18), NSAVM
      integer, allocatable, save :: NTYPE(:),NFIRST(:),NLAST(:)
      integer, allocatable, save :: NUMBER(:), NCNTR(:), MLIST(:),
     . ICNTR(:), LPSKIP(:)

!fk maximum number of primitives in one contracted function
      integer, save :: n_prim

      real*8, allocatable, save :: D(:,:), C(:), ETA(:,:), VLIST(:,:),
     .  TLOPOT(:), S(:)


ccc  DATA arrays
      DATA  ITYPE/3HS  ,3HX  ,3HY  ,3HZ  ,3HXX ,3HYY ,3HZZ ,3HXY ,3HXZ ,
     X3HYZ ,3HXXX,3HYYY,3HZZZ,3HXXY,3HXXZ,3HXYY,3HYYZ,3HXZZ,3HYZZ,3HXYZ/

      DATA  NR/   0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,
     X            0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,
     X            0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1 /
c***

      contains

c----------------------------------------------------------------
c     SUBROUTINE POLY(NREC,NRINT,KEYPOLY,KEYORT,NPRALL)
      subroutine set_poly(nocc,nchl,nchld,ns,npt,nexc2)

      implicit none

      integer, intent(in) :: nocc, nchl, nchld, ns
      integer, intent(in) :: npt(max(1,ns)), nexc2(max(1,nchl))

      integer :: icon
      COMMON/IOIND/ICON(10)
      integer :: icon1, icon2, icon10
      EQUIVALENCE (ICON1,ICON(1)),(ICON2,ICON(2))
      EQUIVALENCE (ICON10,ICON(10))
c***

c***
c  Local variables
      integer :: i, j, k, m
      integer :: ii, jj, ll, il, ir
      integer :: NMOF1
!     integer :: MAXTYP, NLIST, MAXRNG
!     integer :: MAXTYP, NLIST
      integer :: NLIST
      integer :: nd_shell, nmoff, ishl
!     integer :: ISYM, NTRN, MIN
      integer :: ISYM, MIN
      integer, allocatable :: JX(:)
      real*8  :: ORBT(18)
      real*8, allocatable :: DTM(:)

      integer :: stat
c***

c***
c     CALL    RDINPT (NBFNX,NGMX,NCMX,NTMX,MAXTYP,NBFNS,NGAUS,NOC,
c    X                NLIST,C,ETA,NUMBER,NCNTR,NTYPE,NFIRST,NLAST,
c    X                MLIST,ICENTR,VLIST,ITYPE,NR,NON,LPSKIP,NSBF,NMOF)
!     call rdinpt(MAXTYP,NLIST)
      call rdinpt(NLIST)
!     call rdinpt

ccc  PREPARATION TO ORTHOGONALITY (originally in sctin0)
      allocate( noh(nbfns) )
cmal
      IF(KEYORT.NE.0)then
c     old orthogonalization order for ivo's only. Does not work for mobsci
cmal
      noh = 0
      DO 401 I=1,NOCC
 401  NOH(I)=I
      NOOR=NOCC
      IF(NCHL.EQ.1) GO TO 430
      DO 431 I=1,NCHLD
      II=NPT(NEXC2(I))
      DO 432 M=1,NOOR
 432  IF(II.EQ.NOH(M)) GO TO 431
      NOOR=NOOR+1
      NOH(NOOR)=II
 431  CONTINUE
 430  CONTINUE
cmal
      END IF
cmal

ccc  Allocate orbitals array
      allocate( d(nbfns,nbfns) )
c***

      read(5,300,iostat=stat) ISYM
  300      FORMAT(30I2)
      call write_stat('ISYM',stat)
      IF(ISYM.EQ.0) NTRN=0
      WRITE(6,301) ISYM
  301      FORMAT(//,'**** ISYM=',I4,'****')

c***
c     CALL INLAB(NBFNS,ISYM)
      call inlab(ISYM)
c***

      CALL NUMR(NBFNS,NSBF,NMOF)

c***
ccc calculate the number of dummy orbitals if necessary
      nd_shell = 0
      if(ndshel .eq. 5) then
        allocate( jx(nbfns) )
        jx = 0
        call d_shells(nd_shell,jx)
      end if
      nmoff = nmof + nd_shell + nlidep
    
      if((nd_shell.gt.0) .and. (nbfns.gt.nmoff)) then
        write(6,"('NO EXTRA ORBITALS WITH NDSHEL=5')")
        write(6,"('NMOF + ND_SHELL + NLIDEP = ',i6)") nmoff
        write(6,"('NBFN = ',i6)") nbfns
        stop'ABORT...'
      end if

      if((nbfns.ne.nmoff) .or. (nsbf.ne.nmoff)) then
        if(nrfock .ne. 0) then
          write(6,"('NO EXTRA ORBITALS WITH NRFOCK<>0')")
          stop'ABORT...'
        end if
      end if

      if(orbfmt .eq. "ALCHEM") then
        do I=1,nmof
          read(5,113,iostat=stat) (ORBT(K),K=1,18)
          call write_stat('first line of orbitals',stat)
          read(5,112,iostat=stat) (D(J,I),J=1,nmoff)
          call write_stat('orbitals coefficients',stat)
        end do
      else 
        LL = nmoff / 5
        IR = MOD(nmoff,5)
        DO J=1, nmof
           JJ = MOD(J,100)
           K = 0
           DO IL=1, LL
              II = MOD(IL,100)
              read(5,10,iostat=stat) JJ, II, (d(K+I,J), I=1,5)
              call write_stat('orbitals coefficients',stat)
              K = K + 5
           ENDDO
           II = MOD(LL+1,100)
           IF(IR.GT.0) read(5,10,iostat=stat) JJ, II, (d(K+I,J), I=1,IR)
           call write_stat('orbitals coefficients',stat)
        ENDDO
      end if

      if( (ndshel.eq.5) .or. (nlidep.gt.0) )
     .  write(6,"(//,'***  A TOTAL OF ',i6,' DUMMY ORBITALS',
     .    ' WILL BE GENERATED ***',/,'NUMBER OF D SHELLS = ',i6,
     .    /,'NUMBER OF LD MOs   = ',i6)") (nd_shell+nlidep),
     .    nd_shell, nlidep

ccc  Generate dummy orbitals related to d shells
      if(ndshel .eq. 5) then
        ishl = 0
        do i = (nmof+1),(nmoff-nlidep)
          ishl = ishl + 1
          do j = 1,nmoff
            d(j,i) = 0d0
          end do
          d(jx(ishl),i) = 1d0
        end do
        deallocate( jx )
      end if

ccc  Generate dummy orbitals related to LD orbitals
      if(nlidep .ne. 0) then
        ishl = 0
        nmof1 = nmoff - nlidep + 1
        do i = nmof1, nmoff
          ishl = ishl + 1
          do j = 1,nmoff
            d(j,i) = 0d0
          end do
            d(ishl,i) = 1d0
        end do
      end if

!fk
      call zero_small_coef(d,nbfns)

ccc Avoid problems in module CIS_TARGET
      nmof_tar = nmof

ccc Avoid problems in GINTS
      nmof = nmoff

      NMOF1=NMOF+1
      IF(NMOF.EQ.NSBF) GO TO 311
c***

c***
      allocate( jx(nbfns), dtm(nbfns) )
c***

      DO 211 I=NMOF1,NSBF
      DO 426 K=1,NSBF
  426 D(K,I)=0.0D0
CMTV  READ(5,511) MIN,(JX(K),DTM(K),K=1,MIN)
      read(5,*,iostat=stat) MIN,(JX(K),DTM(K),K=1,MIN)
      call write_stat('orbitals coefficients',stat)
      DO 611 K=1,MIN
        IF(DTM(K).EQ.0.0D0) DTM(K)=1.0D0
  611 D(JX(K),I)=DTM(K)
  211 CONTINUE
! 511 FORMAT(I2,2X,4(I2,1X,D15.8))

c***
      deallocate( jx, dtm )
c***

  311 CONTINUE

c***
      allocate( eg(nchl) )
      read(5,*,iostat=stat) (EG(I),I=1,NCHL)
      call write_stat('electronic energies',stat)
      WRITE(6,222) (EG(I),I=1,NCHL)
  222 FORMAT(//'*** STATE ENERGY ****'/6(3X,D15.8))

      call end_flag(flag_o2)
c***

  113   FORMAT(18A4)
C=======================================================================
      IF(NPRALL.EQ.0)GO TO 777
      WRITE(6,1234)
 1234 FORMAT(//,'**** TRASFORMATION MATRIX D ****')
C      CALL PRINT(D,NSBFX,NSBFX,NMOF,NMOF)
  777 CONTINUE
C=======================================================================
  112      FORMAT(5D15.8)

c***
cmod  if(nrfock .ne. 0) then
cmod    call ovlgrn_new
cmod    deallocate( icntr, number, ncntr )
cmod    deallocate( mlist, c, none2, none3 ) 
cmod    return
cmod  end if
c***
C=======================================================================
C
C.....CALCULATE  THE  OVERLAP  INTEGRALS  AND  CHECK  SYMMETRY
C
C=======================================================================
c***
c     CALL   GINTS (NTYPE,NR,NFIRST,NLAST,ETA,NBFNX,NTMX,
c    X    NGMX,NBFNS,NSBF,NMOF,NTRN,KEYORT,NPRALL,KEYPOLY)
!     call gints(NTRN)
      call gints
!fk
      open(unit=n_final_orb,file='ORBS.LIS',form='formatted',
     .     iostat=stat )
      call open_stat('ORBS.LIS',stat)
      call write_gms_orbitals(d,nbfns,nmof_tar,n_final_orb)
      close(n_final_orb)

 10   format(I2,1X,I2,5E15.8)
912   FORMAT(' GINTS - EVALUATE OVERLAP INTEGRALS',//)

      end subroutine set_poly
c----------------------------------------------------------------

c----------------------------------------------------------------
      subroutine zero_small_coef(d,nbfns)

      implicit none

      real*8, dimension(nbfns,nbfns), intent(inout) :: d
      integer, intent(in) :: nbfns

      integer, dimension(nbfns) :: table_ao
      integer, dimension(nbfns) :: table_mo
      integer :: nunit, i, j, k
      character(len=*), parameter :: file_table_ao = 'table_ao'
      logical :: ex
      integer :: stat

      inquire(file=file_table_ao,exist=ex)
      if(ex) then
         write(*,"(2a)") file_table_ao, 'was found'
         open(newunit=nunit,file=file_table_ao,form="formatted",
     .        iostat=stat)
         call open_stat(file_table_ao,stat)
         read(nunit,*,iostat=stat) table_ao 
         call write_stat('table_ao',stat)
         close(nunit)
         table_ao = table_ao / abs(table_ao)
         table_mo = table_ao( (maxloc(d,dim=1)) )
         do i=1,nbfns
           k = table_mo(i)
           write(*,"(2(i0,x))") i, k
           do j=1,nbfns
             if(table_ao(j).ne.k) then
               write(*,"(2(i0,x),1e20.12)") j, i, d(j,i)
               d(j,i) = 0d0
             end if
           end do
           write(*,*)
         end do
      end if

      end subroutine zero_small_coef
c----------------------------------------------------------------

c----------------------------------------------------------------
c     SUBROUTINE POLY(NREC,NRINT,KEYPOLY,KEYORT,NPRALL)
      subroutine poly(ns)

!     use ampoly, only: STR0, STR1, STR2, STR3, STR4, STR5, STR6,
!    .                  STR7, STR8,STR9, STR10,STR11,STR12, cpycmi,
!    .                  generf, savrge
      use ampoly, only: cpycmi, generf, savrge

      implicit none

      integer, intent(in) :: ns

c***
c These common blocks are kept for use in ampoly
      integer :: NITAPE, LSTNAM, NOTAPE, INTNAM,NCTAPE
      COMMON /NAMTAP/  NITAPE, LSTNAM, NOTAPE, INTNAM,NCTAPE

      integer :: ICNT,JCNT,KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,
     1LS,IF,JF,KF,LF,M,I,J,K,L
cccc  COMMON /SPECS/ ICNT,JCNT,KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,
cccc 1LS,IF,JF,KF,LF,M,I,J,K,L

cccc  real*8  :: STR0, STR1, STR2, STR3, STR4, STR5, STR6, STR7, STR8,
cccc .           STR9, STR10,STR11,STR12
cccc  COMMON/STORE/STR0(280),STR1(280),STR2(280),STR3(280),STR4(280),
cccc 1STR5(280),STR6(280),STR7(280),STR8(280),STR9(280),STR10(280),
cccc 2STR11(280),STR12(280)

      integer :: icon
      COMMON/IOIND/ICON(10)
      integer :: icon1, icon2, icon10
      EQUIVALENCE (ICON1,ICON(1)),(ICON2,ICON(2))
      EQUIVALENCE (ICON10,ICON(10))
c***

c***
c  Local variables
c
      logical :: useless = .false.
      integer :: NREC
!     integer :: MAXTYP, NLIST, MAXRNG
!     integer :: MAXTYP, MAXRNG
      integer :: MAXRNG
!     integer :: NEWLBL
!     integer :: ISYM, NTRN, MIN, NINMAX
!     integer :: NTRN, NINMAX
      integer :: NINMAX

ccc     IF(NRINT.NE.0) RETURN
c***
C=======================================================================
C
C.....CALCULATE THE POTENTIAL ENERGY INTEGRALS
C
C=======================================================================
c***
c     CALL   VINTS (NON,VLIST,LPSKIP,NTYPE,NR,NFIRST,NLAST,ETA,
c    X  NBFNX,NCMX,NTMX,NINMAX,NGMX,NBFNS,NSBF,NMOF,NTRN,NREC,NPRALL)
      nrec = nwint
!     call vints (NINMAX,NTRN,NREC)
!     call vints(NTRN,NREC)
      call vints(NREC)
c***

c***
c     CALL   TINTS (NTYPE,NR,NFIRST,NLAST,ETA,NBFNX,NTMX,
c    X              NINMAX,NGMX,NBFNS,NSBF,NMOF,NTRN,NREC,NPRALL)
!     call tints(NINMAX,NTRN,NREC)
!     call tints(NTRN,NREC)
      call tints(NREC)

ccc Avoid useless calculation of 2-e integrals
      useless = ( (ns.eq.0) .and. (nrfock.ne.0) .and. (nroute.ge.0) )
      if( useless ) then
        write(6,"(/,/,20x,'**** WARNING  ****',/,5x,
     .  '(NS.EQ.0) .AND. (NRFOCK.NE.0)',/,5x,
     .  'NO NEED TO CALCULATE 2-E INTEGRALS',/,5x,
     .  'NTWOLEL SWITCHED TO 3',/,/)")
        ntwoel = 3 
        return
      end if
c***

C=======================================================================
C
C     IF  ( ICON(1) .EQ. 1 )   GO TO 876
C.....CALCULATE THE 2-ELECTRON INTEGRALS
C=======================================================================
      IF(ICON2.GE.3) GO TO 570
      I=4
570   IF(ICON1.LE.2) GO TO 580
C=======================================================================
C
C.....COPY 2-ELECTRON INTEGRALS
C
C=======================================================================
cccc  CALL CPYCMI(NBFNS,NINMAX)
      call cpycmi(NBFNS,NINMAX,ICNT,JCNT,KCNT,LCNT,ITYP,JTYP,
     .           KTYP,LTYP,IS,JS,KS,LS,IF,JF,KF,LF,M,I,J,K,L)
C=======================================================================
C
C     IF(ICON1.EQ.3) GO TO 876
C
C=======================================================================
580   CONTINUE
C=======================================================================
C
C...   GENERATE F INTEGRAL TABLES
C
C=======================================================================
      CALL GENERF(MAXTYP,MAXRNG)
CW
      WRITE (6,420) MAXTYP,MAXRNG
C=======================================================================
C
C...COMPUTE AND STORE PRE-EXPONENTIAL FACTOR FOR ALL I,J INDEX PAIRS.
C
C=======================================================================
CW
C=======================================================================
C
C     WRITE  (6,410)
C
C=======================================================================
c***
c     CALL    SAVRGE (NGAUS,ETA,S,NGMX,NSAVMX)
      nsavm = ngaus*(ngaus+1)/2
      allocate( s(nsavm) )
      CALL    SAVRGE (NGAUS,ETA,S,ngaus,nsavm)
c***
CW
C=======================================================================
C
C     WRITE (6,860)
C
C=======================================================================
c***
c     CALL    MINTS (NLIST,MLIST,NCNTR,NTYPE,ETA,NFIRST,NLAST,S,
c    X                     NR,NBFNX,NTMX,NINMAX,NGMX,NSAVMX,NEWLBL,
c    X                     NBFNS,NSBF,NMOF,NTRN,NREC)
!     call mints(NLIST,S,NINMAX,NSAVM,NEWLBL,NTRN,NREC)
      call mints(S,NSAVM)
c***

      if( allocated(s) ) deallocate( s )

ccc   RETURN
CW
ccc 959  NNNA=0
ccc      WRITE(6,960) MLAB(I),NLAB(I),NNNA
ccc      RETURN

7787  FORMAT('1TAPE USED FOR LABELS - ',18A4//)
913   FORMAT(' TINTS - EVALUATE KINETIC ENERGY INTEGRALS',//)
933   FORMAT(' VINTS - EVALUATE POTENTIAL ENERGY INTEGRALS',//)
400   FORMAT(' CPYCMI - COPY 2-ELECTRON INTEGRALS'//)
420   FORMAT(//,' GENERF - GENERATE F-INTEGRAL TABLES',5X,'MAXTYP =',
     2 I3,5X,'MAXRNG =',I5,//)
410   FORMAT(' SAVRGE - COMPUTE PRE-EXPONENTIAL FACTORS',//)
860   FORMAT(' MINTS - EVALUATE 2-ELECTRON INTEGRALS',//)
950   FORMAT(//,10X,'**  LABELS  NBFNS  =',I4,',  DOES  NOT  AGREE WITH
     2 INTEGRALS  NBFNS  =',I4,'  **',/)
960   FORMAT(//,10X,'** EXPECTING  ',A4,'  OR  ',A4,',  FOUND  ',A4,
     2 '  **',/)
 11   format(I2,1X,I2,5(1pe15.8))

      end subroutine poly
c----------------------------------------------------------------

c----------------------------------------------------------------
c     SUBROUTINE   RDINPT (NBFXX,NGXX,NCXX,NTXX,MAXTYP,NBFNS,NGAUS,NOC,
c    $                    NLIST,C,ETA,NUMBER,NCNTR,NTYPE,NFIRST,NLAST,
c    $                    MLIST,ICNTR,VLIST,ITYPE,NR,NON,LPSKIP,
c    $                    NSBF,NMOF)
!     subroutine rdinpt(MAXTYP,NLIST)
      subroutine rdinpt(NLIST)
!     subroutine rdinpt

      use ampoly, only: ovlap

      implicit none

C=======================================================================
C
C.....ICON DEFINITIONS
C
C     ICON(1) - CALCULATE
C         0 = 1E AND 2E
C         1 = 1E ONLY
C         2 = 1E AND CHANGE SOME 2E (MLIST)
C         3 = 1E AND COPY 2E
C         4 = 1E, COPY 2E AND ADD BFNS
C         5 = 1E, COPY 2E AND RESTART 2E
C
C     ICON(2) - TAPE IN
C         0 = NONE
C         1 = POLYATOM
C         2 = POLYIJLK
C         3 = POLYPAIR
C         4 = POLYPAIR + IJLK
C
C     ICON(3) - NORMALIZE
C         0 = YES
C         1 = NO
C
C     ICON(4) - CHECK SYMMETRY
C         0 = YES
C         1 = NO
C
C=======================================================================
c     IMPLICIT REAL*8(A-H,O-Z)
C=======================================================================
C
C     VALID COMBINATIONS OF ICON(1) AND ICON(2)
C
C                    ICON(2)
C             0    1    2    3    4
C     ICON(1)
C        0    X    X    X
C        1    X
C        2         X    X
C        3         X    X    X    X
C        4         X    X
C        5         X
C
C=======================================================================
C=======================================================================
c     INCLUDE 'STRUC.BFN'
C=======================================================================

!     integer, intent(out) :: MAXTYP,NLIST
      integer, intent(out) :: NLIST

c*** 
c  This common block is kept for use in ampoly
c
      integer :: icon, icon1, icon2, icon3, icon4, icon5,
     .                 icon6, icon7, icon8, icon9, icon10
      COMMON/IOIND/ICON(10)
      EQUIVALENCE (ICON1,ICON(1)),(ICON2,ICON(2)),(ICON3,ICON(3)),
     2 (ICON4,ICON(4)),(ICON9,ICON(9)),(ICON10,ICON(10))
c***

c***
c  Local variables
c
      integer :: i, j, k, l, m, n, is, if, ii, jj, kk, jo, ja, jb, jt, 
     .           mj, mm, ij
      integer :: isf, isave, ityp, KCNT, KTYP, INC, NONM1, IP1, NAC
      integer :: kont1, kont2, kont3, ntxx
      integer :: iblnk, ijlk, jkpr, ierr, ncon(30)
      real*8  :: t, t1, t2, t3, soo, gii, prtint, rij, rija
      real*8  :: blank,expnt, zval, zro(3), wzero(3)
      real*8  :: CORE(2),ALFCOR(2),SIGMA(3,3),AION(2,3,3)

      DATA NCON/1,1,1,0,0,1,0,0,0,0,0,1,1,0,0,0,1,1,1,1,0,1,1,0,0,0,1,0,
     * 0,0/
      DATA IBLNK/4H    /
      DATA BLANK/8H        /
      DATA IJLK,JKPR/4HIJLK,4HJKPR/
      DATA ZRO/8H ZERO CO,8HEF SET T,8HO ONE   /

      integer :: stat
c***

      NLIST=0
      ICANON =  2
      IERR=0
C=======================================================================
C
C  READ AND PRINT THE PROBLEM LABEL.
C
C=======================================================================
CR
c***
      call ini_flag(flag_o1)
c***
      read(5,930,iostat=stat) (ILBL(II),II=1,18)
      call write_stat('title',stat)
C=======================================================================
C
C  READ AND PRINT THE CONTROL OPTIONS.
C
C=======================================================================
CR
c***
      read  (5,900,iostat=stat)  ICON
      call write_stat('ICON',stat)
c     icon1  = icon(1)
c     icon2  = icon(2)
c     icon3  = icon(3)
c     icon4  = icon(4)
c     icon5  = icon(5)
c     icon6  = icon(6)
c     icon7  = icon(7)
c     icon8  = icon(8)
c     icon9  = icon(9)
c     icon10 = icon(10)
      ntxx = ntmx
c***
      IF(ICON1.EQ.2) ICON9=1
      IF(ICON1.LT.0.OR.ICON1.GT.5) GO TO 710
      IF(ICON2.LT.0.OR.ICON2.GT.4) GO TO 710
      IF(NCON(5*ICON1+ICON2+1).EQ.0) GO TO 710
      ICON10=ICON2
      IF(ICON1.GE.4) ICON10=0
      ILBL(18)=IBLNK
      IF(ICON1.NE.4.AND.ICON2.EQ.2) ILBL(18)=IJLK
      IF(ICON2.GE.3) ILBL(18)=JKPR
CW
      WRITE(6,931) (ILBL(II),II=1,18)
CW
      WRITE (6,901)  ICON
      GO TO 670
710   IERR=IERR+1
CW
      WRITE(6,931) (ILBL(II),II=1,18)
CW
      WRITE (6,901)  ICON
CW
      WRITE(6,935)
C=======================================================================
C
C  READ AND PRINT THE CENTER COORDINATES.
C
C=======================================================================
CW
670   WRITE(6,671)
CW
ccc   WRITE(6,672) NBFXX,NGXX
CW
ccc   WRITE(6,673) NCXX
CW
      WRITE(6,674)
  671 FORMAT(//5X,20HPROGRAM LIMITATIONS     )
  672 FORMAT(/5X,25HMAX NO BASIS FUNCTIONS=   ,I4
     X       /5X,28HMAX NO GAUSSIAN PRIMATIVES=   ,I4)
  673 FORMAT(/5X,15HMAX NO CENTERS=  ,I4)
  674 FORMAT(/5X,  'S,P,D,F TYPE GAUSSIANS ONLY')
CR
cmb
      read(5,*,iostat=stat) IPST, LGUIMA
      call write_stat('IPST, LGUIMA',stat)
      WRITE(6,990)IPST
      IF(IPST.EQ.1)WRITE(6,991)
cmb

      read (5,918,iostat=stat) NON, NAC
      call write_stat('NON, NAC',stat)
CW
      WRITE (6,904)  NON, NAC
      NOC=NON+NAC

c***
c  Allocate PP arrays
c
      if( non .eq. 0 ) then
        write(6,"('NON = 0 (no atomic centers)')") 
        stop'ABORT...'
      end if
      ncpp = non

      allocate( LPSEUD(NCPP) )
      if(ipst.eq.1) then
        allocate( BLKC1(NCPP), BLKC2(2*NCPP) )
        allocate( BLKC3(2*NCPP), BLKI1(9*NCPP) )
        allocate( BLKI2(18*NCPP) )
        do i = 1,ncpp
          blkc1(i) = 0d0
          blkc2(i) = 0d0
          blkc3(i) = 0d0
          blki1(i) = 0d0
          blki2(i) = 0d0
        end do
        do i = ncpp+1,ncpp*2
          blkc2(i) = 0d0
          blkc3(i) = 0d0
          blki1(i) = 0d0
          blki2(i) = 0d0
        end do
        do i = ncpp*2+1,ncpp*9
          blki1(i) = 0d0
          blki2(i) = 0d0
        end do
        do i = ncpp*9+1,ncpp*18
          blki2(i) = 0d0
        end do
      end if
c***

ccc
c     IF(NON.LE.NCXX.AND.NOC.LE.NCXX) GO TO 60
c       IERR =  IERR + 1
CW
c       WRITE(6,940)
c       NON=MIN0(NON,NCXX)
c       NAC=MIN0(NAC,NCXX-NON)
CW
ccc

c***
c  Allocate variables related to the number of centers
c
      allocate( icntr(noc), vlist(noc,4) )
      allocate( lpskip(noc) )
c***
   60 WRITE(6,905)
      KONT1=0
      KONT2=0
      KONT3=0
      DO  70  I=1,NON
CR
cmb        READ  (5,906)  ICNTR(I), (VLIST(I,J),J=1,4),TLOPOT(I)
        read(5,906,iostat=stat) ICNTR(I), (VLIST(I,J),J=1,4),LPSEUD(I)
        call write_stat('line with atom, coordinates and charge',stat)
CW
cmb        WRITE (6,907)  ICNTR(I), (VLIST(I,J),J=1,4),TLOPOT(I)
        WRITE (6,907)  ICNTR(I), (VLIST(I,J),J=1,4),LPSEUD(I)
cmb
Cmb####################PSEPAR######################
      IF(LPSEUD(I).EQ.1)THEN
c***
      if(ipst.ne.1) then
        write(6,"('PP INCONSISTENCY FOR ATOM ',i3)"), I
        stop'ABORT...'
      end if
c***
      WRITE(6,995)
       read(5,*,iostat=stat)(CORE(M),ALFCOR(M),M=1,2),ZVAL
        call write_stat('pseudopotential parameters',stat)
        DO L=1,3
          read(5,*,iostat=stat)(SIGMA(J,L),J=1,3)
          call write_stat('pseudopotential parameters',stat)
          read(5,*,iostat=stat)((AION(N,J,L),J=1,3),N=1,2)
          call write_stat('pseudopotential parameters',stat)
        END DO
          BLKC1(I)=ZVAL
        DO M=1,2
          KONT1=KONT1+1
          BLKC2(KONT1)=CORE(M)
          BLKC3(KONT1)=ALFCOR(M)
        END DO
        DO L=1,3
        DO J=1,3
          KONT2=KONT2+1
          BLKI1(KONT2)=SIGMA(J,L)
        END DO
         DO N=1,2
          DO J=1,3
           KONT3=KONT3+1
           BLKI2(KONT3)=AION(N,J,L)
          END DO
         END DO
        END DO
      END IF
Cmb####################PSEPAR######################
cmb
   70 CONTINUE
      IF  ( NAC .LE. 0 )   GO TO  90
      K = NON + 1
      L=NON+NAC
CW
      WRITE(6,908)
      DO  80  I=K,L
CW
CR
        read(5,906,iostat=stat)  ICNTR(I),(VLIST(I,J),J=1,4)
        call write_stat('line with atom, coordinates and charge',stat)
CW
CW
        WRITE (6,907)  ICNTR(I), (VLIST(I,J),J=1,4)
   80   CONTINUE
C=======================================================================
C
C  READ AND CHECK THE BASIS FUNCTIONS.
C
C=======================================================================
CR
ccc  90  READ (5,918)NGAUS,NBFNS,NSBF,NMOF
  90  read(5,*,iostat=stat) NGAUS,NBFNS,NSBF,NMOF
      call write_stat('ngaus, nbfns, nsbf, nmof',stat)
      NRESGAU=NGAUS
CW
      WRITE (6,909)  NGAUS, NBFNS,NSBF,NMOF
      IF(NSBF.EQ.0) NSBF=NBFNS
      IF(NMOF.EQ.0) NMOF=NBFNS
ccc   IF  (NGAUS .LE. NGXX )   GO TO  100
ccc     IERR = IERR + 1
CW
ccc     WRITE (6,941)
ccc     NGAUS =  NGXX
c100  IF  ( NBFNS .LE. NBFXX )   GO TO  110
ccc     IERR =  IERR + 1
CW
ccc     WRITE (6,942)
ccc     NBFNS =  NBFXX
CR

c***
c  Allocate variables related to the number of basis functions
c  and primitive gaussians
c
      allocate( number(nbfns), nfirst(nbfns), nlast(nbfns)  )
      allocate( ncntr(nbfns), ntype(nbfns), mlist(nbfns) )
      allocate( c(ngaus), eta(ngaus,5), iresum(ngaus,2) )
      mlist = 0
c***

!110  READ  (5,915,iostat=stat)  (NUMBER(I),I=1,NBFNS)
 110  read(5,*,iostat=stat)  (NUMBER(I),I=1,NBFNS)
      call write_stat('atomic orbitals indexes',stat)
!fk
      n_prim = maxval(number)

      IF  ( NBFNS .NE. NGAUS )   GO TO  114
      DO  112  I=1,NBFNS
 112    NUMBER(I) =  1
 114  NFIRST(1) = 1
      IF(NUMBER(1).EQ.0) NUMBER(1)=1
      NLAST(1) =  NUMBER(1)
      DO  120  I=2,NBFNS
        NFIRST(I) =  NLAST(I-1) + 1
      IF(NUMBER(I).EQ.0) NUMBER(I)=1
  120   NLAST(I)  =  NLAST(I-1) + NUMBER(I)
      IF  ( NLAST(NBFNS) .EQ. NGAUS )  GO TO  130
        IERR =  IERR + 1
CW
        WRITE (6,943)
  130 MAXTYP =  1
      I =  0
      DO  300  JO=1,NBFNS
        ISAVE =  0
        ISF =  NUMBER(JO)
        DO  290  K=1,ISF
          I =  I + 1
          IF  ( ISAVE .NE. 0 )  GO TO  160
CR
         read(5,910,iostat=stat)  KCNT, KTYP, INC, EXPNT, C(I)
         call write_stat('list of basis set',stat)
         IF  ( INC .EQ. 0 )   GO TO  170
         IF  ( INC .GT. 0 .AND. INC .LT. JO )  GO TO   140
            IERR =  IERR + 1
CW
            WRITE (6,944)  JO, KCNT, KTYP, INC
            GO TO  300
  140     IF  ( NUMBER(INC) .EQ. ISF ) GO TO 150
            IERR =  IERR + 1
CW
            WRITE (6,945)  JO, INC
            ISF =  NUMBER(INC)
 150      II =  NFIRST(INC)
          ISAVE =  1
 160      C(I) =  C(II)
          EXPNT =  ETA(II,4)
          II =  II + 1
 170      IF  ( EXPNT .NE. 0.  D0)  GO TO  180
            IERR = IERR + 1
CW
            WRITE (6,946)
  180 DO 6055 MM=1,3
 6055 WZERO(MM)=BLANK
      IF(C(I).NE.0.D0) GO TO 190
            C(I) =  1.D0
      DO 6056 MM=1,3
 6056 WZERO(MM)=ZRO(MM)
 190      IF  ( K .EQ. 1  .OR.  IERR .NE. 0 )   GO TO  230
          IF  ( KCNT .EQ. ICNTR(JA) )   GO TO  210
            IERR = IERR + 1
CW
            WRITE (6,948)  JO, K
 210      IF  ( KTYP .EQ. ITYPE(JB) )   GO TO  270
            IERR = IERR + 1
CW
            WRITE (6,949)  JO, K
         GO TO  270
 230      DO  240  JT=1,NOC
            JA =  JT
            IF  ( KCNT .EQ. ICNTR(JA) )   GO TO  250
 240        CONTINUE
            IERR =  IERR + 1
CW
            WRITE (6,950)  JO
 250      DO  260  JT=1,NTXX
            JB =  JT
            IF  ( KTYP .EQ. ITYPE(JB) )   GO TO  270
 260        CONTINUE
            IERR =  IERR + 1
            WRITE (6,951)  JO
CW
 270      NCNTR(JO) =  JA
          NTYPE(JO) =  JB
CW
         IRESUM(I,1)=KCNT
         IRESUM(I,2)=KTYP
      WRITE(6,1911) I,JO,K,KCNT,KTYP,EXPNT,C(I),(WZERO(MJ),MJ=1,3)
          DO  280  M=1,3
 280        ETA(I,M) =  VLIST(JA,M)
          ETA(I,4) =  EXPNT
          IF  ( JB .GT. MAXTYP )   MAXTYP =  JB
 290      CONTINUE
 300    CONTINUE
CR
      read(5,912,iostat=stat)  ACRCY, SCALE
      call write_stat('ACRCY, SCALE',stat)
      IF  ( ACRCY .EQ. 0. D0)  ACRCY =  1.0D-10
      IF  ( SCALE .EQ. 0. D0)  SCALE =  1.0D0
      SCALE =  SCALE*ACRCY
CW
      WRITE (6,913)  ACRCY, SCALE
      IF  ( ICON(9) .NE. 1 )   GO TO  320
CR
      read(5,918,iostat=stat)  NLIST
      call write_stat('NLIST',stat)
CR
      read(5,918,iostat=stat)  (MLIST(I),I=1,NLIST)
      call write_stat('MLIST',stat)
CW
      WRITE (6,914)  NLIST, (MLIST(I),I=1,NLIST)
      DO  310  I=1,NLIST
C=======================================================================
C
C     ERROR IN POLYATOMIN NEXT STATEMENT-FIXED 7/31/69-WYH
C     (NBFN) REPLACED BY (NBFNS)
C
C=======================================================================
        IF  ( MLIST(I) .GT. 0 .AND. MLIST(I) .LE. NBFNS)   GO TO  310
        IERR =  IERR + 1
CW
        WRITE (6,954)
 310    CONTINUE
320   IF  ( IERR .EQ. 0 )   GO TO  330
C=======================================================================
C
C     CALL LOCPOT(NOC,VLIST,LPSKIP)
C
C=======================================================================
      DO 434 KK=1,NOC
  434      LPSKIP(KK)=1
CW
        WRITE (6,952)  IERR
      STOP
C=======================================================================
C
C  ARE THE BASIS FUNCTIONS IN STANDARD ORDER.
C
C=======================================================================
 330  DO  340  JO=2,NBFNS
        IF  ( NTYPE(JO) .GE. NTYPE(JO-1) )  GO TO  340
        ICANON =  1
CW
        WRITE (6,922)
        GO TO  400
 340    CONTINUE
CW
      WRITE (6,923)
C=======================================================================
C
C400  CALL LOCPOT(NOC,VLIST,LPSKIP)
C
C=======================================================================
  400      DO 433 KK=1,NOC
  433      LPSKIP(KK)=1
C=======================================================================
C
C  NORMALIZE THE PRIMITIVE FUNCTIONS
C
C=======================================================================
      DO  420  I=1,NBFNS
        ITYP =  NTYPE(I)
        L =  NR(ITYP,1)
        M =  NR(ITYP,2)
        N =  NR(ITYP,3)
        IS =  NFIRST(I)
        IF =  NLAST(I)
        DO  410  II=IS,IF
          T =  0.5D0/ETA(II,4)
          SOO =  PITERN*T**1.5D0
          T1 =  OVLAP(L,L,0.0D0,0.0D0,T)
          T2 =  OVLAP(M,M,0.0D0,0.0D0,T)
          T3 =  OVLAP(N,N,0.0D0,0.0D0,T)
          GII =  SOO*T1*T2*T3
 410      ETA(II,5) =  1.0D0/DSQRT(GII)
 420    CONTINUE
C=======================================================================
C
C     IF  ( ICON(3) .EQ. 1 )   GO TO  550
C  RENORMALIZE THE BASIS FUNCTIONS.
C
C=======================================================================
CW
      WRITE (6,916)
      DO  540  I=1,NBFNS
        ITYP =  NTYPE(I)
        L =  NR(ITYP,1)
        M =  NR(ITYP,2)
        N =  NR(ITYP,3)
        IS =  NFIRST(I)
        IF =  NLAST(I)
        PRTINT =  0.D0
        DO  520  II=IS,IF
          DO  510  JJ=IS,IF
            T =  1.0D0/(ETA(II,4) + ETA(JJ,4))
            SOO =  PITERN*(T**1.5D0)*ETA(II,5)*ETA(JJ,5)
          T1 =  OVLAP(L,L,0.0D0,0.0D0,T)
          T2 =  OVLAP(M,M,0.0D0,0.0D0,T)
          T3 =  OVLAP(N,N,0.0D0,0.0D0,T)
 510        PRTINT =  PRTINT + C(II)*C(JJ)*SOO*T1*T2*T3
 520      CONTINUE
        PRTINT =  1.0D0/DSQRT(PRTINT)
        DO  530  K=IS,IF
          C(K) =  C(K)*PRTINT
        ETA(K,5)=ETA(K,5)*C(K)
          IJ =  K - IS + 1
CW
          WRITE (6,917)  K,I,IJ,NCNTR(I),NTYPE(I),NR(ITYP,1),NR(ITYP,2)
     $                  ,NR(ITYP,3),ETA(K,4),C(K),ETA(K,5)
 530      CONTINUE
 540    CONTINUE
C=======================================================================
C
C550  DO  560  K=1,NGAUS
C560    ETA(K,5) =  ETA(K,5)*C(K)
C  CALCULATE THE NUCLEAR REPULSION ENERGY.
C
C=======================================================================
      ENERGY =  0.D0
      IF  ( NOC .LE. 1 )   GO TO  630
CW
      WRITE (6,919)
      NONM1 =  NON - 1
      DO  620  I=1,NONM1
        IP1 =  I + 1
        DO  610  J=IP1,NON
          RIJ =DSQRT( (VLIST(I,1) - VLIST(J,1))**2 + (VLIST(I,2)
     X               - VLIST(J,2))**2 + (VLIST(I,3) - VLIST(J,3))**2 )
          RIJA =  RIJ*0.529167D0
CW
          WRITE (6,920)  ICNTR(I), ICNTR(J), RIJ, RIJA
 610      ENERGY =  ENERGY + VLIST(I,4)*VLIST(J,4)/RIJ
 620    CONTINUE
CW
 630  WRITE (6,921)  ENERGY

c***
      deallocate( c )
      RETURN
 999  write(6,"('END OF FILE IN UNIT 5 (RDINPT)')")
      stop'ABORT...'
c***
 900  FORMAT ( 10I5 )
 901  FORMAT ( / 3X, 32HPROGRAM  CONTROL  OPTIONS  ...  ,10I5   )
 904  FORMAT ( / 3X, 20HNUMBER  OF  NUCLEI =, I5, 15X,
     X    33HNUMBER  OF  ADDITIONAL  CENTERS =, I5  )
 905  FORMAT ( / 10X, 26H* *  NUCLEAR  CENTERS  * * //3X, 6HCENTER, 18X,
     X    12HCOORDINATES, 17X, 6HCHARGE,5X,'LPSEUD',/)
cmb  906 FORMAT(A4,6X,4F12.8,2X,A8)
  906 FORMAT(A4,6X,4F12.8,2X,I5)
cmb 907  FORMAT (3X, A4, 6X, 3F12.8, 6X, F4.1 ,5X,A8)
 907  FORMAT (3X, A4, 6X, 3F12.8, 6X, F7.4 ,5X,I5)
  990 FORMAT(//10X,'*************** IPST = ',I1,' ***************',/)
  991 FORMAT(/10X,'**** THIS IS A PSEUDOPOTENTIAL CALCULATION ****',//)
  995 FORMAT(/10X,'********* PSEUDOPOTENTIAL  PARAMETERS *********',//)
  992 FORMAT(//10X,4F10.4,1X,F4.1)
  993 FORMAT(//10X,3F7.2)
  994 FORMAT(//,2(3F15.4,//))
 908  FORMAT ( / 10X, 29H* *  ADDITIONAL  CENTERS  * * / 3X, 6HCENTER,
     X    18X, 12HCOORDINATES  / )
 909  FORMAT ( // 10X, 44H* *  GAUSSIAN  FUNCTION  SPECIFICATIONS  * *
     X    //3X,31HNUMBER OF PRIMITIVE GAUSSIANS =,I5 / 3X,
     X    31HNUMBER  OF  BASIS  FUNCTIONS  =, I5 //
     X    '    NUMBER OF BASIS F. FOR SCATTERING=',I5//
     X    '    NUMBER OF BASIS F. FOR MO=        ',I5//
     X    3X,8HGAUSSIAN, 3X,
     X    8HFUNCTION, 3X, 9HCOMPONENT, 3X, 6HCENTER, 4X, 4HTYPE, 6X,
     X    8HEXPONENT, 6X, 11HCOEFFICIENT   )
 910  FORMAT ( A4, 6X, A4, I3, 3X, 2D15.8 )
 911  FORMAT (3(3X,I5,3X), 4X,A4,5X,A4,2F15.7 )
 1911 FORMAT( 3(3X,I5,3X), 4X,A4,5X,A4,2F15.7,2X,3A8)
  912 FORMAT(2D15.8)
 913  FORMAT ( // 3X,39HDO NOT CALCULATE PRIMITIVE TWO-ELECTRON   ,
     X    39H INTEGRALS WHOSE PREFACTOR IS LESS THAN,E15.5,//,3X,
     X    53HDO NOT WRITE BASIS FUNCTION INTEGRALS ON TAPE IF THEY ,
     X    14H ARE LESS THAN,E15.5,/ )
 914  FORMAT ( // 3X,35HRECALCULATE ONLY THOSE TWO-ELECTRON
     X 'INTEGRALS WHICH INVOLVE THE   '  ,I5,17H  BASIS FUNCTIONS /
     X    (10X,24I4) )
 915  FORMAT ( 36I2 )
 916  FORMAT (1H1,10X,  44H* *  RENORMALIZE  THE  BASIS  FUNCTIONS  * *
     X    // 3X, 8HGAUSSIAN, 3X, 8HFUNCTION, 3X, 9HCOMPONENT, 3X,
     X    6HCENTER, 4X, 4HTYPE, 5X, 1HL, 5X, 1HM, 5X, 1HN, 6X,
     X    8HEXPONENT, 6X, 11HCOEFFICIENT  ,6X,'   ETA(K,5)')
 917  FORMAT ( 3X,I5,6X,I5,6X,I4,5X,I6,4X,I6,1X,3I6,2X,3D15.8 )
 918  FORMAT ( 24I3 )
 919  FORMAT ( // 10X, 36HINTERNUCLEAR DISTANCES FROM GEOMETRY  //
     X    8X,8H CENTERS,11X, 4HA.U., 10X, 2HA. )
 920  FORMAT ( 5X, A4, 3H - , A4, 2F14.6 )
 921  FORMAT ( // 3X, 27HNUCLEAR REPULSION ENERGY = , F14.8, 6H  A.U.)
 922  FORMAT ( /3X,46HTHE BASIS FUNCTIONS ARE NOT IN STANDARD ORDER )
 923  FORMAT ( /3X,49HTHE BASIS FUNCTIONS ARE LISTED IN STANDARD ORDER )
  930 FORMAT(18A4)
  931 FORMAT(1H1//5X,18A4 //)
935   FORMAT(//10X,'**  INCOMPATIBLE  ICON(1)  AND  ICON(2)  PARAMETERS
     2 **',/)
 940  FORMAT ( // 10X,26H**  TOO  MANY  CENTERS  ** / )
 941  FORMAT ( // 10X,29H**  TOO MANY  PRIMITIVES  ** / )
 942  FORMAT ( // 10X,28H**  TOO  MANY  FUNCTION  **  / )
 943  FORMAT ( // 10X,'**  ERROR  IN  NUMBER  OF  PRIMITIVES  PER
     X FUNCTION  **', / )
 944  FORMAT ( // 10X,39H**  INCORRECT  PREVIOUS  FUNCTION  FOR
     X     9H FUNCTION,I4,6X,6HKCNT= ,A4,8H  KTYP= ,A4,6H  INC=,I4,
     X    4H  ** / )
 945  FORMAT ( // 10X,37H**  NUMBER OF PRIVITIVES IN FUNCTIONS,I5,
     X    4H AND,I5,15H  NOT EQUAL  ** / )
 946  FORMAT ( // 10X,22H**  ZERO  EXPONENT  ** / )
 947  FORMAT (    10X,39H**  ZERO  COEFFICIENT  SET  TO  ONE  **  )
 948  FORMAT ( // 10X,37H**  CENTERS  NOT  SAME FOR  FUNCTION,I5,5X,
     X    10HPRIMITIVE ,I4,4H  ** / )
 949  FORMAT ( // 10X,35H**  TYPES  NOT  SAME  FOR  FUNCTION,I5,5X,
     X    10HPRIMITIVE ,I4,4H  ** / )
 950  FORMAT ( // 10X,36H**  UNDEFINED  CENTER  FOR  FUNCTION,I5,4H  **)
 951  FORMAT ( // 10X,34H**  UNALLOWED  TYPE  FOR  FUNCTION,I5,4H  **/)
952   FORMAT(//10X,'**',I3,'  ERROR(S).  ANOTHER  RUN  FOR  THE  SEUCR
     2MAN  **'/)
 954  FORMAT ( // 10X,38H**  UNDEFINED  FUNCTION  IN  MLIST  ** / )

      end subroutine rdinpt
c----------------------------------------------------------------

c----------------------------------------------------------------
c     SUBROUTINE GINTS (NTYPE,NR,NFIRST,NLAST,ETA,NBFXX,
c    X  NTXX,NGXX,NBFNS,NSBF,NMOF,NTRN,KEYORT,NPRALL,KEYPOLY)
!     subroutine gints(NTRN)
      subroutine gints

      use ampoly, only: ovlap

      implicit none

!     integer :: NTRN
!     integer, intent(in) :: NTRN

c*** 
c  This common block is kept for use in ampoly
c
      integer :: icon
      COMMON/IOIND/ICON(10)
c***

c***
c  Local variables
c
clinux      INTEGER*2 MTRANS,LIST,ITEMP,NOS,NSKIP
c
c     INTEGER*2 MTRANS,LIST,ITEMP
c
clinux
c
      integer :: i, j, m, ip, it, jt, iprdt, itag, mx, ipt
      integer :: ic, ix, jx, nom, index, ityp, jtyp, jj, ij
      integer :: l1, l2, m1, m2, n1, n2, is, if, js, jf, ii
      integer :: k, l, kk, ierr, nbfn, indexa, iflst, NOFILD
      integer :: n, ikeep, jkeep, noory, noor2, nmof1, nsbf1
      real*8  :: a, b, t, p1, p2, p3, ab1, ab2, ab3, distab
      real*8  :: soo, pax, pbx, pay, pby, paz, pbz, t1, t2
      real*8  :: t3, diff, prvint, rawint, x, y
      real*8  :: CHAR(3)

      integer, allocatable :: nskip(:)
      real*8, allocatable :: ovs(:,:), ovl(:), s(:)

      DATA CHAR/1H ,1H+,1H-/
c***

      NBFN = NBFNS
      IERR =  0
      IF  ( ICON(4) .EQ. 1 )   GO TO  10
      WRITE  (6,992)
   10      CONTINUE
      INDEXA=  NBFNS*(NBFNS + 1)/2
c***
c  Allocate relevant variables
c
      allocate( ovl(indexa), s(indexa), nskip(indexa) )
      allocate( ovs(nbfns,nbfns), ovpm(nbfns,nbfns) )
      allocate( list(nbfns,ntrnx), mtrans(nbfns,ntrnx) )
      allocate( nohin(nbfns), nos(nbfns) )
c***
      DO  4  I=1,INDEXA
      OVL(I)=0.0D0
      NSKIP(I)=0
   4   S(I) =  0.D0
      IFLST=0
      LIST(1,6)=0
      DO 25 I=1,NBFN
      LIST(1,1)=I
      DO 25 J=1,I
      LIST(1,2)=J
      NOFILD=1
      IP=NONE3(I)+J
      NOS(1)=IP
        IF(NTRN.EQ.0) GO TO 23
      IF(NSKIP(IP).EQ.1) GO TO 25
          DO  22  M=1,NTRN
            IT =  MTRANS(I,M)
            JT =  MTRANS(J,M)
            IPRDT =  IT*JT
            ITAG =  1
            IF  ( IPRDT )   12,22,13
  12        ITAG =  2
  13        IT =  IABS(IT)
            JT =  IABS(JT)
            IF  ( IT - JT )   14,15,15
  14        MX =  IT
            IT =  JT
            JT =  MX
   15 IPT=NONE3(IT)+JT
      IF(IP-IPT) 17,16,25
   16 IF(IPRDT) 25,22,22
   17 IF(NOFILD-2) 21,18,18
   18 DO 20 IC=2,NOFILD
      IF(IPT-NOS(IC)) 20,22,20
   20 CONTINUE
   21 NOFILD=NOFILD+1
      NOS(NOFILD)=IPT
      LIST(NOFILD,1)=IT
      LIST(NOFILD,2)=JT
      LIST(NOFILD,6)=ITAG
   22 CONTINUE
   23      CONTINUE
      DO 914 M=1,NOFILD
      IX=LIST(M,1)
      JX=LIST(M,2)
      ITAG=LIST(M,6)
      NOM=NOS(M)
      NSKIP(NOM)=1
      IF  ( ICON(4) .EQ. 1 )   GO TO 550
      INDEX=NONE3(IX)+JX
      S(INDEX)=1.D0
      GO TO 403
  550 CONTINUE
      IF (ITAG-1) 403,402,408
  408 OVL(NOM)=-PRVINT
      GO TO 916
  402 OVL(NOM)=PRVINT
      GO TO 916
  403 OVL(NOM)=0.D0
      ITYP=NTYPE(IX)
      JTYP=NTYPE(JX)
      L1=NR(ITYP,1)
      L2=NR(JTYP,1)
      M1=NR(ITYP,2)
      M2=NR(JTYP,2)
      N1=NR(ITYP,3)
      N2=NR(JTYP,3)
      IS=NFIRST(IX)
      IF=NLAST(IX)
      JS=NFIRST(JX)
      JF=NLAST(JX)
      DO 635 II=IS,IF
      A=ETA(II,4)
      DO 1635 JJ=JS,JF
      B=ETA(JJ,4)
      T=1.D0/(A+B)
      P1=(A*ETA(II,1)+B*ETA(JJ,1))*T
      P2=(A*ETA(II,2)+B*ETA(JJ,2))*T
      P3=(A*ETA(II,3)+B*ETA(JJ,3))*T
      AB1=ETA(II,1)-ETA(JJ,1)
      AB2=ETA(II,2)-ETA(JJ,2)
      AB3=ETA(II,3)-ETA(JJ,3)
      DISTAB=AB1*AB1+AB2*AB2+AB3*AB3
      SOO=(PITERN*T **1.5D0)*DEXP(-A*B*DISTAB*T)*ETA(II,5)*ETA(JJ,5)
      PAX=P1-ETA(II,1)
      PBX=P1-ETA(JJ,1)
      PAY=P2-ETA(II,2)
      PBY=P2-ETA(JJ,2)
      PAZ=P3-ETA(II,3)
      PBZ=P3-ETA(JJ,3)
      T1=OVLAP(L1,L2,PAX,PBX,T)
      T2=OVLAP(M1,M2,PAY,PBY,T)
      T3=OVLAP(N1,N2,PAZ,PBZ,T)
 1635 OVL(NOM)=OVL(NOM)+SOO*T1*T2*T3
  635 CONTINUE
      IF  ( ICON(4) .EQ. 1 )   GO TO 510
      IF (ITAG-1) 510,511,512
  511 DIFF=OVL(NOM)-PRVINT
      GO TO 513
  512 DIFF=OVL(NOM)+PRVINT
  513 IF (DABS(DIFF).LT.1.0D-06) GO TO 916
      WRITE  (6,520)  IKEEP, JKEEP, PRVINT, I, J, CHAR(ITAG+1),OVL(NOM)
      IERR =  IERR + 1
      GO TO 916
  510 PRVINT=OVL(NOM)
      IKEEP=I
      JKEEP=J
  916 CONTINUE
  914 CONTINUE
   25      CONTINUE
 915  IF  ( ICON(4) .EQ. 1 )   RETURN
      DO 581 I=1,NBFNS
      DO 580 J=1,I
      INDEX=NONE3(I)+J
      IF (S(INDEX).NE.0.D0)GO TO 580
      RAWINT=0.D0
      ITYP=NTYPE(I)
      JTYP=NTYPE(J)
      L1=NR(ITYP,1)
      L2=NR(JTYP,1)
      M1=NR(ITYP,2)
      M2=NR(JTYP,2)
      N1=NR(ITYP,3)
      N2=NR(JTYP,3)
      IS=NFIRST(I)
      IF=NLAST(I)
      JS=NFIRST(J)
      JF=NLAST(J)
      DO 536 II=IS,IF
      A=ETA(II,4)
      DO 535 JJ=JS,JF
      B=ETA(JJ,4)
      T=1.D0/(A+B)
      P1=(A*ETA(II,1)+B*ETA(JJ,1))*T
      P2=(A*ETA(II,2)+B*ETA(JJ,2))*T
      P3=(A*ETA(II,3)+B*ETA(JJ,3))*T
      AB1=ETA(II,1)-ETA(JJ,1)
      AB2=ETA(II,2)-ETA(JJ,2)
      AB3=ETA(II,3)-ETA(JJ,3)
      DISTAB=AB1*AB1+AB2*AB2+AB3*AB3
      SOO=(PITERN*T **1.5D0)*DEXP(-A*B*DISTAB*T)*ETA(II,5)*ETA(JJ,5)
      PAX=P1-ETA(II,1)
      PBX=P1-ETA(JJ,1)
      PAY=P2-ETA(II,2)
      PBY=P2-ETA(JJ,2)
      PAZ=P3-ETA(II,3)
      PBZ=P3-ETA(JJ,3)
      T1=OVLAP(L1,L2,PAX,PBX,T)
      T2=OVLAP(M1,M2,PAY,PBY,T)
      T3=OVLAP(N1,N2,PAZ,PBZ,T)
  535 RAWINT=RAWINT+SOO*T1*T2*T3
  536 CONTINUE
      IF (DABS(RAWINT).LT.1.0D-07) GO TO 580
      IERR =  IERR + 1
      WRITE (6,585) I,J,RAWINT
  580 CONTINUE
  581 CONTINUE
      IF  ( IERR .EQ. 0 )   GO TO 702
      WRITE (6,993)   IERR
      STOP
  520 FORMAT (3X, ' SYMMETRY ERROR',5X,'I=',I3,3X,'J=',I3,'PRVINT=',F14.
     X8,5X,2HI=,I3,3X,2HJ=,I3,3X,4HTAG=,A1,3X,9HINTEGRAL=,F14.8)
  585 FORMAT (3X,'ZERO INTEGRAL',  2I4,3X,'ACTUALLY IS  ' ,F14.8)
992   FORMAT(9X,'TEST SYMMETRY'//)
 993  FORMAT ( / 3X, 37H**  GINTS  CANNOT  CONTINUE ,  IERR =,I5,4H  **)
  702      CONTINUE
C
      IJ=0
      DO 700 I=1,NBFNS
      DO 700 J=1,I
      IJ=IJ+1
      OVS(I,J)=OVL(IJ)
  700      OVS(J,I)=OVL(IJ)
C=======================================================================
      IF(NPRALL.EQ.0)GO TO 777
      WRITE(6,778)
  778 FORMAT(/30X,'################ OVS (GINTS)################')
C      CALL PRINT(OVS,NBFNX,NBFNX,NBFN,NBFN)
  777 CONTINUE
C=======================================================================
C BEGINS ORTHOGONALIZATION OF VIRTUAL ORBITALS
C=======================================================================
        DO 950 I=1,NMOF
        DO 960 J=1,I
      IF(J.EQ.I) GO TO 968
      X=0.0D0
      DO 966 K=1,NSBF
      Y=0.0D0
      DO 967 L=1,NSBF
  967      Y=Y+OVS(L,K)*d(L,I)
  966      X=X+Y*d(K,J)
      IF(DABS(X).LT.1.0D-12) GO TO 960
      DO 965 K=1,NSBF
  965      d(K,I)=d(K,I)-X*d(K,J)
      GO TO 960
  968      A=0.0D0
      DO 970      K=1,NSBF
      B=0.0D0
      DO 971 L=1,NSBF
  971      B=B+OVS(L,K)*d(L,I)
  970      A=A+B*d(K,I)
      A=1.0D0/DSQRT(A)
      DO 972 K=1,NSBF
  972      d(K,I)=d(K,I)*A
  960      CONTINUE
  950      CONTINUE
C=======================================================================
C ENDS ORTHOGONALIZATION OF VIRTUAL ORBITALS
C=======================================================================
      NMOF1=NMOF+1
      DO 710 I=1,NMOF
      Y=0.0D0
      DO 711 J=1,NMOF
      X=0.0D0
      DO 712 K=1,NMOF
  712      X=X+OVS(K,J)*d(K,I)
  711      Y=Y+X*d(J,I)
      IF(DABS(Y-1.0D0).LE.1.0D-10) GO TO 714
      Y=1.0D0/DSQRT(Y)
      DO 713 J=1,NMOF
  713      d(J,I)=d(J,I)*Y
  714      CONTINUE
      IF(NMOF.EQ.NSBF) GO TO 710
      DO 715 J=NMOF1,NSBF
  715      d(J,I)=0.0D0
  710      CONTINUE
      NOORY=NOOR
      IF(NMOF.EQ.NSBF) GO TO 800
C=======================================================================
C
C      ORTHOGONALIZATION OF ADDITIONAL SCATTERING BASIS FUNCTIONS TO
C      MOLECULAR ORBITALS.
C
C=======================================================================
        IF(KEYORT.NE.0)THEN
          write(6,"('GINTS: KEYORT<>0 NOT ALLOWED')")
          stop'ABORT...'
        ELSE
        DO 370 I=1,NSBF
 370    NOH(I)=I
        NOOR=NMOF
        NOOR2=NSBF
        ENDIF
      WRITE(6,316) NOOR,NOOR2,(NOH(K),K=1,NOOR2)
  316 FORMAT(//'%%% ORTHOGONALIZATION INFORMATION %%%%%'/
     1      '   NOOR,NOOR2=',2I3/'   NOH=',30I3)
        DO 750 II=NOOR+1,NOOR2
        I=NOH(II)
        DO 760 JJ=1,II
      IF(JJ.EQ.II) GO TO 768
        J=NOH(JJ)
      X=0.0D0
      DO 766 K=1,NSBF
      Y=0.0D0
      DO 767 L=1,NSBF
  767      Y=Y+OVS(L,K)*d(L,I)
  766      X=X+Y*d(K,J)
      IF(DABS(X).LT.1.0D-8) GO TO 760
      DO 765 K=1,NSBF
  765      d(K,I)=d(K,I)-X*d(K,J)
      GO TO 760
  768      A=0.0D0
      DO 770      K=1,NSBF
      B=0.0D0
      DO 771 L=1,NSBF
  771      B=B+OVS(L,K)*d(L,I)
  770      A=A+B*d(K,I)
      A=1.0D0/DSQRT(A)
      DO 772 K=1,NSBF
  772      d(K,I)=d(K,I)*A
  760      CONTINUE
  750      CONTINUE
  800      DO 780 I=1,NSBF
      DO 781 M=1,NSBF
      X=0.0D0
      DO 782 N=1,NSBF
  782      X=X+OVS(N,M)*d(N,I)
  781      S(M)=X
      DO 783 J=1,I
      Y=0.0D0
      DO 784 M=1,NSBF
  784      Y=Y+S(M)*d(M,J)
      OVPM(I,J)=Y
  783      OVPM(J,I)=Y
  780      CONTINUE
          IF(NSBF.EQ.NBFN) GO TO 850
      NSBF1=NSBF+1
      DO 810 I=NSBF1,NBFN
      DO 801 J=1,NSBF
      X=0.0D0
      DO 802 K=1,NSBF
  802      X=X+OVS(K,I)*d(K,J)
      OVPM(J,I)=X
  801      OVPM(I,J)=X
      DO 803 J=NSBF1,I
      OVPM(J,I)=OVS(J,I)
  803      OVPM(I,J)=OVS(J,I)
  810      CONTINUE
  850      WRITE(6,860)
  860 FORMAT(//'### ORTHOGONALIZED BASIS FUNCTIONS ###')
C      CALL PRINT(C,NSBFX,NSBFX,NSBF,NSBF)
      WRITE(6,888)
  888 FORMAT(//'### OVERLAP MATRIX OF NBFN BY NBFN ###')

CBUG
c     IF(KEYPOLY.NE.0)THEN
c       OPEN (UNIT=10,FORM='FORMATTED',FILE='fort.10')
c     CALL WRITED(C,NSBF)
c       CLOSE(10)
c     STOP
c     END IF
CBUG
C      CALL PRINT(OVPM,NBFNX,NBFNX,NBFN,NBFN)

c***
c  ovl should be reallocated before calling ovlgrn
c
      deallocate( ovl )
      allocate( ovl(nbfns*nbfns) )
c***

      CALL OVLGRN(NBFN,NSBF,NMOF,OVL,OVS,NPRALL)
      NOOR=NOORY
      DO 444 K=1,NOOR
      KK=NOH(K)
  444      NOHIN(KK)=K
cmod
c     do i=1,nbfns
c       do j= 1,nbfns
c         write(1,"(i4,i4,1x,1pe15.8)") i,j,ovpm(i,j)
c       end do
c     end do
cmod

c***
c   Deallocate local arrays
c
      deallocate( ovs, ovl, s )
      deallocate( nskip, list, mtrans, nos )
c***

      end subroutine gints
c----------------------------------------------------------------

c----------------------------------------------------------------
c     SUBROUTINE  VINTS (NOC,VLIST,LPSKIP,NTYPE,NR,NFIRST,NLAST,
c    $                  ETA,NBFXX,NCXX,NTXX,NINMAX,NGXX,NBFNS,NSBF,NMOF,
c    $                   NTRN,NREC50,NPRALL)
!     subroutine vints (NINMAX,NTRN,NREC50)
!     subroutine vints(NTRN,NREC50)
      subroutine vints(NREC50)

      use pseudo_pp, only : PSEUDO
      use ampoly, only: f_gam, fmch, dawtab, gfunct, vaaa, vbaa, vbca

      implicit none

!     integer :: NINMAX, NTRN, NREC50
!     integer :: NTRN, NREC50
      integer, intent(in) :: NREC50

c***
c  In the f77 version, NOC is (locally) a dummy name for NON. This would
c  cause problems in this version, since NON and NOC are declared in the
c  module. Use NOC_1 instead of NOC *locally* (only in this subroutine).
c
      integer :: noc_1

c***
c  This common block will be kept for use in ampoly
c
cccc  real*8  :: F0,F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12
cccc  COMMON /GAMMA /  F0,F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12
c***

      integer :: LEXP(3,2)
      real*8  :: ALFA(2),CORE(2),ALFCOR(2),CENTRO(3,2),
     $           SIGMA(3,3),AION(2,3,3)

c***
c This common block is kept for use in ampoly
      integer :: ITYP, JTYP, N, IDUM
      real*8  :: XINT,A,B,AVX,AVY,AVZ,AV,BVX,BVY,BVZ,BV,PCX,PCY,PCZ,
     .           PCSQ, PHASE
      COMMON/BFCOM/XINT,ITYP,JTYP,A,B,N,IDUM,AVX,AVY,AVZ,AV,BVX,BVY,BVZ,
     #  BV,PCX,PCY,PCZ,PCSQ,PHASE
c***

c***
c  Local variables
c
clinux      INTEGER*2 MTRANS,LIST,ITEMP,NOS,NSKIP
c
c     INTEGER*2 MTRANS,LIST,ITEMP,NSKIP
c
clinux
c

      logical :: good = .true.
      integer :: i, j, m, ic, ip, it, jt, imn, ixx, jxx, mx
      integer :: l1, l2, m1, m2, n1, n2, my, mz, is, if, ipt
      integer :: l, js, jf, ii, ji, ll, mm, ix, jy, kz, mxyz
      integer :: k, jj, iprdt, itag, nom, itypt, iitemp
      integer :: kont1, kont2, kont3, klx, klm, kl, numb
      integer :: NBFN, INDEXA, NOFILD, NSBF1, njk, ibhs, l_at
      real*8  :: prvint, ax, ay, az, bx, by, bz, t, t1
      real*8  :: p1, p2, p3, ab1, ab2, ab3, distab, soo
      real*8  :: vsoo, vlp, vnai, zval, vx, vy, vz 
      real*8  :: arg, y, pax, pbx, pay, pby, paz, pbz
      real*8  :: rawint, bnuc, avsq
      real*8  :: x_at, y_at, z_at, zval_at, too_far = 1.d-06
      real*8  :: bvsq, xl, xij, x
!     real*8  :: F(13), G(7,3)
      real*8  :: G(7,3)

      integer, allocatable :: nskip(:)
      real*8, allocatable  :: atr(:), ri(:) 
c***
cccc  EQUIVALENCE (F0,F(1))

      NBFN = NBFNS
      INDEXA=NONE1(NBFN)+NSBF
c***
      noc_1 = NON

ccc  Allocate atr and try to read it from VNUC.FILE
c
      allocate( atr(indexa) )
      open(unit=n_nucpp,file='VNUC.FILE',form='formatted')

      read(n_nucpp,*,end=354,err=355) ibhs
      if(ibhs .ne. ipst) then
        write(6,"(//,'WRONG IPST READING VNUCP.FILE',/,
     .    'EXPECTED IPST = ',i4,/,'FOUND IPST = ',i4)"),
     .     ipst, ibhs
        stop'ABORT...'
      end if

      read(n_nucpp,*,end=354,err=355) numb
      do i = 1,numb
        read(n_nucpp,*,end=354,err=355) 
     .       l_at, zval_at, x_at, y_at, z_at
        good = ( lpseud(i) .eq. l_at)
        if(.not. good) then
          write(6,"(//,'WRONG LPSEUD IN VNUC.FILE',/,
     .     'ATOM ',i4,/,'EXPECTED: ',i4,/,'FOUND: ',i4)") 
     .      i, lpseud(i), l_at
          stop'ABORT'
        end if
        good = ( abs(blkc1(i)-zval_at) .lt. too_far)
        if(.not. good) then
          write(6,"(//,'VALENCE CHARGE MISMATCH IN VNUC.FILE',/,
     .     'ATOM ',i4,/,'EXPECTED: ',1(1pe15.8),/,'FOUND: ',
     .     1(1pe15.8))") i, blkc1(i), zval_at
          stop'ABORT'
        end if
        good = ( (abs(vlist(i,1)-x_at) .lt. too_far) .and. 
     .           (abs(vlist(i,2)-y_at) .lt. too_far) .and. 
     .           (abs(vlist(i,3)-z_at) .lt. too_far) )
        if(.not. good) then
          write(6,"(//,'COORDINATE MISMATCH IN VNUC.FILE',/,
     .     'ATOM ',i4,/,'EXPECTED: ',3(1pe15.8),/,'FOUND: ',
     .     3(1pe15.8))") i, (vlist(i,j),j=1,3), x_at, y_at, z_at
          stop'ABORT'
        end if
      end do

      read(n_nucpp,*,end=354,err=355) numb
      if(numb .ne. indexa) then
        write(6,"(//,'WRONG DIMENSION IN VNUC.FILE',/,
     .    'EXEPECTED = 'i8,/,'FOUND = ',i8)") indexa, numb
        stop'ABORT...'
      end if
      read(n_nucpp,*,end=356) (atr(i), i=1,indexa)
      write(6,"(//'*** NUCLEAR ATTRACTION INTEGRALS SUCCESSFULLY',
     .  ' READ FROM VNUC.FILE ***',//)")

cc Allocate those needed for transformation to MO basis
      allocate( ri(nbfns), atrt(nbfns,nbfns) )
      goto 777

 354  continue
      write(6,"(//,5x,'*** VNUC.FILE NOT FOUND. CALCULATE',
     .   ' NUCLEAR ATTRACTION INTEGRALS ***')")
      goto 357

 355  continue
      write(6,"(//,5x,'*** ERROR READING VNUC.FILE. CALCULATE',
     .   ' NUCLEAR ATTRACTION INTEGRALS ***')")
      goto 357

 356  continue
      write(6,"(//,5x,'*** END OF FILE READING ATR ARRAY FROM',
     .   ' VNUC.FILE')")
      stop'ABORT...'

 357  continue

c  Allocate relevant variables (nuclear attraction integrals)
      allocate( nskip(indexa) )
      allocate( list(nbfns,ntrnx), mtrans(nbfns,ntrnx) )
      allocate( nos(nbfns), ri(nbfns) )
      allocate( atrt(nbfns,nbfns) )
c***
       DO 100 I=1,INDEXA
       NSKIP(I)=0
  100  ATR(I)=0.0D0
      CALL DAWTAB
      LIST(1,6)=0
      DO 25 I=1,NBFN
      LIST(1,1)=I
      IMN=MIN0(I,NSBF)
      DO 25 J=1,IMN
      LIST(1,2)=J
      NOFILD=1
      IP=NONE1(I)+J
      NOS(1)=IP
        IF(NTRN.EQ.0) GO TO 23
      IF(NSKIP(IP).EQ.1) GO TO 25
          DO  22  M=1,NTRN
            IT =  MTRANS(I,M)
            JT =  MTRANS(J,M)
            IPRDT =  IT*JT
            ITAG =  1
            IF  ( IPRDT )   12,22,13
  12        ITAG =  2
  13        IT =  IABS(IT)
            JT =  IABS(JT)
            IF  ( IT - JT )   14,15,15
  14        MX =  IT
            IT =  JT
            JT =  MX
   15 IPT=NONE1(IT)+JT
      IF(IP-IPT) 17,16,25
   16 IF(IPRDT) 25,22,22
   17 IF(NOFILD-2) 21,18,18
   18 DO 20 IC=2,NOFILD
      IF(IPT-NOS(IC)) 20,22,20
   20 CONTINUE
   21 NOFILD=NOFILD+1
      NOS(NOFILD)=IPT
      LIST(NOFILD,1)=IT
      LIST(NOFILD,2)=JT
      LIST(NOFILD,6)=ITAG
   22 CONTINUE
   23       CONTINUE
      DO  916  M=1,NOFILD
      IXX=LIST(M,1)
      JXX=LIST(M,2)
      ITAG=LIST(M,6)
      NOM=NOS(M)
      NSKIP(NOM)=1
      IF (ITAG-1) 403,402,408
  408 ATR(NOM)=-PRVINT
      GO TO 916
  402 ATR(NOM)=PRVINT
      GO TO 916
  403 ATR(NOM)=0.D0
      ITYP=NTYPE(IXX)
      JTYP=NTYPE(JXX)
      IF(ITYP.LE.JTYP) GO TO 404
      ITYPT=ITYP
      IITEMP=IXX
      ITYP=JTYP
      IXX=JXX
      JTYP=ITYPT
      JXX=IITEMP
  404 CONTINUE
      L1=NR(ITYP,1)
      L2=NR(JTYP,1)
      M1=NR(ITYP,2)
      M2=NR(JTYP,2)
      N1=NR(ITYP,3)
      N2=NR(JTYP,3)
cmb
      IF(IPST.EQ.0)GO TO 301
      LEXP(1,1)=L1
      LEXP(1,2)=L2
      LEXP(2,1)=M1
      LEXP(2,2)=M2
      LEXP(3,1)=N1
      LEXP(3,2)=N2
  301 CONTINUE
cmb
      MX=L1+L2+1
      MY=M1+M2+1
      MZ=N1+N2+1
      IS=NFIRST(IXX)
      IF=NLAST(IXX)
      JS=NFIRST(JXX)
      JF=NLAST(JXX)
      DO 635 II=IS,IF
      A=ETA(II,4)
      AX=ETA(II,1)
      AY=ETA(II,2)
      AZ=ETA(II,3)
      DO 1635 JJ=JS,JF
      B=ETA(JJ,4)
      BX=ETA(JJ,1)
      BY=ETA(JJ,2)
      BZ=ETA(JJ,3)
cmb
      IF(IPST.EQ.1)THEN
      ALFA(1)=A
      ALFA(2)=B
      END IF
cmb
      T1=A+B
      T=1.D0/T1
      P1=(A*AX+B*BX)*T
      P2=(A*AY+B*BY)*T
      P3=(A*AZ+B*BZ)*T
      AB1=AX-BX
      AB2=AY-BY
      AB3=AZ-BZ
      DISTAB=AB1*AB1+AB2*AB2+AB3*AB3
      PHASE=DEXP(-A*B*DISTAB*T)
      SOO=4.D0*PI*ETA(II,5)*ETA(JJ,5)
      VSOO=T*PHASE*0.5D0
      VLP=0.0D0
      VNAI=0.D0
      KONT1=0
      KONT2=0
      KONT3=0
c     DO 690 N=1,NOC
      DO 690 N=1,noc_1
cmb
cmb=============================================
      IF(LPSEUD(N).EQ.1)THEN
      ZVAL=BLKC1(N)
      DO IC=1,2
       KONT1=KONT1+1
       CORE(IC)=BLKC2(KONT1)
       ALFCOR(IC)=BLKC3(KONT1)
      END DO
      DO LL=1,3
      DO JI=1,3
       KONT2=KONT2+1
       SIGMA(JI,LL)=BLKI1(KONT2)
      END DO
       DO MM=1,2
        DO JI=1,3
         KONT3=KONT3+1
         AION(MM,JI,LL)=BLKI2(KONT3)
        END DO
       END DO
      END DO
      END IF
cmb=============================================
cmb
      VX=VLIST(N,1)
      VY=VLIST(N,2)
      VZ=VLIST(N,3)
cmb
      CENTRO(1,1)=AX-VX
      CENTRO(1,2)=BX-VX
      CENTRO(2,1)=AY-VY
      CENTRO(2,2)=BY-VY
      CENTRO(3,1)=AZ-VZ
      CENTRO(3,2)=BZ-VZ
cmb
      PCX=P1-VX
      PCY=P2-VY
      PCZ=P3-VZ
      PCSQ=PCX*PCX+PCY*PCY+PCZ*PCZ

      ARG=T1*PCSQ
CMTV  Y=DEXP(-ARG)
CALPHA
      IF(ARG.GE.709.0D0) THEN
                          Y = 0.0D0
                         ELSE
                          Y=DEXP(-ARG)
      END IF
CALPHA
cmod
c     F6=FMCH(6,ARG,Y)
c     ARG=2.D0*ARG
c     F5=(ARG*F6+Y)*X11
c     F4=(ARG*F5+Y)*X9
c     F3=(ARG*F4+Y)*X7
c     F2=(ARG*F3+Y)*X5
c     F1=(ARG*F2+Y)*X3
c     F0= ARG*F1+Y
c
      f_gam(7)=FMCH(6,ARG,Y)
      ARG=2.D0*ARG
      f_gam(6)=(ARG*f_gam(7)+Y)*X11
      f_gam(5)=(ARG*f_gam(6)+Y)*X9
      f_gam(4)=(ARG*f_gam(5)+Y)*X7
      f_gam(3)=(ARG*f_gam(4)+Y)*X5
      f_gam(2)=(ARG*f_gam(3)+Y)*X3
      f_gam(1)= ARG*f_gam(2)+Y
cmod

      PAX=P1-AX
      PBX=P1-BX

      CALL  GFUNCT(L1,L2,PAX,PBX,PCX,T,G,1)
      PAY=P2-AY
      PBY=P2-BY
      CALL  GFUNCT(M1,M2,PAY,PBY,PCY,T,G,2)
      PAZ=P3-AZ
      PBZ=P3-BZ
      CALL  GFUNCT(N1,N2,PAZ,PBZ,PCZ,T,G,3)
      RAWINT=0.D0
      DO 506 IX=1,MX
      DO 507 JY=1,MY
      DO 508 KZ=1,MZ
      MXYZ=IX+JY+KZ-2

cmod
c     write(93,"(i8,2x,1pe15.8)")mxyz,f_gam(mxyz)
cmod

cmod
c 508 RAWINT=RAWINT+G(IX,1)*G(JY,2)*G(KZ,3)*F(MXYZ)
  508 RAWINT=RAWINT+G(IX,1)*G(JY,2)*G(KZ,3)*f_gam(MXYZ)
cmod
  507 CONTINUE
  506 CONTINUE
cmb=============================================
      BNUC=-RAWINT*VLIST(N,4)*VSOO
      IF(LPSEUD(N).EQ.0)GO TO 303

      BNUC=PSEUDO(LGUIMA,ZVAL,CORE,ALFCOR,AION,SIGMA,
     $            CENTRO,LEXP,ALFA)

  303 CONTINUE
cmb=============================================
      XINT=0.0D0
      IF(LPSKIP(N)) 600,600,699
cmb!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  600 IF((DISTAB+PCSQ).GT.1.0D-16)  GO TO 602
  601 CALL VAAA
      VLP=VLP+XINT
      GO TO 699
  602 AVX=AX-VX
      AVY=AY-VY
      AVZ=AZ-VZ
      AVSQ=AVX*AVX+AVY*AVY+AVZ*AVZ
      BVX=BX-VX
      BVY=BY-VY
      BVZ=BZ-VZ
      BVSQ=BVX*BVX+BVY*BVY+BVZ*BVZ
      IF(AVSQ) 603,603,604
  603 BV=DSQRT(BVSQ)
      CALL VBAA(BVX,BVY,BVZ,BV,B,JTYP,A,ITYP)
      VLP=VLP+XINT
      GO TO 699
  604 IF(BVSQ) 605,605,606
  605 AV=DSQRT(AVSQ)
      CALL VBAA(AVX,AVY,AVZ,AV,A,ITYP,B,JTYP)
      VLP=VLP+XINT
      GO TO 699
  606 AV=DSQRT(AVSQ)
      BV=DSQRT(BVSQ)
      CALL VBCA
      VLP=VLP+XINT
cmb!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cmb=============================================
  699 CONTINUE
  690 VNAI=VNAI+BNUC

 1635 ATR(NOM)=ATR(NOM)+(VNAI+VLP)*SOO

cmb=============================================
  635 CONTINUE
      PRVINT=ATR(NOM)
  916 CONTINUE
   25 CONTINUE
  915 CONTINUE
cmod
!     open(unit=1,file='atr.lis',form='formatted')
!     do i=1,indexa
!       write(1,"(i8,1x,1pe17.10)") i,atr(i)
!     end do
!     close(1)
cmod
C=======================================================================
      IF(NPRALL.EQ.0)GO TO 777
       WRITE(6,7777)
 7777 FORMAT(//'&&&& ATR (VINTS) &&&&')
c      CALL PRINL(ATR,NINDEX,INDEXA)
  777 CONTINUE
C=======================================================================
C      TRANSFORMATION FOR NUCLEAR ATTRACTION ENERGY
C       DO 711 I=1,NBFN
C      DO 711 J=1,NSBF
C      IJ=NONE1(MAX0(I,J))+MIN0(I,J)
C      XYZ(I,J)=ATR(IJ)
C 711  CONTINUE
C      WRITE(6, 722)
C 722  FORMAT(//' ATR')
C      CALL PRINT(XYZ,40,40,NBFN,NSBF)
C
C=======================================================================
      NSBF1=NSBF+1
      DO 520 I=1,NSBF
      DO 530 L=1,NSBF
      XL=0.0D0
      DO 540 K=1,NSBF
      KLX=MAX0(K,L)
      KLM=MIN0(K,L)
      KL=NONE1(KLX)+KLM
  540 XL=XL+d(K,I)*ATR(KL)
  530 RI(L)=XL
      DO 560 J=1,I
      XIJ=0.0D0
      DO 570 L=1,NSBF
  570 XIJ=XIJ+d(L,J)*RI(L)
      ATRT(I,J)=XIJ
  560 ATRT(J,I)=XIJ
  520 CONTINUE
      IF(NSBF.EQ.NBFN) GO TO 521
      NSBF1=NSBF+1
      DO 522 J=NSBF1,NBFN
      DO 522 I=1,NSBF
      X=0.0D0
      DO 523 K=1,NSBF
      NJK=NONE1(J)+K
  523 X=X+ATR(NJK)*d(K,I)
      ATRT(I,J)=X
  522 ATRT(J,I)=X
  521 CONTINUE
c***
c  Deallocate those no longer needed
c
      if( allocated(nskip) )  deallocate( nskip )
      if( allocated(list) )   deallocate( list )
      if( allocated(mtrans) ) deallocate( mtrans )
      if( allocated(nos) )    deallocate( nos )
!fk   deallocate( ri, atr )
c***
cmod
c     open(unit=1,file='atrt.lis',form='formatted')
c     do i=1,nbfns
c       do j= 1,nbfns
c         write(1,"(i4,i4,1x,1pe15.8)") i,j,atrt(i,j)
c       end do
c     end do
c     close(1)
cmod
      IF(NREC50.EQ.0) GO TO 920
      WRITE(9) NSBF,NBFN
      DO 921 I=1,NSBF
  921 WRITE(9) (ATRT(J,I),J=1,NBFN)
  920      CONTINUE
C=======================================================================
      IF(NPRALL.EQ.0)GO TO 888
       WRITE(6,8888)
 8888 FORMAT(//'&&&& ATRT (VINTS) &&&&')
C       CALL PRINT(ATRT,NBFNX,NBFNX,NBFN,NSBF)
  888 CONTINUE
C=======================================================================
ccc   RETURN
ccc   END

!fk
      open(unit=99,file='VNUC.FILE',form='formatted')
      write(99,*) ipst
      write(99,*) noc_1
      do n = 1,noc_1
        write(99,*) lpseud(n), blkc1(n), (vlist(n,i), i=1,3)
      end do
      write(99,*) indexa
      write(99,*) (atr(i), i=1,indexa)
      close(99)

      deallocate( ri, atr )

      end subroutine vints 
c----------------------------------------------------------------

c----------------------------------------------------------------
c    SUBROUTINE TINTS (NTYPE,NR,NFIRST,NLAST,ETA,NBFXX,
c    $      NTXX,NINMAX,NGXX,NBFNS,NSBF,NMOF,NTRN,NREC40,NPRALL)
!     subroutine tints(NTRN,NREC40)
      subroutine tints(NREC40)

      use ampoly, only: ovlap

      implicit none

!     integer :: NINMAX,NTRN,NREC40
!     integer :: NTRN,NREC40
      integer, intent(in) :: NREC40

c
clinux      INTEGER*2 MTRANS,LIST,ITEMP,NOS,NSKIP
c
c     INTEGER*2 MTRANS,LIST,ITEMP,NSKIP
c
clinux
c

c***
c  Local variables
c
      integer :: i, j, m, imn, ip, it, jt, mx, l1, l2, m1, m2
      integer :: jj, n1, n2, is, if, js, jf, ii, k, l, kl, njk
      integer :: ic, ipt, ix, jx, nom, ityp, jtyp
      integer :: NBFN, INDEXA, NOFILD, NSBF1, iprdt, itag
      real*8  :: a, b, t, p1, p2, p3, ab1, ab2, ab3, distab, soo
      real*8  :: pax, pbx, pay, pby, paz, pbz, t1, t2, t3, s1, xl
      real*8  :: x, xij, s2, s3, part, ttke, prvint

      integer, allocatable :: nskip(:)
      real*8, allocatable :: tke(:), ri(:) 
c***

C=======================================================================
C
C      TKE(NSBF,NSBF) LINEARIZED WITH NONE1
C
C=======================================================================
      NBFN=NBFNS
      INDEXA=NONE1(NBFN)+NSBF
c***
c  Allocate relevant variables
      allocate( nskip(indexa), tke(indexa) )
      allocate( list(nbfns,ntrnx), mtrans(nbfns,ntrnx) )
      allocate( nos(nbfns), ri(nbfns) )
      allocate( tket(nbfns,nbfns) )
c***
      DO 100 I=1,INDEXA
      TKE(I)=0.0D0
  100  NSKIP(I)=0
      LIST(1,6)=0
      DO 25 I=1,NBFNS
      LIST(1,1)=I
      IMN=MIN0(I,NSBF)
      DO 25 J=1,IMN
      LIST(1,2)=J
      NOFILD=1
      IP=NONE1(I)+J
      NOS(1)=IP
        IF(NTRN.EQ.0) GO TO 23
      IF(NSKIP(IP).EQ.1) GO TO 25
          DO  22  M=1,NTRN
            IT =  MTRANS(I,M)
            JT =  MTRANS(J,M)
            IPRDT =  IT*JT
            ITAG =  1
            IF  ( IPRDT )   12,22,13
  12        ITAG =  2
  13        IT =  IABS(IT)
            JT =  IABS(JT)
            IF  ( IT - JT )   14,15,15
  14        MX =  IT
            IT =  JT
            JT =  MX
   15 IPT=NONE1(IT)+JT
      IF(IP-IPT) 17,16,25
   16 IF(IPRDT) 25,22,22
   17 IF(NOFILD-2) 21,18,18
   18 DO 20 IC=2,NOFILD
      IF(IPT-NOS(IC)) 20,22,20
   20 CONTINUE
   21 NOFILD=NOFILD+1
      NOS(NOFILD)=IPT
      LIST(NOFILD,1)=IT
      LIST(NOFILD,2)=JT
      LIST(NOFILD,6)=ITAG
   22 CONTINUE
   23       CONTINUE
      DO  916  M=1,NOFILD
      IX=LIST(M,1)
      JX=LIST(M,2)
      ITAG=LIST(M,6)
      NOM=NOS(M)
      NSKIP(NOM)=1
      IF (ITAG-1) 403,402,408
  408 TKE(NOM)=-PRVINT
      GO TO 916
  402 TKE(NOM)=PRVINT
      GO TO 916
  403 TKE(NOM)=0.D0
      ITYP=NTYPE(IX)
      JTYP=NTYPE(JX)
      L1=NR(ITYP,1)
      L2=NR(JTYP,1)
      M1=NR(ITYP,2)
      M2=NR(JTYP,2)
      N1=NR(ITYP,3)
      N2=NR(JTYP,3)
      IS=NFIRST(IX)
      IF=NLAST(IX)
      JS=NFIRST(JX)
      JF=NLAST(JX)
      DO 635 II=IS,IF
      A=ETA(II,4)
      DO 1635 JJ=JS,JF
      B=ETA(JJ,4)
      T=1.D0/(A+B)
      P1=(A*ETA(II,1)+B*ETA(JJ,1))*T
      P2=(A*ETA(II,2)+B*ETA(JJ,2))*T
      P3=(A*ETA(II,3)+B*ETA(JJ,3))*T
      AB1=ETA(II,1)-ETA(JJ,1)
      AB2=ETA(II,2)-ETA(JJ,2)
      AB3=ETA(II,3)-ETA(JJ,3)
      DISTAB=AB1*AB1+AB2*AB2+AB3*AB3
      SOO=(PITERN*T**1.5D0)*DEXP(-A*B*DISTAB*T)*ETA(II,5)*ETA(JJ,5)
      PAX=P1-ETA(II,1)
      PBX=P1-ETA(JJ,1)
      PAY=P2-ETA(II,2)
      PBY=P2-ETA(JJ,2)
      PAZ=P3-ETA(II,3)
      PBZ=P3-ETA(JJ,3)
      T1=OVLAP(L1,L2,PAX,PBX,T)
      T2=OVLAP(M1,M2,PAY,PBY,T)
      T3=OVLAP(N1,N2,PAZ,PBZ,T)
      S1=OVLAP(L2+2,L1,PBX,PAX,T)
      S2=OVLAP (M2+2,M1,PBY,PAY,T)
      S3=OVLAP (N2+2,N1,PBZ,PAZ,T)
      PART=2*(L2+M2+N2)+3
      TTKE=B*(PART*T1*T2*T3-2.D0*B*(S1*T2*T3+T1*S2*T3+T1*T2*S3))
      IF (L2-1) 190,190,191
  191 PART=(L2*(L2-1))/2
      S1=OVLAP(L1,L2-2,PAX,PBX,T)
      TTKE=TTKE-PART*S1*T2*T3
  190 IF (M2-1) 192,192,193
  193 PART=(M2*(M2-1))/2
      S2=OVLAP(M1,M2-2,PAY,PBY,T)
      TTKE=TTKE-PART*T1*S2*T3
  192 IF (N2-1) 194,194,195
  195 PART=(N2*(N2-1))/2
      S3=OVLAP(N1,N2-2,PAZ,PBZ,T)
      TTKE=TTKE-PART*T1*T2*S3
  194 CONTINUE
 1635 TKE(NOM)=TKE(NOM)+SOO*TTKE
  635 CONTINUE
      PRVINT=TKE(NOM)
  916 CONTINUE
   25      CONTINUE
 915  CONTINUE
C=======================================================================
      IF(NPRALL.EQ.0)GO TO 777
       WRITE(6,7777)
 7777 FORMAT(//'&&&& TKE (TINTS) &&&&')
c      CALL PRINL(TKE,NINDEX,INDEXA)
  777 CONTINUE
C=======================================================================
C       DO 711 I=1,NBFN
C       DO 711 J=1,NSBF
C       IJ=NONE1(MAX0(I,J))+MIN0(I,J)
C       XYZ(I,J)=TKE(IJ)
C 711  CONTINUE
C       WRITE(6,722)
C 722   FORMAT(//'  TKE')
C      CALL PRINT(XYZ,40,40,NBFN,NSBF)
C      TRANSFORMATION FOR KINETIC ENERGY
C
C=======================================================================
      DO 520 I=1,NSBF
      DO 550 L=1,NSBF
      XL=0.0D0
      DO 560 K=1,NSBF
      KL=NONE1(MAX0(K,L))+MIN0(K,L)
  560      XL=XL+d(K,I)*TKE(KL)
  550      RI(L)=XL
      DO 580 J=1,I
      XIJ=0.0D0
      DO 590 L=1,NSBF
  590      XIJ=XIJ+d(L,J)*RI(L)
           TKET(I,J)=XIJ
  580      TKET(J,I)=XIJ
  520   CONTINUE
      IF(NSBF.EQ.NBFN) GO TO 521
      NSBF1=NSBF+1
      DO 522 J=NSBF1,NBFN
      DO 522 I=1,NSBF
      X=0.0D0
      DO 523 K=1,NSBF
      NJK=NONE1(J)+K
  523 X=X+TKE(NJK)*d(K,I)
      TKET(I,J)=X
  522 TKET(J,I)=X
  521 CONTINUE

c***
c  Deallocate those no longer needed
      deallocate( nskip, tke, list )
      deallocate( mtrans, nos, ri )
c***
cmod
c     open(unit=1,file='tket.lis',form='formatted')
c     do i=1,nbfns
c       do j= 1,nbfns
c         write(1,"(i4,i4,1x,1pe15.8)") i,j,tket(i,j)
c       end do
c     end do
c     close(1)
cmod

      IF(NREC40.EQ.0) GO TO 610
      DO 611 J=1,NSBF
  611 WRITE(9) (TKET(I,J),I=1,NBFNS)
  610      CONTINUE
C=======================================================================
      IF(NPRALL.EQ.0)GO TO 888
       WRITE(6,8888)
 8888 FORMAT(//'&&&& TKET (TINTS) &&&&')
C       CALL PRINT(TKET,NBFNX,NBFNX,NSBF,NSBF)
  888 CONTINUE
C=======================================================================
ccc   RETURN
ccc   END

      end subroutine tints
c----------------------------------------------------------------

c----------------------------------------------------------------
c     SUBROUTINE   MINTS (NLIST,MLIST,NCNTR,NTYPE,ETA,NFIRST,NLAST,S,
c    $                          NR,NBFXX,NTXX,NIMAX0,NGXX,
c    $                    NSAVXX,NEWLBL,NBFNS,NSBF,NMOF,NTRN,NREC30)
!     subroutine mints(NLIST,S,NIMAX0,NSAVXX,NEWLBL,NTRN,NREC30)
      subroutine mints(S,NSAVXX)

!$    use omp_lib, only: OMP_get_thread_num, OMP_get_num_threads,
!$   .                   OMP_get_wtime
      use ampoly, only: spones, spints, spdone, spdint, spdfnt

      implicit none

!     integer :: NLIST,NIMAX0,NSAVXX,NEWLBL,NTRN,NREC30
!     integer :: NSAVXX
!     real*8  :: S(NSAVXX)
      integer, intent(in) :: NSAVXX
      real*8, intent(in)  :: S(NSAVXX)

c***
c This common block is kept for use in ampoly
c
      integer :: ICNT,JCNT,KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,
     1LS,IF,JF,KF,LF,I,J,K,L
cccc  COMMON /SPECS/ ICNT,JCNT,KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,
cccc $LS,IF,JF,KF,LF,NINT,IY,JY,KY,LY
c***

c***
c  Local variables
c
      logical :: single_center
!     integer :: ij, in, iy, jy, ky, ly, imn, ij2, lmx, igg
      integer :: ij, in, iy, jy, ky, ly, imn, ij2, lmx
      integer :: ii, iii, jj, ll, kk, kl0, ij0
!     integer :: kl, jl, intap, nbfn2, jl_d, kl_d
      integer :: kl, nbfn2, jl_d, kl_d
      integer :: NINT, NBFN, NINMAX
!     integer :: len_b, nom_new, x_block, n_block, i_block = 1
      integer :: len_b, nom_new, n_block, i_block = 1
      integer :: n_threads

!     integer(i_large) :: indexa, ip, nom, ijkl0, ijkl
      integer(i_large) :: indexa, nom, ijkl0, ijkl
      integer(i_large), allocatable :: n_posit(:)

      real*8 :: x, t_start, t_end
      real*8,  allocatable :: rsl(:), rin(:,:), rin2(:), 
     .                        rin3(:), rrin(:)
      integer :: alloc_stat, stat
c***

C=======================================================================
C
C      ERP(NBFN,NSBF/NSBF,NSBF)
C
C=======================================================================
      NINT=1
      NINMAX=1
      NBFN=NBFNS
      IN=NONE1(NBFN)+NSBF
      INDEXA=NTWO(IN)+NONE1(NSBF)+NSBF
      rpmolx = indexa
      nbfn2 = nbfn * (nbfn+1) / 2
      n_threads = 0

c***
c  Allocate relevant variables: private arrays will be allocated within
c  the parallel region
c
      allocate( erp(indexa), stat=alloc_stat )
      call allocation_stat('erp',alloc_stat)
c***

      allocate( list(1,4), rsl(1) )

      DO 538 I=1,NBFN
      LIST(1,1)=I
      IMN=MIN0(I,NSBF)
      DO 536 J=1,IMN
      LIST(1,2)=J
      IJ2=NONE1(I)+J
      DO 534 K=1,IMN
      LIST(1,3)=K
      IF(I-K) 512,511,512
  511 LMX=J
      GO TO 513
  512 LMX=K
  513 DO 532 L=1,LMX
      LIST(1,4)=L

  514 nom=NTWO(IJ2)+NONE1(K)+L
      IY=LIST(1,1)
      JY=LIST(1,2)
      KY=LIST(1,3)
      LY=LIST(1,4)
   60 RSL(1)=0.0D0

 120    IF  ( NTYPE(JY) - NTYPE(IY) )   140,140,130
  130 III=IY
          IY =  JY
      JY=III
 140    IF  ( NTYPE(LY) - NTYPE(KY) )   160,160,150
  150 III=KY
          KY =  LY
      LY=III
 160  IF  ( NTYPE(KY) - NTYPE(IY) )   180,165,170
 165  IF  ( NTYPE(LY) - NTYPE(JY) )   180,180,170
  170 III=IY
          IY =  KY
      KY=III
      III=JY
          JY =  LY
      LY=III
  180 ICNT=NCNTR(IY)
      JCNT=NCNTR(JY)
      KCNT=NCNTR(KY)
      LCNT=NCNTR(LY)
      ITYP=NTYPE(IY)
      JTYP=NTYPE(JY)
      KTYP=NTYPE(KY)
      LTYP=NTYPE(LY)
      IS=NFIRST(IY)
      JS=NFIRST(JY)
      KS=NFIRST(KY)
      LS=NFIRST(LY)
      IF=NLAST(IY)
      JF=NLAST(JY)
      KF=NLAST(KY)
      LF=NLAST(LY)

      single_center = ( (ICNT.eq.JCNT) .and. (KCNT.eq.LCNT) 
     .                                 .and. (ICNT.eq.KCNT) )

      If( ITYP .le. 4 ) Then

        if( single_center ) then 
          call spones(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .                KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .                KF,LF,NINT,IY,JY,KY,LY)
        else
          call spints(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .                KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .                KF,LF,NINT,IY,JY,KY,LY)
        end if

      Else If( ITYP .le. 10 ) Then

        if( single_center ) then
          call spdone(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .                KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .                KF,LF,NINT,IY,JY,KY,LY)
        else
          call spdint(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .                KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .                KF,LF,NINT,IY,JY,KY,LY)
        end if
      
      Else If( ITYP .gt. 10 ) Then

        if( single_center ) then
          call spdone(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .                KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .                KF,LF,NINT,IY,JY,KY,LY)
        else
          call spdfnt(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .                KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .                KF,LF,NINT,IY,JY,KY,LY)
        end if

      End If

      ERP(NOM)=RSL(1)

  532      CONTINUE
  534      CONTINUE
  536      CONTINUE
  538      CONTINUE

cmod
c     open(unit=99,file='erp.lis',form='formatted')
c     do i=1,indexa
c       write(99,"(i8,2x,1pe15.8)") i,erp(i)
c     end do
c     close(99)
cmod


c
ccc  Deallocate integral arrays and allocate transformation arrays
c
      deallocate( list, rsl )

      if( ntwoel .eq. 3 ) return

      If( ntwoel .eq. 0 ) Then ! parallel transformation

      allocate( erpt(indexa), stat=alloc_stat )
      call allocation_stat('erpt',alloc_stat)

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(erp, erpt, d, nsbf, 
!$OMP.   none1, ntwo, nbfn, nbfn2)

      allocate( rin(nbfn,nbfn2), stat=alloc_stat )
      call allocation_stat('rin',alloc_stat)
      allocate( rin2(nbfn2), rin3(nbfn), stat=alloc_stat )
      call allocation_stat('rin2, rin3',alloc_stat)


!$OMP MASTER
!$      n_threads = OMP_get_num_threads()
!$      write(6,"(/,10x,'START PARALLEL AO-TO-MO BASIS',
!$   .                                ' TRANSFORMATION')")
!$      write(6,"(10x,'USING ',i4,' THREADS')") n_threads
!$      write(6,"(10x,'MEMORY OVERHEAD = ',1pe13.6,' GB')") 
!$   .      dble((nbfn2*nbfn+nbfn2+nbfn)*(n_threads-1))*8.0/1024.0**3 
!$    t_start  = OMP_get_wtime()
!$OMP END MASTER


!$OMP DO
        DO 800 I=1,NSBF
ccc
        jl_d = 0
ccc
        DO 801 LL=1,NSBF
        DO 801 KK=1,LL
ccc
        jl_d = jl_d + 1
ccc
        KL0=NONE1(LL)+KK
        DO 801 JJ=1,NSBF
        X=0.0D0
        DO 802 II=1,NSBF
        IJ0=NONE1(MAX0(II,JJ))+MIN0(II,JJ)
        IJKL0=NTWO(MAX0(IJ0,KL0))+MIN0(IJ0,KL0)
  802   X=X+ERP(IJKL0)*d(II,I)

c***
cc      RIN(JJ,LL,KK)=X
cc801   RIN(JJ,KK,LL)=X
c
        rin(jj,jl_d)=x
  801   continue
c***

        DO 810 J=1,I
        IJ=NONE1(I)+J
        DO 811 LL=1,NSBF
        DO 811 KK=1,LL
        X=0.0D0

c***
cc      DO 812 JJ=1,NSBF
cc812   X=X+RIN(JJ,KK,LL)*d(JJ,J)
cc      RIN2(LL,KK)=X
cc811   RIN2(KK,LL)=X
c
        jl_d = n_get(kk,ll)
        do 812 jj=1,nsbf
        x=x+rin(jj,jl_d)*d(jj,j)
  812   continue
        rin2(jl_d)=x
  811   continue
c***

        DO 820 K=1,I
        LMX=K
        IF(I.EQ.K)  LMX=J
  853   DO 821 LL=1,NSBF
        X=0.0D0

c***
cc      DO 822 KK=1,NSBF
cc822   X=X+RIN2(KK,LL)*d(KK,K)
        do 822 kk=1,nsbf
        jl_d=n_get(kk,ll)
        x=x+rin2(jl_d)*d(kk,k)
  822   continue
c***

  821   RIN3(LL)=X
        DO 830 L=1,LMX
        KL=NONE1(K)+L
        IJKL=NTWO(IJ)+KL
        X=0.0D0
        DO 831 LL=1,NSBF
  831   X=X+RIN3(LL)*d(LL,L)
        ERPT(IJKL)=X
  830   CONTINUE
  820   CONTINUE
  810   CONTINUE
  800   CONTINUE
!$OMP END DO

c***
c  Deallocate transformation arrays
c
      deallocate( rin, rin2, rin3 )
c***

!$OMP MASTER
!$    t_end  = OMP_get_wtime()
!$      write(6,"(10x,'END OF MO CALCULATIONS. TIME SPENT(SECS) = ',
!$   .        1pe13.6  )") t_end - t_start
!$OMP END MASTER

!$OMP END PARALLEL

      Else If( ntwoel. eq. 1 ) Then

      allocate( erpt(indexa), stat=alloc_stat )
      call allocation_stat('erpt',alloc_stat)
      allocate( rin(nbfn,nbfn2), stat=alloc_stat )
      call allocation_stat('rin',alloc_stat)
      allocate( rin2(nbfn2), rin3(nbfn), stat=alloc_stat )
      call allocation_stat('rin2, rin3',alloc_stat)


      write(6,"(/,10x,'START SERIAL AO-TO-MO BASIS TRANSFORMATION')")

        DO 700 I=1,NSBF
ccc
        jl_d = 0
ccc
        DO 701 LL=1,NSBF
        DO 701 KK=1,LL
ccc
        jl_d = jl_d + 1
ccc
        KL0=NONE1(LL)+KK
        DO 701 JJ=1,NSBF
        X=0.0D0
        DO 702 II=1,NSBF
        IJ0=NONE1(MAX0(II,JJ))+MIN0(II,JJ)
        IJKL0=NTWO(MAX0(IJ0,KL0))+MIN0(IJ0,KL0)
  702   X=X+ERP(IJKL0)*d(II,I)

c***
cc      RIN(JJ,LL,KK)=X
cc801   RIN(JJ,KK,LL)=X
c
        rin(jj,jl_d)=x
  701   continue
c***

        DO 710 J=1,I
        IJ=NONE1(I)+J
        DO 711 LL=1,NSBF
        DO 711 KK=1,LL
        X=0.0D0

c***
cc      DO 812 JJ=1,NSBF
cc812   X=X+RIN(JJ,KK,LL)*d(JJ,J)
cc      RIN2(LL,KK)=X
cc811   RIN2(KK,LL)=X
c
        jl_d = n_get(kk,ll)
        do 712 jj=1,nsbf
        x=x+rin(jj,jl_d)*d(jj,j)
  712   continue
        rin2(jl_d)=x
  711   continue
c***

        DO 720 K=1,I
        LMX=K
        IF(I.EQ.K)  LMX=J
  753   DO 721 LL=1,NSBF
        X=0.0D0

c***
cc      DO 822 KK=1,NSBF
cc822   X=X+RIN2(KK,LL)*d(KK,K)
        do 722 kk=1,nsbf
        jl_d=n_get(kk,ll)
        x=x+rin2(jl_d)*d(kk,k)
  722   continue
c***

  721   RIN3(LL)=X
        DO 730 L=1,LMX
        KL=NONE1(K)+L
        IJKL=NTWO(IJ)+KL
        X=0.0D0
        DO 731 LL=1,NSBF
  731   X=X+RIN3(LL)*d(LL,L)
        ERPT(IJKL)=X
  730   CONTINUE
  720   CONTINUE
  710   CONTINUE
  700   CONTINUE

c***
c  Deallocate transformation arrays
c
      deallocate( rin, rin2, rin3 )
c***

      Else If( ntwoel. eq. 2 ) Then

      write(6,"(/,10x,'START SERIAL AO-TO-MO BASIS TRANSFORMATION')")
      write(6,"(10x,  'SAVE MEMORY')")

      len_b  = nbfn2
cc LEN_B is the block length, N_BLOCK is the number of blocks
cc and N_POSIT(:) is the first element of each block
      n_block = indexa / len_b
      if(mod(indexa,len_b) .ne. 0) n_block = n_block + 1
      allocate( n_posit(n_block) )
      n_posit(i_block) = 1

cc Open 2-electron file
      open(unit=n_erpmo,file='ERPT.SCR',form='unformatted',iostat=stat)
      call open_stat('ERPT.SCR',stat)

cc Allocate the relevant arrays
      allocate( rrin(nbfn), stat=alloc_stat )
      call allocation_stat('rin',alloc_stat)
      allocate( rin2(nbfn2), rin3(nbfn), stat=alloc_stat )
      call allocation_stat('rin2, rin3',alloc_stat)
      allocate( erpt(len_b), stat=alloc_stat )
      call allocation_stat('erpt',alloc_stat)

ccc Shell I
      DO 600 I=1,NSBF
        open(unit=n_scrat, file='TWOEL.SCR', form='unformatted',
     .       iostat=stat)
        call open_stat('ERPT.SCR',stat)
        kl_d = 0
        DO 601 LL=1,NSBF
        DO 601 KK=1,LL
        kl_d = kl_d + 1
        KL0=NONE1(LL)+KK
        DO 602 JJ=1,NSBF
        X=0.0D0
        DO 603 II=1,NSBF
        IJ0=NONE1(MAX0(II,JJ))+MIN0(II,JJ)
        IJKL0=NTWO(MAX0(IJ0,KL0))+MIN0(IJ0,KL0)
  603   X=X+ERP(IJKL0)*d(II,I)

  602   rrin(jj)=x
        call wr_rin(rrin,nbfns,kl_d)
  601   continue

ccc Shell J
        DO 610 J=1,I
          rewind(n_scrat)
          IJ=NONE1(I)+J
          kl_d = 0
          DO 611 LL=1,NSBF
          DO 611 KK=1,LL
          X=0.0D0
          kl_d = kl_d + 1
          call rd_rin(rrin,nbfns,kl_d)
          do 612 jj=1,nsbf
          x=x+rrin(jj)*d(jj,j)
  612     continue
          rin2(kl_d)=x
  611     continue

ccc Shell K
          DO 620 K=1,I
            LMX=K
            IF(I.EQ.K)  LMX=J
  653       DO 621 LL=1,NSBF
            X=0.0D0
            do 622 kk=1,nsbf
            kl_d=n_get(kk,ll)
            x=x+rin2(kl_d)*d(kk,k)
  622       continue

  621       RIN3(LL)=X

ccc Shell M
            DO 630 L=1,LMX
              KL=NONE1(K)+L

c***
ccc           IJKL=NTWO(IJ)+KL
              nom_new = NTWO(IJ)+KL
              ijkl = mod(nom_new,len_b)
              if(ijkl .eq. 0) ijkl= len_b
c***

              X=0.0D0
              DO 631 LL=1,NSBF
  631         X=X+RIN3(LL)*d(LL,L)
              ERPT(IJKL)=X

c***
              if( (mod(nom_new,len_b).eq.0) .or.
     .            (nom_new.eq.indexa) ) then
                call wr_erp(erpt, len_b, n_posit(i_block))
                if(i_block .lt. n_block) then
                  i_block = i_block + 1
                  n_posit(i_block) = nom_new + 1
                end if
              end if
c***

  630       CONTINUE
  620     CONTINUE
  610   CONTINUE
        close(n_scrat)
  600 CONTINUE

c***
c  Deallocate transformation arrays and AO integrals and load full ERPT
c
      deallocate( rrin, rin2, rin3 )
      deallocate( erp )

      deallocate( erpt )
      allocate( erpt(indexa), stat=alloc_stat )
      call allocation_stat('erpt',alloc_stat)
      do i = 1,n_block
        if( (i.eq.n_block) .and. (mod(indexa,len_b).ne.0) )
     .    len_b = mod(indexa,len_b)
        call rd_erp(erpt(n_posit(i)), len_b, i, n_posit(i))
      end do
c***

      deallocate( n_posit )

      End If

cmod
c     open(unit=1,file='erpt.lis',form='formatted')
c     do i=1,indexa
c       write(1,"(i4,2x,1pe15.8)") i,erpt(i)
c     end do
c     close(1)
cmod

c***
c  Deallocate transformation AO integrals
c
      if( allocated(erp) ) deallocate( erp )
c***

      end subroutine mints
c----------------------------------------------------------------

c----------------------------------------------------------------
      subroutine wr_erp(x, len, n1)

      implicit none

      integer :: len
      integer(i_large) :: n1
      real*8  :: x(len)

      write(n_erpmo) n1, x
      
      end subroutine wr_erp
c----------------------------------------------------------------

c----------------------------------------------------------------
      subroutine rd_erp(x, len, ib, n1)

      implicit none

      integer :: len, ib
      integer(i_large) :: n1
      real*8  :: x(len)

c***
c  Local variables
c
      integer :: i
      integer(i_large) :: nb
c***
      rewind(n_erpmo)
      do i = 1,ib
        read(n_erpmo) nb, x
      end do

      if(nb .ne. n1) then
        write(6,"('mismatch in RD_ERP')") 
        stop'ABORT...'
      end if
      
      end subroutine rd_erp
c----------------------------------------------------------------

c----------------------------------------------------------------
      subroutine wr_rin(x,len,i)

      implicit none

      integer :: len ,i
      real*8  :: x(len)

      write(n_scrat) i, x

      end subroutine wr_rin
c----------------------------------------------------------------

c----------------------------------------------------------------
      subroutine rd_rin(x,len,i)

      implicit none

      integer :: len ,i
      real*8  :: x(len)

c***
c  Local variable
c
      integer :: n
c***

      read(n_scrat) n, x

      if(n .ne. i) then
        write(6,"('mismatch in RD_RIN')")
        stop 'ABORT'
      end if

      end subroutine rd_rin
c----------------------------------------------------------------

c----------------------------------------------------------------
      function erpt_fun(i,j,k,l)

      implicit none

      integer, intent(in) :: i, j, k, l

c***
c  Local variables
c
!     logical, save :: first_call = .true. 
      integer :: nbfn2, jl_d, ii, jj, ll, kk
      integer :: mi, mj, mk, ml, mm1, mm2
      integer :: ij0, kl0
!     integer(i_large) :: ijkl0, indexa
      integer(i_large) :: ijkl0
      real*8  :: x, erpt_fun
      real*8,  allocatable :: rin(:,:), rin2(:), rin3(:)
c***

      nbfn2 = nbfns * (nbfns+1) / 2
      allocate( rin(nbfns,nbfn2) )
      allocate( rin2(nbfn2), rin3(nbfns) )

      mi = max(i,j)
      mj = min(i,j)
      mk = max(k,l)
      ml = min(k,l)
      mm1 = none1(mi)+mj
      mm2 = none1(mk)+ml
      if(mm2 .gt. mm1) then
        mm1 = mi
        mi  = mk
        mk  = mm1
        mm2 = mj
        mj  = ml
        ml  = mm2
      end if

cmod
c     print*,'mi,mj,mk,ml', mi,mj,mk,ml
c     print*,'allocated(erp)', allocated(erp)
c     print*,'size(erp)', size(erp)
c     print*,'allocated(none1)', allocated(none1)
c     print*,'allocated(ntwo)', allocated(ntwo)
c     print*,'allocated(d)', allocated(d)
c     stop
cmod

!$OMP   PARALLEL DEFAULT(PRIVATE) SHARED(nsbf, none1,
!$OMP.  ntwo, mi, mj, mk, ml, d, rin, rin2, rin3, erp)

      jl_d = 0
!$OMP DO
      DO 801 LL=1,NSBF
      DO 801 KK=1,LL
        jl_d = jl_d + 1
        KL0=NONE1(LL)+KK
        DO 801 JJ=1,NSBF
        X=0.0D0
        DO 802 II=1,NSBF
        IJ0=NONE1(MAX0(II,JJ))+MIN0(II,JJ)
        IJKL0=NTWO(MAX0(IJ0,KL0))+MIN0(IJ0,KL0)
  802   X=X+erp(ijkl0)*d(II,mi)
        rin(jj,jl_d)=x
  801 continue
!$OMP END DO

cmod
c     print*,'passou loop1'
c     stop
cmod

!$OMP DO
      DO 811 LL=1,NSBF
      DO 811 KK=1,LL
        X=0.0D0  
        jl_d = n_get(kk,ll)
        do 812 jj=1,nsbf
          x=x+rin(jj,jl_d)*d(jj,mj)
  812   continue
        rin2(jl_d)=x
  811 continue
!$OMP END DO

!$OMP DO
      DO 821 LL=1,NSBF
        X=0.0D0
        do 822 kk=1,nsbf
          jl_d=n_get(kk,ll)
          x=x+rin2(jl_d)*d(kk,mk)
  822   continue
  821 RIN3(LL)=X
!$OMP END DO
!$OMP END PARALLEL

      X=0.0D0
      DO 831 LL=1,NSBF
  831   X=X+RIN3(LL)*d(LL,ml)

      erpt_fun=X

      deallocate( rin, rin2, rin3 )

      end function erpt_fun
c----------------------------------------------------------------

c----------------------------------------------------------------
      function erp_fun(mi,mj,mk,ml)

      use ampoly, only: spones, spints, spdone, spdint, spdfnt

      implicit none

!     integer :: mi, mj, mk, ml
      integer, intent(in) :: mi, mj, mk, ml

c***
c This common block is kept for use in ampoly
c
      integer :: ICNT,JCNT,KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,
     1LS,IF,JF,KF,LF,I,J,K,L
cccc  COMMON /SPECS/ ICNT,JCNT,KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,
cccc $LS,IF,JF,KF,LF,NINT,IY,JY,KY,LY
c***

c***
c  Local variables
c
      integer :: mm1, mm2, iy, jy, ky, ly, igg, iii
      integer :: nint, jcan, ninmax, nsavxx
      real*8  :: prvint, rsl(1)
      real*8  :: erp_fun
c***

      NINT=1
      NINMAX=1
      JCAN=ICANON-1
      nsavxx = nsavm

      i = max(mi,mj)
      j = min(mi,mj)
      k = max(mk,ml)
      l = min(mk,ml)

      mm1 = none1(i)+j
      mm2 = none1(k)+l
      if(mm2 .gt. mm1) then
        mm1 = i
        i   = k
        k   = mm1
        mm2 = j
        j   = l
        l   = mm2
      end if

      IY = mi
      JY = mj
      KY = mk
      LY = ml
      IGG=0

55    IF(IGG-1) 60,56,58
   56      erp_fun=PRVINT
      GO TO 411
   58      erp_fun=-PRVINT
      GO TO 411
   60      RSL(1)=0.0D0
      IF(JCAN) 120,120,180
 120    IF  ( NTYPE(JY) - NTYPE(IY) )   140,140,130
  130 III=IY
          IY =  JY
      JY=III
 140    IF  ( NTYPE(LY) - NTYPE(KY) )   160,160,150
  150 III=KY
          KY =  LY
      LY=III
 160  IF  ( NTYPE(KY) - NTYPE(IY) )   180,165,170
 165  IF  ( NTYPE(LY) - NTYPE(JY) )   180,180,170
  170 III=IY
          IY =  KY
      KY=III
      III=JY
          JY =  LY
      LY=III
  180 ICNT=NCNTR(IY)
      JCNT=NCNTR(JY)
      KCNT=NCNTR(KY)
      LCNT=NCNTR(LY)
      ITYP=NTYPE(IY)
      JTYP=NTYPE(JY)
      KTYP=NTYPE(KY)
      LTYP=NTYPE(LY)
      IS=NFIRST(IY)
      JS=NFIRST(JY)
      KS=NFIRST(KY)
      LS=NFIRST(LY)
      IF=NLAST(IY)
      JF=NLAST(JY)
      KF=NLAST(KY)
      LF=NLAST(LY)
      IF  ( ITYP - 4 )   320,320,321
  320 IF (ICNT-JCNT) 350,351,350
  351 IF (KCNT-LCNT) 350,352,350
  352 IF (ICNT-KCNT) 350,353,350

c***
ccc 353  CALL    SPONES  (ETA,   RSL,S,NR,NTXX,NGXX,NINMAX,NSAVXX)
ccc 353  call spones(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX)
 353  call spones(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .            KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .            KF,LF,NINT,IY,JY,KY,LY)
c***

      erp_fun=RSL(1)
      GO TO 250

c***
ccc350  CALL    SPINTS  (ETA,   RSL,S,NR,NTXX,NGXX,NINMAX,NSAVXX)
c350  call spints(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX)
 350  call spints(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .            KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .            KF,LF,NINT,IY,JY,KY,LY)
c***

      erp_fun=RSL(1)
      GO TO 250
 321  IF  ( ITYP - 10 )   420,420,421
  421 IF (ICNT-JCNT) 450,451,450
  451 IF (KCNT-LCNT) 450,452,450
  452 IF (ICNT-KCNT) 450,453,450

c***
c453  CALL    SPDONE  (ETA,   RSL,S,NR,NTXX,NGXX,NINMAX,NSAVXX)
c453  call spdone(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX)
 453  call spdone(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .            KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .            KF,LF,NINT,IY,JY,KY,LY)
c***

      erp_fun=RSL(1)
      GO TO 250

c***
c450  CALL    SPDFNT  (ETA,   RSL,S,NR,NTXX,NGXX,NINMAX,NSAVXX)
c450  call spdfnt(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX)
 450  call spdfnt(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .            KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .            KF,LF,NINT,IY,JY,KY,LY)
c***

      erp_fun=RSL(1)
      GO TO 250
  420 IF (ICNT-JCNT) 360,361,360
  361 IF (KCNT-LCNT) 360,362,360
  362 IF (ICNT-KCNT) 360,363,360

c***
c363  CALL    SPDONE  (ETA,   RSL,S,NR,NTXX,NGXX,NINMAX,NSAVXX)
c363  call spdone(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX)
 363  call spdone(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .            KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .            KF,LF,NINT,IY,JY,KY,LY)
c***

      erp_fun=RSL(1)
      GO TO 250

c***
c360  CALL    SPDINT  (ETA,   RSL,S,NR,NTXX,NGXX,NINMAX,NSAVXX)
c360  call spdint(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX)
 360  call spdint(ETA,RSL,S,NR,ntmx,ngaus,NINMAX,NSAVXX,ICNT,JCNT,
     .            KCNT,LCNT,ITYP,JTYP,KTYP,LTYP,IS,JS,KS,LS,IF,JF,
     .            KF,LF,NINT,IY,JY,KY,LY)
c***

      erp_fun=RSL(1)
  250   PRVINT =    erp_fun 
  411      CONTINUE

      end function erp_fun
c----------------------------------------------------------------

c----------------------------------------------------------------
      integer function n_get(ii,jj)

      implicit none

      integer, intent(in) :: ii, jj

c***
c  Local variables
c
      integer :: i, j
c***

ccc Assumes (i .ge. j), i.e., do i=1,n / do j=1,i
      i = max(ii,jj)
      j = min(ii,jj)

      n_get = i*(i-1)/2 + j

      end function
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE OVLGRN(NBFN,NSBF,NMOF,OV,UNIT,NPRALL)
c     IMPLICIT REAL*8(A-H,O-Z)
C=======================================================================
c     INCLUDE 'STRUC.BFN'
C=======================================================================
c     COMMON /TRANSF/D(NSBFX,NSBFX)
c     COMMON /OVRPG/OVPM(NBFNX,NBFNX)
C=======================================================================
c     DIMENSION OV(NBFN,NBFN),UNIT(NBFN,NBFN)

c***
      implicit none

      integer :: NBFN,NSBF,NMOF,NPRALL
      real*8  :: OV(NBFN,NBFN),UNIT(NBFN,NBFN)

      integer :: i, j
      integer :: ier, nsbf1
      real*8  :: eps
c***

C=======================================================================
C
C      A MATRIX WHICH TRANSFORMS GAUSSIAN BASIS GREEN'S FUNCTION TO
C      ITS MOLECULAR ORBITAL REPRESENTATION IS GENERATED.
C
C=======================================================================
      IF(NBFN.EQ.NSBF.AND.NSBF.EQ.NMOF) GO TO 500
      DO 100 I=1,NBFN
      DO 100 J=1,NBFN
  100      OV(J,I)=OVPM(J,I)
      DO 300 I=1,NSBF
      DO 300 J=1,NSBF
  300      UNIT(J,I)=D(I,J)
      IF(NSBF.EQ.NBFN) GO TO 400
      NSBF1=NSBF+1
      DO 350  I=1,NBFN
      DO 350  J=NSBF1,NBFN
      UNIT(J,I)=0.0D0
  350      UNIT(I,J)=0.0D0
      DO 360 I=NSBF1,NBFN
  360      UNIT(I,I)=1.0D0
  400      CONTINUE
      EPS=1.0D-15
      CALL GELG(UNIT,OV,NBFN,NBFN,EPS,IER)
      WRITE(6,200) IER
  200 FORMAT(//'   IER(OVLGRN)=',I4)
        DO 250 I=1,NBFN
      DO 250 J=1,NBFN
  250      OVPM(J,I)=UNIT(I,J)
C=======================================================================
      IF(NPRALL.EQ.0)GO TO 777
      WRITE(6,600)
  600 FORMAT(//'*** TRANSFORMATION MATRIX FOR THE GREEN FUNCTION ***
     1(OVLGRN)')
C      CALL PRINT(OVPM,NBFNX,NBFNX,NBFN,NBFN)
  777 CONTINUE
C=======================================================================
      RETURN
  500   DO 550 I=1,NMOF
      DO 550 J=1,NMOF
  550      OVPM(J,I)=D(J,I)

      end subroutine OVLGRN
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE GELG(R,A,M,N,EPS,IER)

c***
c     IMPLICIT REAL*8(A-H,O-Z)
c     DIMENSION A(1),R(1)

      implicit none

      integer :: M,N,IER
      real*8  :: R(1),A(1),EPS

ccc Local variables

      integer :: mm, nm, l, lst, k, i, j, ll
      integer :: lend, ii, ist
      real*8  :: piv, tb, tol, pivi
c***

      IF(M) 23,23,1
    1      IER=0
      PIV=0.0D0
      MM=M*M
      NM=N*M
      DO 3 L=1,MM
      TB=DABS(A(L))
      IF(TB-PIV) 3,3,2
    2      PIV=TB
      I=L
    3   CONTINUE
      TOL=EPS*PIV
      LST=1
      DO 17 K=1,M
      IF(PIV) 23,23,4
    4      IF(IER) 7,5,7
    5      IF(PIV-TOL) 6,6,7
    6       IER=K-1
    7      PIVI=1.0D0/A(I)
      J=(I-1)/M
      I=I-J*M-K
      J=J+1-K
      DO 8 L=K,NM,M
      LL=L+I
      TB=PIVI*R(LL)
      R(LL)=R(L)
    8      R(L)=TB
      IF(K-M) 9,18,18
    9      LEND=LST+M-K
      IF(J) 12,12,10
   10      II=J*M
      DO 11 L=LST,LEND
      TB=A(L)
      LL=L+II
      A(L)=A(LL)
   11      A(LL)=TB
   12       DO 13 L=LST,MM,M
      LL=L+I
      TB=PIVI*A(LL)
      A(LL)=A(L)
   13      A(L)=TB
      A(LST)=J
      PIV=0.0D0
      LST=LST+1
      J=0
      DO 16 II=LST,LEND
      PIVI=-A(II)
      IST=II+M
      J=J+1
      DO 15 L=IST,MM,M
      LL=L-J
      A(L)=A(L)+PIVI*A(LL)
      TB=DABS(A(L))
      IF(TB-PIV) 15,15,14
   14      PIV=TB
      I=L
   15 CONTINUE
      DO 16 L=K,NM,M
      LL=L+J
   16      R(LL)=R(LL)+PIVI*R(L)
   17      LST=LST+M
   18      IF(M-1) 23,22,19
   19      IST=MM+M
      LST=M+1
      DO 21 I=2,M
      II=LST-I
      IST=IST-LST
      L=IST-M
      L=A(L)+0.5D0
      DO 21 J=II,NM,M
      TB=R(J)
      LL=J
      DO 20 K=IST,MM,M
      LL=LL+1
   20      TB=TB-A(K)*R(LL)
      K=J+L
      R(J)=R(K)
   21      R(K)=TB
   22      RETURN
   23      IER=-1

      end subroutine GELG
c----------------------------------------------------------------

c----------------------------------------------------------------
c     SUBROUTINE INLAB(NBFNS,ISYM)
      subroutine inlab(isym)

      implicit none

      integer, intent(in) :: isym

CR
      READ  (5,913)  JCON
CW
      WRITE(6,600)
CW
      WRITE(6,620) JCON
      IF(ISYM.EQ.0) RETURN

      write(6,"('ERROR IN SUBROUTINE INLAB')")
      write(6,"('(NSBF.GT.NMOF) NOT ALLOWED')")
      stop'ABORT...'

  600 FORMAT(1H1//5X,26HLABEL GENERATION PROGRAM           )
  620 FORMAT(/5X,8HOPTIONS ,10I3)
 913  FORMAT ( 10I5 )

      end subroutine inlab
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE NUMR(NA,NB,NC)

      implicit none

!     logical :: skip = .false.
      integer, intent(in) :: NA, NB, NC
      integer :: I, II

      integer(i_large) :: ibig, init, last

c***
c Allocate arrays for numbering 1-e and 2-e integrals 
c
      allocate( none1(na), none2(na), none3(na) )
!fk   skip = ((ntwoel.eq.3) .and. (nrfock.ne.0))
!     if( .not. skip ) 
!    .   allocate( ntwo(nbfns*(nbfns+1)/2 + (nbfns-nbfns)*nbfns) )
         allocate( ntwo(nbfns*(nbfns+1)/2 + (nbfns-nbfns)*nbfns) )
c***

C=======================================================================
C
C      NUMBERING FOR
C      (NA,NB)=NONE1      NA.GE.NB
C      (NB,NC)=NONE2      NB.GE.NC
C      (NA,NA)=NONE3
C      (NA,NB/NB,NC)=NTWO
C
C=======================================================================
      DO 100 I=1,NC
      II=I*(I-1)/2
      NONE1(I)=II
      NONE2(I)=II
  100      NONE3(I)=II
      IF(NC.EQ.NA) GO TO 330
      DO 120 I=NC+1,NA
  120      NONE3(I)=I*(I-1)/2
  330      CONTINUE
      IF(NC.EQ.NB) GO TO 300
      DO 150 I=NC+1,NB
      NONE1(I)=I*(I-1)/2
  150      NONE2(I)=NONE2(I-1)+NC
  300      IF(NB.EQ.NA) GO TO 400
      DO 200 I=NB+1,NA
  200      NONE1(I)=NONE1(I-1)+NB
  400      LAST=NONE1(NA)+NB
      INIT=NONE1(NB)+NB

!     if( skip ) return

      DO 500 ibig=1,INIT
  500   NTWO(ibig)=ibig*(ibig-1)/2
      IF(INIT.EQ.LAST) RETURN
      DO 550 ibig=INIT+1,LAST
  550      NTWO(ibig)=NTWO(ibig-1)+INIT

      end subroutine numr
c----------------------------------------------------------------

c----------------------------------------------------------------
      subroutine d_shells(nd_shell,jx)

      implicit none

      integer, intent(inout) :: nd_shell, jx(nbfns) 

c***
cc Local variables
c
      integer :: i, l, m, n, ityp
      logical :: dxx_typ, dyy_typ, dzz_typ
      logical :: dxy_typ, dxz_typ, dyz_typ
c***

ccc  Only d shells from molecular orbitals are taken into account
ccc  (i.e., additional scattering orbitals not taken into account)

!     If(orbfmt .eq. " GAMESS") Then
      If(orbfmt .eq. "GAMESS") Then
      do i = 1,nbfns
        ityp = ntype(i)
        l = nr(ityp,1)
        m = nr(ityp,2)
        n = nr(ityp,3)

        dxx_typ = ( (l.eq.2) .and. (m.eq.0)
     .                       .and. (n.eq.0) )
        if( dxx_typ ) then
          if(i .gt. (nbfns-5) ) then
            write(6,103) i
            stop'ABORT...'
          end if
          ityp = ntype(i+1)
          l = nr(ityp,1)
          m = nr(ityp,2)
          n = nr(ityp,3)
          dyy_typ = ( (l.eq.0) .and. (m.eq.2)
     .               .and. (n.eq.0) )
          ityp = ntype(i+2)
          l = nr(ityp,1)
          m = nr(ityp,2)
          n = nr(ityp,3)
          dzz_typ = ( (l.eq.0) .and. (m.eq.0)
     .               .and. (n.eq.2) )
          ityp = ntype(i+3)
          l = nr(ityp,1)
          m = nr(ityp,2)
          n = nr(ityp,3)
          dxy_typ = ( (l.eq.1) .and. (m.eq.1)
     .               .and. (n.eq.0) )
          ityp = ntype(i+4)
          l = nr(ityp,1)
          m = nr(ityp,2)
          n = nr(ityp,3)
          dxz_typ = ( (l.eq.1) .and. (m.eq.0)
     .               .and. (n.eq.1) )
          ityp = ntype(i+5)
          l = nr(ityp,1)
          m = nr(ityp,2)
          n = nr(ityp,3)
          dyz_typ = ( (l.eq.0) .and. (m.eq.1)
     .               .and. (n.eq.1) )

          if(dxx_typ .and. dyy_typ .and. dzz_typ .and.
     .       dxy_typ .and. dxz_typ .and. dyz_typ ) then 
            nd_shell = nd_shell + 1 
            if(nd_shell .gt. nbfns) then
              write(6,"('TOO MANY D SHELLS')")
              stop'ABORT...'
            end if
            jx(nd_shell) = i
          else
            write(6,102)
            stop
          end if
        end if
      end do
      Else
      do i = 1,nmof
        ityp = ntype(i)
        l = nr(ityp,1)
        m = nr(ityp,2)
        n = nr(ityp,3)

        dxx_typ = ( (l.eq.2) .and. (m.eq.0)
     .                       .and. (n.eq.0) )
        if( dxx_typ ) then
            nd_shell = nd_shell + 1 
            if(nd_shell .gt. nbfns) then
              write(6,"('TOO MANY D SHELLS')")
              stop'ABORT...'
            end if
            jx(nd_shell) = i
        end if
      end do
      End If

 102  format(/,20x,'**** ERROR ****',/,'D-SHELL NOT IN',
     . ' GAMESS ORDER. CHECK INPUT. ABORT...')
 103  format(/,20x,'**** ERROR ****',/,'CONTRACTED ATOMIC',
     . ' ORBITAL NUMBER ',i4,/,'IS OF DXX TYPE. ABORT...',/)

      end subroutine d_shells
c----------------------------------------------------------------

c----------------------------------------------------------------
      subroutine wrdown_orb

      implicit none

      integer :: stat

      open(unit=n_orbit,file='ORB.FILE',form='unformatted',iostat=stat)
      call open_stat('ORB.FILE',stat)

c***
c  NBLOCK5
      write(n_orbit) nbfns,ntmx,ngaus,noc,NON,nbfns,NSBF,NMOF
c  NBLOCK3
      write(n_orbit) ntype, nfirst, nlast
c  NBLOCK4
      write(n_orbit) nr
c  NONE1
      write(n_orbit) none1
c  NBLOCK8
      write(n_orbit) noor, noh
c  NGAUS, IRESUM
      write(n_orbit) ngaus
      write(n_orbit) iresum
c  BLOCK1
      write(n_orbit) vlist
c  BLOCK2
      write(n_orbit) eta
!fk
      write(n_orbit) n_prim
c  D, OVPM
      call wr2(D,nbfns,nbfns,nbfns,nbfns,n_orbit)   
      call wr2(OVPM,nbfns,nbfns,nbfns,nbfns,n_orbit)
c***
      close(n_orbit)

!     if(ntwoel .ne. 3) deallocate( ntype, nfirst, nlast, eta )
!     deallocate( noh, nohin, iresum, vlist )

!     if ( allocated(ncntr) )  then
!       if(ntwoel .ne. 3) deallocate( ncntr )
!     end if
!     if ( allocated(icntr) )  deallocate( icntr )
!     if ( allocated(number) ) deallocate( number )
!     if ( allocated(mlist) )  deallocate( mlist )
!     if ( allocated(none2) )  deallocate( none2 )
!     if ( allocated(none3) )  deallocate( none3 )
!     if ( (nrfock .eq. 0) .and. (ntwoel.ne.3) )     deallocate( d )

      end subroutine wrdown_orb
c----------------------------------------------------------------

c----------------------------------------------------------------
      subroutine wrdown_pse

      implicit none

      integer :: stat

      open(unit=n_pseud,file='PSE.FILE',form='unformatted',iostat=stat)
      call open_stat('PSE.FILE',stat)

      WRITE(n_pseud)IPST, ncpp
      IF(IPST.EQ.1)THEN
      WRITE(n_pseud) LPSEUD
      WRITE(n_pseud) BLKC1
      WRITE(n_pseud) BLKC2
      WRITE(n_pseud) BLKC3
      WRITE(n_pseud) BLKI1
      WRITE(n_pseud) BLKI2
      END IF

      close(n_pseud)

!     if( allocated(lpseud) ) deallocate( lpseud )
!     if( allocated(blkc1) )  deallocate( blkc1 )
!     if( allocated(blkc2) )  deallocate( blkc2 )
!     if( allocated(blkc3) )  deallocate( blkc3 )
!     if( allocated(blki1) )  deallocate( blki1 )
!     if( allocated(blki2) )  deallocate( blki2 )

      end subroutine wrdown_pse
c----------------------------------------------------------------

C----------------------------------------------------------------
      SUBROUTINE WR2(A,NDIM1,NDIM2,N1,N2,itap)

      implicit none

      integer :: NDIM1, NDIM2, N1, N2, itap
      real*8  :: A(NDIM1,NDIM2)

      integer :: i, j

      DO 10 J=1,N2
   10 WRITE(itap)(A(I,J),I=1,N1)

      end subroutine wr2
C----------------------------------------------------------------

c----------------------------------------------------------------
      subroutine ovlgrn_new

      implicit none

      integer :: i, j

      allocate( ovpm(nbfns,nbfns) )

      DO 550 I=1,nbfns
      DO 550 J=1,nbfns
 550      OVPM(J,I)=D(J,I)

      end subroutine ovlgrn_new 
c----------------------------------------------------------------

!----------------------------------------------------------------
      subroutine write_gms_orbitals(orb,nao,nmo,n_unit)

      implicit none

      integer, intent(in) :: nao, nmo
      real*8, intent(inout)  :: orb(nao,nmo)
      integer, intent(in) :: n_unit

      character(len=*), parameter :: fmtgms  = '(i2,1x,i2,5e15.8)'
      integer :: i, j, k, ll, ir, jj, ii, il
      real*8 :: small = 1d-13

      do i = 1,nmo
        do j = 1,nao
          if( abs(orb(j,i)) .lt. small ) orb(j,i) = 0d0
        end do
      end do

      write(n_unit,"(a5)") ' $VEC'
        ll = nao / 5
        ir = mod(nao,5)
        do j=1, nmo
           jj = mod(j,100)
           k = 0
           do il=1,ll
              ii = mod(il,100)
              write(n_unit,fmtgms) jj, ii, (orb(k+i,j), i=1,5)
              k = k + 5
           end do
           ii = mod(ll+1,100)
           if(ir.gt.0) write(n_unit,fmtgms) jj, ii, (orb(k+i,j), i=1,ir)
        end do
      write(n_unit,"(a5)") ' $END'

      end subroutine write_gms_orbitals
!----------------------------------------------------------------

      end module poly_gaus
!----------------------------------------------------------------
