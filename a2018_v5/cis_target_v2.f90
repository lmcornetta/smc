module cis_target

use param
use poly_gaus, only: nbfns, nmof_tar

implicit none

integer, save :: nsa_tar
integer, allocatable, save :: nht_tar(:), npt_tar(:)

!***
! NTRN = number of symmetry planes
! NRR  = number of irreducible representations
!
integer, save :: ntrn, nrr
!***

!***
!   common /ph/ns,nht(nsx),npt(nsx),nhth(nsx),nbur
!   common /phnovo/nbu(nsx),npa(nsx),nsas(8)
!
integer, save :: nsas(8)
integer, allocatable, save :: nbu(:), npa(:)
!***

!***
! Relevant for TARGET CIS and IVO runs:
! Z_SING stores the target singlet eigenvectors
! Z_TRIP stores the target triplet eigenvectors
! D_SING stores the target singlet eigenvalues
! D_TRIP stores the target triplet eigenvalues
!
real*8, allocatable, save :: z_sing(:,:), z_trip(:,:)
real*8, allocatable, save :: d_sing(:),   d_trip(:)
!
! NOCC_MVO is the number of doubly occupied orbitals in MVO runs
!
integer, save :: nocc_mvo
!***

contains

!----------------------------------------------------------------
subroutine target_cis

!***
!  This is the driver for CIS target state calculations. A number of
!  parameters have been redefined (number of hole and particle orbitals, etc.)
!  for use in this type of run. This procedure avoids problems with the
!  corresponding variables used for scattering calculations.
!***

implicit none

!***
! Local variables
!
integer :: irr, nposi
integer :: stat
!***

write(6,"(//,10x,'*** CIS TARGET RUN ***')") 

!   Read target active space
call sctin_cis_target

!   Prepare symmetry adapted active space
call targ_symmetry

open( unit=n_scrci, file='CIS.SCR', form='unformatted', iostat=stat )
call open_stat('CIS.SCR',stat)

!***
!  Calculate CIS states in ALVO subroutine
!  The loop (IRR) runs over symmetry blocks

nposi = 0
do irr = 1,nrr
  if(nsas(irr) .ne. 0)then
    call alvo(nsas(irr), nposi, irr)
    nposi = nposi + nsas(irr)
  endif
end do
!***

!***
! Gather and reorder eigenvalues and eigenvectors
!
call cis_finish
!***

end subroutine target_cis
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine target_ivo

!***
!  This is the driver for target IVOs calculations. A number of
!  parameters have been redefined (number of hole and particle orbitals, etc.)
!  for use in this type of run. This procedure avoids problems with the
!  corresponding variables used for scattering calculations.
!***

implicit none

!***
!  Local variables
!
integer :: irr, nposi
integer :: stat
!***

write(6,"(//,10x,'*** TARGET IVO RUN ***')")

! Define IVO active space (hole and virtual orbitals)
call sctin_ivo_target

! Prepare symmetry adapted active space
call targ_symmetry

open( unit=n_scrci, file='IVO.SCR', form='unformatted', iostat=stat )
call open_stat('IVO.SCR',stat)

!***
!  Calculate single-hole CIS states in ALVO subroutine.
!  The loop (IRR) runs over symmetry blocks
!
nposi = 0
do irr = 1,nrr
  if(nsas(irr) .ne. 0)then
    call alvo(nsas(irr), nposi, irr)
    nposi = nposi + nsas(irr)
  endif
end do
!***

!***
! Gather eigenvalues and eigenvectors (no reordering)
call cis_finish
!***

!***
! Generate IVOs
call gen_ivos
!***

end subroutine target_ivo
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine target_mvo

!***
!  This is the driver for target MVOs calculations.
!***

implicit none

!***
!  Local variables
!***

write(6,"(//,10x,'*** TARGET MVO RUN ***')")

!   Define MVO active space (core charge)
call sctin_mvo_target

!   Calculate the MVOs
call gen_mvos

end subroutine target_mvo
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine sctin_cis_target

use poly_conf, only: ns

!***
!  This subroutine was adapted from the SCTIN0 version of the CORINGA code
!  (module poly-2-fsci.f). It is assumed that a WAVEFUN BLOCK was read
!  in SCTIN_NEW subroutine (POLY_CONF module).
!***

implicit none

!***
!  Local variables
!
integer :: i, n10, i1, i2, j
integer :: stat
!***

!***
!  Warn the user the (N+1)-particle active spae is not relevant for this run.
!  Some data is read in SCTIN_NEW for use in POLY (orthogonalization procedure),
!  but only the target active space is relevant for this run.
!
if( ns .ne. 0 ) then
  write(6,"(//,10x,'*** WARNING: (NS .GT. 0) ***')")
  write(6,"(10x,'THE (N+1)-PARTICLE ACTIVE SPACE IS NOT RELEVANT FOR TARGET CIS RUNS')")
end if
!***

!***
! Input CISTARG block. The number of hole-particle pairs for the cis target
! active space (this is the same as NSA in a scattering calculation with
! electronic excitation)

call ini_flag(flag_t1)

call read_int_flag('NPAIRS =',nsa_tar)

!***
! Allocate and read hole and particle orbitals for CIS target active space.
! NHT_NSA and NPT_NSA are the hole and particle orbitals for the target CIS.
!
allocate( nht_tar(nsa_tar), npt_tar(nsa_tar) )

call check_flag('HOLES','yes')
read(5,*,iostat=stat) ( nht_tar(i), i=1,nsa_tar)
call write_stat('list of holes',stat)

call check_flag('PARTS','yes')
read(5,*,iostat=stat) ( npt_tar(i), i=1,nsa_tar)
call write_stat('list of holes',stat)

write(6,"(//,'*** NUMBER OF TARGET EXCITATIONS IN CIS CALCULATION = ',i6,' ***')") nsa_tar
write(6,"(/'*** HOLE ORBITALS ***'/)")
n10 = nsa_tar / 10
do i = 1,n10
  i1 = (i-1)*10 + 1
  i2 = i1 + 9
  write(6,"(10(i6))")(nht_tar(j), j = i1,i2 ) 
end do
if( mod(nsa_tar,10) .ne. 0 ) write(6,"(10(i6))")(nht_tar(j), j = i2+1,nsa_tar )

write(6,"(/,'*** PARTICLE ORBITALS ***',/)")
do i = 1,n10
  i1 = (i-1)*10 + 1
  i2 = i1 + 9
  write(6,"(10(i6))")(npt_tar(j), j = i1,i2 ) 
end do
if( mod(nsa_tar,10) .ne. 0 ) write(6,"(10(i6))")(npt_tar(j), j = i2+1,nsa_tar )

return

end subroutine sctin_cis_target
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine sctin_ivo_target

use poly_conf, only: ns, nocc

implicit none

!***
! Local variables
!
integer :: nivo, i, j, i1, i2, n10
!***

!***
!  Warn the user the (N+1)-particle active spae is not relevant for this run.
!  Some data is read in SCTIN_NEW for use in POLY (orthogonalization procedure),
!  but only the target active space is relevant for this run.
!
if( ns .ne. 0 ) then
  write(6,"(//,10x,'*** WARNING: (NS .GT. 0) ***')")
  write(6,"(10x,'THE (N+1)-PARTICLE ACTIVE SPACE IS NOT RELEVANT FOR TARGET IVO RUNS')")
end if
!***

!***
! Input IVORBTL block. The number of hole-particle pairs for the cis target
! active space (this is the same as NSA in a scattering calculation with
! electronic excitation)

call ini_flag(flag_i1)

nsa_tar = nmof_tar - nocc
allocate( nht_tar(nsa_tar), npt_tar(nsa_tar) )

call read_int_flag('ORBIVO =',nivo)

if( nivo .lt. 0 ) then
  write(6,"('ORBIVO CANNOT BE .LE. 0')")
  stop'ABORT...'
end if
if( nivo .gt. nocc ) then
  write(6,"('ORBIVO CANNOT BE .GT. NOCC')")
  stop'ABORT...'
end if

nht_tar = nivo

i1 = 0
do i = nocc+1,nmof_tar
  i1 = i1 + 1
  npt_tar(i1) = i
end do
!***

!***
! Print out...
!
write(6,"(//,10x, 'GENERATE IVOs FROM THE HOLE ORBITAL ORBIVO = ',i4)") nivo
write(6,"(10x, 'GENERATE IVOs USING THE FOLLOWING PARTICLE ORBITALS:')")
if( nsa_tar .le. 10) then
  write(6,"(10(i5))") (npt_tar(j), j = 1,nsa_tar )
  return
end if

n10 = nsa_tar / 10
do i = 1,n10
  i1 = (i-1)*10 + 1
  i2 = i1 + 9
  write(6,"(10(i5))") (npt_tar(j), j = i1,i2 )
end do
if( mod(nsa_tar,10) .ne. 0 ) write(6,"(10(i5))") (npt_tar(j), j = i2+1,nsa_tar )

return

end subroutine sctin_ivo_target
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine sctin_mvo_target

use poly_conf, only: ns, nocc, nelc

implicit none

!***
! Local variables
!
integer :: mvoq
integer :: i
integer :: stat
!***

!***
!  Warn the user the (N+1)-particle active space is not relevant for this run.
!  Some data is read in SCTIN_NEW for use in POLY (orthogonalization procedure),
!  but only the target active space is relevant for this run.
!
if( ns .ne. 0 ) then
  write(6,"(//,10x,'*** WARNING: (NS .GT. 0) ***')")
  write(6,"(10x,'THE (N+1)-PARTICLE ACTIVE SPACE IS NOT RELEVANT FOR TARGET MVO RUNS')")
end if
!***

!***
! Input MVORBTL block. The number of hole-particle pairs for the cis target
! active space (this is the same as NSA in a scattering calculation with
! electronic excitation)

call ini_flag(flag_m1)

call read_int_flag('MVOCHR =',mvoq)

if( mod(mvoq,2) .ne. 0 ) then
  write(6,"('MVOCHR MUST BE EVEN')")
  stop'ABORT...'
end if
if( mvoq .le. 0 ) then
  write(6,"('MVOCHR CANNOT BE .LE. 0')")
  stop'ABORT...'
end if
if( mvoq .gt. nelc ) then
  write(6,"('MVOCHR CANNOT EXCEED THE NUMBER OF ELECTRONS')")
  stop'ABORT...'
end if

nocc_mvo = nocc - mvoq/2

call check_flag('HOLES','yes')
allocate( nht_tar(mvoq/2) )
read(5,*,iostat=stat)( nht_tar(i), i=1,mvoq/2)
call write_stat('mvo holes',stat)

call end_flag(flag_m2)

return

end subroutine sctin_mvo_target
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine targ_symmetry

use poly_conf, only: nocc
!***
!  Works out the target symmetry with respect to reflections on Cartesian 
!  planes. Adapted from the SIMETRIA subroutine in CORINGA code.
!***

implicit none

!=======================================================================
!     NTRN = No. of symmetry operations
!     NRR = No. of irreducible representations. If NTRN=1 NRR=2.
!     If NTRN=2 NRR=4. If NTRN=3 NRR=8
!     MA(I,M) = Orbital transformation (I=1,NBFN) under the symmetry operation (M=1,NTRN)
!     NPAD(NRR,I) = Result of symmetry operations (M=1,NTRN) on the 
!     irreducible representations. 
!=======================================================================

!***
! Local variables
!
integer :: i, j, k, m, na, irr, ncont1, ncont2, isinal
integer :: npad(8,3)
integer, allocatable :: ma(:,:)

data  npad/+1,-1,+1,-1,+1,-1,+1,-1,+1,+1,-1,-1,+1,+1,-1,-1,+1,+1,+1,+1,-1,-1,-1,-1/
character(len=*), parameter :: symm_tab_file = 'symm_tab'
integer, allocatable :: new_ma(:,:)
integer :: stat

call read_int_flag('NPLANE =',ntrn)

if( (ntrn.lt.0) .or. (ntrn.gt.3) ) then
  write(6,"('Invalid NPLANE (only 0,1,2,3 allowed) = ',i6,/,'Abort...')")
  stop
end if

!***
! Define the number of irreducible representations (NRR) 
! from the number of symmetry planes (NTRN).
!
nrr = 2 * ntrn
if( ntrn .eq. 0 ) nrr = 1
if( ntrn .eq. 3 ) nrr = nrr + 2
!***

!***
! Allocate CIS variables
!
allocate( nbu(nsa_tar), npa(nsa_tar) )
!***

nsas = 0
if(ntrn.eq.0)then
  do i=1,nsa_tar
    nbu(i)=nht_tar(i)
    npa(i)=npt_tar(i)
  end do
  nsas(1)=nsa_tar

  if( nroute .eq. -1 ) then
    call end_flag(flag_t2)
    write(6,"(//,15x,'*** A SINGLE SYMMETRY BLOCK IN THE CIS ACTIVE SPACE ***',//)")
  else if( nroute .eq. -2 ) then
    call end_flag(flag_i2)
    write(6,"(//,15x,'*** A SINGLE SYMMETRY BLOCK IN THE IVO ACTIVE SPACE ***',//)")
  end if

  return
end if

!***
! Allocate orbital transformation array
!
allocate ( ma(nmof_tar,ntrn) )

!***

!***
! Read up symmetry information

do m = 1,ntrn
  read(5,*,iostat=stat) (ma(i,m),i=1,nmof_tar)
  call write_stat('orbitals symmetries',stat)
end do

if( nroute .eq. -1 ) then
  call end_flag(flag_t2)
else if( nroute .eq. -2 ) then
  call end_flag(flag_i2)
end if

do m = 1,ntrn
  do i = 1,nmof_tar
    ma(i,m)=ma(i,m)/abs(ma(i,m))
  end do
end do

!
!  Build up the active space within symmetry blocks
!
ncont1 = 0
do irr = 1,nrr
  ncont2 = 0
  do i = 1,nsa_tar
    na = nht_tar(i)
    k  = npt_tar(i)
    do j = 1,ntrn
      isinal = ma(na,j)*ma(k,j)
      if(isinal .ne. npad(irr,j)) go to 103
    end do
    ncont1 = ncont1 + 1
    ncont2 = ncont2 + 1
    nbu(ncont1) = na
    npa(ncont1) = k
  103     continue
  end do
  nsas(irr) = ncont2
end do
!***

if( nroute.eq.-2 ) then
  allocate ( new_ma(nmof_tar,ntrn) )
! fk
! Compute new_ma
  do i=1,nocc
    new_ma(i,:) = i * ma(i,:)
  end do
  i=nocc+1
  do irr=1,nrr
    do m=1,nsas(irr)
!   new_ma(i,:) = i * npad(irr,:)
      new_ma(i,:) = i * npad(irr,:) * ma(nht_tar(1),:)
      i=i+1
    end do
  end do

! Print new symmetry table 
  open(n_symmt,file=symm_tab_file,form='formatted')
  do m=1,ntrn
    write(n_symmt,"(*(x,i0))") new_ma(:,m)
  end do
  close(n_symmt)
  write(6,'(/2(a)/)') 'NEW TABLE OF SYMMETRIES WRITTEN TO FILE ', symm_tab_file

  deallocate( new_ma )
end if

deallocate( ma )

!***
!  Print out...
!
write(6,"(//,15x,'*** SYMMETRY ADAPTED TARGET CIS ACTIVE SPACE ***')")
write(6,"(5x,'HOLE ORBITALS')")
write(6,130) (nbu(i),  i = 1,nsa_tar)
write(6,"(5x,'PARTCLES ORBITALS')")
write(6,130) (npa(i),  i = 1,nsa_tar)
write(6,"(5x,'NUMBER OF HOLE_PARTICLE PAIRS IN  EACH SYMMETRY')")
write(6,130) (nsas(i), i = 1,nrr)

write(6,"(//,15x,'*** PAIR SYMMETRIES (NPAD)  ***')")
do j = 1,ntrn
  write(6,130) (npad(irr,j),irr=1,nrr)
end do
!***

 130  format(20i4)
! 140  format(10('i= ',i3,'par ',2i3))

end subroutine targ_symmetry
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine alvo(ndim, nposi, iblk)

implicit none

integer, intent(in) :: ndim, nposi, iblk

!***
!  Local variables
!
integer :: i, k, ii, kk, na, nb, ispin, nconti, ncontk
real*8  :: resh
real*8, allocatable, save :: halvo(:,:)
integer :: alloc_stat
!***

!***
!  Allocate target Hamiltonian (within a give symmetry block)
!
allocate( halvo(ndim,ndim), stat=alloc_stat )
call allocation_stat('halvo',alloc_stat)
halvo = 0d0
!***

write(6,"(//,10x,'*** SYMMETRY BLOCK = ',I3,' ***')") iblk

Do ispin = 1,2
!
!   ISPIN = 1,2 means singlets and triplets, respectively
!
  do ii = 1,ndim
    nconti = nposi + ii
    NA = nbu(nconti)
    I =  npa(nconti)
    do KK = 1,II
      ncontk = nposi + kk
      NB = nbu(ncontk)
      K =  npa(ncontk)
      CALL ALVOC(NA,I,NB,K,RESH,ISPIN)
      HALVO(KK,II) = RESH
      HALVO(II,KK) = RESH
    end do
  end do
!
!   For each symmetry block, diagonalize in singlets and then
!   triplets space. The eigenvalues and eigenvectors are stored
!   block by block (first singlets then triplets within blocks).
!
  call halvo_diag(halvo, ndim, ispin)

End Do

deallocate( halvo ) 

end subroutine alvo
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine alvoc(na, i, nb, k, res, n)

use poly_conf, only: eg0, rpmol, f

! =======================================================================
!
!     IN THIS SUBROUTINE IS CALCULATED THE FOLLOWING MATRIX ELEMENT
!
!
!  <NA-I> (J) HN <NB-K> (L)
!
!     WHERE:
!     <NA-I> <J> ==>> MEANS THAT J IS NOT ANTISYMMETRIZED WITH <NA-I>
!
!     <NB-K> (L) ==>> MEANS THAT L IS ANTISYMMETRIZED WITH <NB-K>
!
!     HN = TARGET HAMILTONIAN => SUM i h(i) + SUM i e j  1/rij
!
! =======================================================================

implicit none

integer, intent(in) :: na, i, nb, k, n
real*8, intent(out)  :: res

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
! =======================================================================
!
!    <NA-I>  HN  <NA-K>
!
! ======================================================================

IA=IABS(I)
KA=IABS(K)
NAA=IABS(NA)
!
IF(N.EQ.3) GO TO 20
IF(N.EQ.2) GO TO 10
RES = F(IA,KA)-RPMOL(IA,KA,NAA,NAA)+2*RPMOL(IA,NAA,KA,NAA)
RETURN
!
   10 RES = F(IA,KA)-RPMOL(IA,KA,NAA,NAA)
RETURN
!
   20 RES = F(IA,KA)-RPMOL(IA,KA,NAA,NAA)
RETURN
!
! =======================================================================
!
!    <NA-I>  HN  <NA-I>
!
! =======================================================================
!
  360 NAA=IABS(NA)
IA=IABS(I)
!
!mod
!     print*,'allocated(f)',allocated(f)
!     print*,'naa,ia,',naa,ia
!     print*,'eg0',eg0
!     print*,'RPMOL(IA,IA,NAA,NAA)',RPMOL(IA,IA,NAA,NAA)
!mod

IF(N.EQ.3) GO TO 40
IF(N.EQ.2) GO TO 30
RES = eg0-F(NAA,NAA)+F(IA,IA)-RPMOL(IA,IA,NAA,NAA) &
          +2*RPMOL(IA,NAA,IA,NAA)
RETURN
!
   30 RES = eg0-F(NAA,NAA)+F(IA,IA)-RPMOL(IA,IA,NAA,NAA)
RETURN
!
   40 RES = eg0-F(NAA,NAA)+F(IA,IA)-RPMOL(IA,IA,NAA,NAA)
RETURN
!
  440 IF(I.EQ.K) GO TO 370
!
! =======================================================================
!
!    <NA-I>  HN  <NB-K>
!
! =======================================================================
!
IA=IABS(I)
KA=IABS(K)
NBA=IABS(NB)
NAA=IABS(NA)
!
IF(N.EQ.3) GO TO 60
IF(N.EQ.2) GO TO 50
RES = 2*RPMOL(NAA,IA,NBA,KA)-RPMOL(NAA,NBA,IA,KA)
RETURN
!
   50 RES =  -RPMOL(NAA,NBA,IA,KA)
RETURN
!
   60 RES =  -RPMOL(NAA,NBA,IA,KA)
RETURN
!
! =======================================================================
!
!    <NA-I>  HN  <NB-I>
!
! =======================================================================
!
  370 NAA=IABS(NA)
NBA=IABS(NB)
IA=IABS(I)
!
IF(N.EQ.3) GO TO 80
IF(N.EQ.2) GO TO 70
RES = -F(NBA,NAA)-RPMOL(NBA,NAA,IA,IA)+2*RPMOL(NBA,IA,NAA,IA)
RETURN
!
   70 RES = -F(NBA,NAA)-RPMOL(NBA,NAA,IA,IA)
RETURN
!
   80 RES = -F(NBA,NAA)-RPMOL(NBA,NAA,IA,IA)

end subroutine alvoc
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine halvo_diag(h, n, ispin)

use poly_conf, only: eg0

implicit none

integer, intent(in) :: n, ispin
real*8, intent(inout)  :: h(n,n)

!***
!  Local variables
!
integer :: i, j, i1, i2, n10
integer :: lwork, info, n_eig
real*8 :: w_work, e_low=0d0, e_upp=0d0
character(len=7) :: spin
real*8, allocatable :: e(:), z(:,:)
real*8, allocatable :: work(:), iwork(:), ifail(:)
real*8, parameter :: abstol = 1.d-08
integer :: alloc_stat
!***

!***
!c!   SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
!c!  $ ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
!c!  $ IFAIL, INFO )
!***

!
!   E, Z are the eigenvalues and eigenvectors, respectively.
!   RANGE = 'A' calculates the full eigenvalue spectrum
allocate( e(n), z(n,n), stat=alloc_stat)
call allocation_stat('z(n,n) (in halvo_diag subroutine)',alloc_stat)
!
!   These are related to workspace
allocate( iwork(5*n), ifail(n) )

!
!   First call (LWORK=-1) to optimize work space
lwork = -1
call dsyevx('V', 'A', 'U', n, h, n, e_low, e_upp, 0, 0, abstol, n_eig, e, z, n, w_work, lwork, iwork, ifail, info)
!mod  write(6,"(/,5x,'DSYEVX FIRST CALL: LWORK = ',i8)") int(w_work)

!
!   Allocate workspace and diagonalize
lwork = int( w_work )
allocate( work(lwork), stat=alloc_stat )
call allocation_stat('work (in halvo_diag subroutine)',alloc_stat)
call dsyevx('V', 'A', 'U', n, h, n, e_low, e_upp, 0, 0, abstol, n_eig, e, z, n, work, lwork, iwork, ifail, info)
!
!   Success check
if( info .lt. 0 ) then
  write(6,"(/,5x,'INFO = ',i8,/,5x, 'if INFO = -i, the i-th argument had an illegal value',/,&
        &'Abort...')")
  stop
else if( info .gt. 0) then
  write(6,"(/,5x,'INFO = ',i8,/,5x, 'if INFO = i, then i eigenvectors failed to converge',/,&
        &'Abort...')")
  stop
!mod  else
!mod    write(6,"(/,5x,'INFO = 0, SUCCESSFUL DIAGONALIZATION')")
end if

!***
! Return from here in MVO runs
if( nroute .eq. -3 ) then
  h = z
  write(6,"(//,'MVO ENERGIES IN HARTREE:')")
  n10 = n / 10
  do i = 1,n10
    i1 = (i-1)*10 + 1
    i2 = i1 + 9
    write(6,"(10(1pe13.6,2x))") (e(j), j = i1,i2 )
  end do
  if( mod(n,10) .ne. 0 ) write(6,"(10(1pe13.6,2x))") (e(j), j = i2+1,n )
  return
end if
!***

!***
! Print out...
if( ispin .eq. 1 ) write(spin,"(a7)") 'SINGLET'
if( ispin .eq. 2 ) write(spin,"(a7)") 'TRIPLET'
  
write(6,"(5x,a7,' STATES EXCITATION ENERGIES IN EV:')") spin
n10 = n / 10
if(n10.ge.1) then
  do i = 1,n10
    i1 = (i-1)*10 + 1
    i2 = i1 + 9
    write(6,"(10(1pe13.6,2x))") ((e(j)-eg0)*toHartree, j = i1,i2 )
  end do
else
  write(6,"(10(1pe13.6,2x))") ((e(j)-eg0)*toHartree, j = 1,n )
end if

if( mod(n,10) .ne. 0 ) write(6,"(10(1pe13.6,2x))") ((e(j)-eg0)*toHartree, j = i2+1,n )
!***

!***
! Save to file...
write(n_scrci) n
write(n_scrci) e
write(n_scrci) z
!***

!mod
!     do i = 1,n
! if(ispin.eq.1) write(6,"('singleto, e = ',1pe13.6)")e(i)
! if(ispin.eq.2) write(6,"('tripleto, e = ',1pe13.6)")e(i)
! do j = 1,n
!   write(6,"('z(',i2,',',i2,') = ',1pe13.6)")j,i,z(j,i)
! end do
! write(6,"(' ')")
!     end do
!mod

end subroutine halvo_diag
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine cis_finish

use poly_conf, only: eg0

implicit none

!***
!  Local variables
!
integer :: irr, i, j, is, js, pos, nn, ndim_sum
real*8, allocatable :: z_scr(:,:,:), d_scr(:,:)
integer :: stat
!***

allocate( z_sing(nsa_tar,nsa_tar), z_trip(nsa_tar,nsa_tar) )
allocate( d_sing(nsa_tar), d_trip(nsa_tar) )

z_sing = 0d0
d_sing = 0d0
z_trip = 0d0
d_trip = 0d0

ndim_sum = 0
rewind(n_scrci)

pos = 0
do irr = 1,nrr
  if( nsas(irr) .ne. 0 ) then
    read(n_scrci) nn
    ndim_sum = ndim_sum + nn
    if(nn .ne. nsas(irr)) stop'(nn .ne. nsas(irr))'
    if(ndim_sum .gt. nsa_tar) stop'(ndim_sum .gt. nsa_tar)'
    allocate( z_scr(nn,nn,2), d_scr(nn,2) )
    read(n_scrci) d_scr(:,1)
    read(n_scrci) z_scr(:,:,1)
    read(n_scrci) nn
    read(n_scrci) d_scr(:,2)
    read(n_scrci) z_scr(:,:,2)
!
!   Build block-diagonal matrix
!
    do is = 1,nsas(irr)
      i = pos + is
      d_sing(i) = d_scr(is,1)
      d_trip(i) = d_scr(is,2)
      do js = 1,nsas(irr)
        j = pos + js
        z_sing(j,i) = z_scr(js,is,1)
        z_trip(j,i) = z_scr(js,is,2)
      end do
    end do
    pos = pos + nsas(irr)
    deallocate( z_scr, d_scr )
  end if
end do

close(n_scrci)

!
!   If this is not a CIS target run, return
!
if( nroute .lt. -1) return

!
!   Reorder: first singlets then triplets
!
if( ntrn .gt. 0 ) then
  call reord(d_sing, z_sing, nsa_tar)
  call reord(d_trip, z_trip, nsa_tar)
end if

!
!   Save: first singlets then triplets
!
open( unit=n_hci_1, file='HCI1.FILE', form='unformatted', iostat=stat )
call open_stat('HCI1.FILE',stat)
write(n_hci_1) nsa_tar
write(n_hci_1) d_sing
write(n_hci_1) z_sing
close(n_hci_1)
 
open( unit=n_hci_2, file='HCI2.FILE', form='unformatted', iostat=stat )
call open_stat('HCI2.FILE',stat)
write(n_hci_2) nsa_tar
write(n_hci_2) d_trip
write(n_hci_2) z_trip
close(n_hci_2)

!***
!  Print out...
!
write(6,"(//,10x,'*** TRIPLET EXCITATION ENERGIES IN ASCENDING ORDER (EV) ***',/)")
do j=1,nsa_tar
  write(6,"(f18.10)") (d_trip(j)-eg0)*toHartree
end do
!     n10 = nsa_tar / 10
!     do i = 1,n10
! i1 = (i-1)*20 + 1
! i2 = i1 + 19
! write(6,"(10(1pe13.6,2x))") ((d_trip(j)-eg0)*toHartree, j = i1,i2 )
!     end do
!     if( mod(nsa_tar,10) .ne. 0 ) then
! if( n10 .eq. 0 ) i2 = 0
! write(6,"(10(1pe13.6,2x))") ((d_trip(j)-eg0)*toHartree, j = i2+1,nsa_tar )
!     end if
!
write(6,"(//,10x,'*** SINGLET EXCITATION ENERGIES IN ASCENDING ORDER (EV) ***',/)")
do j=1,nsa_tar
  write(6,"(f18.10)") (d_sing(j)-eg0)*toHartree
end do
!     n10 = nsa_tar / 10
!     do i = 1,n10
! i1 = (i-1)*20 + 1
! i2 = i1 + 19
! write(6,"(10(1pe13.6,2x))") ((d_sing(j)-eg0)*toHartree, j = i1,i2 )
!     end do
!     if( mod(nsa_tar,10) .ne. 0 ) then
! if( n10 .eq. 0 ) i2 = 0
! write(6,"(10(1pe13.6,2x))") ((d_sing(j)-eg0)*toHartree, j = i2+1,nsa_tar )
!     end if
!***

deallocate( d_sing, d_trip, z_sing, z_trip )

end subroutine cis_finish
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine gen_ivos

use poly_gaus, only: d
use poly_conf, only: nocc

!***
!  Local variables
!
integer :: i, j, k, kk
real*8  :: ivoc
real*8, allocatable :: ivorbs(:,:)
integer :: stat
!***

!     allocate( ivorbs(nbfns,nsa_tar) )
allocate( ivorbs(nbfns,nbfns) )

! Complete for the occupied orbitals
do i = 1,nocc
  do j = 1,nbfns
    ivorbs(j,i) = d(j,i)
  end do
end do

!***
! Obtain IVOs (IVORBS) from VOs (D) and single-hole CIS eigenstates (Z)
! First the singlets
!
do i = 1,nsa_tar
  do j = 1,nbfns
    ivoc = 0d0
    do k = 1,nsa_tar
      kk = npa(k)
      ivoc = ivoc + d(j,kk) * z_sing(k,i)
    end do
    ivorbs(j,i+nocc) = ivoc
  end do
end do

open( unit=n_ivorb_s, file='IVOS_SING.LIS', form='formatted', iostat=stat )
call open_stat('IVO_SING.LIS',stat)
call write_orbitals(ivorbs,nbfns,nocc+nsa_tar,n_ivorb_s)
close(n_ivorb_s)

do i = 1,nsa_tar
  do j = 1,nbfns
    ivoc = 0d0
    do k = 1,nsa_tar
      kk = npa(k)
      ivoc = ivoc + d(j,kk) * z_trip(k,i)
    end do
    ivorbs(j,i+nocc) = ivoc
  end do
end do

open( unit=n_ivorb_t, file='IVOS_TRIP.LIS', form='formatted', iostat=stat )
call open_stat('IVO_TRIP.LIS',stat)
call write_orbitals(ivorbs,nbfns,nocc+nsa_tar,n_ivorb_t)
close(n_ivorb_t)

!***

!mod
!     do i = 1,nsa_tar
! write(6,"('i = ',i3)") i
! write(6,"(20(1pe13.6,1x))") (z_trip(j,i), j=1,nsa_tar)
!     end do
!     write(6,"('ivos')")
!     do i = 1,nsa_tar
! write(6,"('i = ',i3)") i
! write(6,"(20(1pe13.6,1x))") (ivorbs(j,i), j=1,nbfns)
!     end do
!mod

deallocate( ivorbs )

end subroutine gen_ivos
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine gen_mvos

use poly_gaus, only: tket, atrt, nmof_tar, d
use poly_conf, only: nocc, rpmol

implicit none

!***
! Local variables
!
integer :: i, j, k, m, ii, jj, kk, nmin
integer :: nbfn, nsbf, insbf, ndim
real*8  :: xm, ym
real*8, allocatable :: mvos(:,:)
integer :: stat
!***

nbfn = nmof_tar
nsbf = nmof_tar
ndim = nmof_tar - nocc
nmin = nocc + 1

!***
! MVO Fock matrix in z_sing array
!
allocate( z_sing(ndim,ndim), mvos(nbfns,nmof_tar) )

ii = 0
do I=nmin,NBFN
  ii = ii + 1
  INSBF=MIN0(NSBF,I)
  jj = 0
  do J=nmin,INSBF
    jj = jj + 1
    XM=0.0D0
    YM=0.0D0
    do M=1,nocc_mvo
      if( any(nht_tar.eq.m) ) cycle
      XM=XM+RPMOL(I,J,M,M)
      YM=YM+RPMOL(I,M,J,M)
    end do
    YM=2.0d0*XM-YM
    z_sing(jj,ii)=YM+TKET(I,J)+ATRT(I,J)
    z_sing(ii,jj)=z_sing(jj,ii)
  end do
end do
!***

!***
! Diagonalize the Fock matrix and generate MVOs
! 
call halvo_diag(z_sing, ndim, 1)

do i = 1,ndim
  do j = 1,nbfns
    xm = 0d0
    do k = 1,ndim
      kk = k + nocc
      xm = xm + d(j,kk) * z_sing(k,i)
    end do
    mvos(j,i+nocc) = xm
  end do
end do
!***

! Complete for the occupied orbitals
do i = 1,nocc
  do j = 1,nbfns
    mvos(j,i) = d(j,i)
  end do
end do

!***
! Save MVOs
open( unit=n_mvorb, file='MVOS.LIS', form='formatted', iostat=stat )
call open_stat('MVOS.LIS',stat)
call write_orbitals(mvos,nbfns,nmof_tar,n_mvorb)
close(n_mvorb)
!***

end subroutine gen_mvos
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine reord(a, b, n)

implicit none

integer, intent(in) :: n
real*8, intent(inout)  :: a(n), b(n,n)
!***
!  Local variables
!
integer :: i, j, tot, ilow
real*8  :: upp, upper
integer, allocatable :: iaux(:)
real*8,  allocatable :: aux(:), baux(:,:)
logical, allocatable :: found(:)
!***

allocate( iaux(n), aux(n), found(n), baux(n,n) )

!
!   Find the highest
!
upp = a(1)
do i = 2, n
  upp = max(upp,a(i))
end do 
    
!
!   and reorder
!
tot = 0
found = .false.
upper = upp
do
  do i = 1,n
    if( .not. found(i) ) then
      if( a(i) .le. upp ) then
        upp  = a(i)
        ilow = i
      end if
    end if
  end do
  tot = tot + 1
  iaux(tot) = ilow
  found(ilow) = .true.

  if( tot .eq. n ) exit
  upp = upper
end do

deallocate( found )

!mod  do i = 1,n
!mod    write(6,"(i3,2x,1pe13.6)") i,a(i) 
!mod    write(6,"(10(1pe13.6,2x))") (b(j,i), j=1,n) 
!mod  end do

do i = 1,n
  aux(i) = a(iaux(i))
  do j = 1,n
    baux(j,i) = b(j,iaux(i))
  end do
end do
a = aux
b = baux

!mod  do i = 1,n
!mod    write(6,"(i3,2x,1pe13.6)") i,a(i) 
!mod    write(6,"(10(1pe13.6,2x))") (b(j,i), j=1,n) 
!mod  end do

deallocate( aux, baux )


end subroutine reord
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine write_orbitals(orb,nao,nmo,n_unit)

implicit none

integer, intent(in) :: nao, nmo
real*8, intent(inout)  :: orb(nao,nmo)
integer, intent(in) :: n_unit

!***
! Local variables
!
integer :: i, j, k, ll, ir, jj, ii, il
character(len=18) orbt
!***

!***
!  Zero out small coefficients
!
do i = 1,nmo
  do j = 1,nao
    if( abs(orb(j,i)) .lt. 1d-13 ) orb(j,i) = 0d0
  end do
end do
!***

if(orbfmt .eq. "ALCHEM") then
  write(orbt,"(a18)") '   ORBITAL ='
  do i=1,nmo
    write(n_unit,"(a18,i3)") orbt, i
    write(n_unit,112) (orb(J,I),J=1,nao)
  end do
else
  write(n_unit,"(a5)") ' $VEC'
  LL = nao / 5
  IR = MOD(nao,5)
  DO J=1, nmo
     JJ = MOD(J,100)
     K = 0
     DO IL=1, LL
  II = MOD(IL,100)
  write(n_unit,10) JJ, II, (orb(K+I,J), I=1, 5)
  K = K + 5
     ENDDO
     II = MOD(LL+1,100)
     IF(IR.GT.0) write(n_unit,10) JJ, II, (orb(K+I,J), I=1, IR)
  ENDDO
  write(n_unit,"(a5)") ' $END'
end if

!  10   format(I2,1X,I2,5E15.8)
   10   format(I2,1X,I2,1p5e15.8)
  112   FORMAT(5D15.8)

end subroutine write_orbitals
!----------------------------------------------------------------

end module cis_target
!----------------------------------------------------------------
