module onk3dk

use paramB
   
use precision

implicit none

!***
!
!  Modifications introduced by Fábris Kossoski in the begining of 2017.
!
!  The module precision-f90.f (same as used in part A) is now used in
!  part B. This is required for calculations which involve many
!  configurations, as 'ndd2' could not fit into a double precision
!  variable.
!
!  The elmbmod subroutine, which evaluates the [(-P+1/(N+1))E] term, 
!  now runs in part A, since the collision energy E appears as a simple
!  multiplicative factor. Part A now writes in WFN.FILE the usual
!  'eaqqd' (the energy-independent part) and also 'eaqqd2',
!  which is the term linear in energy.
!
!  onk and 3dk runs now contain a loop over energy points.
!
!  Within a single SVD decomposition, it is now possible to evaluate
!  the scattering amplitudes for a given combination of discarded
!  singular vectors.
!
!  The flag 'IN_MEM =' reads either 'OFF' or 'DEN', which is relevant
!  for 3dk runs, and is related to the use of memory. Now one can
!  make a much more intensive use of RAM memory, thus avoiding repetitive 
!  and slow acess to the disk. The user should evaluate the best use
!  of available resources before submiting a production job. Particular 
!  attention should be given to the number of energy points of a run.
!  The following strategies have been implemented:
!  a) 'off in mem': off-shell tapes are read first and kept in
!  'offd_all', and afterwards the denominator matrix 'raqqd_tri' is 
!  constructed one energy point at a time, making use of 'off_all'. 
!  b) 'den in mem': the energy dependent 'eaqqd1' denominator matrix
!  is kept in memory, and for each reading of an off-shell tape the 
!  whole array is updated. 
!  The program accesses the vkgkv routines with termination 
!  '_offinmem' for strategy a) and with '_deninmem' for strategy b).
!
!  As a rule of thumb, strategy 'off in mem' is better whenever the
!  number of energy points is greater than the number of off-shell
!  tapes, which in practice is the case for elastic calculations.
!  For multichannel calculations with several open channels, the 'den in
!  mem' strategy might be best option.
!
!  The previous elmnm1 subroutine has been dismembered into
!  elmnm1_s (for singlets) and elmnm1_t (for triplets), which 
!  dismissed if statements and uncessary variables, but more
!  importantly, allowed for a more case-oriented approach.
!
!  The old elmnm1_off has been totally redesigned. The code is much more 
!  vectorized, as represented by the 'cr(nsa)', 'crr(nsa,2)', 
!  'crr1(nsa,3)' and 'crr2(nsa,3)' arrays. The old pvgg2, pvgs2 
!  and pvsg2 parts are now computed locally, and only the pvss2 part 
!  requires external calls. Also, pvss2 has been parallelized with 
!  OpenMP, while for the cheaper pvsg, etc. parts, parellization doesn't 
!  help. Few modifications were introduced when quartet configurations 
!  are present, and in this case the code still counts on the old pvss2 
!  subroutine.
!
!  The content of the old crpmol subroutine (called in fock but mostly
!  in elmnm1) has been incorporated into the calling part.
!
!  inv_luf_kmt, inv_svd_kmt and t_matrix have been generalized for the
!  case of multichannel calculations.
!
!  inv_luf_kmt subroutine has been simplified. Now it calls the 
!  LU decomposition routine (dgetrf) and then the linear system solver
!  (zgetrf).
!
!  inv_svd and inv_svd_kmt now employ the faster ?gesdd, instead of the
!  older ?gesvd Lapack routines.
!
!  The subroutine invmod has been removed, as the inv_luf_kmt and
!  inv_luf alternatives always present better performance.
!
!  cb is kept in memory, which avoids unecessary disk access in
!  crossec
!
!  New variables have been introduced or had its dimensions
!  augmented: 'cb', 'kmat', 'offd_all', 'raqqd_tri', 'raqqd1_tri',
!  'raqqd2_tri', 'iaqqd_tri', and the analogous ones for quartets.
!
!  hemdois is no longer called in off_driver and vkgkv subrotines,
!  since the values of cvec array from nop+1 up to nop2 are not
!  required in the following smn and smn_off routines.
!
!  smn and smn_off now make use of the real and imaginary components
!  cvecr and cveci, which introduce the following advantadges in the 
!  calculation of eaqqd1: 
!  (i) the sum is performed in the fastest index, 
!  (ii) when accounting for the weightings 'rtww3', there are fewer
!  mulltiplications by a factor of 'nd/2',
!  (iii) aimag is called by a factor 'nd' fewer times than the old
!  dconjg function.
!
!  Several declared and unused variables have been removed, as well as
!  some use/only unnecessary statements.
!
!  The log files are much more concise, as uncessary information for the
!  user is no longer printed.
!
!***

!***
! This version was modified to improve the interface with the MPLANE
! module. The latter was modified to allow for parallelization with
! openmp and split in two parts, namely for electronic and nuclear
! < K G_1 |V| G_2 G_3 > integrals.
! 
! Subroutine PLANE was modified for paralellization eith openmp. At this
! point, only the calculation and transformation (AO basis -> MO basis)
! of electronic integrals has been parallelized.
!***

!***
!  This module was modified by MAPL to debug the multichannel calculations.
!  The modifications can be tracked by the cmal comments.
!***

!***
!  Modified to invert the denominator in two steps: the real part of the
!  standard denominator matrix (NDxND) followed by a complex NOP2xNOP2
!  matrix. This is controlled by the flag invtyp.eq."LUF_KMT" or
!  invtyp.eq."SVD_KMT"
!***

!***
!  Modified to facilitate the use of different angular quadratures in
!  the different radial points of the off-shell integration. In OFF
!  runs, the WROFF subroutine saves the quadrature type (GAUSS/LEBEDEV)
!  and the number of qudrature points in every off-shell tape, if saved
!  in the NDNO form (nothing changes if saved in the NDND form). In
!  3DK runs, the radial integration employs the SMN_OFF subroutine in
!  case the off-shell tapes were saved in the NDNO form (nothing changes
!  in the NDND form). SMN_OFF reads in the quadrature type and the
!  number of points and carries out the angular integration (therefore
!  independently for each off-shell tape). The angular quadrature defined 
!  in the input file for 3DK runs (NOTHT, NOPHI, NLEB, etc.) are used only 
!  for the on-shell matrix elements.
!***

!***
!  Modified to employ either Gauss-Legendre or Gauss-Laguerre
!  quadratures for the radial off-shell integration above KMAX.
!  The inner radial integration is always performed with Gauss-Legendre
!  quadratures. This is controlled by the KRADTYP flag ("LEGEN" or
!  "LAGUE").
!***

!***
! Numerator matrix elements
!     COMMON /WORK/CVEC1(NDDX,NOP2X),CVEC2(NDDX,NOP2X),WORK2(7)
complex*16, allocatable, dimension(:,:), save   :: cvec1, cvec2
complex*16, allocatable, dimension(:,:,:), save :: cvec1_off, cvec2_off
integer, allocatable, dimension(:,:), save :: nmf

! Related to the effective number of open channels, relevant in off-shell
! calculations
logical :: do_ground, do_loop_sing, do_loop_trip
integer :: nchlg_eff, nchls_eff, nchlt_eff, nchl_eff
integer, allocatable, dimension(:) :: chl_eff

!     COMMON /MOLTWO/cprt(NSBFX,NONE1X),CPAT(NSBFX)
complex*16, allocatable, dimension(:,:), save :: cprt
complex*16, allocatable, dimension(:), save   :: cpat

!     COMMON /PLFOCK/CG(NSBFX),CG2(NSBFX)
complex*16, allocatable, dimension(:), save :: cg, cg2
!***

!     COMMON /SYMCHK/ischk(ischkX)
integer, allocatable, dimension(:), save :: ischk

!     DIMENSION CAQQD(NDDX,NDDX),CAQQQ(NDQX,NDQX)
  
!***
!  These are related to invtyp.eq."LUF_KMT" or invtyp.eq."SVD_KMT"
logical, save :: inv_real
logical, save :: cb_inmem = .true.
complex*16, allocatable, dimension(:,:,:,:), save :: cb
complex*16, allocatable, dimension(:,:,:), save   :: kmat
 
real*8, allocatable, dimension(:,:), save :: offd_all, offq_all

! raqqd1_tri: denominator matrix from part A (independent of energy)
! raqqd2_tri: denominator matrix from part A (linear in energy)
! raqqd_tri: full real denominator matrix in triangular format
! raqqd: full real denominator matrix in square format
! caqqd: full complex denominator matrix in square format
real*8, allocatable, dimension(:), save   :: raqqd1_tri, raqqq1_tri, raqqd2_tri, raqqq2_tri
real*8, allocatable, dimension(:,:), save :: raqqd_tri,  raqqq_tri, raqqd, raqqq
complex*16, allocatable, dimension(:,:), save :: caqqd, caqqq

contains

!----------------------------------------------------------------
subroutine onk_driver

use legacy,  only: ndd, ndq, nchls, npct, rdup_orb, nchl, keyq
use quadsym, only: nop, nop2, pw, ao2mo_tran 
use mplane1, only: setfn, abint

implicit none

integer :: enen, i, kk, ich, ntape
real*8  :: p(3) 
integer :: alloc_stat

!  Read up orbital information relevant for numerator matrix elements
call rdup_orb(1)

!  Load a few arrays and constants
call setfn

!  Evaluate which atomic to molecular orbital transformations are actually required
call ao2mo_tran

!  Allocate symmetry check array
allocate( ischk(nop*npct) )

!  Allocate cvec1
allocate( cvec1(ndd,nop2), stat=alloc_stat )
call allocation_stat('cvec1 (onk_driver subroutine)',alloc_stat)
! Allocate cvec2 in case spin quartets are taken into account
if( keyq .ne. 0 ) then
  allocate( cvec2(ndq,nop2), stat=alloc_stat )
  call allocation_stat('cvec2 (onk_driver subroutine)',alloc_stat)
end if

kk = 0

call abint

do enen=1,nener

!------------------------------------------------------------------------
!     Run over singlet open channels
!------------------------------------------------------------------------
  do ich=1,nchls+1
    ischk = 0
    do i=1,nop
      if( ischk(i) .eq. 0 ) then
        kk=kk+1
        p(:) = pw(:,i,ich,enen)
        if(keyplm.le.0) call plane(p, nplrtp(ich,enen))
        call fockpl(nplmtp(ich,enen))
        call elmnm1_s(i,ich)
      end if
    end do
    ntape = numtap(ich,enen)
    print*,'WRNUM',' ich=', ich, 'ener =', elenev(enen)
    call wrnum(cvec1,ndd,nop,ntape)
  end do

!------------------------------------------------------------------------
!     Run over triplet open channels
!------------------------------------------------------------------------
  do ich=nchls+2,nchl
    ischk = 0
    do i=1,nop
      if( ischk(i) .eq. 0 ) then
        kk=kk+1
        p(:) = pw(:,i,ich,enen)
        if(keyplm.le.0) call plane(p, nplrtp(ich,enen))
        call fockpl(nplmtp(ich,enen))
        call elmnm1_t(i,ich)
      end if
    end do
    ntape = numtap(ich,enen)
    print*,'WRNUM',' ich=', ich, 'ener =', elenev(enen)
    call wrnum(cvec1,ndd,nop,ntape)
    if( keyq .ne. 0 ) call wrnum(cvec2,ndq,nop,ntape)
  end do

end do ! enen=1,nener

deallocate( ischk )
deallocate( cvec1 )
if( keyq .ne. 0 ) deallocate( cvec2 )

write(6,'(//a,i0)') '** NUMBER OF PLANE WAVES COMPUTED = ', kk

end subroutine onk_driver
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine d3k_driver

use legacy,  only: ndd, ndd2, ndq2, keyq, nchl, free_wfn
use quadsym, only: free_inv 

implicit none

integer(i_large) :: nm 
integer :: enen
integer :: alloc_stat

! Free some memory
call free_wfn
call free_inv

if(nchout.ne.0) then
  nout1 = nchout
  nout2 = nchout
else
  nout1 = 1
  nout2 = nchl
end if

call set_nmf

inv_real = .false.
if( (invtyp.eq."LUF_KMT") .or. (invtyp.eq."SVD_KMT") ) inv_real = .true.
    
! Skip denominator calculation if read from file
!
if( keyden .eq. 1 ) then
  write(6,"(/,5x,'*** SKIP DENOMINATOR CALCULATION ***')")
  write(*,'(/a,i0/)') 'INV FOR NDD = ', ndd
  call allocate_denominator

  call allocate_cb

  do enen=1,nener
    call read_denominator(enen)
    call solve_linear_system(enen)
    call crosec(enen)
  end do 
  return

end if ! keyden.eq.1

! Read denominator matrix elements calculated in part A
call rdup_den

if( off_inmem ) then

  allocate( raqqd_tri(ndd2,1), stat=alloc_stat )
  call allocation_stat('raqqd_tri',alloc_stat)
  if( keyq.ne.0 ) then
    allocate( raqqq_tri(ndq2,1) )
    call allocation_stat('raqqq_tri',alloc_stat)
  end if

  call vkgkv3_off_offinmem

else

  allocate( raqqd_tri(ndd2,nener), stat=alloc_stat )
  call allocation_stat('raqqd_tri',alloc_stat)
  !$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(enen, nm)
  !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nener,ndd2,raqqd_tri,raqqd1_tri,raqqd2_tri,ene)
  do enen=1,nener
    do nm=1,ndd2
      raqqd_tri(nm,enen) = raqqd1_tri(nm) + raqqd2_tri(nm)*ene(enen)
    end do
  end do
  !$OMP  END PARALLEL DO
  deallocate( raqqd1_tri, raqqd2_tri )
  if( keyq .ne. 0 ) then
    allocate( raqqq_tri(ndq2,nener), stat=alloc_stat )
    call allocation_stat('raqqq_tri1',alloc_stat)
    !$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(enen, nm)
    !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nener,ndq2,raqqq_tri,raqqq1_tri,raqqq2_tri,ene)
    do enen=1,nener
      do nm=1,ndq2
        raqqq_tri(nm,enen) = raqqq1_tri(nm) + raqqq2_tri(nm)*ene(enen)
      end do
    end do
    !$OMP  END PARALLEL DO
    deallocate( raqqq1_tri, raqqq2_tri )
  end if
  write(*,'(/a/)') 'ELMBMOD 1, 2'

  call vkgkv3_off_deninmem

end if 

write(*,'(/a,i0/)') 'INV FOR NDD = ', ndd

call allocate_denominator

call allocate_cb

do enen=1,nener

  if( off_inmem ) then
    call vkgkv3_on_offinmem(enen)
  else 
    call vkgkv3_on_deninmem(enen)
  end if

  call solve_linear_system(enen)

  call crosec(enen)

end do ! l=1,nener

end subroutine d3k_driver
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine rdup_den

use legacy, only: ndd2, ndq2, keyq

implicit none

integer :: stat, alloc_stat

open( unit=n_apden, file='DEN.FILE', form='unformatted',iostat=stat)
call open_stat('DEN.FILE',stat)

  allocate( raqqd1_tri(ndd2), stat=alloc_stat )
  call allocation_stat('raqqd1_tri',alloc_stat)
  read(n_apden,iostat=stat) raqqd1_tri
  call write_stat('DEN.FILE',stat)
  allocate( raqqd2_tri(ndd2), stat=alloc_stat )
  call allocation_stat('raqqd2_tri',alloc_stat)
  read(n_apden,iostat=stat) raqqd2_tri
  call write_stat('DEN.FILE',stat)

  if( keyq .ne. 0 ) then
    allocate( raqqq1_tri(ndq2), stat=alloc_stat )
    call allocation_stat('raqqq1_tri',alloc_stat)
    read(n_apden,iostat=stat) raqqq1_tri
    call write_stat('DEN.FILE',stat)
    allocate( raqqq2_tri(ndq2), stat=alloc_stat )
    call allocation_stat('raqqq2_tri',alloc_stat)
    read(n_apden,iostat=stat) raqqq2_tri
    call write_stat('DEN.FILE',stat)
  end if

close(n_apden)

end subroutine rdup_den
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine allocate_denominator

use legacy, only: ndd, ndq, keyq

implicit none

integer :: alloc_stat

if( inv_real ) then
  allocate( raqqd(ndd,ndd), stat=alloc_stat )
  call allocation_stat('raqqd',alloc_stat)
  if( keyq.ne.0 ) then
    allocate( raqqq(ndq,ndq), stat=alloc_stat )
    call allocation_stat('raqqq',alloc_stat)
  end if
else
  allocate( caqqd(ndd,ndd), stat=alloc_stat )
  call allocation_stat('caqqd',alloc_stat)
  if( keyq.ne.0 ) then
    allocate( caqqq(ndq,ndq), stat=alloc_stat )
    call allocation_stat('caqqq',alloc_stat)
  end if
end if

end subroutine allocate_denominator
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine read_denominator(enen)

use legacy, only: ndd2, ndq2, keyq

implicit none

integer, intent(in) :: enen

real*8, allocatable, dimension(:) :: iaqqd_tri, iaqqq_tri
integer :: stat

if( inv_real ) then
  allocate( raqqd_tri(ndd2,1) )
  read(dentap(enen),iostat=stat) raqqd_tri
  call write_stat('raqqd_tri',stat)
  call square_real_matrix(raqqd,raqqd_tri(:,1))
  deallocate( raqqd_tri )
  if( keyq.ne.0 ) then
    allocate( raqqq_tri(ndq2,1) )
    read(dentap(enen),iostat=stat) raqqq_tri
    call write_stat('raqqq_tri',stat)
    call square_real_matrix(raqqq,raqqq_tri(:,1))
    deallocate( raqqq_tri )
  end if
else
  allocate( raqqd_tri(ndd2,1) )
  allocate( iaqqd_tri(ndd2) )
  read(dentap(enen),iostat=stat) raqqd_tri
  call write_stat('raqqd_tri',stat)
  read(dentap(enen),iostat=stat) iaqqd_tri
  call write_stat('iaqqd_tri',stat)
  call square_complex_matrix(caqqd,raqqd_tri(:,1),iaqqd_tri)
  deallocate( raqqd_tri )
  deallocate( iaqqd_tri )
  if( keyq.ne.0 ) then
    allocate( raqqq_tri(ndq2,1) )
    allocate( iaqqq_tri(ndq2) ) 
    read(dentap(enen),iostat=stat) raqqq_tri
    call write_stat('raqqq_tri',stat)
    read(dentap(enen),iostat=stat) iaqqq_tri
    call write_stat('iaqqq_tri',stat)
    call square_complex_matrix(caqqq,raqqq_tri(:,1),iaqqq_tri)
    deallocate( raqqq_tri )
    deallocate( iaqqq_tri )
  end if
end if

end subroutine read_denominator
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine allocate_cb

use legacy, only: nchl
use quadsym, only: nop2

implicit none

integer :: alloc_stat

allocate( cb(nop2,nop2,svd_pts,nchl), stat=alloc_stat )
if(alloc_stat.ne.0) then
  write(6,"(2a)") 'Could not allocate variable cb'
  write(6,"(/a)") 'Thus, it will be allocated witohut the nchl dimension, &
                   and the disk will be accessed for each channel'
  allocate( cb(nop2,nop2,svd_pts,1), stat=alloc_stat )
  call allocation_stat('cb',alloc_stat)
  cb_inmem = .false.
end if

end subroutine allocate_cb
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine solve_linear_system(enen)

implicit none

integer, intent(in) :: enen

if( inv_real ) then
! Invert the real part of the denominator matrix
  if( invtyp .eq. 'LUF_KMT' ) then
    call inv_luf_kmt(enen)
  else if( invtyp .eq. 'SVD_KMT' ) then
    call inv_svd_kmt(enen)
  end if
  call t_matrix(enen)
else 
! Invert the complex denominator matrix
  if( invtyp .eq. 'LUF_LPK' ) then
    call inv_luf(enen)
  else if( invtyp .eq. 'SVD_LPK' ) then
    call inv_svd(enen)
  end if 
end if 

end subroutine solve_linear_system
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine off_driver

use legacy,  only: ndd, ndd2, ndq, ndq2, npct, rdup_orb, keyq
use quadsym, only: nok, nop, nop2, pw, xx, xkmax, kradtyp, wk, nokmax, ao2mo_tran 
use mplane1, only: setfn, abint

implicit none

integer :: i, kk, ich, iq, ik, ikk, ntape
integer :: iqmax, nktp(2), ktp(nokmax,2), kget(nokmax,2)
integer :: n_off, i_off, j_off, ich_eff
real*8  :: p(3) 
integer :: alloc_stat

!  Read up orbital information relevant for numerator matrix elements
call rdup_orb(1)

!  Load a few arrays and constants
call setfn

!  Evaluate which atomic to molecular orbital transformations are actually required
call ao2mo_tran

!  Allocate symmetry check array
allocate( ischk(nop*npct) )

if( io_off .ne. 'NDNO' ) then
  call set_nmf
end if

! Evaluate only for the channels where n_pts(ich).ne.0
call evaluate_effective_channels

!n_off = (nok(1)+nok(2)) * nchl
n_off = (nok(1)+nok(2)) * nchl_eff

if( io_off .eq. 'NDND' ) then
  if( off_inmem ) then
    allocate( offd_all(ndd2,n_off), stat=alloc_stat )
    call allocation_stat('offd_all',alloc_stat)
    if( keyq .ne. 0 ) then
      allocate( offq_all(ndq2,n_off), stat=alloc_stat )
      call allocation_stat('offq_all',alloc_stat)
    end if
  else
    allocate( offd_all(ndd2,1), stat=alloc_stat )
    call allocation_stat('offd_all',alloc_stat)
    if( keyq .ne. 0 ) then
      allocate( offq_all(ndq2,1), stat=alloc_stat )
      call allocation_stat('offq_all',alloc_stat)
    end if
  end if 
end if 

call abint

allocate( cvec1_off(ndd,nop2,nchl_eff), stat=alloc_stat )
call allocation_stat('cvec1_off',alloc_stat)
if( keyq .ne. 0 ) then
  allocate( cvec2_off(ndq,nop2,nchl_eff), stat=alloc_stat )
  call allocation_stat('cvec2_off',alloc_stat)
end if

kk = 0
pw = pw * xkmax / wk(1,1)
call get_tapes(chl_eff(1), iqmax, nktp, ktp, kget, nokmax)
do iq=1,iqmax
  do ik=1,nktp(iq)
    ikk   = kget(ik,iq)
    ischk = 0
    do i=1,nop
      if( ischk(i) .eq. 0 ) then
        kk=kk+1
        if(iq.eq.1) then
          p(:) = pw(:,i,1,1) * xx(ikk,iq)
        else
          if( kradtyp .eq. "LEGEN" ) then
            p(:) = pw(:,i,1,1) / xx(ikk,iq)
          else
            p(:) = pw(:,i,1,1) * xx(ikk,iq)
          end if
        end if
        if(keyplm.le.0) call plane(p, nplrtp(1,1))
        call fockpl(nplmtp(1,1))
        call elmnm1_off(i)
      end if
    end do 
    do ich_eff=1,nchl_eff
      ich = chl_eff(ich_eff)
      call get_tapes(ich, iqmax, nktp, ktp, kget, nokmax)
      ntape = ktp(ik,iq)
      j_off = (ich_eff-1)*(nok(1)+nok(2))
      if( io_off .eq. 'NDNO' ) then
        print*,'WRNUM',' ik=', ik, ' ikk=', ikk, ' iq=', iq, ' ich=', ich
        call wroff(cvec1_off(:,:,ich_eff),ndd,nop,ntape)
        if(keyq.ne.0) call wroff(cvec2_off(:,:,ich_eff),ndq,nop,ntape)
      else
        i_off = 1
        if( off_inmem ) i_off = (iq-1) * nok(1) + j_off + ik
        call smn(cvec1_off(:,:,ich_eff),offd_all(:,i_off),ndd)
        write(ntape) offd_all(:,i_off)
        if(keyq.ne.0) then
          call smn(cvec2_off(:,:,ich_eff),offq_all(:,i_off),ndq)
          write(ntape) offq_all(:,i_off)
        end if
        print *,'WRSMN',' ik=', ik, ' ikk=', ikk, ' iq=', iq, ' ich=', ich
      end if
    end do 
  end do 
end do 

deallocate( cvec1_off )
if( keyq .ne. 0 ) deallocate( cvec2_off )

deallocate( ischk )

write(6,'(//a,i0)') '** NUMBER OF PLANE WAVES COMPUTED = ', kk

end subroutine off_driver
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine evaluate_effective_channels

use legacy, only: nchls, nchl
use quadsym, only: n_pts

implicit none

integer :: ich

allocate( chl_eff(nchl) )
chl_eff = 0
nchl_eff = 0

do ich=1,nchl
  if(n_pts(ich).ne.0) then
    nchl_eff = nchl_eff + 1
    chl_eff(nchl_eff) = ich
  end if
end do

do_ground = .false.
nchlg_eff = 0
if(n_pts(1).ne.0) then
  do_ground = .true.
  nchlg_eff = 1
end if

do_loop_sing = .false.
nchls_eff = 0
do ich=2,nchls+1
  if(n_pts(ich).ne.0) then
    do_loop_sing = .true.
    nchls_eff = nchls_eff + 1
  end if
end do

do_loop_trip = .false.
nchlt_eff = 0
do ich=nchls+2,nchl
  if(n_pts(ich).ne.0) then
    do_loop_trip = .true.
    nchlt_eff = nchlt_eff + 1
  end if
end do

nchl_eff = nchlg_eff + nchls_eff + nchlt_eff

write(*,*)
write(*,'(a,i0,a)') 'Matrix elements will be computed for nchl_eff = ', nchl_eff, ' channels:'
write(*,'(a,i0)') 'Ground state:   ', nchlg_eff
write(*,'(a,i0)') 'Singlet states: ', nchls_eff
write(*,'(a,i0)') 'Triplet states: ', nchlt_eff
write(*,*)

end subroutine evaluate_effective_channels
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine plane(p, ntape)

use legacy, only: nbfns, nbmax, nbmax2, ntype, nfirst, nlast, none1, eta, d, vlist, lpseud, blkc1, blkc2, blkc3, blki1, blki2, non, n_prim 
use cintg_pp, only: cint
use mplane1, only: int2ns2, nqt, cdint, map_ab, map_cd
use mplane2, only: int1nh
use quadsym, only: n_tran1, tran1, n_tran2, tran2, n_tran3, tran3
!$    use omp_lib, only: OMP_get_thread_num, OMP_get_num_threads, OMP_get_wtime

implicit none

!=======================================================================
!     SUBROUTINE TO GIVE INTEGRALS OF WHICH INTEGRAND INCLUDES
!     PLANE WAVE(S)
!=======================================================================
integer, intent(in) :: ntape
real*8, intent(in)  :: p(3)

logical, save :: first_call = .true.
integer :: ntypa, ntypc!, ntypb !, ntypd
real*8  :: xa(3), alpa, xb(3) !, alpb
real*8  :: xc(3), alpc !, xd(3), alpd
integer :: i, m, n, mn, ifi, ila, mfi, mla, nla, nfi 
integer :: j, l, ik, mk, nk, ic, ii, iz, mm, nn, mn0
integer :: i0, n0, m0
integer :: kont1, kont2, kont3, nucl
real*8  :: dcrit, ee, eee, zval
real*8  :: bore(2), RO(2), sigma(3,3), aion(2,3,3)
complex*16 :: cx, cr, cn
complex*16, allocatable, dimension(:,:) :: cpr, cri
complex*16, allocatable, dimension(:) :: cri2, cpa
integer :: i_ab, i_cd, i_ik, i_mk, i_nk, nm_prim
integer :: alloc_stat

!   These are related to OPENMP
!$    integer :: n_id, n_threads
!c!$    real*8  :: t_start, t_end
!c!$    real*8  :: t_start2
!***

dcrit = 1.0d-10

!  Allocate AO-basis numerator, MO-basis numerator and transformation arrays
if( first_call ) then
  allocate( cprt(nbfns,nbmax2), stat=alloc_stat )
  call allocation_stat('cprt',alloc_stat)
  allocate( cpat(nbfns) )
  first_call = .false.
end if

!=======================================================================
!     P=(KX,KY,KZ) PLANE WAVE
!=======================================================================
!     NOTE: THIS PROGRAM USES THE ORIGINAL SUBROUTINES (THANKS ARE
!     DUE TO E.R.DAVIDSON AND T.RESCIGNO) WHICH DO NOT TAKE
!     ACCOUNT OF CONTRACTED BASIS SET EXPLICITLY. THEREFORE,
!     THIS PROGRAM CAN BE MADE FASTER BY THAT MODIFICATION.
!=======================================================================
!     if(keyplr.gt.0)then
!     call crd2(CPR,nbfns,max,nbfns,max,ntape)
!     go to 143
!     endif
!=======================================================================
!ONE PLANE WAVE + THREE GAUSSIANS (CPR)
!=======================================================================

! xd(:)=P(:)
call cdint(p)
!
! Memory requirements for calculation and transformation (complex*16 words):
! Nwords = 2*( nbfns * nbfns*(nbfns+1)/2 ) + nbfns**2 + nbfns
! Nwords = nbfns**3 + nbfns**2 + 2*nbfns
! Nuclear part: Nwords = nbfns**3 + nbfns**2 + 4*nbfns
! Memory(GB) = Nwords * 16 / 1024^3

allocate( cpr(nbfns,nbmax2), stat=alloc_stat )
call allocation_stat('cpr',alloc_stat)

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nfirst, nlast, nbfns, none1, eta, cpr, map_ab, map_cd, n_prim, dcrit)

!c!$OMP. FIRSTPRIVATE(xb)

!c!$OMP MASTER
!c!$    n_threads = OMP_get_num_threads()
!c!$    write(6,*) 'Rodando com n_threads = ',n_threads
!cc   print*,'size(cpr)',size(cpr)
!c!$OMP END MASTER

!     allocate( cri(nbfns,nbfns), cri2(nbfns) )

!c!$    t_start   = OMP_get_wtime()
!c!$    n_id= OMP_get_thread_num()

!$OMP DO
!     ab
do m=1,nbfns
  mfi=nfirst(m)
  mla=nlast(m)
  do n=1,m
    nfi=nfirst(n)
    nla=nlast(n)
    mn=none1(m)+n
!   cd
    do i=1,nbfns
      ifi=nfirst(i)
      ila=nlast(i)
      cx = (0d0,0d0)
      do mk=mfi,mla
        i_mk = mk - mfi + 1
        do nk=nfi,nla
          i_nk = nk - nfi + 1
          ee = eta(mk,5) * eta(nk,5)
          nm_prim = (i_mk-1) * n_prim + i_nk
          i_ab = map_ab(mn,nm_prim)
!         primitives
          do ik=ifi,ila
            i_ik = ik - ifi + 1
            i_cd = map_cd(i,i_ik)
            eee = ee * eta(ik,5)
            call int2ns2(cr,i_ab,i_cd)
            cx = cx + cr * eee
          end do
        end do
      end do
      if(dabs( real(cx)).lt.dcrit) cx = cmplx(0.0d0,aimag(cx))
      if(dabs(aimag(cx)).lt.dcrit) cx = cmplx(real(cx),0.0d0)
      cpr(i,mn) = cx
    end do
  end do
end do
!$OMP END DO
!$OMP END PARALLEL

!   if(keyplr.lt.0)then
!!  call cwrt2(CPR,nbfns,max,nbmax,max,ntape)
!   call cwrt2(CPR,nbfns,nbmax2,nbmax,nbmax2,ntape)
!   endif
!     143 continue
!=======================================================================
!TRANSFORMATION OF CPR
!=======================================================================

!     allocate( cri(nbfns,nbfns), cri2(nbfns) )
!c!$    t_start2  = OMP_get_wtime()
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nbfns, none1, cpr, d, cprt, n_tran1, tran1, n_tran2, tran2, n_tran3, tran3)
allocate( cri(nbfns,nbfns), cri2(nbfns) )
!$OMP DO 
do i0=1,n_tran1
  i = tran1(i0)
  do nn=1,nbfns
    do mm=1,nn
      mn0=none1(nn)+mm
      cri(mm,nn) = sum( cpr(:,mn0)*d(:,i) )
      cri(nn,mm) = cri(mm,nn)
    end do
  end do
  do m0=1,n_tran2(i)
    m = tran2(m0,i0)
    do nn=1,nbfns
      cri2(nn) = sum( cri(:,nn)*d(:,m) )
    end do
    do n0=1,n_tran3(m,i)
      n = tran3(n0,m0,i0)
      mn=none1(m)+n
      cprt(i,mn) = sum( cri2(:)*d(:,n) )
    end do
  end do
end do
!$OMP END DO

!c!$    t_end   = OMP_get_wtime()

!c!$OMP CRITICAL
!cc!$    write(6,*) 'Thread = ', n_id
!cc!$    write(6,*) 'Tempo (segundos) = ', t_end - t_start
!c!$    write(6,*) 'Tempo1(segundos) = ', t_end - t_start
!c!$    write(6,*) 'Tempo2(segundos) = ', t_end - t_start2
!c!$OMP END CRITICAL

!     deallocate( cri, cri2 )
deallocate( cri, cri2 )
!$OMP END PARALLEL
!     deallocate( cri, cri2 )

deallocate( cpr )

!=======================================================================
!NUCLEAR HYBRID INTEGRAL INCLUDING A PLANE WAVE
!=======================================================================
!     if(keyplr.gt.0)then
!     call crd1(cpa,nbmax,ntape)
!     go to 144
!     endif

xb(:) = p(:)

allocate( cpa(nbfns) )

do i=1,nbmax
  ntypa=nqt(ntype(i))
  ifi=nfirst(i)
  xa(1:3)=eta(ifi,1:3)
  ila=nlast(i)
  cr=(0.0D0,0.0D0)
  do ii=ifi,ila
    alpa=eta(ii,4)
    cx=(0.0D0,0.0D0)
    kont1=0
    kont2=0
    kont3=0
    do nucl=1,non
      xc(1:3)=vlist(nucl,1:3)
      if(lpseud(nucl).eq.1)then
!MB---------------------------------------------------------------------
!MB   BEGINS PSEUDOPOTENTIAL CALCULATION. THE INTEGRAL IS CINT
!MB---------------------------------------------------------------------
        zval=blkc1(nucl)
        do ic=1,2
          kont1=kont1+1
          bore(ic)=blkc2(kont1)
          RO(ic)=blkc3(kont1)
        end do
        do l=1,3
          do j=1,3
            kont2=kont2+1
            sigma(j,l)=blki1(kont2)
          end do
          do n=1,2
            do j=1,3
              kont3=kont3+1
              aion(n,j,l)=blki2(kont3)
            end do
          end do
        end do
        iz=idint(vlist(nucl,4))
        cn=cint(zval,bore,RO,aion,sigma,iz,xa,alpa,ntypa,xb,xc)
        cx=cx+cn
!MB---------------------------------------------------------------------
!MB   ENDS PSEUDOPOTENTIAL CALCULATION
!MB---------------------------------------------------------------------
      else
        call int1nh(cn,xa,alpa,ntypa,xb,xc,alpc,ntypc)
        cn=cn*vlist(nucl,4)
        cx=cx-cn
      end if
    end do
    cr=cr+cx*eta(ii,5)
  end do
  cpa(i)=cr
end do
!     if(keyplr.lt.0)then
!     call cwrt1(cpa,nbmax,ntape)
!     endif
!     144 continue
!=======================================================================

!=======================================================================
!TRANSFORMATION OF CPA
!=======================================================================
do i=1,nbfns
  cpat(i) = sum( cpa(:)*d(:,i) )
end do
!=======================================================================

deallocate( cpa )

end subroutine plane
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine fockpl(ntape)

use legacy, only: nocc, nbfns, nbmax, nbmax2, none1

implicit none

integer, intent(in) :: ntape

logical, save :: first_call = .true.
complex*16 :: cx, cy
integer :: j, m

!  Allocate cg and cg2
!
if( first_call ) then
  allocate( cg(nbfns), cg2(nbfns) )
  first_call = .false.
end if

if(keyplm.gt.0) then
  call crd1(cpat,nbmax,ntape)
  call crd2(cprt,nbfns,nbmax2,nbmax,nbmax2,ntape)
endif
do j=1,nbmax
  cx = (0.0d0,0.0d0)
  cy = (0.0d0,0.0d0)
  do m=1,nocc
    cx = cx + cprt(j,none1(m)+m)
    cy = cy + cprt(m,none1(max0(j,m))+min0(j,m))
  end do
  cg2(j) = 2d0*cx + cpat(j)
  cg(j) = cg2(j) - cy
end do
if(keyplm.lt.0) then
  call cwrt1(cpat,nbmax,ntape)
  call cwrt2(cprt,nbfns,nbmax2,nbmax,nbmax2,ntape)
endif

end subroutine fockpl
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine elmnm1_s(ip,ich)

use legacy, only: nctrin, ng, ngt, ns, nht, npt, ng2, ngt2, nsa, hci, numsd, none1

implicit none

integer, intent(in) :: ip, ich 

integer :: na, i, ii, i1
integer :: j, jj, ij1, ij2, l
integer :: ic1, ic1d, ic2, ic2d, ichd
integer :: nb
complex*16 :: c1
complex*16 :: cr(nsa)

!=======================================================================
!  MATRIX ELEMENTS FOR THE NUMERATOR (1). PRELIMINARY.
!=======================================================================

ichd=ich-1
 
ic1=nctrin(1,ich)
ic1d=ic1-1
ic2=nctrin(2,ich)
ic2d=ic2-1
 
if(ich.ge.2) go to 177
!
!=======================================================================
!<0>(I) V <0><PLANE>
!=======================================================================
do ii=1,ng
  i=ngt(ii)
  cvec1(ii,ip)= cg(i)
end do

if(ns.eq.0) go to 950
!
!=======================================================================
!<NA-I>(J) V <0><PLANE>
!=======================================================================
do ii=1,ns
  na = nht(ii)
  i  = npt(ii)
  do jj=1,ng2(ii)
    ij1=numsd(jj,ii,1)
    ij2=numsd(jj,ii,2)
    if( (ij1.ne.0).and.(ij2.ne.0) ) then
      J= ngt2(jj,ii)
      c1 = cprt(i,none1(j)+na)
      cvec1(ij1,ip)= ( - c1 + 2d0 * cprt(j,none1(i)+na) ) * rrt2
      cvec1(ij2,ip)= c1 * rrt6t3
    else if( (ij1.ne.0).and.(ij2.eq.0) ) then
      J= ngt2(jj,ii)
      cvec1(ij1,ip)= ( - cprt(i,none1(j)+na) + 2d0 * cprt(j,none1(i)+na) ) * rrt2
    else if( (ij1.eq.0).and.(ij2.ne.0) ) then
      J= ngt2(jj,ii)
      cvec1(ij2,ip)= cprt(i,none1(j)+na) * rrt6t3
    end if
  end do
end do

go to 950
!
  177 continue
!
!=======================================================================
!  CASE OF ich.ne.1;  <0>(I) QV <NA-K>(PLANE)
!=======================================================================
!-----------------------------------------------------------------------
!     Dublet subspace run
!-----------------------------------------------------------------------
!=======================================================================
!<0>(I) V <NA-K><PLANE>
!=======================================================================

do i=1,ng
  ii=ngt(i)
  do l=1,nsa
    na = nht(l)
    i1 = npt(l)
    cr(l) = 2d0 * cprt(ii,none1(i1)+na) - cprt(na,none1(max0(i1,ii))+min0(i1,ii))
    if(i1.eq.ii) cr(l) = cr(l) - cg(na)
  end do
  cvec1(i,ip)= sum( hci(:,ichd,1) * cr(:) ) * rrt2
end do

!=======================================================================
!<NB-I>(J) V <NA-K><PLANE>
!=======================================================================
 
!c!$OMP   PARALLEL DO DEFAULT(PRIVATE)
!c       SHARED(ns,nht,npt,ng2,numsd,ngt2,cvec1,hci,ip,ichd,rt3d2)
!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ii, nb, i, jj, ij1, ij2, j, cr)
do ii=1,ns
  nb=nht(ii)
  i=npt(ii)
  do jj=1,ng2(ii)
    ij1=numsd(jj,ii,1)
    ij2=numsd(jj,ii,2)
    if( (ij1.ne.0).and.(ij2.ne.0) ) then
      j = ngt2(jj,ii)
      call pvss2_b4_s1(nb,i,j,cr)
      cvec1(ij1,ip)= sum( hci(:,ichd,1)*cr(:) ) * 0.5d0
      call pvss2_b4_s2(nb,i,j,cr)
      cvec1(ij2,ip)= sum( hci(:,ichd,1)*cr(:) ) * rt3d2
    else if( (ij1.ne.0).and.(ij2.eq.0) ) then
      j = ngt2(jj,ii)
      call pvss2_b4_s1(nb,i,j,cr)
      cvec1(ij1,ip)= sum( hci(:,ichd,1)*cr(:) ) * 0.5d0
    else if( (ij1.eq.0).and.(ij2.ne.0) ) then
      j = ngt2(jj,ii)
      call pvss2_b4_s2(nb,i,j,cr)
      cvec1(ij2,ip)= sum( hci(:,ichd,1)*cr(:) ) * rt3d2
    end if
  end do
end do
!$OMP END PARALLEL DO

  950 continue

call cvec_symmetry(ip,ich,ic2)

end subroutine elmnm1_s
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine elmnm1_t(ip,ich)

use legacy, only: nchls, nctrin, ng, ngt, ns, nht, npt, ng2, ngt2, nsa, hci, nsp, numsd, numsq, keyq, none1

implicit none

integer, intent(in) :: ip, ich 

integer :: na, i, ii, i1, na1, na2, na3
integer :: j, jj, ij1, ij2, k1, k2, k3, l
integer :: ip1, ic1, ic1d, ic2, ic2d, ichd
integer :: is, js, nb, nbs, isp, ij3
complex*16 :: c1, c2, c3
complex*16 :: cr(nsa)
complex*16 :: crr1(nsa,3), crr2(nsa,3)

!=======================================================================
!  MATRIX ELEMENTS FOR THE NUMERATOR (1). PRELIMINARY.
!=======================================================================

ichd=ich-nchls-1

ic1=nctrin(1,ich)
ic1d=ic1-1
ic2=nctrin(2,ich)
ic2d=ic2-1
 
!=======================================================================
!CASE OF ich.ne.1;  <0>(I) QV <NA-K>(PLANE)
!=======================================================================
!-----------------------------------------------------------------------
!     Dublet subspace run
!-----------------------------------------------------------------------
!=======================================================================
!<0>(I) V <NA-K><PLANE>
!=======================================================================
do i=1,ng
  ii=ngt(i)
  do l=1,nsa
    na = nht(l)
    i1 = npt(l)
    cr(l) = cprt(na,none1(max0(i1,ii))+min0(i1,ii))
    if(i1.eq.ii) cr(l) = cr(l) + cg(na)
  end do
  cvec1(i,ip)= sum( hci(:,ichd,2) * cr(:) ) * rrt6t3
end do

!c!$OMP   PARALLEL DO DEFAULT(PRIVATE) 
!c!       SHARED(ns,nht,npt,ng2,numsd,ngt2,cvec1,hci,ip,ichd,rt3d2)
!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ii, nb, i, jj, ij1, ij2, j, cr)
do ii=1,ns
  nb=nht(ii)
  i=npt(ii)
  do jj=1,ng2(ii)
    ij1=numsd(jj,ii,1)
    ij2=numsd(jj,ii,2)
    if( (ij1.ne.0).and.(ij2.ne.0) ) then
      j = ngt2(jj,ii)
      call pvss2_b4_t1(nb,i,j,cr)
      cvec1(ij1,ip)= sum( hci(:,ichd,2)*cr(:) ) * rt3d2
      call pvss2_b4_t2(nb,i,j,cr)
      cvec1(ij2,ip)= sum( hci(:,ichd,2)*cr(:) ) * 0.5d0
    else if( (ij1.ne.0).and.(ij2.eq.0) ) then
      j = ngt2(jj,ii)
      call pvss2_b4_t1(nb,i,j,cr)
      cvec1(ij1,ip)= sum( hci(:,ichd,2)*cr(:) ) * rt3d2
    else if( (ij1.eq.0).and.(ij2.ne.0) ) then
      j = ngt2(jj,ii)
      call pvss2_b4_t2(nb,i,j,cr)
      cvec1(ij2,ip)= sum( hci(:,ichd,2)*cr(:) ) * 0.5d0
    end if
  end do
end do
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!     Quartet subspace run
!-----------------------------------------------------------------------
if(keyq.eq.0 .or. ic2.eq.0) go to 950

ip1=1
do II=1,ns
  NB=nht(II)
  I=npt(II)
  do JJ=1,ng2(II)
    ij3=numsq(JJ,II)
    if(ij3.ne.0) then
      J=ngt2(JJ,II)
      do L=1,nsa
        NA1=nht(L)
        K1=npt(L)
        NA2=-nht(L)
        K2=-npt(L)
        NA3=-nht(L)
        K3=npt(L)
        do isp=1,3
          NBS=NB*nsp(isp,1)
          IS=I*nsp(isp,2)
          JS=J*nsp(isp,3)
          call pvss2(NA1,K1,ip1,NBS,IS,JS,c1)
          call pvss2(NA2,K2,ip1,NBS,IS,JS,c2)
          call pvss2(NA3,K3,-ip1,NBS,IS,JS,c3)
          crr1(l,isp)=C1-C2
          crr2(l,isp)=C3
        end do
        cvec2(ij3,ip)= sum( hci(:,ichd,2)*(crr1(:,1)-crr1(:,2)-crr1(:,3)-(crr2(:,1)-crr2(:,2))+crr2(:,3)) ) / 3.0d0
      end do
    end if
  end do
end do

  950 continue

call cvec_symmetry(ip,ich,ic2)

end subroutine elmnm1_t
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine cvec_symmetry(ip,ich,ic2)

use legacy, only: ng, ngt, ns, nht, npt, ng2, ngt2, numsd, numsq, keyq
use quadsym, only: ntrn, ma, map, ngtinv, numinv, ngtin2, machld 

implicit none

integer, intent(in) :: ip, ich, ic2

integer :: m, l, nat, nata
integer :: k, kt
integer :: ktr, ipt, ipta
integer :: kss, no, ichld
integer :: kl, klt, kl2, klt2, kl3, klt3, ks, kta
integer :: kkt, kkta, lt, lta, ls

!=======================================================================
! SYMMETRY CHECK
! PROJECTION OPERATOR P IS ASSUMED TO BE TOTALLY SYMMETRIC.
!=======================================================================
ischk(ip)=1
if(ntrn.eq.0) return
 
do m=1,ntrn

  ktr=1
  ipt=map(ip,m)
  ipta=iabs(ipt)
  if(ischk(ipta).eq.1) cycle
  if(ich.ne.1) then
    ichld=ich-1
    if(machld(ichld,m).lt.0) then
      ktr=-ktr
    elseif(machld(ichld,m).gt.0) then
      ktr=ktr
    else
      return
    end if
  end if
  if(ipt.lt.0) then
    kss=-ktr
  elseif(ipt.gt.0) then
    kss=ktr
  else
    return
  end if

  do k=1,ng
    kt=ma(ngt(k),m)
    if(kt.eq.0) return
    kta=ngtinv(iabs(kt))
!=======================================================================
!  IF kt=0, then kta=ngtinv(0), WHICH GIVES A SUBSCRIPT OUT OF RANGE
!  ERROR...THIS IS NOW AVOIDED BY SKIPPING THIS STATEMENT WHEN
!  kt=0.  TLG/MAL  5/5/83.
!=======================================================================
    if(kt.lt.0) then
      ks=-kss
    elseif(kt.gt.0) then
      ks=kss
    else
      return
    end if
    cvec1(kta,ipta)=cvec1(k,ip)*ks
!roma 31/05/2003  300 if(ic2.ne.0) CVEC2(kta,ipta)=CVEC2(K,IP)*ks
  end do
  if(ns.ne.0) then
    do k=1,ns
      nat=ma(nht(k),m)
      nata=iabs(nat)
      kkt=ma(npt(k),m)
      kkta=iabs(kkt)
      if(nat*kkt.lt.0) then
        ks=-kss
      elseif(nat*kkt.gt.0) then
        ks=kss
      else
        return
      end if
      no=numinv(nata,kkta)
      loop_ng2: do l=1,ng2(k)
        lt=ma(ngt2(l,k),m)
        lta=ngtin2(iabs(lt),no)
        if(lt.lt.0) then
          ls=-ks
        elseif(lt.gt.0) then
          ls=ks
        else
          return
        end if
        kl=numsd(l,k,1)
        if(kl.ne.0) then
          klt=numsd(lta,no,1)
          if(klt.ne.0) then
            cvec1(klt,ipta)=cvec1(kl,ip)*ls
          end if
        end if
        kl2=numsd(l,k,2)
        if(kl2.ne.0) then
          klt2=numsd(lta,no,2)
          if(klt2.ne.0) then
            cvec1(klt2,ipta)=cvec1(kl2,ip)*ls
          end if
        end if
        if(keyq.eq.0 .or. ic2.eq.0) cycle loop_ng2
        kl3=numsq(L,K)
        if(kl3.eq.0) cycle loop_ng2
        klt3=numsq(lta,no)
        if(klt3.eq.0) cycle loop_ng2
        cvec2(klt3,ipta)=cvec2(kl3,ip)*ls
        end do loop_ng2
      end do
    end if
  ischk(ipta)=ischk(ipta)+1

end do

end subroutine cvec_symmetry
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine elmnm1_off(ip)

use legacy, only: ng, ngt, ns, nht, npt, ng2, ngt2, nsa, hci, nsp, numsd, numsq, keyq, none1, nchls, nchl, ndd

implicit none

integer, intent(in) :: ip

integer :: ich_eff, ich, ichd, l, na, i, ii, i1
integer :: j, jj, ij1, ij2, k1, k2, k3, na1, na2, na3
integer :: ip1, is, js, nb, nbs, isp, ij3
complex*16 :: c1, c2, c3 
complex*16 :: cr(nsa,ndd)
complex*16 :: crr1(nsa,3), crr2(nsa,3)

!=======================================================================
! MATRIX ELEMENTS FOR THE NUMERATOR (1). PRELIMINARY.
!=======================================================================
!-----------------------------------------------------------------------
!     Dublet subspace run
!-----------------------------------------------------------------------
!=======================================================================
! <0>(I) V <0><PLANE>
!=======================================================================
if(do_ground) then
  ich_eff = chl_eff(1)

  do ii=1,ng
    i=ngt(ii)
    cvec1_off(ii,ip,ich_eff)= cg(i)
  end do

  if(ns.ne.0) then
!=======================================================================
! <NA-I>(J) V <0><PLANE>
!=======================================================================
    !c!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ii,na,i,jj,ij1,ij2,j)
    do ii=1,ns
      na = nht(ii)
      i  = npt(ii)
      do jj=1,ng2(ii)
        ij1=numsd(jj,ii,1)
        ij2=numsd(jj,ii,2)
        if( (ij1.ne.0).and.(ij2.ne.0) ) then
          j= ngt2(jj,ii)
          c1 = cprt(i,none1(j)+na)
          cvec1_off(ij1,ip,ich_eff)= ( - c1 + 2d0 * cprt(j,none1(i)+na) ) *rrt2
          cvec1_off(ij2,ip,ich_eff)= c1 * rrt6t3
        else if( (ij1.ne.0).and.(ij2.eq.0) ) then
          j= ngt2(jj,ii)
          cvec1_off(ij1,ip,ich_eff)= ( - cprt(i,none1(j)+na) + 2d0 * cprt(j,none1(i)+na) ) * rrt2
        else if( (ij1.eq.0).and.(ij2.ne.0) ) then
          j= ngt2(jj,ii)
          cvec1_off(ij2,ip,ich_eff)= cprt(i,none1(j)+na) * rrt6t3
        end if
      end do
    end do
    !c!$OMP END PARALLEL DO
  end if

end if

if( nchl.eq.1 .and. chl_eff(1).eq.1 ) go to 950

if(do_loop_sing) then

!=======================================================================
! <0>(I) V <NA-K><PLANE>
!=======================================================================
  do i=1,ng
    ii=ngt(i)
     do l=1,nsa
      na = nht(l)
      i1 = npt(l)
      cr(l,i) = 2d0 * cprt(ii,none1(i1)+na) - cprt(na,none1(max0(i1,ii))+min0(i1,ii))
      if(i1.eq.ii) cr(l,i) = cr(l,i) - cg(na)
    end do
  end do

!=======================================================================
! <NB-I>(J) V <NA-K><PLANE>
!=======================================================================

  !c!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ii,nb,i,jj,ij1,ij2,j)
  do ii=1,ns
    nb = nht(ii)
    i = npt(ii)
    do jj=1,ng2(ii)
      ij1 = numsd(jj,ii,1)
      ij2 = numsd(jj,ii,2)
      if( (ij1.ne.0).and.(ij2.ne.0) ) then
        j = ngt2(jj,ii)
        call pvss2_b4_s1(nb,i,j,cr(:,ij1))
        call pvss2_b4_s2(nb,i,j,cr(:,ij2))
      else if( (ij1.ne.0).and.(ij2.eq.0) ) then
        j = ngt2(jj,ii)
        call pvss2_b4_s1(nb,i,j,cr(:,ij1))
      else if( (ij1.eq.0).and.(ij2.ne.0) ) then
        j = ngt2(jj,ii)
        call pvss2_b4_s2(nb,i,j,cr(:,ij2))
      end if
    end do
  end do
  !c!$OMP END PARALLEL DO

  !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ich_eff, ich, ichd, i, ii, jj, ij1, ij2) 
  do ich_eff=nchlg_eff+1,nchlg_eff+nchls_eff
    ich = chl_eff(ich_eff)
    ichd = ich - 1
    do i=1,ng
      cvec1_off(i,ip,ich_eff)= sum(hci(:,ichd,1)*cr(:,i)) * rrt2
    end do
    do ii=1,ns
      do jj=1,ng2(ii)
        ij1 = numsd(jj,ii,1)
        ij2 = numsd(jj,ii,2)
        if( (ij1.ne.0).and.(ij2.ne.0) ) then
          cvec1_off(ij1,ip,ich_eff)=sum(hci(:,ichd,1)*cr(:,ij1))*0.5d0
          cvec1_off(ij2,ip,ich_eff)=sum(hci(:,ichd,1)*cr(:,ij2))*rt3d2
        else if( (ij1.ne.0).and.(ij2.eq.0) ) then
          cvec1_off(ij1,ip,ich_eff)=sum(hci(:,ichd,1)*cr(:,ij1))*0.5d0
        else if( (ij1.eq.0).and.(ij2.ne.0) ) then
          cvec1_off(ij2,ip,ich_eff)=sum(hci(:,ichd,1)*cr(:,ij2))*rt3d2
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

end if


if(do_loop_trip) then

  do i=1,ng
    ii=ngt(i)
    do l=1,nsa
      na = nht(l)
      i1 = npt(l)
      cr(l,i) = cprt(na,none1(max0(i1,ii))+min0(i1,ii))
      if(i1.eq.ii) cr(l,i) = cr(l,i) + cg(na)
    end do
  end do

  !c!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ii,nb,i,jj,ij1,ij2,j)
  do ii=1,ns
    nb = nht(ii)
    i = npt(ii)
    do jj=1,ng2(ii)
      ij1 = numsd(jj,ii,1)
      ij2 = numsd(jj,ii,2)
      if( (ij1.ne.0).and.(ij2.ne.0) ) then
        j = ngt2(jj,ii)
        call pvss2_b4_t1(nb,i,j,cr(:,ij1))
        call pvss2_b4_t2(nb,i,j,cr(:,ij2))
      else if( (ij1.ne.0).and.(ij2.eq.0) ) then
        j = ngt2(jj,ii)
        call pvss2_b4_t1(nb,i,j,cr(:,ij1))
      else if( (ij1.eq.0).and.(ij2.ne.0) ) then
        j = ngt2(jj,ii)
        call pvss2_b4_t2(nb,i,j,cr(:,ij2))
      end if
    end do
  end do
  !c!$OMP END PARALLEL DO

  !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ich_eff, ich, ichd, i, ii, jj, ij1, ij2) 
  do ich_eff=nchlg_eff+nchls_eff+1,nchl_eff
    ich = chl_eff(ich_eff)
    ichd = ich - nchls - 1
    do i=1,ng
      cvec1_off(i,ip,ich_eff)= sum(hci(:,ichd,2)*cr(:,i)) * rrt6t3
    end do
    do ii=1,ns
      do jj=1,ng2(ii)
        ij1 = numsd(jj,ii,1)
        ij2 = numsd(jj,ii,2)
        if( (ij1.ne.0).and.(ij2.ne.0) ) then
          cvec1_off(ij1,ip,ich_eff)=sum(hci(:,ichd,2)*cr(:,ij1))*rt3d2
          cvec1_off(ij2,ip,ich_eff)=sum(hci(:,ichd,2)*cr(:,ij2))*0.5d0
        else if( (ij1.ne.0).and.(ij2.eq.0) ) then
          cvec1_off(ij1,ip,ich_eff)=sum(hci(:,ichd,2)*cr(:,ij1))*rt3d2
        else if( (ij1.eq.0).and.(ij2.ne.0) ) then
          cvec1_off(ij2,ip,ich_eff)=sum(hci(:,ichd,2)*cr(:,ij2))*0.5d0
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

end if

!-----------------------------------------------------------------------
!     Quartet subspace run
!-----------------------------------------------------------------------
if(keyq.ne.0) then
  ip1=1
  do ii=1,ns
    NB=nht(II)
    I=npt(II)
    do jj=1,ng2(ii)
      ij3=numsq(JJ,II)
      if(ij3.ne.0) then
        j=ngt2(jj,ii)
        do l=1,nsa
          NA1=nht(L)
          K1=npt(L)
          NA2=-nht(L)
          K2=-npt(L)
          NA3=-nht(L)
          K3=npt(L)
          do isp=1,3
            NBS=NB*nsp(isp,1)
            IS=I*nsp(isp,2)
            JS=J*nsp(isp,3)
            call pvss2(NA1,K1,ip1,NBS,IS,JS,c1)
            call pvss2(NA2,K2,ip1,NBS,IS,JS,c2)
            call pvss2(NA3,K3,-ip1,NBS,IS,JS,c3)
            crr1(l,isp)=C1+C2
            crr2(l,isp)=C2
          end do
          do ich_eff=nchlg_eff+1,nchlg_eff+nchls_eff
            ich = chl_eff(ich_eff)
            ichd = ich - 1
            cvec2_off(ij3,ip,ich_eff)= sum( hci(:,ichd,1)*(crr1(:,1)-crr1(:,2)-crr1(:,3)-(crr2(:,1)-crr2(:,2))+crr2(:,3)) ) / 3.0d0
          end do
          do isp=1,3
            NBS=NB*nsp(isp,1)
            IS=I*nsp(isp,2)
            JS=J*nsp(isp,3)
            call pvss2(NA1,K1,ip1,NBS,IS,JS,c1)
            call pvss2(NA2,K2,ip1,NBS,IS,JS,c2)
            call pvss2(NA3,K3,-ip1,NBS,IS,JS,c3)
            crr1(l,isp)=C1-C2
            crr2(l,isp)=C3
          end do
          do ich_eff=nchlg_eff+nchls_eff+1,nchl_eff
            ich = chl_eff(ich_eff)
            ichd = ich - nchls - 1
            cvec2_off(ij3,ip,ich_eff)= sum( hci(:,ichd,2)*(crr1(:,1)-crr1(:,2)-crr1(:,3)-(crr2(:,1)-crr2(:,2))+crr2(:,3)) ) / 3.0d0
          end do
        end do
      end if
    end do
  end do
end if

950 continue

call cvec_off_symmetry(ip)

end subroutine elmnm1_off
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine cvec_off_symmetry(ip)

use legacy, only: ng, ngt, ns, nht, npt, ng2, ngt2, numsd, numsq, keyq
use quadsym, only: ntrn, ma, map, ngtinv, numinv, ngtin2, machld 

implicit none

integer, intent(in) :: ip

integer :: ich, ich_eff
integer :: m, l, nat, nata
integer :: k, kt
integer :: ktr, ipt, ipta
integer :: kss, no, ichld
integer :: kl, klt, kl2, klt2, kl3, klt3, ks, kta
integer :: kkt, kkta, lt, lta, ls

!=======================================================================
! SYMMETRY CHECK
! PROJECTION OPERATOR P IS ASSUMED TO BE TOTALLY SYMMETRIC.
!=======================================================================

ischk(ip)=1
if(ntrn.eq.0) return
 
do m=1,ntrn
  ktr=1
  ipt=map(ip,m)
  ipta=iabs(ipt)
  if(ischk(ipta).eq.1) cycle
  if(ipt.lt.0) then
    kss=-ktr
  elseif(ipt.gt.0) then
    kss=ktr
  else 
    return
  end if

  do k=1,ng
    kt=ma(ngt(k),m)
    if(kt.eq.0) return
!=======================================================================
!  IF KT=0, then KTA=ngtinv(0), WHICH GIVES A SUBSCRIPT OUT OF RANGE
!  ERROR...THIS IS NOW AVOIDED BY SKIPPING THIS STATEMENT WHEN
!  KT=0.  TLG/MAL  5/5/83.
!=======================================================================
    kta=ngtinv(iabs(kt))
    if(kt.lt.0) then
      ks=-kss
    elseif(kt.gt.0) then
      ks=kss
    else 
      return
    end if
    loop_chl_0: do ich_eff=1,nchl_eff
      ich = chl_eff(ich_eff)
      if( ich.eq.1 ) then
        cvec1_off(kta,ipta,ich_eff)=cvec1_off(k,ip,ich_eff)*ks
        cycle loop_chl_0
      end if
      ichld=ich-1
      if( machld(ichld,m).ne.0 ) then
        cvec1_off(kta,ipta,ich_eff)=cvec1_off(k,ip,ich_eff)*ks*machld(ichld,m)
      end if
    end do loop_chl_0
!roma 31/05/2003  300 if(ic2.ne.0) CVEC2(KTA,ipta)=CVEC2(K,ip)*KS
  end do

  if(ns.ne.0) then

  do k=1,ns
    nat=ma(nht(k),m)
    nata=iabs(nat)
    kkt=ma(npt(k),m)
    kkta=iabs(kkt)
    if((nat*kkt).lt.0) then
      ks=-kss
    elseif((nat*kkt).gt.0) then
      ks=kss
    else 
      return
    end if
    no=numinv(nata,kkta)
    loop_ng2: do l=1,ng2(k)
      lt=ma(ngt2(l,k),m)
      lta=ngtin2(iabs(lt),no)
      if(lt.lt.0) then
        ls=-ks
      elseif(lt.gt.0) then
        ls=ks
      else 
        return
      end if
      kl=numsd(l,k,1)
      if(kl.ne.0) then
        klt=numsd(lta,no,1)
        if(klt.ne.0) then
          loop_chl_1: do ich_eff=1,nchl_eff
            ich = chl_eff(ich_eff)
            if( ich.eq.1 ) then
              cvec1_off(klt,ipta,ich_eff)=cvec1_off(kl,ip,ich_eff)*ls
              cycle loop_chl_1
            end if
            ichld=ich-1
            if( machld(ichld,m).ne.0 ) then
              cvec1_off(klt,ipta,ich_eff)=cvec1_off(kl,ip,ich_eff)*ls*machld(ichld,m)
            end if
          end do loop_chl_1
        end if
      end if
      kl2=numsd(l,k,2)
      if(kl2.ne.0) then
        klt2=numsd(lta,no,2)
        if(klt2.ne.0) then
          loop_chl_2: do ich_eff=1,nchl_eff
            ich = chl_eff(ich_eff)
            if( ich.eq.1 ) then
              cvec1_off(klt2,ipta,ich_eff)=cvec1_off(kl2,ip,ich_eff)*ls
              cycle loop_chl_2
            end if
            ichld=ich-1
            if( machld(ichld,m).ne.0 ) then
              cvec1_off(klt2,ipta,ich_eff)=cvec1_off(kl2,ip,ich_eff)*ls*machld(ichld,m)
            end if
          end do loop_chl_2
        end if
      end if
      if(keyq.eq.0) cycle loop_ng2
      kl3=numsq(l,k)
      if(kl3.eq.0) cycle loop_ng2
      klt3=numsq(lta,no)
      if(klt3.eq.0) cycle loop_ng2
      loop_chl_3: do ich_eff=1,nchl_eff
        ich = chl_eff(ich_eff)
        if( ich.eq.1 ) then
          cvec2_off(klt3,ipta,ich_eff)=cvec2_off(kl3,ip,ich_eff)*ls
          cycle loop_chl_3
        end if
        ichld=ich-1
        if( machld(ichld,m).ne.0 ) then
          cvec2_off(klt3,ipta,ich_eff)=cvec2_off(kl3,ip,ich_eff)*ls*machld(ichld,m)
        end if
      end do loop_chl_3
    end do loop_ng2
  end do
 
  end if

  ischk(ipta)=ischk(ipta)+1

end do

end subroutine cvec_off_symmetry
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine set_nmf

use legacy, only: ndd

implicit none

integer :: i, j

allocate( nmf(ndd,ndd) )
do i = 1,ndd
  do j = 1,i
    nmf(j,i) = i*(i-1)/2 + j
  end do
end do

end subroutine set_nmf
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine vkgkv3_off_deninmem

use legacy,  only: ndd, ndd2, ndq, ndq2, keyq, nchl, nchls 
use quadsym, only: wk, xx, ww1, nok, kradtyp, xkmax 

implicit none

integer :: ich, ik, ii, l
integer(i_large) :: m 
real*8  :: w, x
real*8  :: xk0(nener), faca1(nener), faca2(nener) 
real*8  :: xf1(nener), xf2(nener), fact2(nener) 
real*8, allocatable, dimension(:) :: offd, offq
integer :: stat, alloc_stat

!***
! The real part (principal value) will be caculated first, and then the residue.
! This makes things a little less confusing regarding the inversion of either
! complex or real denominator matrices
!***

!***
! In this version, NOTHT, NOPHI, NLEB, WW, WW2, WLEB, NOP and NOP2 are
! relevant for the on-shell matrix elements only. In case the off-shell
! matrix elements have been saved in NDNO (CVEC) form, the angular
! integration will be carried out in SMN_OFF with an independent
! quadrature for each radial point.
!***

!***
! The outer radial integration is carried out in either Gauss-Legendre
! or Gauss-Laguerre quadratures.
!***

allocate( offd(ndd2), stat=alloc_stat )
call allocation_stat('offd',alloc_stat)
if( keyq .ne. 0 ) then
  allocate( offq(ndq2), stat=alloc_stat )
  call allocation_stat('offq',alloc_stat)
end If

write(*,"(a)") 'VKGKV3MOD - BEGINS  '

do ich=1,nchl

  do l=1,nener
    xk0(l) = wk(ich,l)
    faca1(l) = xkmax / xk0(l)
    faca2(l) = xk0(l) / xkmax
  end do

  do ik=1,nok(1)
    ii=kkktap(ik,ich,1)
    x=xx(ik,1)
    w=ww1(ik,1)
    xf1(:) = (faca1(:)*x)**2
    if( io_off .eq. 'NDNO' ) then
      call smn_off(offd,ndd,ii)
    else
      read(ii,iostat=stat) offd
      call write_stat('OFF-SHELL',stat)
    end if
    fact2(:) = - xkmax * w * fact * xf1(:) / ( 1.0d0 - xf1(:) )
    !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(l, m)
    !c!$OMP   PARALLEL DO DEFAULT(PRIVATE) SHARED(nener,ndd2,raqqd_tri,fact2,offd)
    do l=1,nener
      do m=1,ndd2
        raqqd_tri(m,l) = raqqd_tri(m,l) + fact2(l) * offd(m)
      end do
    end do
    !$OMP   END PARALLEL DO
    if( (ich.gt.nchls).and.(keyq.ne.0) ) then
      if( io_off .eq. 'NDNO' ) then
        call smn_off(offq,ndq,ii)
      else
        read(ii,iostat=stat) offq
        call write_stat('OFF-SHELL',stat)
      end if
      !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(l, m)
      !c!$OMP   PARALLEL DO DEFAULT(PRIVATE) SHARED(nener,ndq2,raqqq_tri,fact2,offq)
      do l=1,nener
        do m=1,ndq2
          raqqq_tri(m,l) = raqqq_tri(m,l) + fact2(l) * offq(m)
        end do
      end do
      !$OMP   END PARALLEL DO
    end if 
    write(*,"(a,i4,a,i3)") 'Channel ', ich, ', iq = 1, ik = ', ik
  end do ! ik=1,nok(1)

  do ik=1,nok(2)
    ii=kkktap(ik,ich,2)
    x=xx(ik,2)
    w=ww1(ik,2)
    if( kradtyp .eq. "LEGEN" ) then
      xf2(:) = (faca2(:)*x)**2
      fact2(:) = ( xk0(:)**2/xkmax ) * w * fact / xf2(:) / ( 1.0d0 - xf2(:) )
    else
      xf2(:) = (faca1(:)*x)**2
      fact2(:) = xkmax * w * fact * xf2(:) / ( xf2(:) - 1.0d0 )
    end if
    if( io_off .eq. 'NDNO' ) then
      call smn_off(offd,ndd,ii)
    else
      read(ii,iostat=stat) offd
      call write_stat('OFF-SHELL',stat)
    end if
    !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(l, m)
    !c!$OMP   PARALLEL DO DEFAULT(PRIVATE) SHARED(nener,ndd2,raqqd_tri,fact2,offd)
    do l=1,nener
      do m=1,ndd2
        raqqd_tri(m,l) = raqqd_tri(m,l) + fact2(l) * offd(m)
      end do
    end do
    !$OMP   END PARALLEL DO
    if( (ich.gt.nchls).and.(keyq.ne.0) ) then
      if( io_off .eq. 'NDNO' ) then
        call smn_off(offq,ndd,ii)
      else
        read(ii,iostat=stat) offq
        call write_stat('OFF-SHELL',stat)
      end if
      !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(l, m)
      !c!$OMP   PARALLEL DO DEFAULT(PRIVATE) SHARED(nener,ndq2,raqqq_tri,fact2,offq)
      do l=1,nener
        do m=1,ndq2
          raqqq_tri(m,l) = raqqq_tri(m,l) + fact2(l) * offq(m)
        end do
      end do
      !$OMP   END PARALLEL DO
    end if
    write(*,"(a,i4,a,i3)") 'Channel ', ich, ', iq = 2, ik = ', ik
  end do ! ik=1,nok(2)

  write(*,"(a,i4)") 'VKGKV3MOD - channel ', ich
end do ! ich=1,nchl

deallocate( offd )
if( allocated(offq) ) deallocate( offq )

write(*,'(a)') 'vkgkv3_off_deninmem'

end subroutine vkgkv3_off_deninmem
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine vkgkv3_off_offinmem

use legacy,  only: ndd, ndd2, ndq, ndq2, keyq, nchl, nchls
use quadsym, only: nok

implicit none

integer :: ik, ii, ich, n_off, i_off, j_off
integer :: alloc_stat, stat

n_off = (nok(1)+nok(2)) * nchl
allocate( offd_all(ndd2,n_off), stat=alloc_stat )
call allocation_stat('offd_all',alloc_stat)

if( keyq .ne. 0 ) then
  n_off = (nok(1)+nok(2)) * (nchl-nchls)
  allocate( offq_all(ndq2,n_off), stat=alloc_stat )
  call allocation_stat('offq_all',alloc_stat)
end If

write(*,"(a)") 'VKGKV3MOD - BEGINS  '

do ich=1,nchl
  j_off = (ich-1)*(nok(1)+nok(2))

  do ik=1,nok(1)
    ii=kkktap(ik,ich,1)
    i_off = j_off + ik
    if( io_off .eq. 'NDNO' ) then
      call smn_off(offd_all(:,i_off),ndd,ii)
    else
      read(ii,iostat=stat) offd_all(:,i_off)
      call write_stat('OFF-SHELL',stat)
    end if
    if( (ich.gt.nchls).and.(keyq.ne.0) ) then
      if( io_off .eq. 'NDNO' ) then
        call smn_off(offq_all(:,i_off),ndq,ii)
      else
        read(ii,iostat=stat) offq_all(:,i_off)
        call write_stat('OFF-SHELL',stat)
      end if
    end if 
    write(*,"(a,i4,a,i3)") 'Channel ', ich, ', iq = 1, ik = ', ik
  end do 

  do ik=1,nok(2)
    ii=kkktap(ik,ich,2)
    i_off = j_off + nok(1) + ik
    if( io_off .eq. 'NDNO' ) then
      call smn_off(offd_all(:,i_off),ndd,ii)
    else
      read(ii,iostat=stat) offd_all(:,i_off)
      call write_stat('OFF-SHELL',stat)
    end if
    if( (ich.gt.nchls).and.(keyq.ne.0) ) then
      if( io_off .eq. 'NDNO' ) then
        call smn_off(offq_all(:,i_off),ndq,ii)
      else
        read(ii,iostat=stat) offq_all(:,i_off)
        call write_stat('OFF-SHELL',stat)
      end if
    end if
    write(*,"(a,i4,a,i3)") 'Channel ', ich, ', iq = 2, ik = ', ik
  end do 

  write(*,"(a,i4)") 'VKGKV3MOD - channel ', ich
end do 

write(*,'(a)') 'vkgkv3_off_offinmem'

end subroutine vkgkv3_off_offinmem
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine vkgkv3_on_offinmem(enen)

use legacy,  only: ndd, ndd2, ndq, ndq2, keyq, nchl, nchls
use quadsym, only: nop, wk, xx, ww1, nok, kradtyp, xkmax 

implicit none

integer, intent(in) :: enen

integer :: ich, numt, ik, i_off, j_off 
integer(i_large) :: nm
real*8  :: xk0, faca1, faca2
real*8 :: fact2(nok(1)+nok(2)), fact2r, fact2_sum
complex*16, allocatable, dimension(:,:) :: cvec
real*8, allocatable, dimension(:) :: on1d, on2q
real*8, allocatable, dimension(:) :: iaqqd_tri, iaqqq_tri
integer :: stat, alloc_stat

!$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
!c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndd2,raqqd_tri,raqqd1_tri,raqqd2_tri,ene,enen)
do nm=1,ndd2
  raqqd_tri(nm,1) = raqqd1_tri(nm) + raqqd2_tri(nm) * ene(enen)
end do
!$OMP  END PARALLEL DO
if(keyq.ne.0) then
  !$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
  !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndq2,raqqq_tri,raqqq1_tri,raqqq2_tri,ene,enen)
  do nm=1,ndq2
    raqqq_tri(nm,1) = raqqq1_tri(nm) + raqqq2_tri(nm) * ene(enen)
  end do
  !$OMP  END PARALLEL DO
end if
write(*,'(/a/)') 'ELMBMOD 1, 2'

allocate( on1d(ndd2), stat=alloc_stat )
call allocation_stat('on1d',alloc_stat)

if(.not. inv_real ) then
  allocate( iaqqd_tri(ndd2), stat=alloc_stat )
  call allocation_stat('iaqqd_tri',alloc_stat)
  iaqqd_tri = 0d0
  if( keyq.ne.0 ) then
    allocate( iaqqq_tri(ndd2), stat=alloc_stat )
    call allocation_stat('iaqqq_tri',alloc_stat)
    iaqqq_tri = 0d0
  end if
end if

do ich=1,nchl

  xk0 = wk(ich,enen)
  faca1 = xkmax / xk0
  faca2 = xk0 / xkmax

  do ik=1,nok(1)
    fact2(ik) =  ww1(ik,1) / ( 1.0d0 - (faca1 * xx(ik,1))**2 )
  end do
  do ik=1,nok(2)
    if( kradtyp.eq."LEGEN" ) then
      fact2(ik+nok(1)) = - ww1(ik,2) * faca2**2 / ( 1.0d0 - (faca2 * xx(ik,2))**2 )
    else 
      fact2(ik+nok(1)) = - ww1(ik,2) / ( (faca1 * xx(ik,2))**2 - 1.0d0 )
    end if
  end do
  fact2 = fact2 * fact * xkmax

  numt = numtap(ich,enen)
  allocate( cvec(ndd,nop), stat=alloc_stat )
  call allocation_stat('cvec (subroutine vkgkv3_on_offinmem)',alloc_stat)
  call rdnum(cvec,ndd,nop,numt)

! if( quadtyp .eq. "GAUSS" ) call hemdois(cvec1,ndd,nop2,nop)
  call smn(cvec,on1d,ndd)
  deallocate( cvec )

  fact2_sum = sum( fact2(:) )
  !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
  !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndd2,raqqd_tri,fact2_sum,on1d)
  do nm=1,ndd2
    raqqd_tri(nm,1) = raqqd_tri(nm,1) + fact2_sum * on1d(nm)
  end do
  !$OMP   END PARALLEL DO 

  do ik=1,nok(1)
    fact2(ik) = - fact2(ik) * (faca1 * xx(ik,1))**2
  end do
  do ik=1,nok(2)
    if( kradtyp.eq."LEGEN" ) then
      fact2(nok(1)+ik) = - fact2(nok(1)+ik) / (faca2 * xx(ik,2))**2
    else
      fact2(nok(1)+ik) = - fact2(nok(1)+ik) * (faca1 * xx(ik,2))**2
    end if
  end do

  j_off = (ich-1)*(nok(1)+nok(2))
  do ik=1,nok(1)+nok(2)
    i_off = j_off + ik
    !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
    !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndd2,raqqd_tri,fact2,ik,offd_all,i_off)
    do nm=1,ndd2
      raqqd_tri(nm,1) = raqqd_tri(nm,1) + fact2(ik) * offd_all(nm,i_off)
    end do
    !$OMP   END PARALLEL DO 
  end do

  if(.not.inv_real) then
    fact2r = - rfact * wk(ich,enen)
    !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
    !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndd2,iaqqd_tri,fact2r,on1d)
    do nm=1,ndd2
      iaqqd_tri(nm) = iaqqd_tri(nm) + fact2r * on1d(nm)
    end do
    !$OMP   END PARALLEL DO 
  end if

end do ! ich=1,nchl

if( (ich.gt.nchls).and.(keyq.ne.0) ) then

  deallocate( on1d )
  allocate( on2q(ndq2), stat=alloc_stat )
  call allocation_stat('on2q',alloc_stat)

  do ich=nchls+1,nchl

    xk0 = wk(ich,enen)
    faca1 = xkmax / xk0
    faca2 = xk0 / xkmax

    do ik=1,nok(1)
      fact2(ik) =  ww1(ik,1) / ( 1.0d0 - (faca1 * xx(ik,1))**2 )
    end do
    do ik=1,nok(2)
      if( kradtyp.eq."LEGEN" ) then
        fact2(ik+nok(1)) = - ww1(ik,2) * faca2**2 / ( 1.0d0 - (faca2 * xx(ik,2))**2 )
      else 
        fact2(ik+nok(1)) = - ww1(ik,2) / ( (faca1 * xx(ik,2))**2 - 1.0d0 )
      end if
    end do
    fact2 = fact2 * fact * xkmax

    numt = numtap(ich,enen)
    allocate( cvec(ndq,nop), stat=alloc_stat )
    call allocation_stat('cvec (subroutine vkgkv3_on_offinmem)',alloc_stat)
    call rdnum(cvec,ndq,nop,numt)

!   if( quadtyp .eq. "GAUSS" ) call hemdois(cvec2,ndq,nop2,nop)
    call smn(cvec,on2q,ndq)
    deallocate( cvec )

    fact2_sum = sum( fact2(:) )
    !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
    !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndq2,raqqq_tri,fact2_sum,ik,on2q)
    do nm=1,ndq2
      raqqq_tri(nm,1) = raqqq_tri(nm,1) + fact2_sum * on2q(nm)
    end do
    !$OMP   END PARALLEL DO 

    do ik=1,nok(1)
      fact2(ik) = - fact2(ik) * (faca1 * xx(ik,1))**2
    end do
    do ik=1,nok(2)
      if( kradtyp.eq."LEGEN" ) then
        fact2(nok(1)+ik) = - fact2(nok(1)+ik) / (faca2 * xx(ik,2))**2
      else
        fact2(nok(1)+ik) = - fact2(nok(1)+ik) * (faca1 * xx(ik,2))**2
      end if
    end do

    j_off = (ich-1)*(nok(1)+nok(2))
    do ik=1,nok(1)+nok(2)
      i_off = j_off + ik
      !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
      !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndq2,raqqq_tri,fact2,ik,offq_all,i_off)
      do nm=1,ndq2
        raqqq_tri(nm,1) = raqqq_tri(nm,1) + fact2(ik) * offq_all(nm,i_off)
      end do
      !$OMP   END PARALLEL DO 
    end do

    if(.not.inv_real) then
      fact2r = - rfact * wk(ich,enen)
      !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
      !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndq2,iaqqq_tri,fact2r,on2q)
      do nm=1,ndq2
        iaqqq_tri(nm) = iaqqq_tri(nm) + fact2r * on2q(nm)
      end do
      !$OMP   END PARALLEL DO 
    end if

  end do ! i=1,nchls+1,nchl

end if ! key.ne.0

if(inv_real) then

  if( keyden .lt. 0 ) then
    write(dentap(enen),iostat=stat) raqqd_tri(:,1)
  end if
  call square_real_matrix(raqqd,raqqd_tri(:,1))

  if( keyq.ne.0 ) then
    if( keyden .lt. 0 ) then
      write(dentap(enen),iostat=stat) raqqq_tri(:,1)
    end if
    call square_real_matrix(raqqq,raqqq_tri(:,1))
  end if

else

  if( keyden .lt. 0 ) then
    write(dentap(enen),iostat=stat) raqqd_tri(:,1)
    write(dentap(enen),iostat=stat) iaqqd_tri
  end if

  call square_complex_matrix(caqqd,raqqd_tri(:,1),iaqqd_tri)

  if( keyq.ne.0 ) then
    if( keyden .lt. 0 ) then
      write(dentap(enen),iostat=stat) raqqq_tri(:,1)
      write(dentap(enen),iostat=stat) iaqqq_tri
    end if
    call square_complex_matrix(caqqq,raqqq_tri(:,1),iaqqq_tri)
  end if

end if

write(*,*) 'vkgkv3_on_offinmem'

end subroutine vkgkv3_on_offinmem
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine vkgkv3_on_deninmem(enen)

use legacy,  only: ndd, ndd2, ndq, ndq2, keyq, nchl, nchls 
use quadsym, only: nop, wk, xx, ww1, nok, kradtyp, xkmax 

implicit none

integer, intent(in) :: enen

integer :: ich, numt
integer(i_large) :: nm
real*8  :: xk0, faca1, faca2, fact2, fact2r
complex*16, allocatable, dimension(:,:) :: cvec
real*8, allocatable, dimension(:) :: on1d, on2q
real*8, allocatable, dimension(:) :: iaqqd_tri, iaqqq_tri
integer :: stat, alloc_stat

allocate ( on1d(ndd2), stat=alloc_stat )
call allocation_stat('on1d',alloc_stat)

if(.not. inv_real ) then
  allocate( iaqqd_tri(ndd2), stat=alloc_stat )
  call allocation_stat('iaqqd_tri',alloc_stat)
  iaqqd_tri = 0d0
  if( keyq.ne.0 ) then
    allocate( iaqqq_tri(ndd2), stat=alloc_stat )
    call allocation_stat('iaqqq_tri',alloc_stat)
  iaqqq_tri = 0d0
  end if
end if

do ich=1,nchl

  xk0 = wk(ich,enen)
  faca1 = xkmax / xk0
  faca2 = xk0 / xkmax

  fact2 = sum( ww1(1:nok(1),1) / ( 1.0d0 - (faca1 * xx(1:nok(1),1))**2 ) )
  if( kradtyp.eq."LEGEN" ) then
    fact2 = fact2 - sum( ww1(1:nok(2),2) * faca2**2 / ( 1.0d0 - (faca2 * xx(1:nok(2),2))**2 ) )
  else 
    fact2 = fact2 - sum( ww1(1:nok(2),2) / ( (faca1 * xx(1:nok(2),2))**2 - 1.0d0 ) )
  end if
  fact2 = fact2 * fact * xkmax

  numt = numtap(ich,enen)
  allocate( cvec(ndd,nop), stat=alloc_stat )
  call allocation_stat('cvec (subroutine vkgkv3_on_deninmem)',alloc_stat)
  call rdnum(cvec,ndd,nop,numt)

! if( quadtyp .eq. "GAUSS" ) call hemdois(cvec1,ndd,nop2,nop)
  call smn(cvec,on1d,ndd)
  deallocate( cvec )

  !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
  !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndd2,raqqd_tri,enen,fact2,on1d)
  do nm=1,ndd2
    raqqd_tri(nm,enen) = raqqd_tri(nm,enen) + fact2 * on1d(nm)
  end do
  !$OMP   END PARALLEL DO 

  if(.not.inv_real) then
    fact2r = - rfact * wk(ich,enen)
    !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
    !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndd2,iaqqd_tri,fact2r,on1d)
    do nm=1,ndd2
      iaqqd_tri(nm) = iaqqd_tri(nm) + fact2r * on1d(nm)
    end do
    !$OMP   END PARALLEL DO
  end if

end do ! ich=1,nchl


if( (ich.gt.nchls).and.(keyq.ne.0) ) then

  deallocate( on1d )
  allocate( on2q(ndq2), stat=alloc_stat )
  call allocation_stat('on2q',alloc_stat)

  do ich=nchls+1,nchl

    xk0 = wk(ich,enen)
    faca1 = xkmax / xk0
    faca2 = xk0 / xkmax

    fact2 = sum( ww1(1:nok(1),1) / ( 1.0d0 - (faca1 * xx(1:nok(1),1))**2 ) )
    if( kradtyp.eq."LEGEN" ) then
      fact2 = fact2 - sum( ww1(1:nok(2),2) * faca2**2 / ( 1.0d0 - (faca2 * xx(1:nok(2),2))**2 ) )
    else 
      fact2 = fact2 - sum( ww1(1:nok(2),2) / ( (faca1 * xx(1:nok(2),2))**2 - 1.0d0 ) )
    end if
    fact2 = fact2 * fact * xkmax

    numt = numtap(ich,enen)
    allocate( cvec(ndq,nop), stat=alloc_stat )
    call allocation_stat('cvec (subroutine vkgkv3_on_deninmem)',alloc_stat)
    call rdnum(cvec,ndq,nop,numt)

!   if( quadtyp .eq. "GAUSS" ) call hemdois(cvec1,ndd,nop2,nop)
    call smn(cvec,on2q,ndq)
    deallocate( cvec )

    !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
    !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndq2,raqqq_tri,enen,fact2,on2q)
    do nm=1,ndq2
      raqqq_tri(nm,enen) = raqqq_tri(nm,enen) + fact2 * on2q(nm)
    end do
    !$OMP   END PARALLEL DO 

    if(.not.inv_real) then
      fact2r = - rfact * wk(ich,enen)
      !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(nm)
      !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(ndq2,iaqqq_tri,fact2r,on2q)
      do nm=1,ndq2
        iaqqq_tri(nm) = iaqqq_tri(nm) + fact2r * on2q(nm)
      end do
      !$OMP   END PARALLEL DO 
    end if

  end do ! i=1,nchls+1,nchl

end if ! key.ne.0


if(inv_real) then

  if( keyden .lt. 0 ) then
    write(dentap(enen),iostat=stat) raqqd_tri(:,enen)
  end if

  call square_real_matrix(raqqd,raqqd_tri(:,enen))

  if( keyq.ne.0 ) then
    if( keyden .lt. 0 ) then
      write(dentap(enen),iostat=stat) raqqq_tri(:,enen)
    end if

    call square_real_matrix(raqqq,raqqq_tri(:,enen))
  end if

else

  if( keyden .lt. 0 ) then
    write(dentap(enen),iostat=stat) raqqd_tri(:,enen)
    write(dentap(enen),iostat=stat) iaqqd_tri
  end if

  call square_complex_matrix(caqqd,raqqd_tri(:,enen),iaqqd_tri)

  if( keyq.ne.0 ) then
    if( keyden .lt. 0 ) then
      write(dentap(enen),iostat=stat) raqqq_tri(:,enen)
      write(dentap(enen),iostat=stat) iaqqq_tri
    end if

    call square_complex_matrix(caqqq,raqqq_tri(:,enen),iaqqq_tri)
  end if

end if

write(*,*) 'vkgkv3_on_deninmem'

end subroutine vkgkv3_on_deninmem
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine square_real_matrix(x,x_tri)

implicit none

real*8, dimension(:,:), intent(out) :: x
real*8, dimension(:), intent(in) :: x_tri

integer :: n, i, j
integer(i_large) :: ij

n = size(x,1)

!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, ij)
!c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(n,x,x_tri)
do i=1,n
  do j=1,i
    ij = i*(i-1)/2 + j
    x(j,i) = x_tri(ij)
    x(i,j) = x(j,i)
  end do
end do
!$OMP   END PARALLEL DO 

end subroutine square_real_matrix
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine square_complex_matrix(z,x_tri,y_tri)

implicit none

complex*16, dimension(:,:), intent(out) :: z
real*8, dimension(:), intent(in) :: x_tri
real*8, dimension(:), intent(in) :: y_tri

integer :: n, i, j
integer(i_large) :: ij

n = size(z,1)

!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, ij)
!c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(n,z,x_tri,y_tri)
do i=1,n
  do j=1,i
    ij = i*(i-1)/2 + j
    z(j,i) = cmplx( x_tri(ij), y_tri(ij) )
    z(i,j) = z(j,i)
  end do
end do
!$OMP   END PARALLEL DO 

end subroutine square_complex_matrix
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine symmetrize_complex_matrix(x)

complex*16, dimension(:,:), intent(inout) :: x

integer :: n, i, j

n = size(x,1)

!$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j)
!c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(x,n)
do j=2,n
  do i=1,j-1
    x(j,i) = x(i,j)
  end do
end do
!$OMP  END PARALLEL DO 

end subroutine symmetrize_complex_matrix
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine symmetrize_real_matrix(x)

real*8, dimension(:,:), intent(inout) :: x

integer :: n, i, j

n = size(x,1)

!$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j)
!c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(x,n)
do j=2,n
  do i=1,j-1
    x(j,i) = x(i,j)
  end do
end do
!$OMP  END PARALLEL DO 

end subroutine symmetrize_real_matrix
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine inv_luf(enen)

use legacy , only: ndd
use quadsym, only: nop, nop2, quadtyp

implicit none

integer, intent(in) :: enen

integer :: i, j, jj, ich, iich

!   Related to lapack routines
!
integer :: info
integer, allocatable, dimension(:) :: ipiv
integer :: alloc_stat

!   Use of lapack subroutines:
! SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
! M,N   = dimensions of rank-2 array A(M,N). Actial size: M=N=NDD
! A     = caqqd
! LDA   = physical dimension (LDA.ge.N), i.e., NDD or NDDX. A(LDA,N)
! IPIV  = information on row exchange ( IPIV(NX,NX) ), where NX=NDD or NDDX
! INFO  = returns zero if successfull
!
! SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
! TRANS = 'N' (always)
! N     = NDD (acutal dimension of A)
! NRHS  = NOP2 (number of columns of B)
! A     = CAQQD
! LDA   = NDD or NDDX (physical dimension of A)
! IPIV  = same as above
! B     = CVEC1
! LDB   = NDD or NDDX
! INFO  = same as above

if( allocated(cvec1) ) deallocate( cvec1 )
if( allocated(cvec2) ) deallocate( cvec2 )

!   LU factorization 

allocate( ipiv(ndd) )
call zgetrf(ndd,ndd,caqqd,ndd,ipiv,info)

if( info .ne. 0 ) then
  write(6,"(/,5x,'ZGETRS: FACTOR U IS SINGULAR.')")
  write(6,"('Abort...')")
  stop
end if

write(6,"(/,5x,'INFO = 0. COMPUTE INVERSE')")

allocate(cvec1(ndd,nop2), stat=alloc_stat)
call allocation_stat('cvec1 (inv_luf subroutine)',alloc_stat)

rewind numtap(1,enen)
call rdnum(cvec1,ndd,nop,numtap(1,enen))

if( quadtyp .eq. "GAUSS" ) call hemdois(cvec1,ndd,nop2,nop)

call zgetrs('N',ndd,nop2,caqqd,ndd,ipiv,cvec1,ndd,info)

deallocate( ipiv )

allocate(cvec2(ndd,nop), stat=alloc_stat)
call allocation_stat('cvec2 (inv_luf subroutine)',alloc_stat)

!***
! Generate amplitudes

do ich=nout1,nout2

  iich = 1
  if( cb_inmem ) iich = ich
  rewind numtap(ich,enen)
  call rdnum(cvec2,ndd,nop,numtap(ich,enen))

! if( quadtyp .eq. "GAUSS" ) call hemdois(cvec2,ndd,nop2,nop)

  !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jj)
  !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nop2,nop,cb,iich,cvec2,cvec1)
  do i=1,nop2
    do j=1,nop
      jj = j + nop
      cb(j,i,1,iich)  = sum( conjg(cvec2(:,j)) * cvec1(:,i) ) 
      cb(jj,i,1,iich) = sum( cvec2(:,j) * cvec1(:,i) ) 
    end do
  end do
  !$OMP   END PARALLEL DO
  write(ncbtap(1,ich,enen)) cb(:,:,1,iich)

end do

end subroutine inv_luf
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine inv_svd(enen)

use legacy , only: ndd
use quadsym, only: nop, nop2, quadtyp

implicit none

integer, intent(in) :: enen

integer :: i, j, k, l, m, ip, jj, ich, iich, nds 
complex*16, allocatable, dimension(:,:,:) :: cvec_svd, cvec1k

!   Related to lapack routines
!
integer :: info
integer(i_large) :: lwork
real*8,  allocatable, dimension(:) :: s, rwork
integer, allocatable, dimension(:) :: iwork
integer(i_large) :: lrwork
complex*16 :: u0(1,1), vt0(1,1), w_buff
complex*16, allocatable, dimension(:) :: work
complex*16, allocatable, dimension(:,:) :: vt
integer :: alloc_stat

! Free some memory (workspace required for svd) 
!
if( allocated(cvec1) ) deallocate( cvec1 )
if( allocated(cvec2) ) deallocate( cvec2 )

!   Singular value decomposition (SVD) of matrix A and its inverse A^-1
!  
!    A(i,j)    = sum( U(i,:) * VT(j,:) * S(:) )
!    A^-1(i,j) = sum( dconjg(U(i,:)) * dconjg(VT(:,j) ) / S(:) ) 
 
allocate( s(ndd), iwork(8*ndd) )
!lrwork = 5*ndd*ndd + 7*ndd
lrwork = 6*ndd*ndd 
allocate( rwork(lrwork), stat=alloc_stat )
call allocation_stat('rwork (inv_svd subroutine)',alloc_stat)
!
!   First call: get optimal LWORK (workspace allocation)
lwork = -1
call zgesdd('N',ndd,ndd,caqqd,ndd,s,u0,1,vt0,1,w_buff,lwork,rwork,iwork,info)
if( info .ne. 0 ) then
  write(6,"(/,'  FIRST CALL OF ZGESVD: INFO = ',i3)"), info
  write(6,"('Abort...')")
  stop
end if

lwork = 2*ndd*ndd + 3*ndd 
!write(*,*) 'lwork = ', lwork
!lwork = int( real(w_buff), i_large )
write(6,"(/,' ALLOCATED WORKSPACE FOR SVD: LWORK =',i11 )") lwork
allocate( work(lwork), vt(ndd,ndd), stat=alloc_stat )
call allocation_stat('work, vt (inv_svd subroutine)',alloc_stat)

!   Left eigenvectors overwritten on CAQQD. 
!   Right eigenvectors stored in VT.
!
call zgesdd('O',ndd,ndd,caqqd,ndd,s,u0,1,vt,ndd,work,lwork,rwork,iwork,info)

if( info .ne. 0 ) then
  write(6,"(/,' SECOND CALL OF ZGESVD: INFO = ',i3)"), info
  write(6,"('Abort...')")
  stop
end if

write(6,"(/,' SUCCESSFULL SVD DECOMPOSITION ' /)") 

nds = ndd - svd_max

if( nds .lt. ndd ) then
  do i = nds+1, ndd
    write(6,"(a,f11.8,a,i4,a,1pe13.6)") ' AT E = ', elenev(enen), ' eV, REMOVED SINGULAR VALUE # ', ndd-i+1, ' = ',  s(i)
  end do
end if
write(6,"(a,f11.8,a,1pe13.6)") ' AT E = ', elenev(enen), ' eV, CONDITION NUMBER = ', s(1) / s(ndd)

deallocate( rwork, work )
if(allocated(iwork))  deallocate( iwork )

!    Zero out singular eigenvalues and perform the first part of
!    the x = (A^-1) * B transformation, namely (VT^H) * B, where
!    B is the entrance-channel numerator matrix (CVEC1) and VT is 
!    the right-eigenvector SVD matrix (VT). The result is stored
!    in CVEC2 (used as workspace).
!
! if( nds .lt. ndd ) s(nds+1:ndd) = 0d0

allocate( cvec2(ndd,nop2), stat=alloc_stat )
call allocation_stat('cvec2 (inv_svd subroutine)',alloc_stat)
allocate( cvec1(ndd,nop2), stat=alloc_stat )
call allocation_stat('cvec1 (inv_svd subroutine)',alloc_stat)
rewind numtap(1,enen)
call rdnum(cvec1,ndd,nop,numtap(1,enen))
if( quadtyp .eq. "GAUSS" ) call hemdois(cvec1,ndd,nop2,nop)
  
vt = conjg( transpose(vt) )
s = 1d0 / s

!$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(i, m)
!c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nop2,ndd,cvec2,vt,cvec1,s)
do i = 1,nop2
  do m = 1,ndd
   cvec2(m,i) = sum( vt(:,m) * cvec1(:,i) ) * s(m)
  end do
end do
!$OMP  END PARALLEL DO

deallocate( vt )
deallocate ( s )

!   Now complete the linear system solution, X = U^H * CVEC2, where 
!   CVEC2 is the partial result generated above, U^H is the hermitean
!   conjugate of the left-eigenvector matrix (U) and X is the 
!   solution (stored in CVEC1).

 caqqd = conjg( transpose(caqqd) )

!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ip, j)
!c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nop2,ndd,cvec1,caqqd,nds,cvec2)
do ip = 1,nop2
  do j = 1,ndd
    cvec1(j,ip) = sum( caqqd(1:nds,j) * cvec2(1:nds,ip) )
  end do
end do
!$OMP   END PARALLEL DO

if( svd_pts.gt.1 ) then

  allocate( cvec1k(ndd,nop2,svd_max), stat=alloc_stat )
  call allocation_stat('cvec1d (inv_svd subroutine)',alloc_stat)
  !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ip, j, k)
  !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(svd_max,nop2,ndd,cvec1k,caqqd,cvec2)
  do k=1,svd_max
    do ip = 1,nop2
      do j = 1,ndd
        cvec1k(j,ip,k) = caqqd(ndd+1-k,j) * cvec2(ndd+1-k,ip)
      end do
    end do
  end do
  !$OMP   END PARALLEL DO

  allocate( cvec_svd(ndd,nop2,svd_pts), stat=alloc_stat )
  call allocation_stat('cvec_svd (inv_svd subroutine)',alloc_stat)
 
  !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(l, ip, j)
  !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(svd_pts,nop2,ndd,cvec_svd,cvec1,cvec1k,svd_map)
  do l = 1, svd_pts
    do ip = 1,nop2
      do j = 1,ndd
        cvec_svd(j,ip,l) = cvec1(j,ip) + sum( cvec1k(j,ip,:) * svd_map(:,l) )
      end do
    end do
  end do
  !$OMP   END PARALLEL DO

  deallocate( cvec1k )

else if( svd_pts.eq.1 ) then

  allocate( cvec_svd(ndd,nop2,1), stat=alloc_stat )
  call allocation_stat('cvec_svd (inv_svd subroutine)',alloc_stat)
  cvec_svd(:,:,1) = cvec1(:,:)

end if

do ich=nout1,nout2

  iich = 1
  if( cb_inmem ) iich = ich
  rewind numtap(ich,enen)
  call rdnum(cvec2,ndd,nop,numtap(ich,enen))
  if( quadtyp .eq. "GAUSS" ) call hemdois(cvec2,ndd,nop2,nop)

  !$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l,jj)
  !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(svd_pts,nop2,nop,cb,iich,cvec2,cvec_svd)
  do l = 1, svd_pts
    do i=1,nop2
      do j=1,nop
        jj = j + nop
        cb(j,i,l,iich) = sum( cvec2(:,jj) *cvec_svd(:,i,l) ) 
        cb(jj,i,l,iich)= sum( cvec2(:,j)  *cvec_svd(:,i,l) ) 
      end do
   end do
  end do
  !$OMP END PARALLEL DO

  do l = 1, svd_pts
    write(ncbtap(l,ich,enen)) cb(:,:,l,iich)
  end do

end do

end subroutine inv_svd
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine inv_luf_kmt(enen)

use legacy , only: ndd, nchl 
use quadsym, only: nop, nop2, quadtyp, ww3, wleb

implicit none

integer, intent(in) :: enen

integer :: l, m, mpt, npt, mptex 
integer :: nout, ich1, ich2, nop3
integer :: n1, n2, n3

!   Related to lapack routines
!
integer :: info
integer, allocatable, dimension(:) :: ipiv
integer :: alloc_stat

!***
! The lapack subroutine DGETRF will be used to LU-decompose the
! real part of the denominator matrix. From the the decomposition,
! the linear problem Rx=B will be solved, where R is the real part
! of the denominator and B is essentially a numerator tape. Note
! that R is real but B is complex, so the problem is not 'purely'
! real or 'purely' complex.
!*** 

!nout = nout2 - nout1 + 1
nout = nchl

if( allocated(cvec1) ) deallocate( cvec1 )
if( allocated(cvec2) ) deallocate( cvec2 )

allocate( ipiv(ndd) )
call dgetrf(ndd, ndd, raqqd, ndd, ipiv, info)

if( info .ne. 0 ) then
  write(6,"(/,5x,'DGETRF: FACTOR U IS SINGULAR.')")
  write(6,"('Abort...')")
  stop
end if

write(6,"(5x,'INFO = 0. COMPUTE K-MATRIX')")

nop3 = nop2 * nout
allocate( kmat(nop3,nop3,1), stat=alloc_stat )
call allocation_stat('kmat (inv_luf_kmt subroutine)',alloc_stat)

allocate( cvec1(ndd,nop2), stat=alloc_stat )
call allocation_stat('cvec1 (inv_luf_kmt subroutine)',alloc_stat)
allocate( cvec2(ndd,nop2), stat=alloc_stat )
call allocation_stat('cvec2 (inv_luf_kmt subroutine)',alloc_stat)


do ich1=1,nchl

  rewind numtap(ich1,enen)
  call rdnum(cvec1,ndd,nop,numtap(ich1,enen))
  if( quadtyp .eq. "GAUSS" ) call hemdois(cvec1,ndd,nop2,nop)

  call zgetrs('N',ndd,nop2,dcmplx(raqqd,0d0),ndd,ipiv,cvec1,ndd,info)

! Generate amplitudes
  n1 = nop2*(ich1-1)

!  Now form the K matrix
  do ich2=1,nout

    rewind numtap(ich2,enen)
    call rdnum(cvec2,ndd,nop,numtap(ich2,enen))
    if( quadtyp .eq. "GAUSS" ) call hemdois(cvec2,ndd,nop2,nop)

    n2 = nop2*(ich2-1) 

    If( quadtyp .eq. "GAUSS" ) Then

      !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(npt, mpt, mptex, n3)
      !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nop2,n1,nop,kmat,n2,ww3,cvec2,cvec1)
      do npt = 1,nop2
        n3 = n1 + npt
        do mpt = 1,nop
          mptex = mpt + nop
          kmat(n2+mpt,n3,1)   = ww3(mpt)* sum( cvec2(:,mptex) * cvec1(:,npt) )
          kmat(n2+mptex,n3,1) = ww3(mpt)* sum( cvec2(:,mpt)   * cvec1(:,npt) )
        end do
      end do
      !$OMP   END PARALLEL DO
    Else
      !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(l, m)
      !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nop2,kmat,n2,n1,wleb,cvec2,cvec1)
      do l = 1,nop2
        do m = 1,nop2
          kmat(n2+m,n1+l,1) = wleb(m) * sum( conjg(cvec2(:,m)) * cvec1(:,l) )
        end do
      end do
      !$OMP   END PARALLEL DO
    End If

!     write(96,"('KMAT should be Hermitian'))") 
!     do l = 1, nop2
!       do m = 1, nop2
!         write(96,"(2i4,2x,2(1e13.6,1x))") l,m,kmat(l,m)
!       end do
!     end do

  end do ! ich2=1,nchl

end do ! ich1=1,nchl

end subroutine inv_luf_kmt
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine inv_svd_kmt(enen)

use legacy , only: ndd, nchl 
use quadsym, only: nop, nop2, quadtyp, ww3, wleb

implicit none

integer, intent(in) :: enen

integer :: i, j, l, m, ip, nds, npt, mpt, mptex, k, k2
integer :: nout, ich1, ich2, nop3 
integer :: n1, n2, n3
complex*16, allocatable, dimension(:,:,:) :: cvec_svd, cvec1k

!
!   These are related to lapack subroutines
!
integer :: info
integer(i_large) :: lwork
real*8  :: u(1,1), vt0(1,1), w_buff
real*8,  allocatable, dimension(:) :: s, work
real*8,  allocatable, dimension(:,:) :: vt
integer, allocatable, dimension(:) :: iwork
integer :: alloc_stat

!nout = nout2 - nout1 + 1
nout = nchl

! Free some memory (workspace required for svd)
!
if( allocated(cvec1) ) deallocate( cvec1 )
if( allocated(cvec2) ) deallocate( cvec2 )

!   Singular value decomposition (SVD) of matrix A and its inverse A^-1
!  
!    A(i,j)    = sum( U(i,:) * VT(j,:) * S(:) )
!    A^-1(i,j) = sum( dconjg(U(i,:)) * dconjg(VT(:,j) ) / S(:) )
 
allocate( s(ndd) )
allocate( iwork(8*ndd) )
 
!   First call: get optimal LWORK (workspace allocation)
lwork = -1
call dgesdd('N',ndd,ndd,raqqd,ndd,s,u,1,vt0,1,w_buff,lwork,iwork,info)
!info = 0
if( info .ne. 0 ) then
  write(6,"(/,'  FIRST CALL OF DGESVD: INFO = ',i3)"), info
  write(6,"('Abort...')")
  stop
end if

lwork = 6 * ndd * ndd
!write(*,*) 'lwork = ', lwork
!lwork = int( w_buff, i_large )
write(6,"(/,' ALLOCATED WORKSPACE FOR SVD: LWORK =',i11 )") lwork
allocate( work(lwork), vt(ndd,ndd), stat=alloc_stat )
call allocation_stat('work, vt (inv_svd_kmt subroutine)',alloc_stat)

!   Left eigenvectors overwritten on RAQQ.
!   Right eigenvectors stored in VT.
!
call dgesdd('O',ndd,ndd,raqqd,ndd,s,u,1,vt,ndd,work,lwork,iwork,info)

if( info .ne. 0 ) then
  write(6,"(/,' SECOND CALL OF DGESVD: INFO = ',i3)"), info
  write(6,"('Abort...')")
  stop
end if

write(6,"(/,' SUCCESSFULL SVD DECOMPOSITION ' /)")

nds = ndd - svd_max

if( nds .lt. ndd ) then
  do i = nds+1, ndd
    write(6,"(a,f11.8,a,i4,a,1pe13.6)") ' AT E = ', elenev(enen), ' eV, REMOVED SINGULAR VALUE # ', ndd-i+1, ' = ',  s(i)
  end do
end if
write(6,"(a,f11.8,a,1pe13.6)") ' AT E = ', elenev(enen), ' eV, CONDITION NUMBER = ', s(1) / s(ndd)

deallocate( work )
deallocate( iwork )

!    Zero out singular eigenvalues and perform the first part of
!    the x = (A^-1) * B transformation, namely (VT^H) * B, where
!    B is the entrance-channel numerator matrix (CVEC1) and VT is
!    the right-eigenvector SVD matrix (VT). The result is stored
!    in CVEC2 (used as workspace).

nop3 = nop2 * nout
allocate( kmat(nop3,nop3,svd_pts), stat=alloc_stat )
call allocation_stat('kmat (inv_svd_kmt subroutine)',alloc_stat)

allocate( cvec1(ndd,nop2), stat=alloc_stat )
call allocation_stat('cvec1 (inv_svd_kmt subroutine)',alloc_stat)
allocate( cvec2(ndd,nop2), stat=alloc_stat )
call allocation_stat('cvec2 (inv_svd_kmt subroutine)',alloc_stat)
if(.not.allocated(cvec_svd)) then
  allocate( cvec_svd(ndd,nop2,svd_pts), stat=alloc_stat)
  call allocation_stat('cvec_svd (inv_svd_kmt subroutine)',alloc_stat)
end if
if(.not.allocated(cvec1k)) then
  allocate( cvec1k(ndd,nop2,svd_max), stat=alloc_stat )
  call allocation_stat('cvec1k (inv_svd_kmt subroutine)',alloc_stat)
end if

s = 1d0/s

do ich1=1,nchl

  rewind numtap(ich1,enen)
  call rdnum(cvec1,ndd,nop,numtap(ich1,enen))
  if( quadtyp .eq. "GAUSS" ) call hemdois(cvec1,ndd,nop2,nop)

  !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(i, m)
  !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nop2,ndd,cvec2,vt,cvec1,s)
  do i = 1,nop2
    do m = 1,ndd
      cvec2(m,i) = sum( vt(m,:) * cvec1(:,i) ) * s(m)
    end do
  end do
  !$OMP   END PARALLEL DO

  !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ip, j)
  !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nop2,ndd,cvec1,raqqd,nds,cvec2)
  do ip = 1,nop2
    do j = 1,ndd
      cvec1(j,ip) = sum( raqqd(j,1:nds)*cvec2(1:nds,ip) )
    end do
  end do
  !$OMP   END PARALLEL DO

  if( svd_pts.gt.1 ) then
    !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ip, j, k, k2)
    !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(svd_max,ndd,nop2,cvec1k,raqqd,cvec2)
    do k=1,svd_max
      k2 = ndd + 1 - k
      do ip = 1,nop2
        do j = 1,ndd
          cvec1k(j,ip,k) = raqqd(j,k2) * cvec2(k2,ip)
        end do
      end do
    end do
    !$OMP   END PARALLEL DO

    !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(ip, j, l)
    !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(svd_pts,nop2,ndd,cvec_svd,cvec1,cvec1k,svd_map)
    do l = 1, svd_pts
      do ip = 1,nop2
        do j = 1,ndd
          cvec_svd(j,ip,l) = cvec1(j,ip) + sum( cvec1k(j,ip,:) * svd_map(:,l) )
        end do
      end do
    end do
    !$OMP   END PARALLEL DO

  else if( svd_pts.eq.1 ) then
    cvec_svd(:,:,1) = cvec1(:,:)
  end if

!   Now complete the linear system solution, X = U^H * CVEC2, where
!   CVEC2 is the partial result generated above, U^H is the hermitean
!   conjugate of the left-eigenvector matrix (U) and X is the
!   solution (stored in CVEC1).

!  Now form the K matrix
!
  n1 = nop2*(ich1-1)

  do ich2=1,nchl

    rewind numtap(ich2,enen)
    call rdnum(cvec2,ndd,nop,numtap(ich2,enen))
    if( quadtyp .eq. "GAUSS" ) call hemdois(cvec2,ndd,nop2,nop)

    n2 = nop2*(ich2-1)

    If( quadtyp .eq. "GAUSS" ) Then
      !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(i, npt, mpt, mptex, n3)
      !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(svd_pts,nop2,n1,nop,kmat,n2,ww3,cvec2,cvec_svd)
      do i=1,svd_pts
        do npt = 1,nop2
          n3 = n1 + npt
          do mpt = 1,nop
            mptex = mpt + nop
            kmat(n2+mpt,n3,i)   =  ww3(mpt) * sum( cvec2(:,mptex) * cvec_svd(:,npt,i) )
            kmat(n2+mptex,n3,i) =  ww3(mpt) * sum( cvec2(:,mpt)   * cvec_svd(:,npt,i) )
          end do
        end do
      end do
      !$OMP   END PARALLEL DO
    Else
      !$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(i, l, m)
      !c!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(svd_pts,nop2,kmat,n2,n1,wleb,cvec_svd,cvec1)
      do i=1,svd_pts
        do l = 1,nop2
          do m = 1,nop2
            kmat(n2+m,n1+l,i) = wleb(m) * sum( conjg(cvec_svd(:,m,i)) * cvec1(:,l) )
          end do
        end do
      end do
      !$OMP   END PARALLEL DO
    End If

  end do ! ich2

end do ! ich1

end subroutine inv_svd_kmt
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine t_matrix(enen)

use legacy,  only: nchl 
use quadsym, only: nop, nop2, quadtyp, wk, wleb, ww3

implicit none

integer, intent(in) :: enen

integer :: i, j, ich, iich, l, m, mpt, mptex, npt, j1, l2
integer :: nout, ich1, ich2, n1, n2, nop3
real*8, allocatable, dimension(:) :: ww3i
complex*16, allocatable, dimension(:,:,:) :: jost_svd
complex*16 :: ci 
 
!   Related to lapack routines
integer :: info
integer, allocatable, dimension(:) :: ipiv
integer :: alloc_stat

!nout = nout2 - nout1 + 1
nout = nchl

nop3 = nop2 * nout

allocate( jost_svd(nop3,nop3,svd_pts), stat=alloc_stat )
call allocation_stat('jost_svd (t_matrix subroutine)',alloc_stat)

! fpii  = 0.25d0 / pi

! (1-iK) is stored in JOST

!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ich1, n1, ci, ich2, n2, j, j1, l2)
!c!$OMP  PARALLEL DO DEFAULT(PRIVATE) 
!c!      SHARED(svd_pts,nchl,nop2,wk,enen,fpii,jost_svd,kmat)
do i=1,svd_pts
  do ich1=1,nchl
    n1 = nop2*(ich1-1)
    ci = (0d0,-1d0) * wk(ich1,enen) * fpii
    do ich2=1,nchl
      n2 = nop2*(ich2-1)
      do j=1,nop2
        j1 = j + n1
        do l=1,nop2
          l2 = l + n2
          jost_svd(l2,j1,i) = kmat(l2,j1,i) * ci 
        end do 
      end do 
    end do 
  end do 
  do j=1,nop3
    jost_svd(j,j,i) = jost_svd(j,j,i) + 1d0
  end do 
end do 
!$OMP   END PARALLEL DO

allocate( ipiv(nop3) )
do i=1,svd_pts
! LU decomposition of JOST
  call zgetrf(nop3,nop3,jost_svd(:,:,i),nop3,ipiv,info)
  if( info .ne. 0 ) then
    write(6,"(/,5x,'ZGETRS: FACTOR U IS SINGULAR IN JOST.')")
    write(6,"('Abort...')")
    stop
  end if
  write(6,"(5x,'INFO = 0. COMPUTE JOST INVERSE')")
  call zgetrs('N',nop3,nop3,jost_svd(:,:,i),nop3,ipiv,kmat(:,:,i),nop3,info)
end do

deallocate( ipiv, jost_svd )

!   Dividing by the quadrature weights is not very elegant but
!   retrieves the original definition of the scattering amplitude
!   thus facilitating the interface with part C.
!
allocate(ww3i(nop2))
If( quadtyp .eq. "GAUSS" ) Then
  ww3i = 1d0 / ww3
  do ich=nout1,nout2
    iich = 1
    if( cb_inmem ) iich = ich
    n2 = nop2*(ich-1)
    do i=1,svd_pts
      do npt = 1,nop2
        do mpt = 1,nop
          mptex = mpt + nop
          cb(mpt,npt,i,iich)   = kmat(n2+mpt,npt,i)   * ww3i(mpt)
          cb(mptex,npt,i,iich) = kmat(n2+mptex,npt,i) * ww3i(mpt)
        end do
      end do
      write(ncbtap(i,ich,enen)) cb(:,:,i,iich)
    end do
  end do
Else
  ww3i = 1d0 / wleb
  do ich=nout1,nout2
    iich = 1
    if( cb_inmem ) iich = ich
    n2 = nop2*(ich-1)+1
    do i=1,svd_pts
      do l = 1,nop2
        do m = 1,nop2
          cb(m,l,i,iich) = kmat(n2+m,l,i) * ww3i(m)
        end do
      end do
      write(ncbtap(i,ich,enen)) cb(:,:,i,iich)
    end do
  end do
End If

!do ich=nout1,nout2
!  do i=1,svd_pts
!    write(ncbtap(i,ich,enen)) cb(:,:,i,ich)
!  end do
!end do

deallocate( ww3i )
deallocate( kmat )

end subroutine t_matrix
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine hemdois(cvec,nd,nop2,nop)

implicit none

integer, intent(in) :: nd, nop2, nop
complex*16, intent(inout) :: cvec(nd,nop2)

integer :: i, j, m

!$OMP   PARALLEL DO DEFAULT(PRIVATE) SHARED(nop,nd,cvec)
!c!$OMP   PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, m)
do i=1,nop
  j=i+nop
  do m=1,nd
    cvec(m,j) = conjg(cvec(m,i))
  end do
end do
!$OMP   END PARALLEL DO 

end subroutine hemdois
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine smn(cvec,eaqq,nd)

use quadsym, only: nop, rtww3

implicit none

integer, intent(in) :: nd
real*8, intent(out)  :: eaqq( nd * (nd+1) / 2 )
complex*16, intent(in) :: cvec(nd,nop)

integer(i_large) :: i, m, n, nm, nd2
real*8 :: cvecr(nop,nd), cveci(nop,nd)
real*8, allocatable, dimension(:) :: iaqq

!$OMP   PARALLEL DO DEFAULT(PRIVATE) SHARED(nd, cvec, cvecr, cveci, rtww3)
do i=1,nd
  cvecr(:,i) =  real(cvec(i,:)) * rtww3(:)
  cveci(:,i) = aimag(cvec(i,:)) * rtww3(:)
end do
!$OMP   END PARALLEL DO

nd2 = nd*(nd+1)/2
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nd, nmf, cvecr, eaqq)
do m = 1,nd
  do n = 1,m
!   nm = m*(m-1)/2 + n
    nm = nmf(n,m)
    eaqq(nm) = sum( cvecr(:,n)*cvecr(:,m) )
  end do
end do
!$OMP  END PARALLEL DO
allocate( iaqq(nd2) )
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nd, nmf, cveci, iaqq)
do m = 1,nd
  do n = 1,m
!   nm = m*(m-1)/2 + n
    nm = nmf(n,m)
    iaqq(nm) = sum( cveci(:,n)*cveci(:,m) )
  end do
end do
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nd2, eaqq, iaqq)
do nm = 1,nd2
  eaqq(nm) = eaqq(nm) + iaqq(nm)
end do
!$OMP  END PARALLEL DO
deallocate( iaqq )

end subroutine smn
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine smn_off(eaqq,nd,itape)

use mathlib,   only : gauss, lebedev_driver
use gaussquad, only : cpquad

implicit none

integer, intent(in) :: nd, itape
real*8, intent(out)  :: eaqq( nd * (nd+1) / 2 )

!  Be aware NOTHT, NOPHI, NOP, NOP2, WW and WW2 are declared LOCALLY.
!  The quadrature variables declared in the QUADSYM module
!  are only used for the on-shell matrix elements.
!
character(len=6) :: qty
complex*16, allocatable, dimension(:,:) :: cvec
integer :: notht, nophi, nop, nop2
integer :: i, l, m, n, np, nt, inl
integer(i_large) :: nm, nd2
real*8, allocatable, dimension(:) :: tht, ww, ww2, phi, wleb, rtww3
real*8, allocatable, dimension(:,:) :: cvecr, cveci !, eaqq_sq
real*8, allocatable, dimension(:) :: iaqq
integer :: stat, alloc_stat

nd2 = nd*(nd+1)/2

!   Identify the quadrature type, allocate and read CVEC
!
read(itape,iostat=stat) qty

if( qty .eq. "GAUSS" ) then

  read(itape,iostat=stat) notht, nophi
  call write_stat('notht, nophi (smn_off subroutine)',stat)
  nop  = notht * nophi
! nop2 = nop * 2
  allocate( cvec(nd,nop), stat=alloc_stat )
  call allocation_stat('cvec (smn_off subroutine)',alloc_stat)
  call rdnum(cvec,nd,nop,itape)
! call hemdois(cvec, nd, nop2, nop)

!   Load the appropriate quadrature weights
  allocate( tht(notht), ww(notht) )
  call cpquad(notht, 0.d0, "Legendre", ww, tht)
  allocate( phi(nophi), ww2(nophi) )
  call cpquad(nophi, 0.d0, "Legendre", ww2, phi)
  ww2 = ww2 * pi / 2d0

!   Business as usual
  allocate(rtww3(nop))
  do np=1,nophi
    inl = (np-1) * notht
    do nt=1,notht
      l = inl + nt
      rtww3(l) = dsqrt( 2d0 * ww(nt) * ww2(np) )
    end do
  end do

else ! quadtyp.eq."LEBED "

  read(itape,iostat=stat) nop
  call write_stat('nop (smn_off subroutine)',stat)
  nop2 = nop
  allocate( cvec(nd,nop), stat=alloc_stat )
  call allocation_stat('cvec (smn_off subroutine)',alloc_stat)
  call rdnum(cvec,nd,nop,itape)

!   Load the appropriate quadrature weights
  allocate( tht(nop2), phi(nop2), ww(nop2), wleb(nop2) )
  call lebedev_driver(tht, phi, ww, wleb, nop2)
! Check whether we have to multiply by 4pi
!  wleb = wleb * 4d0 * pi

!   Business as usual
  allocate( rtww3(nop) )
  rtww3 = dsqrt( 2d0 * wleb )

end if ! quadtyp

allocate( cvecr(nop,nd), cveci(nop,nd), stat=alloc_stat )
call allocation_stat('cvecr, cveci (smn_off subroutine)',alloc_stat)
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nd, cvec, cvecr, cveci, rtww3)
do i=1,nd
  cvecr(:,i) =  real(cvec(i,:)) * rtww3(:)
  cveci(:,i) = aimag(cvec(i,:)) * rtww3(:)
end do
!$OMP END PARALLEL DO
deallocate( cvec )

!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nd, nmf, cvecr, eaqq)
do m = 1,nd
  do n = 1,m
    nm = nmf(n,m)
    eaqq(nm) = sum( cvecr(:,n)*cvecr(:,m) )
  end do
end do
!$OMP  END PARALLEL DO
allocate( iaqq(nd2) )
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nd, nmf, cveci, iaqq)
do m = 1,nd
  do n = 1,m
    nm = nmf(n,m)
    iaqq(nm) = sum( cveci(:,n)*cveci(:,m) )
  end do
end do
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(nd2, eaqq, iaqq)
do nm = 1,nd2
  eaqq(nm) = eaqq(nm) + iaqq(nm)
end do
!$OMP  END PARALLEL DO
deallocate( iaqq )
deallocate( cvecr, cveci )

end subroutine smn_off
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine crosec(enen)

use legacy,  only: nchl, nchls
use quadsym, only: wk, nop, quadtyp, wleb, ww3

implicit none

integer, intent(in) :: enen

integer :: ich, jch, itp, jtp, itpex, jtpex, l, state, iich
real*8  :: wki, wkj, wkij, opth, react, sgt4(nchl)

write(*,'(a)') 'CROSECMOD'

!=======================================================================
!     ich : exit channel
!     jch : entrance channel
!=======================================================================

jch=1
wkj = wk(jch,enen)

if( quadtyp .eq. "GAUSS" ) then 

do ich=nout1,nout2

  if( cb_inmem ) iich = ich
  wki = wk(ich,enen)
  wkij = wki/wkj
  do l=1,svd_pts
    sgt4(ich) = 0.0d0
    if( .not.cb_inmem ) then
      read(ncbtap(l,ich,enen)) cb(:,:,l,1)
      iich = 1
    end if
!=======================================================================
!   BODY FRAME: INTEGRAL CROSS SECTION FROM LINEAR MOMENTUM REPRESENTATION
!=======================================================================
    do jtp=1,nop
      jtpex=jtp+nop
      do itp=1,nop
        itpex=itp+nop
        sgt4(ich)=sgt4(ich)+ &
        (real(cb(itp,  jtp,  l,iich))**2 + aimag(cb(itp,  jtp,  l,iich))**2 + &
        &real(cb(itpex,jtp,  l,iich))**2 + aimag(cb(itpex,jtp,  l,iich))**2 + &
        &real(cb(itp  ,jtpex,l,iich))**2 + aimag(cb(itp,  jtpex,l,iich))**2 + &
        &real(cb(itpex,jtpex,l,iich))**2 + aimag(cb(itpex,jtpex,l,iich))**2 ) * &
          ww3(itp)*ww3(jtp)
      end do
    end do

    sgt4(ich) = sgt4(ich) * fpii * wkij
    if(ich.eq.1) then
      opth = 0.0d0
      do itp=1,nop
        itpex=itp+nop
        opth = opth + (aimag(cb(itp,itp,l,iich)) + aimag(cb(itpex,itpex,l,iich))) * ww3(itp)
      end do
      opth = opth / wk(1,enen)
      react = opth - sgt4(1)
      write(6,670) elenev(enen), l-1, opth, elenev(enen), l-1, react
    end if
    if(ich.eq.1) then
      write(6,'("### ENERGY ", f11.8, " SVD ", i2.2, " FROM GROUND    TO GROUND     ", " ### CROSS=", &
      &1Pe15.8, 1X)') elenev(enen), l-1, sgt4(ich)
    else if(ich.le.nchls+1) then
      state = ich - 1
      write(6,'("### ENERGY ", f11.8, " SVD ", i2.2, " FROM GROUND    TO SINGLET ", i3, " ### CROSS=", &
      &1Pe15.8, 1X)') elenev(enen), l-1, state, sgt4(ich)
    else if(ich.gt.nchls+1) then
      state = ich - nchls - 1
      write(6,'("### ENERGY ", f11.8, " SVD ", i2.2, " FROM GROUND    TO TRIPLET ", i3, " ### CROSS=", &
      &1Pe15.8, 1X)') elenev(enen), l-1, state, sgt4(ich)
    end if
  end do ! l=1,svd_pts
end do ! ich=nout1,nout2

else ! quadtyp .eq. "LEBED"

!   Lebedev-Laikov quadrature. At this point, NOP=NOP2=NOLEB.
!
do ich=nout1,nout2

  if( cb_inmem ) iich = ich
  wki = wk(ich,enen)
  wkij = wki/wkj
  do l=1,svd_pts
    sgt4(ich) = 0.0d0
    if( .not.cb_inmem ) then
      read(ncbtap(l,ich,enen)) cb(:,:,l,1)
      iich = 1
    end if
!=======================================================================
!   BODY FRAME: INTEGRAL CROSS SECTION FROM LINEAR MOMENTUM REPRESENTATION
!=======================================================================
    do itp = 1,nop
      do jtp = 1,nop
        sgt4(ich) = sgt4(ich) + (real(cb(itp,jtp,l,iich))**2 + aimag(cb(itp,jtp,l,iich))**2) * wleb(itp) * wleb(jtp)
      end do
    end do

    sgt4(ich) = sgt4(ich) * fpii * wkij
    if(ich.eq.1) then
      opth = 0d0
      do itp = 1,nop
        opth = opth + aimag(cb(itp,itp,l,iich)) * wleb(itp)
      end do
      opth = opth / wk(1,enen)
      react = opth - sgt4(1)
      write(6,670) elenev(enen), l-1, opth, elenev(enen), l-1, react
    end if
    if(ich.eq.1) then
      write(6,'("### ENERGY ", f11.8, " SVD ", i2.2, " FROM GROUND    TO GROUND     ", " ### CROSS=", &
       &1Pe15.8, 1X)') elenev(enen), l-1, sgt4(ich)
    else if(ich.le.nchls+1) then
      state = ich - 1
      write(6,'("### ENERGY ", f11.8, " SVD ", i2.2, " FROM GROUND    TO SINGLET ", i3, " ### CROSS=", &
       &1Pe15.8, 1X)') elenev(enen), l-1, state, sgt4(ich)
    else if(ich.gt.nchls+1) then
      state = ich - nchls - 1
      write(6,'("### ENERGY ", f11.8, " SVD ", i2.2, " FROM GROUND    TO TRIPLET ", i3, " ### CROSS=", &
       &1Pe15.8, 1X)') elenev(enen), l-1, state, sgt4(ich)
    end if
  end do ! l=1,svd_pts
end do ! ich=nout1,nout2

end if ! quadtyp

  670 format('### ENERGY ', f11.8,  ' SVD ', i2.2, ' OPTICAL THEOREM CROSS SECTION ###  OPTH=', e15.8, / &
        '### ENERGY ', f11.8, ' SVD ', i2.2, ' INELASTIC CROSS SECTION       ### REACT=', e15.8)

write(*,'(a)') 

end subroutine crosec
!----------------------------------------------------------------

!----------------------------------------------------------------
double complex function crpmol(I,J,K)

use legacy, only: none1

implicit none

integer, intent(in) :: i, j, k

integer :: ma, mi, jk

!=======================================================================
!     <PLANE WAVE(1) I(1) / J(2)K(2)>
!=======================================================================

ma=max0(j,k)
mi=min0(j,k)
jk=none1(ma)+mi
crpmol=cprt(i,jk)

end function crpmol
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine pvss2(na,i,jp,nb,k,l,cr)

implicit none

integer, intent(in) :: na, i, jp, nb, k, l
complex*16, intent(out) :: cr

!=======================================================================
!     <NB-K>(L) V <NA-I><JP>
!     SEE PVGG.
!     JP IS THE PLANE-WAVE
!=======================================================================

integer :: ia, la, ka, naa, nba 

CR=(0.0D0,0.0D0)
if(K.eq.L.OR.NA.eq.0) return
if(NA.eq.NB) go to 100
if(I.eq.K) go to 150
if(I.eq.L) go to 160
return
!=======================================================================
!     <NB-I>(L) V <NA-I><JP>
!=======================================================================
  150 continue
LA=iabs(L)
NBA=iabs(NB)
NAA=iabs(NA)
if(NB*NA.gt.0) CR=-crpmol(LA,NBA,NAA)
if(NB*L.gt.0) CR=CR+crpmol(NAA,NBA,LA)
go to 800
!=======================================================================
!     <NB-K>(I) V <NA-I><JP>
!=======================================================================
  160 continue
KA=iabs(K)
NAA=iabs(NA)
NBA=iabs(NB)
if(NB*NA.gt.0) CR=crpmol(KA,NBA,NAA)
if(NB*K.gt.0) CR=CR-crpmol(NAA,NBA,KA)
go to 800
!=======================================================================
!     <NA-K>(L) V <NA-I><JP>
!=======================================================================
  100 IA=iabs(I)
NAA=iabs(NA)
KA=iabs(K)
LA=iabs(L)
if(I.ne.K) go to 400
CR=cg(LA)-crpmol(LA,NAA,NAA)+crpmol(LA,IA,IA)
if(I*L.gt.0) CR=CR-crpmol(IA,IA,LA)
if(L*NA.gt.0) CR=CR+crpmol(NAA,LA,NAA)
go to 800
  400 if(I.ne.L) go to 450
CR=-(cg(KA)-crpmol(KA,NAA,NAA)+crpmol(KA,IA,IA))
if(I*K.gt.0) CR=CR+crpmol(IA,IA,KA)
if(K*NA.gt.0) CR=CR-crpmol(NAA,KA,NAA)
go to 800
  450 continue
if(I*K.gt.0) CR=CR+crpmol(LA,IA,KA)
if(I*L.gt.0) CR=CR-crpmol(KA,IA,LA)
  800 continue

end subroutine pvss2
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine pvss2_b4_s1(nb,k,l,cr)

use legacy, only: nht, npt, nsa , none1

implicit none

integer, intent(in) :: nb, k, l
complex*16, intent(out) :: cr(nsa)

integer :: ll, na, i

do ll=1,nsa
  na = nht(ll)
  i = npt(ll)
  if(na.eq.nb) then
    if(i.eq.k) then
      if(i.eq.l) then ! i.eq.k , i.eq.l
        cr(ll) = cg(i) - cprt(i,none1(na)+na) + cprt(i,none1(i)+i) - cprt(na,none1(i)+na)
      else ! i.eq.k , i.ne.l
        cr(ll) = cprt(na,none1(l)+na) - cprt(i,none1(max0(i,l))+min0(i,l)) + &
                 2d0 * ( cg(l) - cprt(l,none1(na)+na) + cprt(l,none1(i)+i) )
      end if 
    else ! i.ne.k
      if(i.eq.l) then ! i.ne.k , i.eq.l
        cr(ll) = - ( cg(k) - cprt(k,none1(na)+na) + cprt(k,none1(i)+i) ) + &
                 2d0 * ( - cprt(na,none1(k)+na) + cprt(i,none1(max0(i,k))+min0(i,k)) )
      else ! i.ne.k , i.ne.l
        cr(ll) = 2d0 * cprt(l,none1(max0(i,k))+min0(i,k)) - cprt(k,none1(max0(i,l))+min0(i,l))
      end if
    end if
  else ! na.ne.nb
    if(i.eq.k) then
      if(i.eq.l) then ! i.eq.k , i.eq.l
        cr(ll) = - ( cprt(i,none1(max0(nb,na))+min0(nb,na)) + cprt(na,none1(i)+nb) )
      else ! i.eq.k , i.ne.l
        cr(ll) = cprt(na,none1(l)+nb) - 2d0 * cprt(l,none1(max0(nb,na))+min0(nb,na))
      end if
    else ! i.ne.k
      if(i.eq.l) then ! i.ne.k , i.eq.l
        cr(ll) = cprt(k,none1(max0(nb,na))+min0(nb,na)) - 2d0 * cprt(na,none1(k)+nb)
!
      else ! i.ne.k , i.ne.l
        cr(ll) = 0d0
      end if
    end if
  end if
end do

end subroutine pvss2_b4_s1
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine pvss2_b4_s2(nb,k,l,cr)

use legacy, only: nht, npt, nsa , none1

implicit none

integer, intent(in) :: nb, k, l
complex*16, intent(out) :: cr(nsa)

integer :: ll, na, i

do ll=1,nsa
  na = nht(ll)
  i = npt(ll)
  if(na.eq.nb) then
    if(i.eq.k) then
      if(i.eq.l) then ! i.eq.k , i.eq.l
        cr(ll) = cg(i) - cprt(i,none1(na)+na) + cprt(i,none1(i)+i) - cprt(na,none1(i)+na)
      else ! i.eq.k , i.ne.l
        cr(ll) = - cprt(na,none1(l)+na) + cprt(i,none1(max0(i,l))+min0(i,l))
      end if 
    else ! i.ne.k
      if(i.eq.l) then ! i.ne.k , i.eq.l
        cr(ll) = cg(k) - cprt(k,none1(na)+na) + cprt(k,none1(i)+i) 
      else ! i.ne.k , i.ne.l
        cr(ll) = cprt(k,none1(max0(i,l))+min0(i,l))
      end if
    end if
  else ! na.ne.nb
    if(i.eq.k) then
      if(i.eq.l) then ! i.eq.k , i.eq.l
        cr(ll) = - ( cprt(i,none1(max0(nb,na))+min0(nb,na)) + cprt(na,none1(i)+nb) )
      else ! i.eq.k , i.ne.l
        cr(ll) = - cprt(na,none1(l)+nb)
      end if
    else ! i.ne.k
      if(i.eq.l) then ! i.ne.k , i.eq.l
        cr(ll) = - cprt(k,none1(max0(nb,na))+min0(nb,na))
!
      else ! i.ne.k , i.ne.l
        cr(ll) = 0d0
      end if
    end if
  end if
end do

end subroutine pvss2_b4_s2
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine pvss2_b4_t1(nb,k,l,cr)

use legacy, only: nht, npt, nsa, none1

implicit none

integer, intent(in) :: nb, k, l
complex*16, intent(out) :: cr(nsa)

integer :: ll, na, i

do ll=1,nsa
  na = nht(ll)
  i = npt(ll)
  if(na.eq.nb) then
    if(i.eq.k) then
      if(i.eq.l) then ! i.eq.k , i.eq.l
        cr(ll) =  cg(i) - cprt(i,none1(na)+na) + cprt(i,none1(i)+i) + cprt(na,none1(i)+na)
      else ! i.eq.k , i.ne.l
        cr(ll) = cprt(i,none1(max0(i,l))+min0(i,l)) - cprt(na,none1(l)+na)
      end if 
    else ! i.ne.k
      if(i.eq.l) then ! i.ne.k , i.eq.l
        cr(ll) = 2d0 * cprt(na,none1(k)+na) + cg(k) - cprt(k,none1(na)+na) + cprt(k,none1(i)+i)
      else ! i.ne.k , i.ne.l
        cr(ll) = cprt(k,none1(max0(i,l))+min0(i,l))
      end if
    end if
  else ! na.ne.nb
    if(i.eq.k) then
      if(i.eq.l) then ! i.eq.k , i.eq.l
        cr(ll) = cprt(na,none1(k)+nb) - cprt(k,none1(max0(nb,na))+min0(nb,na))
      else ! i.eq.k , i.ne.l
        cr(ll) = - cprt(na,none1(l)+nb)
      end if
    else ! i.ne.k
      if(i.eq.l) then ! i.ne.k , i.eq.l
        cr(ll) = 2d0 * cprt(na,none1(k)+nb) - cprt(k,none1(max0(nb,na))+min0(nb,na))
!
      else ! i.ne.k , i.ne.l
        cr(ll) = 0d0
      end if
    end if
  end if
end do

end subroutine pvss2_b4_t1
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine pvss2_b4_t2(nb,k,l,cr)

use legacy, only: nht, npt, nsa, none1

implicit none

integer, intent(in) :: nb, k, l
complex*16, intent(out) :: cr(nsa)

integer :: ll, na, i

do ll=1,nsa
  na = nht(ll)
  i = npt(ll)
  if(na.eq.nb) then
    if(i.eq.k) then
      if(i.eq.l) then ! i.eq.k , i.eq.l
        cr(ll) = 3d0 * ( cg(i) - cprt(i,none1(na)+na) + cprt(i,none1(i)+i) + cprt(na,none1(i)+na) )
      else ! i.eq.k , i.ne.l
        cr(ll) = 2d0 * ( cg(l) - cprt(l,none1(na)+na) + cprt(l,none1(i)+i) ) + &
                 3d0 * cprt(na,none1(l)+na) + cprt(i,none1(max0(i,l))+min0(i,l))
      end if 
    else ! i.ne.k
      if(i.eq.l) then ! i.ne.k , i.eq.l
        cr(ll) = 2d0 * cprt(i,none1(max0(i,k))+min0(i,k)) + cg(k) - cprt(k,none1(na)+na) + cprt(k,none1(i)+i)
      else ! i.ne.k , i.ne.l
        cr(ll) = 2d0 * cprt(l,none1(max0(i,k))+min0(i,k)) + cprt(k,none1(max0(i,l))+min0(i,l))
      end if
    end if
  else ! na.ne.nb
    if(i.eq.k) then
      if(i.eq.l) then ! i.eq.k , i.eq.l
        cr(ll) = 3d0 * (cprt(na,none1(k)+nb) - cprt(k,none1(max0(nb,na))+min0(nb,na)) )
      else ! i.eq.k , i.ne.l
        cr(ll) = - 2d0 * cprt(l,none1(max0(nb,na))+min0(nb,na)) + 3d0 * cprt(na,none1(l)+nb)
      end if
    else ! i.ne.k
      if(i.eq.l) then ! i.ne.k , i.eq.l
        cr(ll) = - cprt(k,none1(max0(nb,na))+min0(nb,na))
!
      else ! i.ne.k , i.ne.l
        cr(ll) = 0d0
      end if
    end if
  end if
end do

end subroutine pvss2_b4_t2
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine get_tapes(ich, iqmax, nktp, ktp, kget, nokmax)

use quadsym, only: nok, kkkkey

implicit none

integer, intent(in) :: ich, nokmax
integer, intent(out) :: iqmax, nktp(2), kget(nokmax,2), ktp(nokmax,2)

integer :: iq, ik

do iq = 1,2
  nktp(iq) = 0
  do ik = 1,nok(iq)
    if( kkkkey(ik,ich,iq) .ne. 0 ) then
      iqmax    = iq
      nktp(iq) = nktp(iq) + 1
      ktp( nktp(iq), iq ) = kkktap(ik,ich,iq)
      kget( nktp(iq), iq) = ik
    end if
  end do
end do
     
end subroutine get_tapes
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine wrnum(cvec,n1,n2,itape)

implicit none

integer, intent(in) :: n1, n2, itape
complex*16, intent(in) :: cvec(n1,n2)

write(itape) cvec

end subroutine wrnum
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine wroff(cvec, n1, n2, itape) 

use quadsym, only: quadtyp, notht, nophi

implicit none

integer, intent(in) :: n1, n2, itape
complex*16, intent(in) :: cvec(n1,n2)

!***
! Save the quadrature type and the number of points to facilitate
! independent angular integration for different radial points in
! the off-shell matrix elements
!
write(itape) quadtyp
if( quadtyp .eq. "GAUSS" ) then
  write(itape) notht, nophi
else 
  write(itape) n2
end if
write(itape) cvec

end subroutine wroff
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine rdnum(cvec,n1,n2,itape)

implicit none

integer, intent(in) :: n1, n2, itape
complex*16, intent(out) :: cvec(n1,n2)

integer :: stat

read(itape,iostat=stat) cvec
call write_stat('ON-SHELL',stat)

end subroutine rdnum
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine wrsmn(EAQQ,N1,ITAPE)

implicit none

integer, intent(in) :: n1, itape
real*8, intent(in)  :: eaqq(n1*(n1+1)/2)

write(itape) eaqq

end subroutine wrsmn
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine rdsmn(eaqq,n1,itape)

implicit none

integer, intent(in) :: n1, itape
real*8, intent(out)  :: eaqq(n1,n1)

integer :: i, j, stat

do j=1,n1
  read(itape,iostat=stat) (eaqq(i,j),i=1,j)
  call write_stat('ON-SHELL',stat)
end do

end subroutine rdsmn
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine cwrt1(c,n,ntape)

implicit none

integer, intent(in) :: n, ntape
complex*16, intent(in) :: c(n)

write(ntape) c

end subroutine cwrt1
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine cwrt2(C,N1X,N2X,N1,N2,ntape)

implicit none

integer, intent(in) :: N1X, N2X, N1, N2, ntape
complex*16, intent(in) :: C(N1X,N2X)

integer :: i, j

do J=1,N2
  write(ntape)(C(I,J),I=1,N1)
end do

end subroutine cwrt2
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine crd1(C,N,ntape)

implicit none

integer, intent(in) :: N, ntape
complex*16, intent(out) :: C(N)

integer :: stat

read(ntape,iostat=stat) C
call write_stat('crd1',stat)

end subroutine crd1
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine crd2(C,N1X,N2X,N1,N2,ntape)

implicit none

integer, intent(in) :: N1X, N2X, N1, N2, ntape
complex*16, intent(out) ::  C(N1X,N2X)

integer :: i, j, stat

do J=1,N2
  read(ntape,iostat=stat) (C(I,J),I=1,N1)
  call write_stat('crd2',stat)
end do

end subroutine crd2
!----------------------------------------------------------------

end module onk3dk
!----------------------------------------------------------------
