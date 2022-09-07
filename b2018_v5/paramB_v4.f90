module paramB

implicit none

!***
!
!  Modifications introduced by Fábris Kossoski in the begining of 2017.
!
!  Some unimportant variables have been discarded ('keystop', 'keykkk',
!  'linv', 'numkey' and 'ires').
!
!  In subroutine setup_run the dimension of variables related to tape
!  numbers and files names has been increased, allowing for several
!  energy points and different discarded singular vectors. The name of 
!  the amplitude files have been modified acordingly, and now it 
!  specifies the index for the discarded combination of singular
!  vectors.
! 
!  The program now allows for the computation of scattering amplitudes
!  for different singular vector 
!
!  Hartree to eV convertion factor modified to the more accurate value
!  of 27.21138386.
!
!***

!***
!  This module was modified by MAPL to debug the multichannel calculations.
!  The modifications can be tracked by the cmal comments.
!***

!***
!  Modified to allow for the diagonalization of the K-matrix (as opposed to
!  the T-matrix of the former implementations). The diagonalization can be
!  performed with either the LUF or SVD algorithms employing the flags
!  LUF_KMT and SVD_KMT. In both cases, LAPACK subroutines are used.
!***

!    Units
integer, save :: n_pseud = 2
integer, save :: n_hci_1 = 3
integer, save :: n_hci_2 = 4
integer, save :: n_orbit = 7
integer, save :: n_wvfun = 8
integer, save :: n_apden = 9
integer, save :: n_label = 10

!    I/O and run control
integer, save :: keyplr = 0
integer, save :: keyplm = 0
integer, save :: keysol = 0
integer, save :: nchout = 0
integer, save :: keyden 
integer, save :: ndebug
integer, save :: nout1, nout2

character(len=16), save :: typrun
character(len=16), save :: io_den, io_off
character(len=16), save :: invtyp

character(len=13), save :: off_name
character(len=25), allocatable, dimension(:), save :: den_name
character(len=25), allocatable, dimension(:,:), save :: num_name
character(len=35), allocatable, dimension(:,:,:), save :: amp_name
integer, allocatable, dimension(:), save :: dentap
integer, allocatable, dimension(:,:), save :: nplrtp, nplmtp
integer, allocatable, dimension(:,:), save :: numtap
integer, allocatable, dimension(:,:,:), save :: ncbtap
integer, allocatable, dimension(:,:,:), save :: kkktap

!    Related to SVD
integer, save :: svd_max = 0
integer, save :: svd_pts = 1
integer, allocatable, dimension(:,:), save :: svd_map

!    INPUT DATA
character(len=*), parameter :: flag_c1 = '%INI_CONTROL_BLOCK'
character(len=*), parameter :: flag_c2 = '%END_CONTROL_BLOCK'
character(len=*), parameter :: flag_c3 = '%INI_SYMETRY_BLOCK'
character(len=*), parameter :: flag_c4 = '%END_SYMETRY_BLOCK'
character(len=*), parameter :: flag_c5 = '%INI_QUADTRE_BLOCK'
character(len=*), parameter :: flag_c6 = '%END_QUADTRE_BLOCK'

!  Constants
real*8, parameter :: pi     = 3.141592653589793238462643d0
real*8, parameter :: rt2    = 1.4142135623730950488d0 ! sqrt(2)
real*8, parameter :: rt3    = 1.7320508075688772935d0 ! sqrt(3)
real*8, parameter :: rt6    = 2.449489742783177881d0 ! sqrt(6)
real*8, parameter :: zero   = 0.0d0
real*8, parameter :: one    = 1.0d0
real*8, parameter :: rrt2   = 1.0d0/rt2 ! 1/sqrt(2)
real*8, parameter :: rrt3   = 1.0d0/rt3 ! 1/sqrt(3)
real*8, parameter :: rrt6t3 = 3.0d0/rt6 ! 3/sqrt(6)
real*8, parameter :: rt3d2  = rt3/2.0d0 ! sqrt(3)/2
real*8, parameter :: toHartree = 27.21138386d0

real*8, parameter :: fact = - 0.5d0 / pi / pi
real*8, parameter :: rfact = 0.25d0 / pi
real*8, parameter :: fpii  = 0.25d0 / pi

!   Number and list of energy points
integer, save :: nener
real*8, allocatable, dimension(:), save :: elenev, ene

!   Parsed energy label
integer, allocatable, dimension(:), save :: len_lbl
character(len=20), allocatable, dimension(:), save :: ene_label

logical, save :: off_inmem

contains

!----------------------------------------------------------------
subroutine read_flags_b

implicit none

logical :: good = .true.
integer :: stat
integer :: i, j
character(len=16) :: svd_flag

!  Read run control flags
!
call ini_flag(flag_c1)

call read_int_flag('NDEBUG =', ndebug)

if(ndebug .ne. 0) then 
  read(5,*,iostat=stat) nchout
  call write_stat('debug options',stat)
  read(5,*,iostat=stat) keyplr, keyplm, keysol 
  call write_stat('debug options',stat)
end if

write(6,"(//,'*** NDEBUG = ',i2,' ***')") ndebug
write(6,"('NCHOUT = ', i2)") nchout
write(6,"('KEYPLR = ', i2)") keyplr
write(6,"('KEYPLM = ', i2)") keyplm
write(6,"('KEYSOL = ', i2)") keysol

call read_char_flag('RUNTYP =', typrun)
typrun = trim(adjustl(typrun))

good = ( (typrun.eq."ONK").or.(typrun.eq."OFF").or.(typrun.eq."3DK") )
if(.not. good) then
  write(6,"('UNKNOWN TYPE OF RUN: RUNTYP = ',a7)") typrun
  write(6,"('IT SHOULD BE ONE OF: ONK, OFF, 3DK')")
  stop'ABORT...'
end if

keyden = 0
!     if(typrun.eq." ONK") then
! keyden  = 0
! keykkk  = 0
! keystop = 1 
! linv    = 1
!     else if(typrun.eq." OFF") then
! keyden  = 0
! keykkk  = 1
! keystop = 1 
! linv    = 1
!     else if(typrun.eq." 3DK") then
! keyden  = 0
! keykkk  = 1
! keystop = 0 
! linv    = 1
!     end if
!

write(6,"(//,'*** RUNTYP = ',a4,' ***')") typrun
if(typrun.eq."ONK")  write(6,"('COMPUTE NUMERATOR MATRIX')")
if(typrun.eq."OFF")  write(6,"('COMPUTE TAPES FOR OFF-SHELL INTEGRATION')")
if(typrun.eq."3DK")  write(6,"('COMPUTE CROSS SECTION (OFF-SHELL)')")
!***

!***
!  Save/Read denominator (relevant for typrun=3DK)
!
call read_char_flag('KEYDEN =', io_den)
io_den = trim(adjustl(io_den))

good = ( (io_den.eq."SAVE") .or. (io_den.eq."READ") .or. (io_den.eq."NONE") )
if(.not. good) then
  write(6,"('UNKNOWN KEYDEN = ',a7)") io_den
  write(6,"('IT SHOULD BE ONE OF: SAVE, READ, NONEK')")
  stop'ABORT...'
end if

if( io_den.eq."SAVE") keyden = -1
if( io_den.eq."READ") keyden =  1

if(typrun.eq."3DK") then 
  write(6,"(//,'*** KEYDEN = ',a5,' ***')") io_den
  if(io_den.eq."SAVE")  write(6,"('DENOMINATOR MATRIX WILL BE SAVED')")
  if(io_den.eq."READ")  write(6,"('DENOMINATOR MATRIX WILL BE READ FROM FILE')")
  if(io_den.eq."NONE")  write(6,"('DO NOT READ/WRITE DENOMINATOR MATRIX')")
end if
!***

!***
!  Save off-shell tapes with dimension NDxND or NDxNOP
!
call read_char_flag('WRTOFF =', io_off)
io_off = trim(adjustl(io_off))

good = ( (io_off.eq."NDND") .or. (io_off.eq."NDNO") )
if(.not. good) then
  write(6,"('UNKNOWN WRTOFF = ',a7)") io_off
  write(6,"('IT SHOULD BE ONE OF: NDND, NDNO')")
  stop'ABORT...'
end if

if(typrun.ne."ONK") then 
  write(6,"(//,'*** WRTOFF = ',a5,' ***')") io_off
  if(io_off.eq."NDND") then
    write(6,"('OFF-SHELL TAPES I/O IN NDxND FORMAT')")
  else
    write(6,"('OFF-SHELL TAPES I/O IN NDxNOP FORMAT')")
  end if
end if
!***

!***
!  Choose denominator inversion algorithm
!
call read_char_flag('INVERT =', invtyp)
invtyp=trim(adjustl(invtyp))

if( invtyp.eq.'DEFAULT' ) stop 'invmod routine has been disabled'

good = ( (invtyp.eq."LUF_LPK") .or.  (invtyp.eq."SVD_LPK") &
       .or. (invtyp.eq."LUF_KMT") .or. (invtyp.eq."SVD_KMT") )

if(.not. good) then
  write(6,"('UNKNOWN INVERT = ',a7)") invtyp
  write(6,"('IT SHOULD BE ONE OF: LUF_LPK, SVD_LPK, LUF_KMT, SVD_KMT')")
  stop'ABORT...'
end if

if( (invtyp .eq. "SVD_LPK") .or. (invtyp .eq. "SVD_KMT") ) then

  call read_int_flag('SVDMAX =', svd_max)
  if( svd_max .eq. 0 ) then
    write(6,"('YOU SHOULD USE LU DECOMPOSITION INSTEAD')")
    stop'ABORT...'
  end if
  call read_char_flag('SVDDEL =', svd_flag)
  svd_flag = trim(adjustl(svd_flag))
  good = ((svd_flag.eq."ALL") .or. (svd_flag.eq."SEL"))
  if( .not. good ) then
    write(6,"('UNKNOWN SVDDEL = ',a4)") svd_flag
    write(6,"('IT SHOULD BE ONE OF: ALL, SEL')")
    stop'ABORT...'
  end if

  if( svd_flag.eq."ALL") then
    svd_pts = svd_max + 1
    allocate( svd_map(svd_max,svd_pts) )
    svd_map = 1
    do j=2,svd_pts
      svd_map(1:j-1,j) = 0
    end do
  else
    call read_int_flag('SVDPTS =', svd_pts)
    allocate( svd_map(svd_max,svd_pts) )
    do j=1,svd_pts
      read(5,*,iostat=stat) ( svd_map(i,j), i=1,svd_max )
      call write_stat('svd_map',stat)
    end do
  end if

end if

if(typrun.eq."3DK") then 
  write(6,"(//,'*** INVERT = ',a8,' ***')") invtyp
  if( (invtyp.eq."DEFAULT") )  write(6,"('INVERSION WITH CINVRS ROUTINE')")
  if( (invtyp.eq."LUF_LPK") )  write(6,"('INVERSION WITH LUF LAPACK ROUTINE')")
  if( (invtyp.eq."SVD_LPK") )  write(6,"('INVERSION WITH SVD LAPACK ROUTINE')")
  if( (invtyp.eq."LUF_KMT") )  write(6,"('K-MATRIX INVERSION WITH LUF LAPACK ROUTINE')")
  if( (invtyp.eq."SVD_KMT") )  write(6,"('K-MATRIX INVERSION WITH SVD LAPACK ROUTINE')")
  if(  (invtyp.eq."SVD_LPK") .or. (invtyp.eq."SVD_KMT") ) then
    write(6,"('NUMBER OF MAXIMUM SINGULAR VECTOR TO BE REMOVED: SVDMAX = ',i6)") svd_max
    write(6,"('NUMBER OF COMBINATIONS OF SINGULAR VECTORS TO BE REMOVED: SVDPTS = ',i6)") svd_pts
    write(6,*) 
    write(6,*) 'MAPPING OF SVD AND COMBINATION OF SINGULAR VECTORS REMOVED:'
    do j=1,svd_pts
      write(6,"(a3,i3,a3,*(i6))") 'SVD',j,' - ', ( svd_map(i,j), i=1,svd_max )
    end do
  end if
end if

call in_mem

!return

end subroutine read_flags_b
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine read_ene(eg1)

implicit none

real*8, intent(in) :: eg1

integer :: i
integer :: stat

call read_int_flag('NOENER =', nener)

write(6,'(/a,i0,a)') '*** NONER = ', nener,' ***'
write(6,'(a,i0,a/)') 'CALCULATIONS WILL RUN FOR ', nener, ' ENERGY POINTS'

! Allocate variables on number of energy points
allocate( elenev(nener), ene(nener) )

call check_flag('ENERGY =','no')
read(5,*,iostat=stat) (elenev(i),i=1,nener)
call write_stat('list of energies', stat)

call end_flag(flag_c2)

ene(:) = eg1 + elenev(:) / toHartree

allocate( len_lbl(nener), ene_label(nener) )
call parse_energy(elenev)

write(6,"('INCIDENT ELECTRON ENERGIES = ', 1pe13.6,' eV')") (elenev(i),i=1,nener)

end subroutine read_ene
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine parse_energy(e)

implicit none

real*8, intent(in) :: e(nener)
!***
!  The energy is assumed smaller than 9999 eV
!  Precision is limited to 10**-8
!***
integer :: i, j 
integer :: ie, lab(5,nener)

if( any(e.gt.9999.d0) ) then
  write(6,"('ENERGIES ABOVE 9999eV CANNOT BE PARSED.')")
  stop 'ABORT...'
end if

do j=1,nener

  lab(1,j) = int(e(j)) /  10000
  lab(2,j) = int(e(j)) /   1000
  lab(3,j) = int(e(j)) /    100
  lab(4,j) = int(e(j)) /     10
  lab(5,j) = int(e(j))   

  i = 0 
  do
    i = i + 1
    if( lab(i,j) .gt. 0 ) exit
  if( i .eq. 5 ) exit
  end do

  len_lbl(j) = 6 - i

  ie = anint( 100000000 * ( e(j) - int(e(j)) ) )

  if( len_lbl(j) .eq. 1 ) then
    write(ene_label(j),"(i1,a,i8.8)") int(e(j)), '.', ie
  else if( len_lbl(j) .eq. 2 ) then 
    write(ene_label(j),"(i2,a,i8.8)") int(e(j)), '.', ie
  else if( len_lbl(j) .eq. 3 ) then 
    write(ene_label(j),"(i3,a,i8.8)") int(e(j)), '.', ie
  else if( len_lbl(j) .eq. 4 ) then 
    write(ene_label(j),"(i4,a,i8.8)") int(e(j)), '.', ie
  else
    write(ene_label(j),"(i5,a,i8.8)") int(e(j)), '.', ie
  end if

  len_lbl(j) = len_lbl(j) + 9

end do

end subroutine parse_energy
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine setup_run(nchl,nchls,nok,kkkkey,nokmax)

implicit none

integer, intent(in) :: nchl, nchls, nok(2), nokmax
integer, intent(in) :: kkkkey(nokmax,nchl,2)

integer :: i, j, k, n_off, jj, nchls1
integer :: stat

nchls1=nchls+1

!   These were used nowhere
!     allocate( numkey(nchl), ires(nchl) )
!     ires   = 1
!     numkey = 1
!     if(typrun .eq. " 3DK") numkey = 0

allocate( nplrtp(nchl,nener), nplmtp(nchl,nener) )
allocate( numtap(nchl,nener), ncbtap(svd_pts,nchl,nener) )
allocate( num_name(nchl,nener), amp_name(svd_pts,nchl,nener) )

do i = 1,nchl
  do k = 1,nener
    nplrtp(i,k) = n_label + i + k * nchl
    nplmtp(i,k) = n_label + i + k * nchl + nchl * nener
    numtap(i,k) = n_label + i + k * nchl + 2 * nchl * nener
    do j = 1,svd_pts
      ncbtap(j,i,k) = n_label + i + k * nchl + j * nener * nchl + 2 * nchl * nener
    end do

    if(i.eq.1)then
      jj=i
      if( typrun .eq. "ONK" ) then
        call assign_label("NUMCHG",num_name(i,k),jj,k)
        open(unit=numtap(i,k), file=num_name(i,k), form='unformatted',iostat=stat)
        call open_stat(num_name(i,k),stat)
      end if
      if( typrun .eq. "3DK" ) then
        call assign_label("NUMCHG",num_name(i,k),jj,k)
        call check_file(num_name(i,k))
        open(unit=numtap(i,k), file=num_name(i,k), form='unformatted',iostat=stat)
        call open_stat(num_name(i,k),stat)
        do j=1,svd_pts
          call assign_label_svd("AMPCHG",amp_name(j,i,k),jj,j-1,k)
          open(unit=ncbtap(j,i,k), file=amp_name(j,i,k), form='unformatted',iostat=stat)
          call open_stat(amp_name(j,i,k),stat)
        end do
      end if
    endif

    if(i.gt.1.and.i.le.nchls1)then
      jj=i-1
      if( typrun .eq. "ONK" ) then
        call assign_label("NUMCHS",num_name(i,k),jj,k)
        open(unit=numtap(i,k), file=num_name(i,k), form='unformatted',iostat=stat)
        call open_stat(num_name(i,k),stat)
      end if
      if( typrun .eq. "3DK" ) then
        call assign_label("NUMCHS",num_name(i,k),jj,k)
        call check_file(num_name(i,k))
        open(unit=numtap(i,k), file=num_name(i,k), form='unformatted',iostat=stat)
        call open_stat(num_name(i,k),stat)
        do j=1,svd_pts
          call assign_label_svd("AMPCHS",amp_name(j,i,k),jj,j-1,k)
          open(unit=ncbtap(j,i,k), file=amp_name(j,i,k), form='unformatted',iostat=stat)
          call open_stat(amp_name(j,i,k),stat)
        end do
      end if
    endif

    if(i.gt.nchls1)then
      jj=i-nchls1
      if( typrun .eq. "ONK" ) then
       call assign_label("NUMCHT",num_name(i,k),jj,k)
       open(unit=numtap(i,k), file=num_name(i,k), form='unformatted',iostat=stat)
       call open_stat(num_name(i,k),stat)
      end if
      if( typrun .eq. "3DK" ) then
        call assign_label("NUMCHT",num_name(i,k),jj,k)
        call check_file(num_name(i,k))
        open(unit=numtap(i,k), file=num_name(i,k), form='unformatted',iostat=stat)
        call open_stat(num_name(i,k),stat)
        do j=1,svd_pts
          call assign_label_svd("AMPCHT",amp_name(j,i,k),jj,j-1,k)
          open(unit=ncbtap(j,i,k), file=amp_name(j,i,k), form='unformatted',iostat=stat)
          call open_stat(amp_name(j,i,k),stat)
        end do
      end if
    endif

  end do ! k=1,nener

!   n_label = n_label + 1
!   The indexing now uses the dummy indexes i, j, k

end do ! i=1,nchl


if(typrun .eq. "ONK") return

!     n_label = 1 + 9 + nchl + nener * nchl + svd_pts * nener * nchl + 2 * nchl * nener
n_label = 1 + 10 + nchl + (nener * nchl) + (svd_pts * nener * nchl) + (2 * nchl * nener)

allocate( kkktap(nokmax,nchl,2) )
kkktap = 0

do k = 1,2
  do j = 1, nchl
    do i = 1, nok(k)

      if(j.eq.1)then
        jj=j
        if( kkkkey(i,j,k) .ne. 0 ) then
          n_off = i
          if( k .eq. 2 ) n_off = n_off + nok(1)
          kkktap(i,j,k) = n_label
          write(off_name,"(a6,i3.3,a1,i3.3)") 'OFFCHG',jj,'_',n_off
          if(typrun .eq. "3DK") call check_file(off_name)
          open(unit=n_label,file=off_name,form='unformatted',iostat=stat)
          call open_stat(off_name,stat)
          n_label = n_label + 1
        end if
      endif

      if(j.gt.1.and.j.le.nchls1) then
        jj=j-1
        if( kkkkey(i,j,k) .ne. 0 ) then
          n_off = i
          if( k .eq. 2 ) n_off = n_off + nok(1)
          kkktap(i,j,k) = n_label
          write(off_name,"(a6,i3.3,a1,i3.3)") 'OFFCHS',jj,'_',n_off
          if(typrun .eq. "3DK") call check_file(off_name)
          open(unit=n_label,file=off_name,form='unformatted',iostat=stat)
          call open_stat(off_name,stat)
          n_label = n_label + 1
        end if
      endif

      if(j.gt.nchls1)then
        jj=j-nchls1
        if( kkkkey(i,j,k) .ne. 0 ) then
          n_off = i
          if( k .eq. 2 ) n_off = n_off + nok(1)
          kkktap(i,j,k) = n_label
          write(off_name,"(a6,i3.3,a1,i3.3)") 'OFFCHT',jj,'_',n_off
          if(typrun .eq. "3DK") call check_file(off_name)
          open(unit=n_label,file=off_name,form='unformatted',iostat=stat)
          call open_stat(off_name,stat)
          n_label = n_label + 1
        end if
      endif

    end do ! i=1,nok(k)
  end do ! j=1,nchl
end do ! k=1,2

if( (typrun .eq. "3DK") .and. (keyden .ne. 0) )  then
  allocate( den_name(nener), dentap(nener) )
  do k=1,nener
    dentap(k) = n_label
    n_label = n_label + 1
    call assign_label("DEN_CH",den_name(k),0,k)
    open(unit=dentap(k), file=den_name(k), form='unformatted',iostat=stat)
    call open_stat(den_name(k),stat)
  end do
end if
  
end subroutine setup_run
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine assign_label(str,filename,ich,enen)

implicit none

integer, intent(in) :: ich
character(len=6), intent(in)   :: str
character(len=25), intent(out) :: filename
integer, intent(in) :: enen
 
if( len_lbl(enen) .eq. 10 ) then
  write(filename,"(a6,i3.3,a1,a10,a2)") str, ich, "_", ene_label(enen), "eV"
else if( len_lbl(enen) .eq. 11 ) then
  write(filename,"(a6,i3.3,a1,a11,a2)") str, ich, "_", ene_label(enen), "eV"
else if( len_lbl(enen) .eq. 12 ) then
  write(filename,"(a6,i3.3,a1,a12,a2)") str, ich, "_", ene_label(enen), "eV"
else if( len_lbl(enen) .eq. 13 ) then
  write(filename,"(a6,i3.3,a1,a13,a2)") str, ich, "_", ene_label(enen), "eV"
else
  write(6,"('Invalid LEN_LBL in ASSIGN_LABEL')")
  stop'Abort...'
end if

end subroutine assign_label
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine assign_label_svd(str,filename,ich,svd,enen)

implicit none

integer, intent(in) :: ich
character(len=6), intent(in)    :: str
character(len=35), intent(out)   :: filename
integer, intent(in) :: svd
integer, intent(in) :: enen

if( len_lbl(enen) .eq. 10 ) then
  write(filename,"(a6,i3.3,a1,a3,i3.3,a1,a10,a2)") str, ich, "_", "SVD", svd, "_", ene_label(enen), "eV"
else if( len_lbl(enen) .eq. 11 ) then
  write(filename,"(a6,i3.3,a1,a3,i3.3,a1,a11,a2)") str, ich, "_", "SVD", svd, "_", ene_label(enen), "eV"
else if( len_lbl(enen) .eq. 12 ) then
  write(filename,"(a6,i3.3,a1,a3,i3.3,a1,a12,a2)") str, ich, "_", "SVD", svd, "_", ene_label(enen), "eV"
else if( len_lbl(enen) .eq. 13 ) then
  write(filename,"(a6,i3.3,a1,a3,i3.3,a1,a13,a2)") str, ich, "_", "SVD", svd, "_", ene_label(enen), "eV"
else
  write(6,"('Invalid LEN_LBL in ASSIGN_LABEL_SVD')")
  stop'Abort...'
end if

end subroutine assign_label_svd
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine assign_den(str,filename,enen)

implicit none

character(len=3), intent(in)  :: str
character(len=19), intent(out) :: filename
integer, intent(in) :: enen

if( len_lbl(enen) .eq. 10 ) then
  write(filename,"(a3,a1,a10,a2)") str, "_", ene_label(enen), "eV"
else if( len_lbl(enen) .eq. 11 ) then
  write(filename,"(a3,a1,a11,a2)") str, "_", ene_label(enen), "eV"
else if( len_lbl(enen) .eq. 12 ) then
  write(filename,"(a3,a1,a12,a2)") str, "_", ene_label(enen), "eV"
else if( len_lbl(enen) .eq. 13 ) then
  write(filename,"(a3,a1,a13,a2)") str, "_", ene_label(enen), "eV"
else
  write(6,"('Invalid LEN_LBL in ASSIGN_LABEL_DEN')")
  stop'Abort...'
end if

end subroutine assign_den
!----------------------------------------------------------------

!----------------------------------------------------------------
!subroutine start_param

!implicit none
 
!     integer :: nelc
! XNFC  = 1.0d0/DFLOAT(NELC+1)
! XNFC2 = 1.0d0/DSQRT(DFLOAT(NELC+1))
! RRT2=1.0d0/DSQRT(2.0d0)
! RRT3=1.0d0/DSQRT(3.0d0)
! rrt6t3 = 3d0/dsqrt(6.0d0)
! rt3d2 = dsqrt(3.0d0)/2.0d0

!end subroutine start_param
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine in_mem

implicit none

character(len=16) :: in_mem_char
logical :: good = .true.

call read_char_flag('IN_MEM =', in_mem_char)
in_mem_char = trim(adjustl(in_mem_char))

good = ( (in_mem_char.eq."OFF").or.(in_mem_char.eq."DEN") )
if(.not. good) then
  write(6,"('UNKNOWN TYPE OF IN_MEM = ',a7)") in_mem_char
  write(6,"('IT SHOULD BE ONE OF: OFF, DEN')")
  stop'ABORT...'
end if

if( in_mem_char.eq."OFF") off_inmem = .true.
if( in_mem_char.eq."DEN") off_inmem = .false.

write(6,"(//,'*** IN_MEM = ',a5,' ***')") in_mem_char
if(in_mem_char.eq."OFF") then
  write(6,"('CALCULATIONS WILL RUN WITH OFF-SHELLS IN MEMORY')")
else
  write(6,"('CALCULATIONS WILL RUN WITH DENOMINATORS IN MEMORY')")
end if
write(6,*)

end subroutine in_mem
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine read_char_flag(flag,var)

implicit none

character(len=*), intent(in)  :: flag
character(len=*), intent(out) :: var

character(len_trim(flag)) :: f_name
integer :: stat
character(len=*), parameter :: flag_invalid  = 'INVALID FLAG NAME READING '

read(5,"(a)",advance='no',iostat=stat) f_name
if( f_name .ne. flag ) then
  write(6,*) flag_invalid, flag
  stop'ABORT...'
end if
!call check_flag(f_name,'no')
read(5,*,iostat=stat) var
call write_stat(f_name,stat)

end subroutine read_char_flag
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine read_int_flag(flag,var)

implicit none

character(len=*), intent(in)  :: flag
integer, intent(out) :: var

character(len_trim(flag)) :: f_name
integer :: stat
character(len=*), parameter :: flag_invalid  = 'INVALID FLAG NAME READING '

read(5,"(a)",advance='no',iostat=stat) f_name
if( f_name .ne. flag ) then
  write(6,*) flag_invalid, flag
  stop'ABORT...'
end if
!call check_flag(f_name,'no')
read(5,*,iostat=stat) var
call write_stat(f_name,stat)

end subroutine read_int_flag
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine read_real_flag(flag,var)

implicit none

character(len=*), intent(in)  :: flag
real*8, intent(out) :: var

character(len_trim(flag)) :: f_name
integer :: stat
character(len=*), parameter :: flag_invalid  = 'INVALID FLAG NAME READING '

read(5,"(a)",advance='no',iostat=stat) f_name
if( f_name .ne. flag ) then
  write(6,*) flag_invalid, flag
  stop'ABORT...'
end if
!call check_flag(f_name,'no')
read(5,*,iostat=stat) var
call write_stat(f_name,stat)

end subroutine read_real_flag
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine ini_flag(flag)

implicit none

character(len=*), intent(in) :: flag

character(len=*), parameter :: flag_notfound = 'DID NOT FIND FLAG '
character(len_trim(flag)) :: f_name
integer :: stat

rewind(5)
do
  read(5,"(a)",iostat=stat) f_name
  if(stat.lt.0) then
    write(*,*) flag_notfound, flag
    stop
  end if
  if( f_name.eq.flag ) exit
end do

end subroutine ini_flag
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine end_flag(flag)

implicit none

character(len=*), intent(in) :: flag

character(len=*), parameter :: flag_notfound = 'DID NOT FIND FLAG '
character(len_trim(flag)) :: f_name

read(5,"(a)") f_name
if( f_name.lt.flag ) then
  write(6,*) flag_notfound, flag
  stop'ABORT...'
end if

end subroutine end_flag
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine check_flag(flag,adv)

implicit none

character(len=*), intent(in)  :: flag
character(len=*), intent(in)  :: adv

character(len_trim(flag)) :: f_name
integer :: stat
character(len=*), parameter :: flag_invalid  = 'INVALID FLAG NAME READING '

read(5,"(a)",advance=adv,iostat=stat) f_name
if( f_name .ne. flag ) then
  write(6,*) flag_invalid, flag
  stop'ABORT...'
end if

end subroutine check_flag
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine write_stat(flag,stat)

implicit none

character(len=*), intent(in) :: flag
integer, intent(in) :: stat

if(stat.eq.0) then
!  write(6,"(2a)") trim(adjustl(flag)), ' found'
else if(stat.lt.0) then
  write(6,"(2a)") trim(adjustl(flag)), ' not found'
  stop
else if(stat.gt.0) then
  write(6,"(2a)") 'Problem reading ', trim(adjustl(flag))
  stop
end if

end subroutine write_stat
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine open_stat(filename,stat)

implicit none

character(len=*), intent(in) :: filename
integer, intent(in) :: stat

if(stat.eq.0) then
!  write(6,"(2a)") 
else if(stat.ne.0) then
  write(6,"(2a)") 'Could not open ', trim(adjustl(filename))
  stop
end if

end subroutine open_stat
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine allocation_stat(var,stat)

implicit none

character(len=*), intent(in) :: var
integer, intent(in) :: stat

if(stat.eq.0) then
!  write(6,"(2a)") 
else if(stat.ne.0) then
  write(6,"(2a)") 'Could not allocate variable ', trim(adjustl(var))
  stop
end if

end subroutine allocation_stat
!----------------------------------------------------------------

!-----------------------------------------------------
subroutine check_file(filename)

implicit none

character(len=*) :: filename
logical :: ex

inquire(file=filename,exist=ex)
if(.not.ex) then
 write(*,'(3a)') 'Stop: file ', trim(adjustl(filename)), ' not found !'
 stop
end if

end subroutine check_file
!-----------------------------------------------------

end module paramB
!----------------------------------------------------------------
