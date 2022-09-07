module quadsym

use paramB

implicit none

!***
!
!  Modifications introduced by Fábris Kossoski in the begining of 2017.
!
!  The new ao2mo_tran subroutine evaluates which are the sets of
!  molecular orbitals for which the 2-electron AO integrals must be 
!  transformed. Since many 2-electron MO integrals were never really used,
!  this initial evaluation implies in a huge time economy in the AO to
!  MO transformation performed in the plane subroutine.
!
!  Array 'ww3(nop)' has been introduced, such that multiplications
!  between 'ww' and 'ww2' (that used to appear in inv_luf_kmt,
!  svd_luf_kmt and crossec subroutines) are now performed here, 
!  only once. The new array 'rtww3' is simpy the usual 
!  'ww3' or 'wleb' multiplied by sqrt(2), which simplifies the product 
!  that gives rise to eaqqd_1d in smn and smn_off subroutines.
!
!  Arrays 'pw' and 'wk' now have one extra dimension, to accomodate the
!  'nener' energy points.
!
!  Hartree to eV convertion factor modified to the more accurate value
!  of 27.21138386.
!
!  'numinv', 'ngtinv' and 'ngtin2' are now integers, rather than reals.
!
!  'nothtx' and 'nophix' have been removed.
!
!   Information on the symmetry of plane waves is no longer printed.
!
!***

!***
!  Modified to work with either Gauss-Legendre or Lebedev-Laikov
!  angular quadratures (according to QUADTYP). In the latter case,
!  NOP = NOP2 = NOLEB, where NOLEB is the number of LL quadrature
!  points. At present only the usual X, Z and XZ symmetry operations
!  are taken into account in LL quadratures.
!***

!***
!  Modified to generate the Gauss-Legendre quadrature points and weights 
!  with the GAUSSQUAD module. The allows for an arbitrary even number
!  quadratures points.
!***

!***
!  Modified to employ either Gauss-Legendre or Gauss-Laguerre
!  quadratures for the radial off-shell integration above KMAX.
!  The inner radial integration is always performed with Gauss-Legendre
!  quadratures. This is controlled by the KRADTYP flag ("LEGEN" or
!  "LAGUE"). An additional keyword has been added to the QUADTRE
!  block to defined the outer quadrature (QUADK).
!***

!***
! These are related to angular quadratures and orbitals
!
!     COMMON /HELP3/NOP,NOTHT,NOPHI,NQUAD(NOTHTX,NOPHIX),NTRNP,
!    3NUMINV(NOCCX,NMOFX),NGTINV(NSBFX),NGTIN2(NSBFX,NSX)

!   Gauss-Legendre:
integer, save :: notht, nophi
integer, allocatable, dimension(:,:), save :: nquad
real*8,  allocatable, dimension(:), save :: tht0, phi0, ww, ww2, ww3

!   Lebedev-Laikov:
integer, save :: noleb
real*8,  allocatable, dimension(:), save :: xleb, yleb, zleb, wleb

! Angular quadrature
integer, save :: nop, nop2, ntrn
real*8, allocatable, dimension(:) :: rtww3
   
!
integer, allocatable, dimension(:), save :: ngtinv
integer, allocatable, dimension(:,:), save :: numinv, ngtin2

character(len=3), save :: sym_lab
character(len=6), save :: quadtyp
character(len=6), save :: kradtyp

!***
! These are related to plane waves
!
!     DIMENSION PW(3,NOPX,NCHLX),WK(NCHLX),YOBI(NOPX,3)
real*8, allocatable, dimension(:,:), save :: wk
real*8, allocatable, dimension(:,:,:,:), save :: pw
!*** 

!**
! These are related to Green's function radial quadrature
!
!     COMMON /KKK/KKKTAP(NOKX,NCHLX,2),XX(NOKX,2),WW1(NOKX,2),NOK(2),
!    1KKKKEY(NOKX,NCHLX,2),XKMAX(NCHLX)
integer, save :: nokmax, noq, nok(2)
integer, allocatable, dimension(:), save :: n_pts
integer, allocatable, dimension(:,:,:), save :: kkkkey
real*8,  allocatable, dimension(:,:), save :: xx, ww1
real*8,  save :: xkmax
!***

!***
! These are related to MO symmetry
!
!     COMMON /SYM/MA(NBFNX,NTRNX),MAP(NOP2X,NTRNX),MACHLD(NCHLDX,NTRNX)
!
integer, allocatable, dimension(:,:), save :: ma, map, machld

!***
! Related to the AO to MO transformation
integer, save :: n_tran1
integer, allocatable, dimension(:), save :: n_tran2
integer, allocatable, dimension(:,:), save :: n_tran3

integer, allocatable, dimension(:), save :: tran1
integer, allocatable, dimension(:,:), save :: tran2
integer, allocatable, dimension(:,:,:), save :: tran3

contains

!----------------------------------------------------------------
subroutine sctin1

use mathlib, only : lebedev_driver
use gaussquad, only : cpquad
use legacy, only : nchl

implicit none

character(len=16) :: kkk_flag
integer :: i, j, k, ik, nk, k8, kk, noksum
integer, allocatable, dimension(:,:) :: k_pts
real*8  :: end
logical :: good
real*8, allocatable, dimension(:,:) :: off_elenev
integer :: stat

call ini_flag(flag_c5)

call read_char_flag('QUADT =',quadtyp)
quadtyp=trim(adjustl(quadtyp))

good = ( (quadtyp.eq."GAUSS").or.(quadtyp.eq."LEBED") )
if( .not. good ) then
  write(6,"('UNKNOWN QUADTYP = ',a4)") quadtyp
  write(6,"('IT SHOULD BE ONE OF: GAUSS, LEBED')")
  stop'ABORT...'
end if

!   Gauss-Legendre angular quadrature
!  
if( quadtyp.eq."GAUSS" ) then

  call read_int_flag('NOTHT =',notht)

  call read_int_flag('NOPHI =',nophi)

  nop = notht*nophi
  nop2 = 2*nop

!  Allocate Gauss-Legendre angle-quadrature arrays
!  
  allocate( tht0(notht), ww(notht) )
  allocate( phi0(notht), ww2(notht) )
  allocate( nquad(notht,nophi ) )
  allocate( ww3(nop) )
  allocate( rtww3(nop) )

  k=0
  do i=1,nophi
    do j=1,notht
      k=k+1
      nquad(j,i)=k
    end do
  end do

  end=pi/2.0d0
  call cpquad(notht, 0.d0, "Legendre", ww, tht0)

  do k=1,notht
    tht0(k)=acos(-tht0(k))
  end do
  if(nophi.ne.1) then
    call cpquad(nophi, 0.d0, "Legendre", ww2, phi0)
    do k=1,nophi
      phi0(k)=end*(phi0(k)+1.0d0)
      ww2(k)=end*ww2(k)
    end do
  else
    read(5,*,iostat=stat) phi0(1)
    call write_stat('phi0',stat)
    ww2(1)=2.0d0*pi
  endif

  k=0
  do i=1,nophi
    do j=1,notht
      k=k+1
      ww3(k) = ww(j) * ww2(i)
      rtww3(k) = dsqrt( 2d0 * ww3(k) )
    end do
  end do

  write(6,"(//,'*** GAUSS-LEGENDRE ANGULAR QUADRATURE ***')")
  write(6,400) notht
  400   FORMAT('# OF QUADRATURE POINTS FOR THETA = ',I3)
  write(6,402) (tht0(i),i=1,notht)
  write(6,401) nophi
  401   FORMAT('# OF QUADRATURE POINTS FOR PHI = ',I3)
  write(6,402) (phi0(i),i=1,nophi)
  402   FORMAT(50(1X,F9.6))

else

!   Lebedev-Laikov angular quadrature
!  
  call read_int_flag('NOPTS =',noleb)

!  Allocate Lebedev-Laikov angle-quadrature arrays
!  
  allocate( xleb(noleb), yleb(noleb) )
  allocate( zleb(noleb), wleb(noleb) )
  allocate( rtww3(noleb) )

  call lebedev_driver(xleb, yleb, zleb, wleb, noleb)
  wleb = wleb * 4d0 * pi
  rtww3 = dsqrt( 2d0 * wleb )

  nop  = noleb
  nop2 = noleb

  write(6,"(//,'*** LEBEDEV-LAIKOV ANGULAR QUADRATURE ***')")
  write(6,"('# OF QUADRATURE POINTS = ',i4)") noleb
  write(6,"(6x,'X',13x,'Y',13x,'Z',13x,'W')") 
  do k = 1,noleb
    write(6,"(4(1pe13.6,1x))") xleb(k), yleb(k), zleb(k), wleb(k)
  end do

end if

if( typrun .eq. "ONK" ) then
! noq = 1
! allocate( xx(1,1) )
! xx(1,1) = 1.0d0
  return
end if

!***
! Read radial quadrature

!   Quadrature type above XKMAX (Gauss-Legendre or Gauss-Laguerre)
call read_char_flag('QUADK =',kradtyp)
kradtyp = trim(adjustl(kradtyp))

good = ( (kradtyp.eq."LEGEN").or.(kradtyp.eq."LAGUE") )

if( .not. good ) then
  write(6,"('UNKNOWN QUADK = ',a4)") kradtyp
  write(6,"('IT SHOULD BE ONE OF: LEGEN, LAGUE')")
  stop'ABORT...'
end if

!   Number of k-points below and above XKMAX
noq = 2

call read_int_flag('KRAD1 =', nok(1))

call read_int_flag('KRAD2 =', nok(2))

! XKMAX is the same for each collision channel
call read_real_flag('XKMAX =', xkmax)

!   Flags to control off-shell calculation (OFF runs)
nokmax = max( nok(1), nok(2) )
allocate( kkkkey(nokmax,nchl,2) )

if( typrun .eq. "OFF" ) then
  call read_char_flag('KCALC =', kkk_flag)
  kkk_flag = trim(adjustl(kkk_flag))
  good = ((kkk_flag.eq."ALL").or.(kkk_flag.eq."SEL"))
  if( .not. good ) then
    write(6,"('UNKNOWN KKK_FLAG = ',a4)") kkk_flag
    write(6,"('IT SHOULD BE ONE OF: ALL, SEL')")
    stop'ABORT...'
  end if
  if(kkk_flag .eq. "ALL") then
    kkkkey = 1
!   Avoids problems in off_driver since N_PTS is noe used there
!   (MAL modification)
    allocate( n_pts(nchl) )
    n_pts = nok(1) + nok(2)
  else
    kkkkey = 0
    allocate( n_pts(nchl) )

    call check_flag('NKPTS =','no')
    read(5,*,iostat=stat) ( n_pts(i), i=1,nchl )
    call write_stat('n_pts',stat)

! N_PTS is the number of radial points calculated in this run for each channel
!  
    k = 0
    do j = 1, nchl
      if( n_pts(j) .gt. k) k = n_pts(j)
    end do

    allocate( k_pts(k,nchl) )
    k_pts = 0
!  
! K_PTS determines which points will be calculated. The k-points are numbered from
! 1 to N = (nok(1)+nok(2)) for each channel (the counting starts from 1 in each channel)
! in consistency with the KKKKEY labeling inherited from the old code.
!  
!mal
    noksum=nok(1)+nok(2)
    do j = 1, nchl
      if(n_pts(j) .ne. 0)then
        if(n_pts(j) .ne. noksum)then
          read(5,*,iostat=stat) ( k_pts(i,j), i = 1,n_pts(j) )
          call write_stat('k_pts',stat)
        else
          do i=1,noksum
            k_pts(i,j)=i
          end do
        endif
      endif
!mal-teste if(n_pts(j).ne.0) read(5,*) ( k_pts(i), i = 1,n_pts(j) )
    end do

! Generate the appropriate KKKKEY positions
!  
    do j = 1,nchl
      do i = 1,n_pts(j)
       if( k_pts(i,j) .ne. 0 ) then
         if( k_pts(i,j) .le. nok(1) ) then
           kkkkey(k_pts(i,j), j, 1) = 1
!mal-teste if(n_pts(j).ne.0) kkkkey(k_pts(i), j, 1) = 1
         else
           ik = k_pts(i,j) - nok(1)
           kkkkey(ik, j, 2) = 1
         end if
       end if
      end do
    end do

  end if

else 
  kkkkey = 1
end if

call end_flag(flag_c6)

!***
!   Allocate radial-quadrature arrays
!
allocate( xx(nokmax,2), ww1(nokmax,2) )

allocate( off_elenev(nokmax,2) )

!   Gauss-Legendre for inner radial integration
!
nk = 1
call cpquad(nok(nk), 0.d0, "Legendre", ww1(:,nk), xx(:,nk))
do ik = 1,nok(nk)
  xx(ik,nk)=0.5d0*(xx(ik,nk)+1.0d0)
  ww1(ik,nk)=0.50d0*ww1(ik,nk)
end do

!   Gauss-Legendre or Gauss-Laguerre for outer radial integration

nk = 2
if( kradtyp .eq. "LEGEN" ) then
  call cpquad(nok(nk), 0.d0, "Legendre", ww1(:,nk), xx(:,nk))
  do ik = 1,nok(nk)
    xx(ik,nk)=0.5d0*(xx(ik,nk)+1.0d0)
    ww1(ik,nk)=0.50d0*ww1(ik,nk)
  end do
else
  call cpquad(nok(nk), 0.d0, "Laguerre", ww1(:,nk), xx(:,nk))
  do ik = 1,nok(nk)
    ww1(ik,nk) = ww1(ik,nk) * exp( xx(ik,nk) )
    xx(ik,nk)  = xx(ik,nk) + 1.d0
  end do
end if

write(6,"(//,10x,'*** GAUSS-LEGENDRE RADIAL QUADRATURE BELOW KMAX ***')")
if( kradtyp .eq. "LEGEN" ) then
  write(6,"(10x,'*** GAUSS-LEGENDRE RADIAL QUADRATURE ABOVE KMAX ***')")
else
  write(6,"(10x,'*** GAUSS-LAGUERRE RADIAL QUADRATURE ABOVE KMAX ***')")
end if
write(6,"('KMAX (au) = ',1pe13.6)") xkmax
write(6,"('EMAX (eV) = ',1pe13.6)") xkmax*xkmax/2.0d0*toHartree

write(6,"('# OF QUADRATURE POINTS BELOW KMAX = ',i4)") nok(1)
k8 = nok(1)/8
do k = 1,k8
  kk = (k-1)*8 + 1  
  write(6,"(8(1pe13.6,1x))") (xx(ik,1), ik=kk,kk+7)
end do
if( mod(nok(1),8) .ne. 0 ) then
  kk = k8*8 + 1  
  write(6,"(8(1pe13.6,1x))") (xx(ik,1), ik=kk,nok(1))
end if

write(6,"('# OF QUADRATURE POINTS ABOVE KMAX = ',i4)") nok(2)
k8 = nok(2)/8
do k = 1,k8
  kk = (k-1)*8 + 1  
  write(6,"(8(1pe13.6,1x))") (xx(ik,2), ik=kk,kk+7)
end do
if( mod(nok(2),8) .ne. 0 ) then
  kk = k8*8 + 1  
  write(6,"(8(1pe13.6,1x))") (xx(ik,2), ik=kk,nok(2))
end if

off_elenev(1:nok(1),1) = xkmax**2 / 2d0 * toHartree * xx(1:nok(1),1)
if( kradtyp .eq. "LEGEN" ) then
  off_elenev(1:nok(2),2) =  xkmax**2 / 2d0 * toHartree / xx(1:nok(2),2)
else
  off_elenev(1:nok(2),2) =  xkmax**2 / 2d0 * toHartree * xx(1:nok(2),2)
end if

write(6,'(a)') 
write(6,'(a)') 'OFF-SHELL ENERGIES (eV): '
!write(6,*) (off_elenev(ik,1), ik=1,nok(1)), (off_elenev(ik,2), ik=nok(2),1,-1)
do ik=1,nok(1)
  write(6,'(a,i4,a,f16.8)') 'OFF ', ik, ' : ', off_elenev(ik,1)
end do
do ik=nok(2),1,-1
  write(6,'(a,i4,a,f16.8)') 'OFF ', nok(1)+ik, ' : ', off_elenev(ik,2)
end do
write(6,'(a)') 

deallocate( off_elenev )

end subroutine sctin1
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine sctin2

use legacy, only: nchl, eg

implicit none

real*8, allocatable, dimension(:) :: deg
integer :: i, j, k, ij
real*8  :: yc, ys, xc, xs, xxc, xxs
integer :: l

allocate( deg(nchl) )

! Allocate arrays related to plane waves
allocate( pw(3,nop,nchl,nener), wk(nchl,nener) )

do l=1,nener
  do i=1,nchl
    deg(i) = ene(l) - eg(i)
    write(6,"(a8,i3,2(a16,f12.8))") 'Channel ',i, ', Energy (eV) = ', elenev(l), ', Deg (eV) =    ', deg(i) * toHartree
    if(deg(i).lt.0.0d0) then
      write(6,"('(DEG.LT.0) IN SUBROUTINE SCTIN2')")
      stop 'ABORT...'
    else
      wk(i,l) = dsqrt(deg(i)*2.0d0)
    end if
  end do
end do
!
!     IF(NOP.NE.1) GO TO 1111
!     RT3=1.732050807D0
!     PW(1,1,1)=WK(1)/RT3
!     PW(2,1,1)=WK(1)/RT3
!     PW(3,1,1)=WK(1)/RT3
!     GO TO 2222
!   1111 CONTINUE
!=======================================================================
!=======================================================================
!DO 200 I=1,NCHL
!DEG(I)=ENE-EG(I)
!
!roma 13/11/2003
!
! Eh preciso comentar a linha que segue para poder rodar o programa bnb.f
! para energias abaixo do threshold de excitacao da molecula em questao
!
!IF(DEG(I).LT.0.0D0) GO TO 999
!
! Eh preciso tambem tomar apenas o valor absoluto da diferenca de energias
! DEG(I) para poder calcular a raiz quadrada abaixo
!
!  200 WK(I)=DSQRT(DEG(I)*2.0D0)
!  200 WK(I)=DSQRT(DABS(DEG(I))*2.0D0)
!
!roma 13/11/2003
!
!IF(NOP.NE.1) GO TO 1111
!RT3=1.732050807D0
!PW(1,1,1)=WK(1)/RT3
!PW(2,1,1)=WK(1)/RT3
!PW(3,1,1)=WK(1)/RT3
!GO TO 2222
! 1111 CONTINUE
!=======================================================================
! DETERMINE PW.
!=======================================================================

if( quadtyp .eq. 'GAUSS' ) then
  do l=1,nener
    do j=1,nophi
      yc=cos(phi0(j))
      ys=sin(phi0(j))
      do i=1,notht
        xc=cos(tht0(i))
        xs=sin(tht0(i))
        ij=nquad(i,j)
        do k=1,nchl
          xxc=xc*wk(k,l)
          xxs=xs*wk(k,l)
          pw(1,ij,k,l)=xxs*yc
          pw(2,ij,k,l)=xxs*ys
          pw(3,ij,k,l)=xxc
        end do
      end do
    end do
  end do
else
  do l=1,nener
    do ij = 1,noleb
      do k = 1,nchl 
        pw(1,ij,k,l) = wk(k,l) * xleb(ij)
        pw(2,ij,k,l) = wk(k,l) * yleb(ij)
        pw(3,ij,k,l) = wk(k,l) * zleb(ij)
      end do
    end do
  end do
end if
!
!   2222 CONTINUE

!***
!  At this point, nop=nop2=noleb, in case quadtyp='LEBED'
!***

write(6,500)
  500 FORMAT(1H1//,'**** INPUT DATA FOR PLANE WAVES ****')

! do l=1,nener
!   write(*,"(a14,f14.8,5x,a24,i5)") "ENERGY (eV) = ",  toHartree * deg(l), "NUMBER OF PLANE WAVES = ", nop
!   write(*,"(/'CHN, WAVE NUMBER'/(6('   (',I2,') ',D15.8)))") (k,wk(k,l),k=1,nchl)
! end do

!return

!999  write(6,"('(DEG.LT.0) IN SUBROUTINE SCTIN2')")
!     stop 'ABORT...'

end subroutine sctin2
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine molsym

use legacy, only: hci, nsa, nchls, nchld, nocc, ng, nbfns, ns, ngt, nht, npt, ng2, ngt2

implicit none

logical :: good, unique
integer :: i, j, l, m, ii, ip, it, jj, na 
integer :: nat, kkt, ipt, isp, itt
integer :: ispin, noth1, nophi1
real*8  :: xek
integer :: stat
integer, allocatable, dimension(:) :: mvth, mvpi

!  First, allocate and load NGTINV, NGTIN2 and NUMINV. Somewhat
!  misplaced, following the old version of the code
!
allocate( ngtinv(nbfns) )
if( ns .gt. 0 ) then
  allocate(numinv(nocc,nbfns), ngtin2(nbfns,ns) )
end if

do i=1,ng
  ii=ngt(i)
  ngtinv(ii)=i
end do
do i=1,ns
  na=nht(i)
  ii=npt(i)
  numinv(na,ii)=i
  do j=1,ng2(i)
    jj=ngt2(j,i)
    ngtin2(jj,i)=j
  end do
end do

!=======================================================================
!     READ SYMMETRY OF MOLECULAR ORBITALS AND ADDITIONAL BASIS.
!     THIS NTRN CAN BE DIFFERENT FROM OLD NTRN, SINCE HEREAFTER
!     SYMMETRY OF PLANE WAVES IS MADE USE OF.
!     TRANSFORMATION MATRIX FOR BASIS ORBITALS. IF I GOES TO J BY M
!     TRANSFORMATION, THEN MA(I,M)=J.
!=======================================================================

! Read symmetry information
!
call ini_flag(flag_c3)

call read_char_flag('SYMOP =',sym_lab)
sym_lab = trim(adjustl(sym_lab))

good = ( (sym_lab.eq.'X').or.(sym_lab.eq.'Z').or.(sym_lab.eq.'XZ').or.(sym_lab.eq.'ZX').or.(sym_lab.eq.'NO') )

if(.not. good) then
  write(6,"('UNKNOWN SYMMETRY FORMAT: SYMOP = ',a7)") sym_lab
  write(6,"('IT SHOULD BE ONE OF: X, Z, XZ, ZX, NO')")
  stop'ABORT...'
end if

ntrn = 0
if( (sym_lab .eq. 'X') .or. (sym_lab .eq. 'Z') )   ntrn = 1
if( (sym_lab .eq. 'XZ') .or. (sym_lab .eq. 'ZX') ) ntrn = 3

!   Allocate symmetry arrays. 
!
if( ntrn .gt. 0 ) then
  allocate( ma(nbfns,ntrn), map(nop2,ntrn) )
  if( quadtyp .eq. 'GAUSS') then
    allocate( mvth(ntrn), mvpi(ntrn) )
    if( sym_lab .eq. 'Z') then
      mvth(1) = 999
      mvpi(1) = 000
    else if( sym_lab .eq. 'X') then
      mvth(1) = 000
      mvpi(1) = 999
    else if( sym_lab .eq. 'XZ') then
      mvth(1) = 000
      mvpi(1) = 999
      mvth(2) = 999
      mvpi(2) = 000
      mvth(3) = 999
      mvpi(3) = 999
    else if( sym_lab .eq. 'ZX') then
      mvth(1) = 999
      mvpi(1) = 000
      mvth(2) = 000
      mvpi(2) = 999
      mvth(3) = 999
      mvpi(3) = 999
    end if
  end if
  if( nchld .gt. 0) allocate( machld(nchld,ntrn) )
  do m=1,ntrn
    read(5,*,iostat=stat) (ma(i,m),i=1,nbfns)
    call write_stat('orbitals symmetries',stat)
  end do
end if

call end_flag(flag_c4)

if( ntrn .eq. 0 ) return
!
!roma 27/09/2005
!
if(nsa.ne.0) then
  do i=1,nchld
    ispin=1
    if(i.gt.nchls) then
      ispin=2
      ii=i-nchls
    else
      ii=i
    endif
    xek=0.0d0
    l=0
    do
      l=l+1
      xek=hci(l,ii,ispin)
      if(dabs(xek).gt.0.1d0) then
        do m=1,ntrn
          nat=ma(nht(l),m)/abs(ma(nht(l),m))
          kkt=ma(npt(l),m)/abs(ma(npt(l),m))
          machld(i,m)=nat*kkt
        end do
        exit
      else
        if(l.eq.nsa) stop 'THERE IS NO CI COEFFICIENT LARGER THAN 0.1d0'
      end if
    end do
  end do
end if ! nsa.ne.0

!roma 27/09/2005DO 110 N=1,NTRN
!roma 27/09/2005  110 READ(5,*) (MACHLD(I,N),I=1,NCHLD)
!
!roma 27/09/2005
!
! 120 FORMAT(20I4)
!
!=======================================================================
!
!     TRANSFORMATION MATRIX FOR OPEN CHANNEL STATES.
!
!=======================================================================
!
!     MVTH(M) : M-TH SYMMETRY OPERATION IN THETA  COORDINATE.
!     999(REFLECTION) OR 0.
!     MVPI(M) : M-TH SYMMETRY OPERATION IN PHI COORDINATE.
!     999 AND ROTATION(SUCH AS 2 OR 3 ETC.) ARE ALLOWED.
!
!=======================================================================

!mod
! No loop abaixo, ISP=1, sempre. Assim:
! mvpi(m) = 000 => ipt = 1,...,nophi (isto e, ipt=ip)
! mvpi(m) = 999 => ipt = nophi,...,1 (isto e, ipt=nophi+1-ip)
! mvth(m) = 000 => itt = 1,...,notht (isto e, itt=it)
! mvth(m) = 999 => itt = notht,...,1 (isto e, itt=notht+1-it)
!
! MAP(I,M) = J => O ponto de quadratura I e levado ao ponto
!de quadratura J pela operacao de simetria M. Vale notar
!que I,J carregam informacao sobre theta e phi simultaneamente
!(NQUAD).
!
! Quando a operacao de simetria nao se aplica, o ponto e levado nele mesmo,
! isto e:
! mvpi(m)=0 => ipt=ip => MAP(I,M) = I
! mvth(m)=0 => itt=it => MAP(I,M) = I
!mod

if( quadtyp .eq. 'GAUSS' ) then
! Symmetry for Gauss-Legendre angular quadrature

  noth1=notht+1
  nophi1=nophi+1
  do m=1,ntrn
    do ip=1,nophi
      if(mvpi(m).eq.999) then
        ipt=nophi1-ip
      else
        ipt=ip+mvpi(m)
        if(ipt.gt.nophi) ipt=nophi-ipt
      end if
      isp=ipt/iabs(ipt)
      do it=1,notht
        itt=noth1-it
        if(mvth(m).ne.999) itt=it
        if(isp.ge.0) then
          map(nquad(it,ip),m)=nquad(itt,iabs(ipt))
          map(nquad(it,ip)+nop,m)=nquad(itt,iabs(ipt))+nop
        else
          itt=noth1-itt
          map(nquad(it,ip),m)=nquad(itt,iabs(ipt))+nop
          map(nquad(it,ip)+nop,m)=nquad(itt,iabs(ipt))
        end if
      end do
    end do
  end do

else
! Symmetry for Lebedev-Laikov angular quadrature

  do i = 1,noleb
    map(i,:) = i
  end do
  do m = 1,ntrn
    do i = 1,noleb-1
      j = i
      do 
        j = j + 1
        if( j .gt. noleb) exit
        call get_sym(unique,xleb(i),yleb(i),zleb(i),xleb(j),yleb(j),zleb(j),m)
        if( .not. unique ) then
          ip = map(i,m)
          map(i,m) = map(j,m)
          map(j,m) = ip
          exit
        end if
      end do
    end do
  end do

end if

!call print_symorbs

if( allocated( mvth) ) deallocate( mvth )
if( allocated( mvpi) ) deallocate( mvpi )

return

end subroutine molsym
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine free_inv

use legacy, only: ns

implicit none

deallocate( ngtinv )
if( ns .gt. 0 ) deallocate( numinv, ngtin2 )

end subroutine free_inv
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine get_sym(unique,x1,y1,z1,x2,y2,z2,m)

implicit none

logical, intent(out) :: unique
integer, intent(in) :: m
real*8, intent(in)  :: x1, y1, z1, x2, y2, z2

logical :: doX, doZ
real*8  :: crit
real*8  :: small = 1.0d-6

If( ntrn .eq. 1 ) Then
  if( sym_lab .eq. 'X' ) crit = abs(x1+x2) + abs(y1-y2) + abs(z1-z2)
  if( sym_lab .eq. 'Z' ) crit = abs(x1-x2) + abs(y1-y2) + abs(z1+z2)
End If

If( (ntrn .eq. 3) .and. (m .eq. 1) ) Then
  doX = ( sym_lab .eq. 'XZ' )
  doZ = ( sym_lab .eq. 'ZX' )
  if( doX ) crit = abs(x1+x2) + abs(y1-y2) + abs(z1-z2)
  if( doZ ) crit = abs(x1-x2) + abs(y1-y2) + abs(z1+z2)
Else If( (ntrn .eq. 3) .and. (m .eq. 2) ) Then
  doX = ( sym_lab .eq. 'ZX' )
  doZ = ( sym_lab .eq. 'XZ' )
  if( doX ) crit = abs(x1+x2) + abs(y1-y2) + abs(z1-z2)
  if( doZ ) crit = abs(x1-x2) + abs(y1-y2) + abs(z1+z2)
Else If( (ntrn .eq. 3) .and. (m .eq. 3) ) Then
  crit = abs(x1+x2) + abs(y1-y2) + abs(z1+z2)
End If  

unique = ( (crit/3d0) .gt. small )

end subroutine get_sym
!----------------------------------------------------------------

!----------------------------------------------------------------
!subroutine print_symorbs

!implicit none

!     WRITE(6,300) NTRN
! 300 FORMAT(//'*** TRANSFORMATION OF PLANE WAVES ***'/
!    1'NTRN=',I3)

!     write(6,"('SYMMETRY OPERATIONS = ',a3)") sym_lab
!     if( ntrn .eq. 1 ) then
! if( sym_lab .eq. 'X' ) call symorbs('X',ma(:,1),nbfns,1)
! if( sym_lab .eq. 'Z' ) call symorbs('Z',ma(:,1),nbfns,1)
!     else
! if( sym_lab .eq. 'XZ' ) then
!   call symorbs('X',ma(:,1),nbfns,1)
!   call symorbs('Z',ma(:,2),nbfns,1)
!   call symorbs('XZ',ma(:,3),nbfns,1)
! else
!   call symorbs('Z',ma(:,1),nbfns,1)
!   call symorbs('X',ma(:,2),nbfns,1)
!   call symorbs('ZX',ma(:,3),nbfns,1)
! end if
!     end if

!     WRITE(6,302)
! 302 FORMAT(//'   MACHLD(OPEN CHANNELS - 1)')
!     DO 308 N=1,NTRN
! 308 WRITE(6,335) (MACHLD(J,N),J=1,NCHLD)

!     WRITE(6,303)
! 303 FORMAT(//'   MAP(PLANE WAVES)')
!     if( ntrn .eq. 1 ) then
! if( sym_lab .eq. 'X' ) call symorbs('X',map(:,1),nop2,2)
! if( sym_lab .eq. 'Z' ) call symorbs('Z',map(:,1),nop2,2)
!     else
! if( sym_lab .eq. 'XZ' ) then
!   call symorbs('X',map(:,1),nop2,2)
!   call symorbs('Z',map(:,2),nop2,2)
!   call symorbs('XZ',map(:,3),nop2,2)
! else
!   call symorbs('Z',map(:,1),nop2,2)
!   call symorbs('X',map(:,2),nop2,2)
!   call symorbs('ZX',map(:,3),nop2,2)
! end if
!     end if
!335  FORMAT(20(1X,I3))

!end subroutine print_symorbs
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine symorbs(lab,m,n,flag)

implicit none

character(len=2), intent(in) :: lab
integer, intent(in) :: flag, n, m(n)

integer :: i, n15, i0, j

if( flag .eq. 1 ) write(6,"(/'MA: THESE TRANSFORM AS ',a3)") lab
if( flag .eq. 2 ) write(6,"(/'MAP: THESE TRANSFORM AS ',a3)") lab
n15 = n/15
do i = 1,n15
  i0 = (i-1)*15 + 1
  write(6,"(15(i5))") ( m(j), j=i0,i0+14 )
end do

if( mod(n,15) .ne. 0 ) then
  i0 = n15*15 + 1
  write(6,"(15(i5))") ( m(j), j=i0,n )
end if
 
end subroutine symorbs
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine ao2mo_tran

use legacy, only: nsa, nchls, nchl, nocc, ng, nbfns, ns, ngt, nht, npt, ng2, ngt2, numsd

implicit none

logical, allocatable, dimension(:) :: l_tran1
logical, allocatable, dimension(:,:) :: l_tran2
logical, allocatable, dimension(:,:,:) :: l_tran3
logical, allocatable, dimension(:) :: l_virt_act
integer, allocatable, dimension(:) :: virt_act
integer :: n_virt_act
integer :: i, j, ii, jj, l,  na
integer :: i1, nb, ij1, ij2, k, j2, i2, k2

allocate( l_tran1(nbfns) )
l_tran1 = .false.
n_tran1 = 0

allocate( l_tran2(nbfns,nbfns) )
allocate( n_tran2(nbfns) )
l_tran2 = .false.
n_tran2 = 0

allocate( l_tran3(nbfns,nbfns,nbfns) )
allocate( n_tran3(nbfns,nbfns) )
l_tran3 = .false.
n_tran3 = 0

!  The (1,2,3) indexes of the tran variables are related to
!  the 1st, 2nd and 3rd transformations from to atomic to the
!  molecular orbital basis. l_tran is true when the transformation 
!  to a given molecular orbital is required.

!  Elastic case:

!  SE terms
l_tran1(1:nocc) = .true.
do i=1,ng
  l_tran1(ngt(i)) = .true.
end do

do j=nocc+1,nbfns
  if( l_tran1(j) ) then
    do i=1,nocc
      l_tran2(i,j) = .true.
      l_tran2(j,i) = .true.
      do k=1,nocc
        l_tran3(k,i,j) = .true.
        l_tran3(k,j,i) = .true.
      end do
    end do
  end if
end do

!  SEP terms
do ii=1,ns
  na = nht(ii)
  i  = npt(ii)
  do jj=1,ng2(ii)
    ij1=numsd(jj,ii,1)
    ij2=numsd(jj,ii,2)
    if( (ij1.ne.0).or.(ij2.ne.0) ) then
      j= ngt2(jj,ii)
      l_tran1(i) = .true.
      l_tran1(j) = .true.
      l_tran2(i,j) = .true.
      l_tran2(j,i) = .true.
      l_tran3(na,i,j) = .true.
      l_tran3(na,j,i) = .true.
    end if
  end do
end do


!  And for more than one open channel:
if( nchl.gt.1 ) then

!     do ich = 2,nchl

!  GG terms
  l_tran1(1:nocc) = .true.
  do i=1,ng
    l_tran1(ngt(i)) = .true.
  end do
  do j=1,nbfns
    if( l_tran1(j) ) then
      do i=1,nocc
        l_tran2(i,j) = .true.
        l_tran2(j,i) = .true.
        do k=1,nocc
          l_tran3(k,i,j) = .true.
          l_tran3(k,j,i) = .true.
        end do
      end do
    end if
  end do

!     end do

!  GS terms, for singlets:
if( nchls.gt.0) then
  do i=1,ng
    ii=ngt(I)
    do l=1,nsa
      na = nht(l)
      i1 = npt(l)
      l_tran1(na) = .true.
      l_tran2(max0(i1,ii),na) = .true.
      l_tran3(min0(i1,ii),max0(i1,ii),na) = .true.
    end do
  end do
end if

!  GS terms, for triplets:
if( nchl.gt.(nchls+1) ) then
  do i=1,ng
    ii=ngt(I)
    do l=1,nsa
      na = nht(l)
      i1 = npt(l)
      l_tran1(na) = .true.
      l_tran2(max0(i1,ii),na) = .true.
      l_tran3(min0(i1,ii),max0(i1,ii),na) = .true.
      l_tran1(ii) = .true.
      l_tran2(i1,ii) = .true.
      l_tran3(na,i1,ii) = .true.
    end do
  end do
end if

!  SS terms:
!  For the moment, this part is not as case specific
!  as for the elastic channel

allocate( l_virt_act(nbfns) )
l_virt_act = .false.
l_virt_act(1:nocc) = .true.
do ii=1,ns
  nb=nht(ii)
  i=npt(ii)
  l_virt_act(i) = .true.
  do jj=1,ng2(ii)
    ij1=numsd(jj,ii,1)
    ij2=numsd(jj,ii,2)
    if( (ij1.ne.0).or.(ij2.ne.0) ) then
      j = ngt2(jj,ii)
      l_virt_act(j) = .true.
    end if
  end do
end do

n_virt_act = 0
do i=1,nbfns
  if( l_virt_act(i) ) n_virt_act = n_virt_act + 1
end do

allocate( virt_act(n_virt_act) )

n_virt_act = 0
do i=1,nbfns
  if( l_virt_act(i) ) then
    n_virt_act = n_virt_act + 1
    virt_act(n_virt_act) = i
  end if
end do

!     do ich = 2, nchl
 do j2=1,n_virt_act
   j = virt_act(j2)
   l_tran1(j) = .true.
   do i2=1,j2
     i = virt_act(i2)
     l_tran2(i,j) = .true.
     l_tran2(j,i) = .true.
     do k2=1,i2
       k = virt_act(k2)
       l_tran3(k,i,j) = .true.
       l_tran3(k,j,i) = .true.
       l_tran3(j,i,k) = .true.
       l_tran3(j,k,i) = .true.
       l_tran3(i,j,k) = .true.
       l_tran3(i,k,j) = .true.
     end do
   end do
 end do
!     end do

end if ! if( nchl.gt.1 ) then

! test
!     l_tran1 = .true.
!     l_tran2 = .true.
!     l_tran3 = .true.

!  n_tran counts the number of transformations
!     do ich = 1, nchl
 j2 = 0
 do j=1,nbfns
   if( l_tran1(j) ) j2 = j2 + 1
   i2 = 0
   do i=1,nbfns
     if( l_tran2(i,j) ) i2 = i2 + 1
     k2 = 0
     do k=1,nbfns
       if( (l_tran3(k,i,j)).and.(k.le.i) ) k2 = k2 + 1
     end do
     n_tran3(i,j) = k2
   end do
   n_tran2(j) = i2
 end do
 n_tran1 = j2
!     end do

allocate( tran1(n_tran1) )
allocate( tran2(maxval(n_tran2),n_tran1) )
allocate( tran3(maxval(n_tran3),maxval(n_tran2),n_tran1) )

!  Finally, the orbitals that show up in the transformations
!  are assigned to tran1, tran2 and tran3, which are the arrays
!  to be employed within the plane subroutine
!     do ich = 1, nchl
 j2 = 0
 do j=1,nbfns
   if( l_tran1(j) ) then
     j2 = j2 + 1
     tran1(j2) = j
   end if
   i2 = 0
   do i=1,nbfns
     if( l_tran2(i,j) ) then
       i2 = i2 + 1
       tran2(i2,j2) = i
     end if
     k2 = 0
     do k=1,nbfns
       if( (l_tran3(k,i,j)).and.(k.le.i) ) then
        k2 = k2 + 1
        tran3(k2,i2,j2) = k
       end if
     end do
   end do
 end do

deallocate( l_tran1, l_tran2, l_tran3 )
if( nchl.gt.1 ) deallocate( l_virt_act, virt_act )

end subroutine ao2mo_tran
!----------------------------------------------------------------

end module quadsym
!----------------------------------------------------------------
