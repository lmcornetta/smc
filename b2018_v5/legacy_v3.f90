module legacy

use paramB

use precision

implicit none

!***
!  Modifications introduced by Fábris Kossoski in the begining of 2017.
!
!  Previous scalar 'ene' has been transformed into an array,
!  which contains all the 'nener' energy points of onk and 3dk runs.
!  For off runs, 'ene' has dimension 1.
!
!  Also, subroutine rdup_wfn now reads matrix 'eaqqd2' (computed in part
!  A), which is the denominator term linear in the collision energy.
!***

!    COMMON /VAI/ IPST,LGUIMA,LPSEUD(30),BLKC1(30),BLKC2(60),BLKC3(60),
!   # BLKI1(270),BLKI2(540)
!  
!   Sendo N a dimensao de LPSEUD (numero de centros), as dimensoes
!   devem seguir: BLKC1(N), BLKC2(2*N), BLKC3(2*N), BLKI1(9*N),
!     BLKI2(18*N)
integer, save :: ncpp, ipst
integer, allocatable, dimension(:), save :: lpseud
real*8,  allocatable, dimension(:), save :: blkc1, blkc2, blkc3, blki1, blki2

!    COMMON /PH/NG,NS,NGT(NGX),NHT(NSX),NPT(NSX),NG2(NSX),NGT2(NGX,NSX)
!   1,NUMSD(NSBFX,NSX,2),NUMSQ(NSBFX,NSX)
integer, save :: ng, ngmax, ns
integer, allocatable, dimension(:), save :: ngt, nht, npt, ng2
integer, allocatable, dimension(:,:), save :: ngt2, numsq
integer, allocatable, dimension(:,:,:), save :: numsd

!    COMMON /CHNLTR/NEXC(NPSCX),NCSPN(NPSCX),NCTR(NPCTX),
!   1NCTRIN(2,NCHLX)
integer, allocatable, dimension(:), save :: nexc, ncspn, nctr
integer, allocatable, dimension(:,:), save :: nctrin

!    COMMON /HAMIL/ENE,EG(NCHLX),DEG(NCHLX),XNFC2
!real*8, allocatable, save :: ene(:), eg(:)
real*8, allocatable, dimension(:), save :: eg

!    COMMON /ORTGNL/NOOR,NOH(NSBFX)
integer, save :: noor
integer, allocatable, dimension(:), save :: noh

!    COMMON /HELP1/NELC,NOCC,NCHL,JSCHW,NWHAQQ,NPEAQQ,NPRALL,ND,
!   1NPSC,NOIN(3,NSX),NPCT,NDD,NDQ,KEYQ,NPHAQQ,NSA,NCHLT,NCHLS
integer ,save :: projectile
integer, save :: nelc, nocc, nchl, jschw, nwhaqq, npeaqq, nprall, nd, npsc, npct, ndd, ndq, keyq, nphaqq, nsa, nchlt, nchls, nchld
integer, allocatable, dimension(:,:), save :: noin
integer(i_large), save :: ndd2, ndq2

!    COMMON /PRESV/EAQQD(NDDX,NDDX),EAQQQ(NDQX,NDQX),
!   1HAQQD(NDOFDX,NDOFDX),HAQQQ(NDOFQX,NDOFQX),
!   2CGRN(NBFNX,NBFNX,NCHLX)
real*8, allocatable, dimension(:,:), save :: eaqqd, eaqqq, haqqd, haqqq
real*8, allocatable, dimension(:,:), save :: eaqqd2, eaqqq2

!    COMMON /SPIN/NSP(3,3),SP(3,3,2),NCSPN2(NCHLDX),MST(2,NPSCX),
!   1NAD(NPSCX)
integer, save :: nsp(3,3)
real*8,  save :: sp(3,3,2)
integer, allocatable, dimension(:), save :: ncspn2, nad
integer, allocatable, dimension(:,:), save :: mst

real*8, allocatable, dimension(:), save  :: egn
real*8, allocatable, dimension(:,:,:), save  :: hci

!   Orbital data
integer, allocatable, dimension(:), save :: ntype, nfirst, nlast
integer, allocatable, dimension(:), save :: none1
integer, allocatable, dimension(:,:), save :: nr, iresum
real*8,  allocatable, dimension(:,:), save :: vlist, eta, d, ovpm

!   Maximum number of primitives in one contracted function
integer, save :: n_prim

!   Other relevant information from part A
integer, save :: nsbf, nbfns, nbfn, nmof
integer, save :: ntmx, ngm, ncm, non
integer, allocatable, dimension(:), save :: nsnglt, ntrplt

integer, save :: nbmax, nbmax2

!   Run label
!real*8, save :: title(18)

contains

!----------------------------------------------------------------
subroutine rdup_pse

implicit none

integer :: stat

open( unit=n_pseud, file='PSE.FILE', form='unformatted', iostat=stat)
call open_stat('PSE.FILE',stat)

read(n_pseud,iostat=stat) ipst, ncpp
call write_stat('PSE.FILE',stat)

if(ipst.eq.1) then

!***
! Allocate PP arrays
!
  allocate( lpseud(ncpp) )
  allocate( blkc1(ncpp), blkc2(2*ncpp) )
  allocate( blkc3(2*ncpp), blki1(9*ncpp) )
  allocate( blki2(18*ncpp) )
!***

  read(n_pseud,iostat=stat) lpseud
  call write_stat('PSE.FILE',stat)
  read(n_pseud,iostat=stat) blkc1
  call write_stat('PSE.FILE',stat)
  read(n_pseud,iostat=stat) blkc2
  call write_stat('PSE.FILE',stat)
  read(n_pseud,iostat=stat) blkc3
  call write_stat('PSE.FILE',stat)
  read(n_pseud,iostat=stat) blki1
  call write_stat('PSE.FILE',stat)
  read(n_pseud,iostat=stat) blki2
  call write_stat('PSE.FILE',stat)
end if

!***
! Avoid problems in all-electron-runs
!
if( ipst .eq. 0 ) then
  allocate( lpseud(ncpp) )
  lpseud = 0
end if
!***

close(n_pseud)

end subroutine rdup_pse
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine rdup_orb(flag)

implicit none

integer, intent(in) :: flag

integer :: n_junk
integer :: stat

if( flag .eq. 0 ) then
  open( unit=n_orbit, file='ORB.FILE', form='unformatted', iostat=stat)
  call open_stat('ORB.FILE',stat)

  read(n_orbit,iostat=stat) nbfns, ntmx, ngm, ncm, non, nbfn, nsbf, nmof
  call write_stat('ORB.FILE',stat)
  return
end if

if( flag .eq. 1 ) then

  allocate( ntype(nbfns), nfirst(nbfns), nlast(nbfns) )
  allocate( nr(ntmx,3), none1(nbfns), noh(nbfns) )
  allocate( iresum(ngm,2), eta(ngm,5), vlist(ncm,4) )
  allocate( d(nbfns,nbfns), ovpm(nbfns,nbfns) )

  read(n_orbit,iostat=stat) ntype, nfirst, nlast
  call write_stat('ORB.FILE',stat)
  read(n_orbit,iostat=stat) nr
  call write_stat('ORB.FILE',stat)
  read(n_orbit,iostat=stat) none1
  call write_stat('ORB.FILE',stat)
  read(n_orbit,iostat=stat) noor, noh
  call write_stat('ORB.FILE',stat)
  read(n_orbit,iostat=stat) n_junk
  call write_stat('ORB.FILE',stat)
  read(n_orbit,iostat=stat) iresum
  call write_stat('ORB.FILE',stat)
  read(n_orbit,iostat=stat) vlist
  call write_stat('ORB.FILE',stat)
  read(n_orbit,iostat=stat) eta
  call write_stat('ORB.FILE',stat)
  read(n_orbit,iostat=stat) n_prim
  call write_stat('ORB.FILE',stat)

  call rd2(d,nbfns,nbfns,nbfns,nbfns,n_orbit)
  call rd2(ovpm,nbfns,nbfns,nbfns,nbfns,n_orbit)

!   No need...
  deallocate( iresum, noh )

  close(n_orbit)

  nbmax = nbfns
  nbmax2 = none1(nbmax)+nbmax

end if

end subroutine rdup_orb
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine rdup_wfn

implicit none

integer :: stat

open( unit=n_wvfun, file='WFN.FILE', form='unformatted', iostat=stat)
call open_stat('WFN.FILE',stat)

read(n_wvfun,iostat=stat) projectile
call write_stat('WFN.FILE',stat)

read(n_wvfun,iostat=stat) nwhaqq, npeaqq, nprall, nphaqq, keyq, jschw
call write_stat('WFN.FILE',stat)
read(n_wvfun,iostat=stat) ng, ngmax, ns, nelc, nocc, nchl, nd, ndd, ndq
call write_stat('WFN.FILE',stat)
read(n_wvfun,iostat=stat) nsa, nchls, nchlt, npct, npsc
call write_stat('WFN.FILE',stat)

nchld = nchl - 1

ndd2 = ndd * (ndd+1) / 2
ndq2 = ndq * (ndq+1) / 2

!read(n_wvfun) title
!call write_stat('WFN.FILE',stat)

allocate( eg(nchl) )
read(n_wvfun,iostat=stat) eg
call write_stat('WFN.FILE',stat)
read(n_wvfun,iostat=stat) sp
call write_stat('WFN.FILE',stat)

allocate( ngt(ng) )
read(n_wvfun,iostat=stat) ngt
call write_stat('WFN.FILE',stat)

if( ns .ge. 1 ) then
  allocate( nsnglt(ns), ntrplt(ns) )
  allocate( nht(ns), npt(ns) )
  allocate( ng2(ns), ngt2(ngmax,ns) )
  read(n_wvfun,iostat=stat) nsnglt, ntrplt
  call write_stat('WFN.FILE',stat)
  read(n_wvfun,iostat=stat) nht, npt
  call write_stat('WFN.FILE',stat)
  read(n_wvfun,iostat=stat) ng2, ngt2
  call write_stat('WFN.FILE',stat)

  allocate( numsd(nbfns,ns,2) )
  read(n_wvfun,iostat=stat) numsd
  call write_stat('WFN.FILE',stat)
  if(keyq .ne. 0) then
    allocate( numsq(nbfns,ns) )
    read(n_wvfun,iostat=stat) numsq
    call write_stat('WFN.FILE',stat)
  end if
end if

allocate( nctrin(2,nchl) )

read(n_wvfun,iostat=stat) nctrin
call write_stat('WFN.FILE',stat)
read(n_wvfun,iostat=stat) nsp
call write_stat('WFN.FILE',stat)

if(nchl .ge. 2) then
  allocate( mst(2,npsc), nad(npsc), ncspn2(nchld) )
  allocate( nexc(npsc), ncspn(npsc), nctr(npct)  )
  read(n_wvfun,iostat=stat) nexc, ncspn, nctr
  call write_stat('WFN.FILE',stat)
  read(n_wvfun,iostat=stat) mst, nad, ncspn2
  call write_stat('WFN.FILE',stat)
end if

close(n_wvfun)

end subroutine rdup_wfn
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine rdup_den

implicit none

integer :: alloc_stat

allocate( eaqqd(ndd,ndd), stat=alloc_stat )
call allocation_stat('eaqqd',alloc_stat)
CALL RD2(eaqqd,ndd,ndd,ndd,ndd,n_apden)

allocate( eaqqd2(ndd,ndd), stat=alloc_stat )
call allocation_stat('eaqqd2',alloc_stat)
CALL RD2(eaqqd2,ndd,ndd,ndd,ndd,n_apden)

if(keyq .ne. 0) then
  allocate( eaqqq(ndq,ndq), stat=alloc_stat )
  call allocation_stat('eaqqq',alloc_stat)
  CALL RD2(eaqqq,ndq,ndq,ndq,ndq,n_apden)
  allocate( eaqqq2(ndq,ndq), stat=alloc_stat )
  call allocation_stat('eaqqq2',alloc_stat)
  CALL RD2(eaqqq2,ndq,ndq,ndq,ndq,n_apden)
end if

end subroutine rdup_den
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine free_wfn

implicit none

deallocate( ngt )

if( ns .ge. 1 ) then
  deallocate( nsnglt, ntrplt )
  deallocate( nht, npt )
  deallocate( ng2, ngt2 )
  deallocate( numsd )
  if(keyq .ne. 0) deallocate( numsq )
end if

if(nchl .ge. 2) then
  deallocate( mst,  nad,   ncspn2 )
  deallocate( nexc, ncspn, nctr   )
end if

end subroutine free_wfn
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine rdup_eaqq(c, n)

implicit none

integer, intent(in) :: n
complex*16, intent(out) :: c(n,n)

integer :: i
real*8, allocatable, dimension(:) :: buffer

allocate( buffer(n) )

do i = 1,n
  call rd22(buffer, n, n_wvfun)
  c(:,i) = buffer(:)
end do

deallocate( buffer )

end subroutine rdup_eaqq
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine rdup_hci

implicit none

integer :: nss
integer :: stat

if( nsa .eq. 0 ) return

allocate( hci(nsa,nsa,2), egn(nsa) )

if(nsa.eq.1) then
  hci(1,1,1)=1.0d0
else
  open(unit=n_hci_1,file='HCI1.FILE',form='unformatted',iostat=stat)
  call open_stat('HCI1.FILE',stat)
  read(n_hci_1,iostat=stat)NSS
  call write_stat('nss in HCI1.FILE',stat)
  if(nss.ne.nsa) stop 'NSA IS DIFFERENT FROM NSS IN HCI1.FILE'
  read(n_hci_1,iostat=stat) egn
  call write_stat('egn in HCI1.FILE',stat)
  read(n_hci_1,iostat=stat) hci(:,:,1)
  call write_stat('hci in HCI1.FILE',stat)
  close(n_hci_1)
endif
 
if(nsa.eq.1) then
  hci(1,1,2)=1.0d0
else
  open(unit=n_hci_2,file='HCI2.FILE',form='unformatted',iostat=stat)
  call open_stat('HCI2.FILE',stat)
  read(n_hci_2,iostat=stat)NSS
  call write_stat('nss in HCI2.FILE',stat)
  if(nss.ne.nsa) stop 'NSA IS DIFFERENT FROM NSS IN HCI2.FILE'
  read(n_hci_2,iostat=stat) egn
  call write_stat('egn in HCI2.FILE',stat)
  read(n_hci_2,iostat=stat) hci(:,:,2)
  call write_stat('hci in HCI2.FILE',stat)
  close(n_hci_2)
endif
 
end subroutine rdup_hci
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine rd2(A,NDIM1,NDIM2,N1,N2,itap)

implicit none

integer, intent(in) :: NDIM1, NDIM2, N1, N2, itap
real*8, intent(out)  :: A(NDIM1,NDIM2)

integer :: i, j

do J=1,N2
  read(itap)(A(I,J),I=1,N1)
end do

end subroutine rd2
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine rd22(a, n, itap)

implicit none

integer, intent(in) :: n, itap
real*8, intent(out)  :: a(n)

read(itap) a

end subroutine rd22
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine rd22_z(a, n, itap)

implicit none

integer, intent(in) :: n, itap
complex*16, intent(out) :: a(n)

read(itap) a

end subroutine rd22_z
!----------------------------------------------------------------

end module legacy
!----------------------------------------------------------------
