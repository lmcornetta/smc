program apart

!***
!  Modified to use the NROUTE flag explained in module param
!***

!***
!  This is the main driver for the code as debugged by MAPL for
!  multuchannel scattering. Major changes were made to the POLY-CONF
!  module.
!***

!***
!  Modified to obtain CIS target states, acording to the former
!  CORINGA code. This option is RELATED to NROUTE=-1 and to the
!  module CIS_TARGET. 
!***

!***
!  Matrix diagonalization now carried out using the LAPACK/BLAS
!  subroutine DSYEVX (no longer EIGRS).
!*** 

!***
!  Modified to obtain target IVOS. This option is RELATED to 
!  NROUTE=-2 and to the module CIS_TARGET. 
!***

!***
!  Modified to obtain target MVOS. This option is RELATED to 
!  NROUTE=-3 and to the module CIS_TARGET. 
!***

use param, only: n_pseud, n_orbit, n_wvfun, n_apden, read_flags, ntwoel, load_pi, nroute, n_hci_1, n_hci_2, &
                 n_ivorb_s, n_ivorb_t, n_scrci, n_mvorb, projectile, open_stat

use poly_gaus, only: atrt, tket, erpt, set_poly, poly, wrdown_orb, wrdown_pse, ntype, ncntr, nfirst, nlast, eta, ntwo, d

use poly_conf, only: nocc, ns, nchl, nchld, npt, nexc2, vvd, vvq, vvdg, vvqg, sctin_new, ivo, fock, &
                     ident_electron, ident_positron, sctin0, apresult, pvresult, & 
                     elmdn0mod_electron, elmdn0_positron, vva00smod_electron, vvamulmod_electron, vva00s_positron, vvamul_positron,& 
                     wrdown_wfn, wrdown_den, read_all, elmbmod, load_spin

use cis_target, only: target_cis, target_ivo, target_mvo

implicit none

!***
!  Local variables
 
!real*8  :: title(18)
!integer :: stat
!***

!***
!  Initialize variables piterm and pitern
call load_pi
!***

!***
!  Read I/O control flags
call read_flags!(title)
!***

!***
! Open files used in part B (scattering calculations)
!if( nroute .ge. 0 ) then
! open( unit=n_pseud, file='PSE.FILE', form='unformatted', iostat=stat)
! call open_stat('PSE.FILE',stat)
! open( unit=n_orbit, file='ORB.FILE', form='unformatted', iostat=stat)
! call open_stat('ORB.FILE',stat)
! open( unit=n_wvfun, file='WFN.FILE', form='unformatted', iostat=stat)
! call open_stat('WFN.FILE',stat)
! open( unit=n_apden, file='DEN.FILE', form='unformatted', iostat=stat)
! call open_stat('DEN.FILE',stat)
!end if
!***

!***
!  Let's get it on

! Read up some information used in poly
call sctin_new
write(6,"(//,'*** SCTIN_NEW DONE ***')")

!fk
! Orbital integrals
call set_poly(nocc, nchl, nchld, ns, npt, nexc2)
write(6,"(//,'*** SET_POLY DONE ***')")

! Read up wave function information
call sctin0
write(6,"(//,'*** SCTIN0 DONE ***')")

! Load numsd and numsq
if( projectile.eq.-1 ) then
  call ident_electron
else
  call ident_positron
end if
write(6,"(//,'*** IDENT DONE ***')")

if( projectile.eq.-1 ) then
! Load spin coefficients
  call load_spin
end if

! Write down orbital-related arrays used in part B and deallocate
if((nroute.eq.0) .or. (nroute.eq.4)) then
  call wrdown_orb
  write(6,"(//,'*** WRDOWN_ORB DONE ***')")
end if

! Write down PP-related arrays used in part B and deallocate
if((nroute.eq.0) .or. (nroute.eq.4)) then
  call wrdown_pse
  write(6,"(//,'*** WRDOWN_PSE DONE ***')")
end if

! Write down wave-function-related arrays used in part B
if((nroute.eq.0) .or. (nroute.eq.4)) then
  call wrdown_wfn!(title)
  write(6,"(//,'*** WRDOWN_WFN DONE ***')")
end if

!fk
call poly(ns)
write(6,"(//,'*** POLY DONE ***')")

! CIS target calculation

If( nroute .lt. 0) Then

  call fock
  write(6,"(//,'*** FOCK DONE ***')")

  if( nroute .eq. -1 ) then
!   open( unit=n_hci_1, file='HCI1.FILE', form='unformatted', iostat=stat )
!   call open_stat('HCI1.FILE',stat)
!   open( unit=n_hci_2, file='HCI2.FILE', form='unformatted', iostat=stat )
!   call open_stat('HCI2.FILE',stat)
!   open( unit=n_scrci, file='CIS.SCR', form='unformatted', iostat=stat )
!   call open_stat('CIS.SCR',stat)
    call target_cis
    write(6,"(//,'*** TARGET_CIS DONE ***')")
  else if( nroute .eq. -2 ) then
!   open( unit=n_scrci, file='IVO.SCR', form='unformatted', iostat=stat )
!   call open_stat('IVO.SCR',stat)
!   open( unit=n_ivorb_s, file='IVOS_SING.LIS', form='formatted', iostat=stat )
!   call open_stat('IVO_SING.LIS',stat)
!   open( unit=n_ivorb_t, file='IVOS_TRIP.LIS', form='formatted', iostat=stat )
!   call open_stat('IVO_TRIP.LIS',stat)
    call target_ivo
    write(6,"(//,'*** TARGET_IVO DONE ***')")
  else if( nroute .eq. -3 ) then
!   open( unit=n_mvorb, file='MVOS.LIS', form='formatted', iostat=stat )
!   call open_stat('MVOS.LIS',stat)
    call target_mvo
    write(6,"(//,'*** TARGET_MVO DONE ***')")
  end if

  stop
End If

! Fock matrix
call fock
write(6,"(//,'*** FOCK DONE ***')")
if( projectile.eq.-1 ) then
  if( (allocated(d)) .and. (ntwoel.ne.3) ) deallocate( d )
end if

! Excited-state energies
call ivo
write(6,"(//,'*** IVO DONE ***')")

if( projectile.eq.-1 ) then
  deallocate( atrt, tket )
end if

! Read up wave function information
!     call sctin0
!     write(6,"(//,'*** SCTIN0 DONE ***')")

! Load numsd and numsq
!     call ident
!     write(6,"(//,'*** IDENT DONE ***')")

! Hamiltonian matrix elements
if((nroute.eq.0) .or. (nroute.eq.1)) then
  if( projectile.eq.-1 ) then
    call elmdn0mod_electron
  else
    call elmdn0_positron
  end if
  write(6,"(//,'*** ELMDN0 DONE ***')")
end if

! Prepare V elements (VVDG and VVQG arrays)
!fk   if((nroute.eq.0) .or. (nroute.eq.2)) then
if((nroute.eq.0) .or. (nroute.eq.2) .or. (nroute.eq.5)) then
  if( projectile.eq.-1 ) then
    call vva00smod_electron
  else
    call vva00s_positron
  end if
  write(6,"(//,'*** VVA00SMOD DONE ***')")
end if

! Prepare V elements (VVD and VVQ arrays) 
if(nchl .ge. 2) then
  if((nroute.eq.0) .or. (nroute.eq.3)) then
    if( projectile.eq.-1 ) then
      call vvamulmod_electron
    else
      call vvamul_positron
    end if
    write(6,"(//,'*** VVAMULMOD DONE ***')")
  end if
end if

! No more calls to rpmol
if( allocated(erpt) ) deallocate( erpt )
if( allocated(ntwo) ) deallocate( ntwo )

if( projectile.eq.-1 ) then

! Nothing else to do if NROUTE = 1,2,3
if((nroute.gt.0) .and. (nroute.lt.4)) goto 999

! Allocate and read auxiliary matrices to assemble denominator
if(nroute .eq. 4) call read_all

! Denominator matrix elements
call apresult
write(6,"(//,'*** APRESULT DONE ***')")
!     deallocate( haqqd )
!     if( allocated( haqqq ) ) deallocate( haqqq )

! Denominator matrix elements
call pvresult
write(6,"(//,'*** PVRESULT DONE ***')")

! No call to subroutine ONK (not implemented in mobsci)
deallocate( vvdg )
if( allocated( vvd ) )   deallocate( vvd )
if( allocated( vvq ) )   deallocate( vvq )
if( allocated( vvqg ) )  deallocate( vvqg )

if(ntwoel .eq. 3) deallocate( ntype, ncntr, nfirst, nlast, eta, d )

!fk Matrix elements where the energy dependence is linear
call elmbmod
write(6,"(//,'*** ELMBMOD DONE ***')")

end if

! Write down energy independent and dependent denominator and deallocate
call wrdown_den
write(6,"(//,'*** WRDOWN_DEN DONE ***')")

 999  continue

write(6,"(//,10x,'*** GOOD LUCK IN PART B ***')")
!***

end program apart
