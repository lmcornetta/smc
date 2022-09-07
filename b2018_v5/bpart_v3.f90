program bpart

use paramB,  only: read_flags_b, typrun, setup_run, read_ene
use legacy,  only: nchl, nchls, eg, rdup_wfn, rdup_orb, rdup_hci, rdup_pse
use quadsym, only: sctin1, sctin2, molsym, nok, kkkkey, nokmax
use onk3dk,  only: onk_driver, off_driver, d3k_driver

implicit none

!  Read I/O control flags
!
call read_flags_b

! Read up PP arrays
!
call rdup_pse

! Read up array dimension parameters
!
call rdup_orb(0)

! Read up wave function arrays
!
call rdup_wfn

! Read up target CI matrix
!
call rdup_hci

! Read up electronic energy
!
call read_ene(eg(1))

!  Set up angular and radial quadratures
!
call sctin1
write(6,"(//,'*** SCTIN1 DONE ***')")

!  Set up run-dependent parameters
!
call setup_run(nchl,nchls,nok,kkkkey,nokmax)
write(6,"(//,'*** SETUP_RUN DONE ***')")

!  Set up plane waves
!
call sctin2
write(6,"(//,'*** SCTIN2 DONE ***')")

!  Orbital symmetry
!
call molsym
write(6,"(//,'*** MOLSYM DONE ***')")

! Call the driver of each type of run
!
if(     typrun .eq. "ONK") then
  call onk_driver
  write(6,"(//,'*** ONK_DRIVER DONE ***')")
else if(typrun .eq. "OFF") then
  call off_driver
  write(6,"(//,'*** OFF_DRIVER DONE ***')")
else if(typrun .eq. "3DK") then
  call d3k_driver
  write(6,"(//,'*** 3DK_DRIVER DONE ***')")
end if

end program bpart
