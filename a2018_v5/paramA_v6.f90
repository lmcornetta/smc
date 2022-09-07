module param

!  Modified to include the NROUTE flag which controls
!  the execution of the code:
!  NROUTE = 0 => usual run
!  NROUTE = 1 => calculate HAQQD and save
!  NROUTE = 2 => calculate VVDG  and save
!  NROUTE = 3 => calculate VVD   and save
!  NROUTE = 4 => read the above matrices and assemble EAQQD
!

!  Modified to read NROUTE, NRFOCK and NTWOEL with explicit
!  reference to the variables.
!
!  Modified to read orbitals in ALCHEMY or GMESS format.
!  Controlled by the flag ORBFMT='ALCHEM' or 'GAMESS'
!
!  Modified to automatically generate dummy orbitals if 
!  spherical d-shells (5 components) are used. Controlled
!  by flag NDSHEL
!
!  Modified to include the flag NDEBUG. If zero, all the
!  old flags (NWINT, NWDATA, ..., KEYORT, KEYQ) are
!  automatically set to zero. In case NDEBUG<>0, these
!  are read from the input file. The only exception is
!  the useful NWHAQQ (save the (N+1)-electron Hamiltonian),
!  which is input independently.
!
!  Modified to include the flag NLIDEP which should equal
!  the number of linera dependent molecular orbitals dropped
!  in GAMESS (see the QMTTOL keyword in GAMESS manual).
!  The code automatically generates dummy orbitals to account
!  for the dimension of the AO space (like it does for ndshel)
!

!  Modified to include solvent effects employing the COSMO
!  PCM approach as implemented in GAMESS. Controlled by the
!  flag SOLVAT.
!

!  Modified to write down the operator 0.5*(PV+VP). This option
!  is related to NROUTE=5.
!

!  Modified to obtain CIS target states, acording to the former
!  CORINGA code. This option is RELATED to NROUTE=-1 and to the
!  module CIS_TARGET.
!

!  Modified to obtain improved virtual orbitals (IVOs), acording 
!  to the former CORINGA code. This option is RELATED to NROUTE=-2 and to the
!  module CIS_TARGET.
!

!  Modified to obtain modified virtual orbitals (MVOs), acording 
!  to the former CORINGA code. This option is RELATED to NROUTE=-3 and to the
!  module CIS_TARGET.
!

!  Modified to use the new NTWOEL options: parallel transformation (0),
!  serial transformation, save memory mode (2, same as before) and new
!  on-demand mode (3, only transform when required).
!

implicit none

!    Units
integer, save :: n_pseud =  2
integer, save :: n_hci_1 =  3
integer, save :: n_hci_2 =  4
integer, save :: n_orbit =  7
integer, save :: n_wvfun =  8
integer, save :: n_apden =  9
integer, save :: n_hamil = 12
integer, save :: n_fockm = 20
integer, save :: n_nucpp = 21
integer, save :: n_scrat = 22
integer, save :: n_erpmo = 23
integer, save :: n_vvd_1 = 24
integer, save :: n_vvd_2 = 25
integer, save :: n_pvvp2 = 26
integer, save :: n_ivorb_s = 50
integer, save :: n_ivorb_t = 51
integer, save :: n_scrci = 52
integer, save :: n_mvorb = 53
integer, save :: n_final_orb = 54
integer, save :: n_symmt = 55

!    I/O control

integer, save :: NWINT=0,NWDATA=0
integer, save :: NPEAQQ=0,NPHAQQ=0,NPRALL=0
integer, save :: KEYORB=0,KEYORT=0,KEYQ=0
integer, save :: NRINT=0, NRDATA=0
integer, save :: NDEBUG,NWHAQQ
integer, save :: NRFOCK,NTWOEL,NROUTE
integer, save :: NDSHEL,NLIDEP
character(7)  :: ORBFMT, SOLVAT

! projectile
integer, save :: projectile

!    INPUT DATA
character(len=*), parameter :: flag_c1 = '%INI_CONTROL_BLOCK'
character(len=*), parameter :: flag_c2 = '%END_CONTROL_BLOCK'
character(len=*), parameter :: flag_o1 = '%INI_ORBITAL_BLOCK'
character(len=*), parameter :: flag_o2 = '%END_ORBITAL_BLOCK'
character(len=*), parameter :: flag_w1 = '%INI_WAVEFUN_BLOCK'
character(len=*), parameter :: flag_w2 = '%END_WAVEFUN_BLOCK'
character(len=*), parameter :: flag_t1 = '%INI_CISTARG_BLOCK'
character(len=*), parameter :: flag_t2 = '%END_CISTARG_BLOCK'
character(len=*), parameter :: flag_i1 = '%INI_IVORBTL_BLOCK'
character(len=*), parameter :: flag_i2 = '%END_IVORBTL_BLOCK'
character(len=*), parameter :: flag_m1 = '%INI_MVORBTL_BLOCK'
character(len=*), parameter :: flag_m2 = '%END_MVORBTL_BLOCK'

!    COMMON /NMBRS /  PI, PITERM, PITERN, ACRCY, SCALE, ICANON
!     COMMON /INC   /  X3,X5,X7,X9,X11,X13,X15,X17,X19,X21,X23,X25
integer, save :: ICANON
real*8, parameter :: PI = 3.14159265358979D0
real*8, save :: PITERM
real*8, save :: PITERN
real*8, save :: ACRCY, SCALE
real*8, parameter :: X3 = 1.D0/3.D0
real*8, parameter :: X5 = 1.D0/5.D0
real*8, parameter :: X7 = 1.D0/7.D0
real*8, parameter :: X9 = 1.D0/9.D0
real*8, parameter :: X11= 1.D0/11.D0
real*8, parameter :: X13= 1.D0/13.D0
real*8, parameter :: X15= 1.D0/15.D0
real*8, parameter :: X17= 1.D0/17.D0
real*8, parameter :: X19= 1.D0/19.D0
real*8, parameter :: X21= 1.D0/21.D0
real*8, parameter :: X23= 1.D0/23.D0
real*8, parameter :: X25= 1.D0/25.D0

!  Other constans used in elmdnomod
real*8, parameter :: RT2  = 1.4142135623730950488D0
real*8, parameter :: RT3  = 1.7320508075688772935D0
real*8, parameter :: RT6  = 2.4494897427831780982D0

real*8, parameter :: RRT2 = 1.0D0/RT2
real*8, parameter :: RRT3 = 1.0D0/RT3
real*8, parameter :: RRT6 = 1.0D0/RT6
real*8, parameter :: RRT3_d2   = 0.5D0/RT3
real*8, parameter :: RRT3_x3d2 = 1.5D0/RT3
!     real*8, parameter :: RRT2 = 0.707106781186547d0
!     real*8, parameter :: RRT3 = 0.577350269189625d0
!     real*8, parameter :: RRT6 = 0.408248290463863d0
!     real*8, parameter :: RRT3_d2   = 0.288675134594812d0
!     real*8, parameter :: RRT3_x3d2 = 0.866025403784438d0
real*8, parameter :: toHartree = 27.21138386d0

!  Discrepancy between input and calculated eg0
real*8, parameter :: discrep = 1.d-04

contains

!----------------------------------------------------------------
subroutine load_pi

implicit none

!
!  IFORT does not like parameter declarations with SQRT intrinsic.
!  GFORTRAN does not like parameter declarations with
!  floating-point exponents.
!
!  This routine initializes PITERN and PITERM in a
!  compiler-independent fashion
!

PITERM=2.D0/PI**0.5D0
PITERN=PI**1.5D0

end subroutine load_pi
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine read_flags!(title)

implicit none

!real*8, intent(out) :: title(18)

!
!  Local variables
!
character(len=50) :: title
integer :: i
logical :: good = .true.
integer :: stat

!
!  PROJEC=+1 : positron
!  PROJEC=-1 : electron
!
!  NDEBUG=0  : set the old flags to zero. Read them otherwise
!  (the useful NWHAQQ is the only exception)
!
!  NWHAQQ<>0 : save (N+1)-electron Hamiltonian
!
!  NRFOCK<>0 : read the Fock matrix (AO basis) from file 
!
!  NTWOEL=0  : old-style SMC run (ERP/ERPT in memory with N**3 transformation)
!  NTWOEL=1  : save memory in 2-e trasformation (ERP and ERPT in memory)
!  NTWOEL=2  : save memory in transformation and do not keep ERP and ERPT
!  simultaneously in memory (intense disk use, N**4)
!  NTWOEL=3  : ERP and ERPT calculated on-the-fly (little memory use, N**3)
!
!  NROUTE=0  : default run (all subroutine calls)
!  NROUTE=1  : calculate HAQQD and save (used with NRFOCK<>0 and NTWOEL=3)
!  NROUTE=2  : calculate VVDG  and save (used with NRFOCK<>0 and NTWOEL=3)
!  NROUTE=3  : calculate VVD   and save (used with NRFOCK<>0 and NTWOEL=3)
!  NROUTE=4  : assemble  EAQQD and save
!  NROUTE=5  : calculate 0.5*(PV+VP) an save.
!  NROUTE=-1 : CIS calculation for excited target states
!  NROUTE=-2 : IVO calculation 
!  NROUTE=-3 : MVO calculation 
!
!  NDSHEL=5  : use 5-component atomic d shells (no need to input dummy orbitals)
!  NDSHEL=6  : use 6-component atomic d shells
!
!  NLIDEP=N  : number of linear dependent MOs dropped in GAMESS 9se the keyword
!  QMTTOL in gamess manual). It should be zero if no LD MOs are found
!  and non-zero otherwise. Dummy orbitals will be automatically 
!  generated in the latter case (i.e., if LD orbitals are found)
!
!  ORBFMT    : orbital format in input files. This character variable shuold be
!  either "GAMESS " or "ALCHEM " for gamess and alchemy formats
!
!  SOLVAT    : solvation model. The only model implemented so far is the COSMO
!  PCM apparoach (COSPCM). It should be set to VACUUM in case no solvent
!  effects are accounted for.
!
!  Read I/O control flags
!
call ini_flag(flag_c1)

!read(5,500,iostat=stat) (title(I),I=1,18)
!  500 FORMAT(18A4)
read(5,'(a)',iostat=stat) title
call write_stat('title',stat)

! projectile
call read_int_flag('PROJEC =',projectile)
if( projectile.ne.-1 .and. projectile.ne.+1 ) then
  write(6,"('PROJEC SHOULD BE -1 (ELECTRON) OR +1 (POSITRON)')")
  stop'ABORT...'
end if

call read_int_flag('NDEBUG =',ndebug)
if(ndebug .ne. 0) then
  read(5,200,iostat=stat) NWINT,NWDATA
  call write_stat('debug options',stat)
  read(5,200,iostat=stat) NRINT,NRDATA
  call write_stat('debug options',stat)
  read(5,200,iostat=stat) NPEAQQ,NPHAQQ,NPRALL
  call write_stat('debug options',stat)
  read(5,200,iostat=stat) KEYORB,KEYORT,KEYQ
  call write_stat('debug options',stat)
end if

call read_int_flag('NWHAQQ =',nwhaqq)

call read_int_flag('NRFOCK =',nrfock)

if( projectile.eq.+1 .and. nrfock.eq.+1 ) then
  write(6,"('WHEN PROJ EQUALS +1 (POSITRON), NRFOCK MUST BE 0')")
  stop'ABORT...'
end if

call read_int_flag('NTWOEL =',ntwoel)

call read_int_flag('NROUTE =',nroute)

call read_int_flag('NDSHEL =',ndshel)

call read_int_flag('NLIDEP =',nlidep)

call read_char_flag('ORBFMT =',orbfmt)
orbfmt = trim(adjustl(orbfmt))

 good = ( (orbfmt.eq."ALCHEM") .or. (orbfmt.eq."GAMESS")  )
if(.not. good) then
  write(6,"('UNKNOWN ORBITAL FORMAT: ORBFMT =',a7)") orbfmt
  write(6,"('IT SHOULD BE ONE OF: ALCHEM, GAMESS')")
  stop'ABORT...'
end if

call read_char_flag('SOLVAT =',solvat)
solvat = trim(adjustl(solvat))

good = ( (solvat.eq."COSPCM").or.(solvat.eq."VACUUM") )
if(.not. good) then
  write(6,"('UNKNOWN SOLVATION MODEL: SOLVAT =',a7)") solvat
  write(6,"('IT SHOULD BE ONE OF: VACUUM, COSPCM')")
  stop'ABORT...'
end if

good = ( (NDSHEL .eq. 5) .or. (NDSHEL .eq. 6) )
if(.not. good) then
  write(6,"('WRONG NUMBER OF D SHELLS: NDSHEL =',i4)") ndshel
  write(6,"('IT SHOULD BE ONE OF: 5, 6')")
  stop'ABORT...'
end if

if(nlidep .lt. 0) then
  write(6,"('NEGATVE NUMBER OF LD MOs. NLIDEP =',i4)") nlidep
  stop'ABORT...'
end if

if( (ntwoel.lt.0) .or. (ntwoel.gt.3) ) then
  write(6,"('BAD NTWOEL. ONLY 0,1,2,3 ARE ACCEPTABLE')")
  stop 'ABORT...'
end if

if( (nroute.lt.-3) .or. (nroute.gt.5) ) then
  write(6,"('BAD NROUTE. ','ONLY -3,-2,-1,0,1,2,3,4,5 ARE ACCEPTABLE')")
  stop 'ABORT...'
end if

call end_flag(flag_c2)

write(6,'(//a//)') '#################################################&
                  &########################################################################'
write(6,'(a)') title
!WRITE(6,510)(title(I),I=1,18)
!  510 FORMAT(//,'#################################################&
!                 &########################################################################'//,18A4)

!roma
  200 FORMAT(15I4)
 
!=======================================================================
! NOTE: IF NRINT.NE.0  NWINT HAS NO EFFECT.
!=======================================================================
WRITE(6,210) NWINT,NWDATA,NWHAQQ,NRINT,NRDATA
  210 FORMAT(//'*** I-O INFORMATION ***'/&
        '    NWINT=',I2,'   NWDATA=',I2,'   NWHAQQ=',I2,/&
       '    NRINT=',I2,'   NRDATA=',I2)
WRITE(6,220) NPEAQQ,NPHAQQ,NPRALL
  220 FORMAT(//'*** PRINT INFORMATION ***'/&
        '    NPEAQQ=',I2,'   NPHAQQ=',I2,'    NPRALL=',I2)
WRITE(6,230) KEYORB,KEYORT,KEYQ
  230 FORMAT(//'*** KEY INFORMATION ***'/&
        '    KEYORB=',I2,'   KEYORT=',I2,'    KEYQ=',I2)

write(6,"(//,'*** NRFOCK = ',i2,' ***')") NRFOCK
if(nrfock .ne. 0 ) write(6,"('FOCK MATRIX WILL BE READ FROM FILE')")

write(6,"(//,'*** NTWOEL = ',i2,' ***')") NTWOEL
if(ntwoel .eq. 0 )  write(6,"('2-E INTGRALS: PARALLEL AO-TO-MO TRANSFORMATION')")
if(ntwoel .eq. 1 )  write(6,"('2-E INTGRALS: SERIAL AO-TO-MO TRANSFORMATION')")
if(ntwoel .eq. 2 )  write(6,"('2-E INTEGRALS MEMORY SAVE MODE. BE PATIENT...')")
if(ntwoel .eq. 3 )  write(6,"('2-E INTEGRALS ON DEMAND. THE AO-TO-MO BASIS TRANSFORMATION WILL BE PERFORMED ONLY WHEN REQUIRED')")

write(6,"(//,'*** NROUTE = ',i2,' ***')") NROUTE
if(nroute .eq. 0) then 
  write(6,"('DEFAULT RUN (ALL SUBROUTINE CALLS)')")
else if(nroute .eq. 1) then 
  write(6,"('CALCULATE ONLY HAQQD AND SAVE')")
else if(nroute .eq. 2) then 
  write(6,"('CALCULATE ONLY VVDG AND SAVE')")
else if(nroute .eq. 3) then 
  write(6,"('CALCULATE ONLY VVD AND SAVE')")
else if(nroute .eq. 4) then 
  write(6,"('ASSEMBLE EAQQD AND SAVE')")
else if(nroute .eq. 5) then 
  write(6,"('CALCULATE 0.5*(PV+VP) AND SAVE')")
else if(nroute .eq. -1) then 
  write(6,"('CIS CALCULATION FOR TARGET EXCITED STATES')")
else if(nroute .eq. -2) then 
  write(6,"('IVO CALCULATION')")
else if(nroute .eq. -3) then 
  write(6,"('MVO CALCULATION')")
end if

if( (nroute.gt.0) .and. (nroute.lt.4) ) then
  good = (ntwoel.eq.3)
  if(.not. good) then 
    write(6,"('THIS NROUTE SHOULD BE USED WITH NTWOEL = 3')")
    stop'ABORT...'
  end if
  good = (nrfock.ne.0)
  if(.not. good) then 
    write(6,"('THIS NROUTE SHOULD BE USED WITH NRFOCK <> 0')")
    stop'ABORT...'
  end if
end if

if(nroute .eq. 1) nwhaqq = 1

if(ndshel .eq. 5 ) then
  write(6,"(//,'*** USE 5-COMPONENT D SHELLS ***')")
else
  write(6,"(//,'*** USE 6-COMPONENT D SHELLS ***')")
end if

if(nlidep .eq. 0 ) then
  write(6,"(//,'*** NO LINEAR DEPENDENT MOs ***')")
else
  write(6,"(//,'*** LINEAR DEPENDENP MOs ***',/)")
  write(6,"(i4,' ORBITALS DROPPED')"), nlidep
end if

if(ORBFMT.eq."ALCHEM")  then
  write(6,"(//,'*** READ ORBITALS IN ALCHEMY FORMAT ***')")
else
  write(6,"(//,'*** READ ORBITALS IN GAMESS FORMAT ***')")
end if

if(SOLVAT.eq."COSPCM")  then
  write(6,"(//,'*** SOLVATION EFFECTS: COSMO PCM MODEL ***')")
else
  write(6,"(//,'*** SOLVATION EFFECTS NOT INCLUDED ***')")
end if

return

end subroutine read_flags
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine read_char_flag(flag,var)

implicit none

character(len=*), intent(in)  :: flag
character(len=*), intent(out) :: var

character(8) :: f_name
integer :: stat
character(len=*), parameter :: flag_invalid  = 'INVALID FLAG NAME READING '

read(5,"(a8)",advance='no',iostat=stat) f_name
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

character(8) :: f_name
integer :: stat
character(len=*), parameter :: flag_invalid  = 'INVALID FLAG NAME READING '

read(5,"(a8)",advance='no',iostat=stat) f_name
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

!----------------------------------------------------------------
end module param
