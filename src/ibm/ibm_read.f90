#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/07/21 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   called once in dns_main.f90 before time integration starts,
!#   reads all relevant ibm parameters and with consistency check 
!#    
!# 
!########################################################################
!# ARGUMENTS 
!#
!# 
!########################################################################
!# REQUIREMENTS
!#
!# 
!########################################################################

subroutine IBM_READ_INI(inifile)

  use TLAB_CONSTANTS, only : efile, lfile
  use TLAB_PROCS,     only : TLAB_STOP, TLAB_WRITE_ASCII
  use IBM_VARS
  
  implicit none

  character*(*), intent(in) :: inifile

  character                 :: bakfile*32, sRes*512

  ! ================================================================== !
  ! initialization 
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'
  call TLAB_WRITE_ASCII(lfile, 'Reading IBM input data from tlab.ini.')

! IBM Status Parameter
  call TLAB_WRITE_ASCII(bakfile, '#')
  call TLAB_WRITE_ASCII(bakfile, '#[IBMParameter]')
  call TLAB_WRITE_ASCII(bakfile, '#Status=<on/off>')
  call SCANINICHAR(bakfile, inifile, 'IBMParameter', 'Status', 'off', sRes)
  if (trim(adjustl(sRes)) == 'off') then; imode_ibm = 0
  else if (trim(adjustl(sRes)) == 'on') then; imode_ibm = 1
  else
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. Wrong IBM Status option.')
      call TLAB_STOP(DNS_ERROR_OPTION)
  end if
  
  ! proceed if the IBM is turned on
  if (imode_ibm == 1) then 
    ! read IBM parameters
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[IBMParameter]')
    call TLAB_WRITE_ASCII(bakfile, '#IBMScalar=<on/off>')
    call TLAB_WRITE_ASCII(bakfile, '#RestartGeometry=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#DataTypeGeometry=<real/int/bit>')
    call TLAB_WRITE_ASCII(bakfile, '#MaxNumberObj=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#FluidPoints=<value>')

    call SCANINICHAR(bakfile, inifile, 'IBMParameter', 'IBMScalar', 'off', sRes)
    if      (TRIM(ADJUSTL(sRes)) == 'off') then; imode_ibm_scal = 0
    else if (TRIM(ADJUSTL(sRes)) == 'on' ) then; imode_ibm_scal = 1
    else
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. Wrong IBMScalar option.')
      call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'IBMParameter', 'RestartGeometry', 'no', sRes)
    if      ( TRIM(ADJUSTL(sRes)) == 'yes' ) then; ibm_restart = .true.
    else if ( TRIM(ADJUSTL(sRes)) == 'no'  ) then; ibm_restart = .false.
    else
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. Wrong IBM Restart option.')
      call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'IBMParameter', 'DataTypeGeometry', 'int', sRes)
    if      ( TRIM(ADJUSTL(sRes)) == 'real' ) then; ibm_io = IBM_IO_REAL
    else if ( TRIM(ADJUSTL(sRes)) == 'int'  ) then; ibm_io = IBM_IO_INT
    else if ( TRIM(ADJUSTL(sRes)) == 'bit'  ) then; ibm_io = IBM_IO_BIT
    else
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. Wrong IBM Data type option.')
      call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINIINT(bakfile, inifile, 'IBMParameter', 'MaxNumberObj', '0', nob_max)

    call SCANINIINT(bakfile, inifile, 'IBMParameter', 'FluidPoints', '3', nflu)

    ! read geometry parameters
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[IBMGeometry]')
    call TLAB_WRITE_ASCII(bakfile, '#Type=<none,XBars,Hill,Valley,Box>')       
    call TLAB_WRITE_ASCII(bakfile, '#Mirrored=<yes/no>')       
    call TLAB_WRITE_ASCII(bakfile, '#Number=<value>') ! max number of elements in one spatial direction
    call TLAB_WRITE_ASCII(bakfile, '#Height=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Width=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Alpha=<value>')
    
    call SCANINICHAR(bakfile, inifile, 'IBMGeometry', 'Type', 'none', sRes)
    if ( TRIM(ADJUSTL(sRes)) == 'none' )  then
      if ( .not. ibm_restart ) then
        call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. No IBM geometry available.')
        call TLAB_STOP(DNS_ERROR_OPTION)
      end if
      continue
    else if ( TRIM(ADJUSTL(sRes)) == 'xbars' ) then; ibm_geo%name = 'xbars'
    else if ( TRIM(ADJUSTL(sRes)) == 'hill'  ) then; ibm_geo%name = 'hill' 
    else if ( TRIM(ADJUSTL(sRes)) == 'valley') then; ibm_geo%name = 'valley'
    else if ( TRIM(ADJUSTL(sRes)) == 'box'   ) then; ibm_geo%name = 'box'
    else
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. Wrong IBMGeometryType option.')
      call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'IBMGeometry', 'Mirrored', 'no', sRes)
    if      ( TRIM(ADJUSTL(sRes)) == 'yes' ) then; ibm_geo%mirrored = .true.
    else if ( TRIM(ADJUSTL(sRes)) == 'no'  ) then; ibm_geo%mirrored = .false.; end if
    call SCANINIINT(bakfile, inifile, 'IBMGeometry', 'Number', '0', ibm_geo%number)
    call SCANINIINT(bakfile, inifile, 'IBMGeometry', 'Height', '0', ibm_geo%height)
    call SCANINIINT(bakfile, inifile, 'IBMGeometry', 'Width',  '0', ibm_geo%width)
    call SCANINIINT(bakfile, inifile, 'IBMGeometry', 'Alpha',  '0', ibm_geo%hill_slope)
  end if

  return
end subroutine IBM_READ_INI

!########################################################################

subroutine IBM_READ_CONSISTENCY_CHECK()

  use TLAB_CONSTANTS, only : efile, MAX_VARS, wi, wp
  use TLAB_VARS,      only : imax, g
  use TLAB_PROCS,     only : TLAB_STOP, TLAB_WRITE_ASCII
  use IBM_VARS
  
  implicit none

  ! ================================================================== !
  ! consistency check of IBM input data
  if ( g(3)%size == 1 ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Not implemented for 2D-Flow.')
    call TLAB_STOP(DNS_ERROR_OPTION)
  end if
  if ( ibm_io == IBM_IO_BIT ) then
    if ( wp == 8 ) then
      if ( mod( imax, 8 ) /= 0 ) then
        call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. IBM_IO bitwise not possible, restriction: mod(imax/8)=0.')
        call TLAB_STOP(DNS_ERROR_OPTION)
      end if
    else
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. IBM_IO bitwise not tested/changed yet for "working precision = single precision".')
      call TLAB_STOP(DNS_ERROR_OPTION)
    end if
  end if
  if ( nob_max <= 0 ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. No objects in flow, set number of objects correctly.')
    call TLAB_STOP(DNS_ERROR_OPTION)
  end if 
  if ( nflu < 2 ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Too less FluidPoints (nflu>=2).')
    call TLAB_STOP(DNS_ERROR_OPTION)
  end if
  if ( ( ibm_geo%name == 'xbars' ) .or. &
       ( ibm_geo%name == 'hill'  ) .or. &
       ( ibm_geo%name == 'valey' ) .or. &
       ( ibm_geo%name == 'box' ) ) then
    if ( ( mod(g(3)%size,2*ibm_geo%number) == 0 ) .and. ( mod(ibm_geo%width,2) /= 0 ) ) then
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Interfaces of bars have to be on gridpoints.')
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Requirenments: mod(kmax_total,(2*nbars))==0 & mod(wbar,2)==0.')
      call TLAB_STOP(DNS_ERROR_UNDEVELOP)
    else if ( ( mod(g(3)%size,2*ibm_geo%number) /= 0 ) .and. & 
              ( mod(real(g(3)%size/(2*ibm_geo%number), wp),0.5) == 0 ) .and. &
              ( mod(ibm_geo%width,2) /= 1) ) then
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Interfaces of bars have to be on gridpoints.')
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. Requirenments: mod(kmax_total/(2*nbars),0.5)==0 & mod(wbar,2)==1.')
      call TLAB_STOP(DNS_ERROR_UNDEVELOP)    
    end if
    if ( ibm_geo%mirrored .and. imode_ibm_scal == 1 ) then
      call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. No IBM for scalars possible with objects on upper domain.')
      call TLAB_STOP(DNS_ERROR_UNDEVELOP)   
    end if
  else if ( ( ibm_geo%name == 'none' ) .and. .not. ibm_restart ) then
    call TLAB_WRITE_ASCII(efile, 'IBM_READ_INI. IBM. No IBM geometry defined.')
    call TLAB_STOP(DNS_ERROR_UNDEVELOP)    
  end if

  return
end subroutine IBM_READ_CONSISTENCY_CHECK

!########################################################################