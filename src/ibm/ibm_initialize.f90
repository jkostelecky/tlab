#include "types.h"
#include "dns_error.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   called once in dns_main.f90 before time integration starts,
!#   all relevant geometry informations needed for the IBM are 
!#   available after calling this routine
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

subroutine IBM_INITIALIZE_GEOMETRY(txc, wrk3d)  
  
  use DNS_IBM
  use TLAB_VARS,      only : isize_field, inb_txc
  use TLAB_VARS,      only : istagger
  use TLAB_CONSTANTS, only : efile
  use TLAB_PROCS
  
  implicit none

  TREAL, dimension(isize_field,inb_txc), intent(inout) :: txc
  TREAL, dimension(isize_field),         intent(inout) :: wrk3d

  target                                               :: txc

  TREAL, dimension(:), pointer                         :: epsi, epsj, epsk
#ifdef IBM_DEBUG
  TREAL, dimension(:), pointer                         :: tmp1, tmp2, tmp3
#endif

  ! ================================================================== !
  ! assigning pointer to scratch
  txc = C_0_R; epsi => txc(:,1); epsj => txc(:,2); epsk => txc(:,3)

  ! options for geometry generation
  !   0. use existing eps field (ibm_restart==.true.)
  !   1. generate eps field outside tlab and use as existing eps field with
  !      (ibm_restart==.true.)
  !   2. write own routine to generate geometry 
  !      (cf. IBM_GENERATE_GEOMETRY_XBARS, ibm_restart==.false.)
  if ( ibm_restart ) then
    call IBM_IO_READ_GEOMETRY(wrk3d)
  else if ( xbars_geo%name == 'xbars' ) then
    call IBM_GENERATE_GEOMETRY_XBARS(wrk3d)
  else
    call TLAB_WRITE_ASCII(efile, 'IBM_GEOMETRY no objects in flow.')
    call TLAB_STOP(DNS_ERROR_IBM_MISS_GEO)
  end if

  ! transpose eps (epsi, epsj, epsk)
  call IBM_GEOMETRY_TRANSPOSE(epsi, epsj, epsk, wrk3d)

  ! generate relevant geometry fields for IBM routines (nobi, nobj, nobk)
  call IBM_GENERATE_GEOMETRY(epsi, epsj, epsk)

  ! verify geometry
  call IBM_VERIFY_GEOMETRY()

  ! generate epsp on pressure mesh
  if ( istagger == 1 ) then
    call IBM_STAGGER_GEOMETRY(eps, epsp, wrk3d)
  end if

  ! check idle procs
  if ( ibm_procs_idle ) then
    call IBM_CHECK_PROCS(epsi, epsj, epsk)
  else
    ims_pro_ibm_x = .true.; ims_pro_ibm_y = .true.; ims_pro_ibm_z = .true.
  end if

#ifdef IBM_DEBUG
  ! io of all geometry fields in debugging mode 
  tmp1 => txc(:,4); tmp2 => txc(:,5); tmp3 => txc(:,6)
  call IBM_GEOMETRY_DEBUG_IO(epsi, epsj, epsk, tmp1, tmp2, tmp3)
  nullify(tmp1, tmp2, tmp3)
#endif

  ! switch to true in routines if IBM is needed
  ibm_burgers = .false. 
  ibm_partial = .false.
  
  ! disassociate pointers
  nullify(epsi, epsj, epsk)

  return
end subroutine IBM_INITIALIZE_GEOMETRY

!########################################################################