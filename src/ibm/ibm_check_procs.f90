#include "types.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   checks IBM status of each proc, whether it is idle (do nothing) 
!#   or active (spline reconstruction) 
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

subroutine IBM_CHECK_PROCS(epsi, epsj, epsk)
  
  use DNS_IBM
  use TLAB_VARS,     only : isize_field
#ifdef USE_MPI 
  use MPI
#ifdef IBM_DEBUG
  use TLAB_MPI_VARS, only : ims_pro
#endif
#endif
  
  implicit none
  
#include "integers.h"

  TREAL, dimension(isize_field), intent(in) :: epsi, epsj, epsk

#ifdef IBM_DEBUG
#ifdef USE_MPI
#else
  TINTEGER, parameter                       :: ims_pro = 0 
#endif
#endif

  ! ================================================================== !
#ifdef USE_MPI 

#ifdef IBM_DEBUG
  if ( ims_pro == 0 ) write(*,*) '======== Check Procs for IBM usage ======================'
#endif

  ! Check in X
  if ( sum(epsi) == 0 ) then
    ims_pro_ibm_x = .false.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' idle   for IBM spline generation in x'
#endif
  else 
    ims_pro_ibm_x = .true.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' active for IBM spline generation in x'
#endif
  end if 

  ! Check in Y
  if ( sum(epsj) == 0 ) then
    ims_pro_ibm_y = .false.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' idle   for IBM spline generation in y' 
#endif
  else 
    ims_pro_ibm_y = .true.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' active for IBM spline generation in y'
#endif
  end if 

  ! Check in Z
  if ( sum(epsk) == 0 ) then
    ims_pro_ibm_z = .false.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' idle   for IBM spline generation in z' 
#endif
  else 
    ims_pro_ibm_z = .true.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' active for IBM spline generation in z'
#endif
  end if 

#else
  ! in case of serial mode: one task with full domain
  ims_pro_ibm_x = .true.; ims_pro_ibm_y = .true.; ims_pro_ibm_z = .true.
#endif

  return
end subroutine IBM_CHECK_PROCS

!########################################################################