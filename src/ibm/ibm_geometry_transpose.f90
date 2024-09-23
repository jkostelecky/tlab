#ifdef USE_MPI 
#include "dns_const_mpi.h"
#endif

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   transposes geometry field eps (epsi, epsj, epsk)
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

subroutine IBM_GEOMETRY_TRANSPOSE(epsi, epsj, epsk, tmp)
  
  use IBM_VARS
  use TLAB_VARS,      only : g, imax, jmax, kmax, isize_field 
  use TLAB_CONSTANTS, only : wi, wp
#ifdef USE_MPI
  use MPI 
  use TLabMPI_VARS,  only : ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  use TLabMPI_VARS,  only : ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  use TLabMPI_VARS,  only : ims_npro_i, ims_npro_k
  use TLabMPI_VARS,  only : ims_size_i 
  use TLabMPI_PROCS
#endif

  implicit none
  
  real(wp), dimension(isize_field), intent(  out) :: epsi, epsj, epsk
  real(wp), dimension(isize_field), intent(inout) :: tmp

#ifdef USE_MPI 
  integer(wi), parameter                          :: idi = TLabMPI_I_PARTIAL 
  integer(wi), parameter                          :: idk = TLabMPI_K_PARTIAL 
#endif
  integer(wi)                                     :: nyz, nxy

  ! ================================================================== !
  ! MPI  and local transposition in x
#ifdef USE_MPI
  if ( ims_npro_i > 1 ) then
    call TLabMPI_TRPF_I(eps, tmp, ims_ds_i(1,idi), ims_dr_i(1,idi), ims_ts_i(1,idi), ims_tr_i(1,idi))
    nyz = ims_size_i(idi)
  else
#endif
  tmp = eps
  nyz = jmax * kmax 
#ifdef USE_MPI
  end if
#endif

#ifdef USE_ESSL
  call DGETMO       (tmp, g(1)%size, g(1)%size, nyz,       epsi, nyz)
#else
  call DNS_TRANSPOSE(tmp, g(1)%size, nyz,       g(1)%size, epsi, nyz)
#endif
  ! -------------------------------------------------------------------
  ! local transposition in y
  nxy = imax * jmax
#ifdef USE_ESSL
  call DGETMO       (eps, nxy, nxy, kmax, epsj, kmax)
#else
  call DNS_TRANSPOSE(eps, nxy, kmax, nxy, epsj, kmax)
#endif
  ! -------------------------------------------------------------------
  ! MPI transposition in z
#ifdef USE_MPI
  if ( ims_npro_k > 1 ) then
    call TLabMPI_TRPF_K(eps, epsk, ims_ds_k(1,idk), ims_dr_k(1,idk), ims_ts_k(1,idk), ims_tr_k(1,idk))
  else
#endif
  epsk = eps
#ifdef USE_MPI
  end if
#endif
  
return
end subroutine IBM_GEOMETRY_TRANSPOSE

!########################################################################