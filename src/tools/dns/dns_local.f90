#include "types.h"

MODULE DNS_LOCAL
  USE DNS_TYPES,  ONLY : filter_dt, grid_dt
  USE DNS_GLOBAL, ONLY : MAX_NSP
#ifdef USE_PSFFT
  USE NB3DFFT,    ONLY : NB3DFFT_SCHEDLTYPE
#endif
  IMPLICIT NONE
  SAVE

! ###################################################################
! Iteration
! ###################################################################
  TINTEGER :: rkm_mode, rkm_substep, rkm_endstep
  TREAL    :: cfla, cfld, cflr, dtime
  TINTEGER :: nitera_first, nitera_last, nitera_save, nitera_stats, nitera_log, nitera_pln
  TINTEGER :: nitera_stats_spa ! Accumulate statistics in spatial mode

  TINTEGER :: idivergence, imode_rhs

! ###################################################################
! Control
! ###################################################################
  TINTEGER :: ilimit_flow, ilimit_scal
  TREAL    :: p_bound_min, p_bound_max, r_bound_min, r_bound_max ! pressure and density
  TREAL    :: s_bound_min(MAX_NSP), s_bound_max(MAX_NSP)         ! scalars
  TREAL    :: d_bound_max                                        ! dilatation

! ###################################################################
! Variable viscosity
! ###################################################################
  LOGICAL :: flag_viscosity
  TREAL   :: visc_stop, visc_time, visc_rate

! ###########################################################
! Filters
! ###########################################################
  TINTEGER :: FilterDomainStep

! ###########################################################
! Immersed Boundary Method (IBM)
! ###########################################################
  TINTEGER, DIMENSION(3)              :: xbars_geo               ! bars in x, xbars_geo(3)=[nbars,hbar,wbar]
  
! ###################################################################
! Output data
! ###################################################################
  TINTEGER, DIMENSION(3)              :: tower_stride           ! Towers
  TINTEGER                            :: tower_mode

! ###################################################################
  TREAL    :: logs_data(20)

! ###################################################################
! NB3DFFT library
! ###################################################################
#ifdef USE_PSFFT
  TYPE(NB3DFFT_SCHEDLTYPE), SAVE :: nbcsetup
#endif

END MODULE DNS_LOCAL
