#include "types.h"
#include "dns_const.h"

PROGRAM VBURGERS

  USE DNS_CONSTANTS
  USE DNS_GLOBAL

  ! IBM
  USE DNS_IBM
#ifdef USE_MPI
  USE DNS_MPI, ONLY: ims_pro
#endif  

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(:,:),   ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:,:), ALLOCATABLE               :: a, b, c
  TREAL, DIMENSION(:,:),   ALLOCATABLE               :: wrk1d, wrk2d
  TREAL, DIMENSION(:),     ALLOCATABLE               :: wrk3d, tmp1

  TINTEGER                                           :: i, j, k,  bcs(2,2), itera
  TREAL                                              :: dummy, error

  ! IBM
  LOGICAL,  PARAMETER                                :: ibm           = .TRUE.
  LOGICAL                                            :: ibm_allocated
  
  ! time messurement
  TINTEGER                                           :: iter 
  TINTEGER, parameter                                :: num_iter = 100
#ifdef USE_MPI
#include "mpif.h"
  TREAL                                              :: t_start, t_end 
#else
  TINTEGER                                           :: c1, c2, c3, c11, c22, c33  
#endif

  ! IBM
#ifdef USE_MPI 
#else
  TINTEGER, PARAMETER                                :: ims_pro       = 0
#endif

! ###################################################################
  CALL DNS_START

  CALL DNS_READ_GLOBAL('dns.ini')
  
  ! IBM
  IF (ibm) THEN
    imode_ibm     = 1                                              ! IBM on
    ibm_allocated = .FALSE.                                        ! not allocated yet
    xbars_geo(1)  = 4; xbars_geo(2) = g(2)%size; xbars_geo(3) = 10 ! geometry description: xbars_geo(3)=[number,height,width]
    ibm_burgers   = .TRUE.                                         ! use IBM in burgers routines
  END IF

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

  isize_wrk3d = isize_txc_field

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d,inb_wrk1d+1))
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))
  ALLOCATE(a(imax,jmax,kmax),b(imax,jmax,kmax),c(imax,jmax,kmax))
  ALLOCATE(tmp1(isize_txc_field),wrk3d(isize_wrk3d))
  
  ! IBM
  IF (ibm) THEN
    CALL IBM_ALLOCATE(ibm_allocated)
    CALL IBM_INITIALIZE_GEOMETRY(tmp1) ! tmp1 needed for transpose calls
  END IF 

  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

  CALL FI_PROFILES_INITIALIZE(wrk1d)

  bcs = 0

! ###################################################################
! Define forcing term
! ###################################################################
  CALL DNS_READ_FIELDS('flow.0', i1, imax,jmax,kmax, i1,i0, isize_wrk3d, a, wrk3d)

! ###################################################################

  IF (ims_pro == 0) WRITE(*,*) '============ENTER BURGERS ROUTINES======================='
  IF (ims_pro == 0) WRITE(*,*) 'Number of Iterations: ', num_iter

#ifdef USE_MPI
  t_start = MPI_WTIME()
#else
  call system_clock(c1,c2,c3)
#endif

  DO iter = 1, num_iter

      ! ###################################################################
      CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), a,b, c, wrk2d,wrk3d)
      DO k = 1,kmax
         DO j = 1,jmax
            DO i = 1,imax
      !           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
               b(i,j,k) = b(i,j,k) *visc *ribackground(j)- a(i,j,k) *c(i,j,k)
            ENDDO
         ENDDO
      ENDDO

      CALL OPR_BURGERS_X(i0,i0, imax,jmax,kmax, bcs, g(1), a,a,a, c, tmp1, wrk2d,wrk3d)
      c = c -b

      error = C_0_R
      dummy = C_0_R
      DO k = 1,kmax
         DO j = 1,jmax
            DO i = 1,imax
               error = error + c(i,j,k)*c(i,j,k)
               dummy = dummy + b(i,j,k)*b(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      ! IF (ims_pro == 0) WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
      !  CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, e, wrk3d)

      ! ###################################################################
      CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), a,b, c, wrk2d,wrk3d)
      DO k = 1,kmax
         DO j = 1,jmax
            DO i = 1,imax
      !           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
               b(i,j,k) = b(i,j,k) *visc *ribackground(j)- a(i,j,k) *c(i,j,k)
            ENDDO
         ENDDO
      ENDDO

      CALL OPR_BURGERS_Y(i0,i0, imax,jmax,kmax, bcs, g(2), a,a,a, c, tmp1, wrk2d,wrk3d)
      c = c -b

      error = C_0_R
      dummy = C_0_R
      DO k = 1,kmax
         DO j = 1,jmax
            DO i = 1,imax
               error = error + c(i,j,k)*c(i,j,k)
               dummy = dummy + b(i,j,k)*b(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      ! IF (ims_pro == 0) WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
      !  CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, c, wrk3d)

      ! ###################################################################
      IF ( g(3)%size .GT. 1 ) THEN

         CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), a,b, c, wrk2d,wrk3d)
         DO k = 1,kmax
            DO j = 1,jmax
               DO i = 1,imax
         !           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
                  b(i,j,k) = b(i,j,k) *visc *ribackground(j)- a(i,j,k) *c(i,j,k)
               ENDDO
            ENDDO
         ENDDO

         CALL OPR_BURGERS_Z(i0,i0, imax,jmax,kmax, bcs, g(3), a,a,a, c, tmp1, wrk2d,wrk3d)
         c = c -b

         error = C_0_R
         dummy = C_0_R
         DO k = 1,kmax
            DO j = 1,jmax
               DO i = 1,imax
                  error = error + c(i,j,k)*c(i,j,k)
                  dummy = dummy + b(i,j,k)*b(i,j,k)
               ENDDO
            ENDDO
         ENDDO
         ! IF (ims_pro == 0) WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
         !  CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, e, wrk3d)
      END IF

  END DO

#ifdef USE_MPI
  t_end = MPI_WTIME()
  write(*,30) ims_pro, t_end - t_start
  30 format(1X, 'MPI_Task: ', I5, '      delta_t (ms) : ', f7.4) 
#else
  call system_clock(c11,c22,c33)
  write(*,30) ims_pro, c11 - c1 
  30 format(1X, 'MPI_Task: ', I5, '      delta_t (ms) : ', I4) 
#endif

  CALL DNS_STOP(0)

END PROGRAM VBURGERS