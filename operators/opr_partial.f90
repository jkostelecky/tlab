#include "types.h"
#include "dns_const.h"

SUBROUTINE OPR_PARTIAL1(nlines, bcs, g, u,result, wrk2d)

  USE DNS_TYPES, ONLY : grid_dt
  
  IMPLICIT NONE

  TINTEGER,                        INTENT(IN)    :: nlines ! # of lines to be solved
  TINTEGER, DIMENSION(2),          INTENT(IN)    :: bcs    ! BCs at xmin (1) and xmax (2):
                                                           !     0 biased, non-zero
                                                           !     1 forced to zero
  TYPE(grid_dt),                   INTENT(IN)    :: g
  TREAL, DIMENSION(nlines*g%size), INTENT(IN)    :: u
  TREAL, DIMENSION(nlines*g%size), INTENT(OUT)   :: result
  TREAL, DIMENSION(nlines),        INTENT(INOUT) :: wrk2d

! -------------------------------------------------------------------
  TINTEGER ip

! ###################################################################
  IF ( g%periodic ) THEN
     SELECT CASE( g%mode_fdm )
        
     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C1N4P_RHS(g%size,nlines, u, result)
        
     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT ) ! Direct = Jacobian because uniform grid
        CALL FDM_C1N6P_RHS(g%size,nlines, u, result)
        
     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8P_RHS(g%size,nlines, u, result)
        
     END SELECT
     
     CALL TRIDPSS(g%size,nlines, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3),g%lu1(1,4),g%lu1(1,5), result,wrk2d)

! -------------------------------------------------------------------
  ELSE
     SELECT CASE( g%mode_fdm )
        
     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C1N4_RHS(g%size,nlines, bcs(1),bcs(2), u, result)

     CASE( FDM_COM6_JACOBIAN )
        CALL FDM_C1N6_RHS(g%size,nlines, bcs(1),bcs(2), u, result)

     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8_RHS(g%size,nlines, bcs(1),bcs(2), u, result)

     CASE( FDM_COM6_DIRECT   ) ! Not yet implemented
        CALL FDM_C1N6_RHS(g%size,nlines, bcs(1),bcs(2), u, result)

     END SELECT
     
     ip = (bcs(1) + bcs(2)*2)*3 
     CALL TRIDSS(g%size,nlines, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3), result)

  ENDIF
  
  RETURN
END SUBROUTINE OPR_PARTIAL1

! ###################################################################
! ###################################################################
SUBROUTINE OPR_PARTIAL2(nlines, bcs, g, u,result, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : grid_dt
  
  IMPLICIT NONE

  TINTEGER,                        INTENT(IN)    :: nlines     ! # of lines to be solved
  TINTEGER, DIMENSION(2,*),        INTENT(IN)    :: bcs    ! BCs at xmin (1,*) and xmax (2,*):
                                                           !     0 biased, non-zero
                                                           !     1 forced to zero
  TYPE(grid_dt),                   INTENT(IN)    :: g
  TREAL, DIMENSION(nlines,g%size), INTENT(IN)    :: u
  TREAL, DIMENSION(nlines,g%size), INTENT(OUT)   :: result
  TREAL, DIMENSION(nlines),        INTENT(INOUT) :: wrk2d
  TREAL, DIMENSION(nlines,g%size), INTENT(INOUT) :: wrk3d      ! First derivative, in case needed

! -------------------------------------------------------------------
  TINTEGER ip
  
! ###################################################################
! Check whether to calculate 1. order derivative
  IF ( .NOT. g%uniform ) THEN
     IF ( g%mode_fdm .eq. FDM_COM4_JACOBIAN .OR. &
          g%mode_fdm .eq. FDM_COM6_JACOBIAN .OR. &
          g%mode_fdm .eq. FDM_COM8_JACOBIAN      ) THEN
        CALL OPR_PARTIAL1(nlines, bcs, g, u,wrk3d, wrk2d)
     ENDIF
  ENDIF
  
! ###################################################################
  IF ( g%periodic ) THEN
     SELECT CASE( g%mode_fdm )
        
     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C2N4P_RHS(g%size,nlines, u, result)
        
     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT ) ! Direct = Jacobian because uniform grid
        CALL FDM_C2N6P_RHS(g%size,nlines, u, result)
        
     CASE( FDM_COM8_JACOBIAN )                  ! Not yet implemented
        CALL FDM_C2N6P_RHS(g%size,nlines, u, result)
        
     END SELECT
     
     CALL TRIDPSS(g%size,nlines, g%lu2(1,1),g%lu2(1,2),g%lu2(1,3),g%lu2(1,4),g%lu2(1,5), result,wrk2d)

! -------------------------------------------------------------------
  ELSE
     SELECT CASE( g%mode_fdm )
        
     CASE( FDM_COM4_JACOBIAN )
        IF ( g%uniform ) THEN
           CALL FDM_C2N4_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
        ELSE ! Not yet implemented
        ENDIF
     CASE( FDM_COM6_JACOBIAN )
        IF ( g%uniform ) THEN
           CALL FDM_C2N6_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
        ELSE
           CALL FDM_C2N6NJ_RHS(g%size,nlines, bcs(1,2),bcs(2,2), g%jac, u, wrk3d, result)
        ENDIF

     CASE( FDM_COM8_JACOBIAN ) ! Not yet implemented; defaulting to 6. order
        IF ( g%uniform ) THEN
           CALL FDM_C2N6_RHS  (g%size,nlines, bcs(1,2),bcs(2,2),        u,        result)
        ELSE
           CALL FDM_C2N6NJ_RHS(g%size,nlines, bcs(1,2),bcs(2,2), g%jac, u, wrk3d, result)
        ENDIF
        
     CASE( FDM_COM6_DIRECT   )
        CALL FDM_C2N6ND_RHS(g%size,nlines, g%lu2(1,4), u, result)

     END SELECT
     
     ip = (bcs(1,2) + bcs(2,2)*2)*3 
     CALL TRIDSS(g%size,nlines, g%lu2(1,ip+1),g%lu2(1,ip+2),g%lu2(1,ip+3), result)

  ENDIF
  
  RETURN
END SUBROUTINE OPR_PARTIAL2

! ###################################################################
! ###################################################################
! Add here OPR_DX(order,...), OPR_DY(order,...), and OPR_DZ(order,...)
