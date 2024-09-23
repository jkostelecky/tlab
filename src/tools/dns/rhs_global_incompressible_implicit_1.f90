#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 2012/07/10 - C. Ansorge
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Derived from RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1
!#
!# Semi-implicit time stepping (implicit for viscous term)
!#
!#
!########################################################################
subroutine RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_1(kex, kim, kco, &
                                                q, hq, u, v, w, h1, h2, h3, s, hs, &
                                                tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9)

#ifdef USE_OPENMP
    use OMP_LIB
#endif
    use TLAB_CONSTANTS
    use TLAB_VARS, only: g
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: isize_field, isize_txc_field, inb_scal, inb_flow
    use TLAB_VARS, only: scal_on
    use TLAB_VARS, only: visc, schmidt, rossby
    use TLAB_VARS, only: buoyancy, coriolis
    use TLAB_ARRAYS, only: wrk2d, wrk3d
    use TLAB_PROCS
    use TIME, only: dte
    use DNS_LOCAL, only: remove_divergence
    use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    use OPR_PARTIAL
    use OPR_ELLIPTIC
    use FI_SOURCES, only: bbackground, FI_BUOYANCY

    implicit none

    real(wp) kex, kim, kco
    real(wp), dimension(isize_field, *) :: q, hq
    real(wp), dimension(isize_field), intent(INOUT) :: u, v, w, h1, h2, h3
    real(wp), dimension(isize_field, inb_scal), intent(INOUT) :: s, hs
    real(wp), dimension(isize_txc_field), intent(OUT) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9

    target tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, h2, u, v, w

! -----------------------------------------------------------------------
    integer(wi) iq, is, ij, k, nxy, ip, ip_b, ip_t
    integer(wi) ibc, bcs(2, 2)
    real(wp) dummy, visc_exp, visc_imp, visc_tot, diff, alpha, beta

    real(wp), dimension(:), pointer :: p_bcs

    integer, parameter :: i0 = 0

#ifdef USE_BLAS
    integer ilen
#endif

! #######################################################################
    nxy = imax*jmax

#ifdef USE_BLAS
    ilen = isize_field
#endif

    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    visc_exp = kex*visc
    visc_imp = kim*visc
    visc_tot = visc

! SAVE old tendencies until end of implicit substep of integration

    tmp4(1:isize_field) = hq(:, 1)
    tmp5(1:isize_field) = hq(:, 2)
    tmp6(1:isize_field) = hq(:, 3)

! ################################################################################
! SAVE VALUES AT BOUNDARIES
! ################################################################################

    ip = 1
    do k = 1, kmax
        ip_b = (k - 1)*imax*jmax + 1; 
        ip_t = ip_b + (jmax - 1)*imax
        if (BcsFlowJmin%type(1) == DNS_BCS_DIRICHLET) BcsFlowJmin%ref(1:imax, k, 1) = u(ip_b:ip_b + imax - 1)
        if (BcsFlowJmin%type(3) == DNS_BCS_DIRICHLET) BcsFlowJmin%ref(1:imax, k, 3) = w(ip_b:ip_b + imax - 1)
        if (BcsFlowJmax%type(1) == DNS_BCS_DIRICHLET) BcsFlowJmax%ref(1:imax, k, 1) = u(ip_t:ip_t + imax - 1)
        if (BcsFlowJmax%type(3) == DNS_BCS_DIRICHLET) BcsFlowJmax%ref(1:imax, k, 3) = w(ip_t:ip_t + imax - 1)
        do is = 1, inb_scal
            if (BcsScalJmin%type(is) == DNS_BCS_DIRICHLET .and. &
                BcsScalJmax%type(is) == DNS_BCS_DIRICHLET) then
                BcsScalJmin%ref(1:imax, k, is) = s(ip_b:ip_b + imax - 1, is)
                BcsScalJmax%ref(1:imax, k, is) = s(ip_t:ip_t + imax - 1, is)
            else  ! Only Dirichlet BCs implemented for scalar
                call TLAB_WRITE_ASCII(efile, 'Only Dirichlet BCs implemented for scalar in implicit mode')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)
            end if

            if (BcsScalJmin%SfcType(is) == DNS_SFC_STATIC .and. &
                BcsScalJmax%SfcType(is) == DNS_SFC_STATIC) then
                ! Nothing to do
            else
                call TLAB_WRITE_ASCII(efile, 'Only static surface implemented in implicit mode')
            end if
        end do
    end do

! ################################################################################
! EXPLICIT PART OF SUBSTEP - change viscosity accordingly
! ################################################################################
    visc = visc_exp

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
    ! h3   contains explicit nonlinear w-tendency (nonlinear operator N(u_n))
    ! tmp9 contains explicit diffusive tendency for w
    if (g(3)%size > 1) then
        !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),w,    tmp2)
        !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),tmp2, tmp1)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), w, tmp1, tmp2)   ! tmp2 used for BCs below
        do ij = 1, isize_field; 
            h3(ij) = -v(ij)*tmp2(ij)
            tmp9(ij) = visc*tmp1(ij)
        end do

        !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),w,    tmp3)
        !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),tmp3, tmp1)
        call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), w, tmp1, tmp3)
        do ij = 1, isize_field; 
            h3(ij) = h3(ij) - u(ij)*tmp3(ij)
            tmp9(ij) = tmp9(ij) + visc*tmp1(ij)
        end do

        !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),w,    tmp3)
        !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),tmp3, tmp1)
        call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), w, tmp1, tmp3)
! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
        if (buoyancy%active(3)) then
            call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, wrk3d, bbackground)
            dummy = buoyancy%vector(3)
            do ij = 1, isize_field
                h3(ij) = h3(ij) + dummy*wrk3d(ij)
            end do
        end if

! -----------------------------------------------------------------------
! Coriolis (so far, rotation only in the Oy direction)
! -----------------------------------------------------------------------
        if (coriolis%type == EQNS_COR_NORMALIZED) then
            dummy = 1.0_wp/rossby
            do ij = 1, isize_field
                h3(ij) = h3(ij) - w(ij)*tmp3(ij) + dummy*(u(ij) - 1.0_wp)
                tmp9(ij) = tmp9(ij) + visc*tmp1(ij)
            end do
! -----------------------------------------------------------------------
        else
            do ij = 1, isize_field
                h3(ij) = h3(ij) - w(ij)*tmp3(ij)
                tmp9(ij) = tmp9(ij) + visc*tmp1(ij)
            end do
        end if
    end if

! #######################################################################
! Diffusion and convection terms in Ox momentum eqn
! #######################################################################
!   h1 contains explicit nonlinear u-tendency
! tmp7 contains explicit diffusive u-tendency
    !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),u,    tmp2)
    !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),tmp2, tmp1)
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), u, tmp1, tmp2)   ! tmp2 used for BCs below
    do ij = 1, isize_field; 
        h1(ij) = -v(ij)*tmp2(ij)
        tmp7(ij) = visc*tmp1(ij)
    end do

    !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),u,    tmp3)
    !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),tmp3, tmp1)
    call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), u, tmp1, tmp3)
    do ij = 1, isize_field; 
        h1(ij) = h1(ij) - u(ij)*tmp3(ij); 
        tmp7(ij) = tmp7(ij) + visc*tmp1(ij)
    end do

    !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),u,    tmp3)
    !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),tmp3, tmp1)
    call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), u, tmp1, tmp3)

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
    if (buoyancy%active(1)) then
        call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, wrk3d, bbackground)
        dummy = buoyancy%vector(1)
        do ij = 1, isize_field
            h1(ij) = h1(ij) + dummy*wrk3d(ij)
        end do
    end if

! -----------------------------------------------------------------------
! Coriolis. So far, rotation only in the Oy direction.
! -----------------------------------------------------------------------

    if (coriolis%type == EQNS_COR_NORMALIZED) then
        dummy = 1.0_wp/rossby
        do ij = 1, isize_field
            h1(ij) = h1(ij) - w(ij)*tmp3(ij) - dummy*w(ij)
            tmp7(ij) = tmp7(ij) + visc*tmp1(ij)
        end do
! -----------------------------------------------------------------------
    else
        do ij = 1, isize_field
            h1(ij) = h1(ij) - w(ij)*tmp3(ij)
            tmp7(ij) = tmp7(ij) + visc*tmp1(ij)
        end do
    end if

! #######################################################################
! Diffusion and convection terms in Oy momentum eqn
! #######################################################################
!   h2 contains explicit nonlinear v-tendency
! tmp8 contains explicit diffusive v-tendency
    !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),v,    tmp2)
    !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),tmp2, tmp1)
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), v, tmp1, tmp2)
    do ij = 1, isize_field
        h2(ij) = -v(ij)*tmp2(ij)
        tmp8(ij) = visc*tmp1(ij)
    end do

    !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),v,    tmp3)
    !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),tmp3, tmp1)
    call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), v, tmp1, tmp3)
    do ij = 1, isize_field
        h2(ij) = h2(ij) - u(ij)*tmp3(ij)
        tmp8(ij) = tmp8(ij) + visc*tmp1(ij)
    end do

    !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),v,    tmp3)
    !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),tmp3, tmp1)
    call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), v, tmp1, tmp3)

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
    if (buoyancy%active(2)) then
        call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, wrk3d, bbackground)
        dummy = buoyancy%vector(2)
        do ij = 1, isize_field
            h2(ij) = h2(ij) - w(ij)*tmp3(ij) + dummy*wrk3d(ij)
            tmp8(ij) = tmp8(ij) + visc*tmp1(ij)
        end do
    else
        do ij = 1, isize_field
            h2(ij) = h2(ij) - w(ij)*tmp3(ij)
            tmp8(ij) = tmp8(ij) + visc*tmp1(ij)
        end do
    end if

! #######################################################################
! Impose buffer zone as relaxation terms (Flow)
! #######################################################################
    if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
        call BOUNDARY_BUFFER_RELAX_FLOW()
    end if

! Explicit part of time integration
! tmp4-6 contain tendencies, h1-3 contain old tendencies
    do ij = 1, isize_field
        tmp1(ij) = u(ij) + dte*(h1(ij) + tmp7(ij) + kco*tmp4(ij)); 
        tmp2(ij) = v(ij) + dte*(h2(ij) + tmp8(ij) + kco*tmp5(ij)); 
        tmp3(ij) = w(ij) + dte*(h3(ij) + tmp9(ij) + kco*tmp6(ij)); 
    end do
! tmp1-tmp3:  u,v,w updated with explicit scheme
! h1-h3:      explicit tendencies from this substep
! u-w:        old velocities
! tmp7-tmp9:  vacant

! ################################################################################
! ADVECTION-DIFFUSION FOR SCALAR
! ################################################################################
! hs  -> explicit tendency without diffusion
! tmp7-> diffusion contribution to explicit tendency
    if (scal_on) then
        do is = 1, inb_scal
            diff = visc_exp/schmidt(is)
            tmp6 = hs(:, is)
            !CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), s(1,is), tmp4, tmp5)
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s(1, is), tmp5)
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp5, tmp4)
            do ij = 1, isize_field; 
                hs(ij, is) = -u(ij)*tmp5(ij)
                tmp7(ij) = diff*tmp4(ij)
            end do

            !CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s(1,is), tmp4, tmp5)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(1, is), tmp5)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp5, tmp4)
            do ij = 1, isize_field; 
                hs(ij, is) = hs(ij, is) - v(ij)*tmp5(ij)
                tmp7(ij) = tmp7(ij) + diff*tmp4(ij)
            end do

            !CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), s(1,is), tmp4, tmp5)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s(1, is), tmp5)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp5, tmp4)

            do ij = 1, isize_field; 
                hs(ij, is) = hs(ij, is) - w(ij)*tmp5(ij)
                tmp7(ij) = tmp7(ij) + diff*tmp4(ij)
            end do

! #######################################################################
! Impose buffer zone as relaxation terms (Scalar #is)
! #######################################################################
            if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
                call BOUNDARY_BUFFER_RELAX_SCAL_I(is, s(1, is), hs(1, is))
            end if

            do ij = 1, isize_field
                s(ij, is) = s(ij, is) + dte*(tmp7(ij) + hs(ij, is) + kco*tmp6(ij)); 
            end do
            diff = visc_imp/schmidt(is)
            alpha = dte*diff
            beta = -1.0_wp/alpha

            ip_b = 1; ip_t = 1 + imax*kmax
            do k = 1, kmax
                ip = (k - 1)*imax*jmax + 1; wrk2d(ip_b:ip_b + imax - 1, 1) = -alpha*BcsScalJmin%ref(1:imax, k, is); ip_b = ip_b + imax
                ip = ip + (jmax - 1)*imax; wrk2d(ip_t:ip_t + imax - 1, 1) = -alpha*BcsScalJmax%ref(1:imax, k, is); ip_t = ip_t + imax
            end do
            ip_b = 1; ip_t = 1 + imax*kmax
            ! call OPR_Helmholtz_FourierXZ_Factorize(imax, jmax, kmax, g, 0, beta, &
            !                        s(1, is), tmp5, tmp6, wrk2d(ip_b, 1), wrk2d(ip_t, 1))
            call OPR_Helmholtz(imax, jmax, kmax, g, 0, beta, s(1, is), tmp5, tmp6, wrk2d(ip_b, 1), wrk2d(ip_t, 1))

            do ij = 1, isize_field
                s(ij, is) = s(ij, is)*beta
            end do

        end do

    end if

! ################################################################################
! IMPLICIT PART OF SUBSTEP - change viscosity accordingly
! ################################################################################
! TO BE CLEANED: in this section wrk2d is used for the BCS of the implicit solver

    visc = visc_imp
    alpha = dte*visc_imp
    beta = -1.0_wp/alpha

    ip_b = 1; ip_t = 1 + imax*kmax
    do k = 1, kmax
        ip = (k - 1)*imax*jmax + 1; wrk2d(ip_b:ip_b + imax - 1, 1) = -alpha*u(ip:ip + imax - 1); ip_b = ip_b + imax
        ip = ip + (jmax - 1)*imax; wrk2d(ip_t:ip_t + imax - 1, 1) = -alpha*u(ip:ip + imax - 1); ip_t = ip_t + imax
    end do

    ip_b = 1; ip_t = 1 + imax*kmax
    ! call OPR_Helmholtz_FourierXZ_Factorize(imax, jmax, kmax, g, 0, beta, &
    !                                        tmp1, tmp5, tmp6, wrk2d(ip_b, 1), wrk2d(ip_t, 1))
    call OPR_Helmholtz(imax, jmax, kmax, g, 0, beta, tmp1, tmp5, tmp6, wrk2d(ip_b, 1), wrk2d(ip_t, 1))

    do k = 1, kmax
        ip = (k - 1)*imax*jmax + 1; wrk2d(ip_b:ip_b + imax - 1, 1) = -alpha*v(ip:ip + imax - 1); ip_b = ip_b + imax
        ip = ip + (jmax - 1)*imax; wrk2d(ip_t:ip_t + imax - 1, 1) = -alpha*v(ip:ip + imax - 1); ip_t = ip_t + imax
    end do
    ip_b = 1; ip_t = 1 + imax*kmax
    ! call OPR_Helmholtz_FourierXZ_Factorize(imax, jmax, kmax, g, 0, beta, &
    !                                        tmp2, tmp5, tmp6, wrk2d(ip_b, 1), wrk2d(ip_t, 1))
    call OPR_Helmholtz(imax, jmax, kmax, g, 0, beta, tmp2, tmp5, tmp6, wrk2d(ip_b, 1), wrk2d(ip_t, 1))

    do k = 1, kmax
        ip = (k - 1)*imax*jmax + 1; wrk2d(ip_b:ip_b + imax - 1, 1) = -alpha*w(ip:ip + imax - 1); ip_b = ip_b + imax
        ip = ip + (jmax - 1)*imax; wrk2d(ip_t:ip_t + imax - 1, 1) = -alpha*w(ip:ip + imax - 1); ip_t = ip_t + imax
    end do
    ip_b = 1; ip_t = 1 + imax*kmax
    ! call OPR_Helmholtz_FourierXZ_Factorize(imax, jmax, kmax, g, 0, beta, &
    !                                        tmp3, tmp5, tmp6, wrk2d(ip_b, 1), wrk2d(ip_t, 1))
    call OPR_Helmholtz(imax, jmax, kmax, g, 0, beta, tmp3, tmp5, tmp6, wrk2d(ip_b, 1), wrk2d(ip_t, 1))
    
    do ij = 1, isize_field
        u(ij) = tmp1(ij)*beta
        v(ij) = tmp2(ij)*beta
        w(ij) = tmp3(ij)*beta
    end do

! ################################################################################
! END OF IMPLICIT PART OF SUBSTEP - set viscosity back to actual value
! ################################################################################
    visc = visc_tot

! #######################################################################
! Pressure term
! #######################################################################
    if (remove_divergence) then ! remove residual divergence
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, tmp3)
    end if

! -----------------------------------------------------------------------
    do ij = 1, isize_field
        tmp1(ij) = (tmp1(ij) + tmp2(ij) + tmp3(ij))/dte ! forcing term in tmp1
    end do

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    ip_b = 1
    ip_t = imax*(jmax - 1) + 1
    do k = 1, kmax
        p_bcs => v(ip_b:); BcsFlowJmin%ref(1:imax, k, 2) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
        p_bcs => v(ip_t:); BcsFlowJmax%ref(1:imax, k, 2) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
    end do

! pressure in tmp1, Oy derivative in tmp3
    ! select case (imode_elliptic)
    ! case (FDM_COM6_JACOBIAN)
    !     call OPR_Poisson_FourierXZ_Factorize(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)
    ! case (FDM_COM4_DIRECT, FDM_COM6_DIRECT)
    !     call OPR_Poisson_FourierXZ_Direct(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)
    ! end select
    call OPR_Poisson(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)

! horizontal derivatives
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp1, tmp4)

! -----------------------------------------------------------------------
! Add pressure gradient
! -----------------------------------------------------------------------
    do ij = 1, isize_field
        u(ij) = u(ij) - dte*tmp2(ij); h1(ij) = h1(ij) - tmp2(ij)
        v(ij) = v(ij) - dte*tmp3(ij); h2(ij) = h2(ij) - tmp3(ij)
        w(ij) = w(ij) - dte*tmp4(ij); h3(ij) = h3(ij) - tmp4(ij)
    end do

! #######################################################################
! Boundary conditions
! #######################################################################
! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
! Default is Dirichlet -> values at boundary are kept const;
! BcsFlowJmin/Jmax for u and w initialized to old BC values above
    BcsFlowJmin%ref(:, :, 2) = 0.0_wp
    BcsFlowJmax%ref(:, :, 2) = 0.0_wp
    do iq = 1, inb_flow
        ibc = 0
        if (BcsFlowJmin%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsFlowJmax%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
! WATCH OUT - THIS IS REALLY GOING TO GIVE NO-FLUX (instead of constant flux)
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), q(1, iq), &
                                        BcsFlowJmin%ref(1, 1, iq), BcsFlowJmax%ref(1, 1, iq), tmp1)
        end if
    end do

    do is = 1, inb_scal
        ibc = 0
        if (BcsScalJmin%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsScalJmax%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
! WATCH OUT - THIS IS REALLY GOING TO GIVE NO-FLUX (instead of constant flux)
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), s(1, is), &
                                        BcsScalJmin%ref(1, 1, is), BcsScalJmax%ref(1, 1, is), tmp1)
        end if
    end do

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin
! -----------------------------------------------------------------------
    ip_b = 1
    do k = 1, kmax
        u(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 1)
        v(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 2)
        w(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 3)
        do is = 1, inb_scal
            s(ip_b:ip_b + imax - 1, is) = BcsScalJmin%ref(1:imax, k, is)
        end do
        ip_b = ip_b + nxy
    end do

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
    ip_t = imax*(jmax - 1) + 1
    do k = 1, kmax
        u(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 1)
        v(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 2)
        w(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 3)
        do is = 1, inb_scal
            s(ip_t:ip_t + imax - 1, is) = BcsScalJmax%ref(1:imax, k, is)
        end do
        ip_t = ip_t + nxy
    end do

! ! -----------------------------------------------------------------------
! ! Preliminaries
! ! -----------------------------------------------------------------------

!   ibc = 0
! ! Default is Dirichlet -> values at boundary are kept const;
! ! bcs_hb and bcs_ht initialized to old BC values above
!   IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
!   IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
!   IF ( ibc .GT. 0 ) THEN
!      ! WATCH OUT - THIS IS REALLY GOING TO GIVE NO-FLUX (instead of constant flux)
!      CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), u, &
!           bcs_hb(1,1,1),bcs_ht(1,1,1), tmp1)
!      CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), w, &
!           bcs_hb(1,1,2),bcs_ht(1,1,2), tmp1)
!   ENDIF

! ! -----------------------------------------------------------------------
! ! Impose bottom BCs at Jmin
! ! -----------------------------------------------------------------------
!   ip_b =                 1
!   DO k = 1,kmax
!      u(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,1)
!      v(ip_b:ip_b+imax-1) = 0.0_wp               ! no penetration
!      w(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,2);
!      IF ( scal_on ) THEN
!         DO is=1,inb_scal
!            s(ip_b:ip_b+imax-1,is) = BcsScalJmin%ref(1:imax,k,is)
!         ENDDO
!      ENDIF
!      ip_b = ip_b + nxy
!   ENDDO

! ! -----------------------------------------------------------------------
! ! Impose top BCs at Jmax
! ! -----------------------------------------------------------------------
!   ip_t = imax*(jmax-1) + 1
!   DO k = 1,kmax
!      u(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,1)
!      v(ip_t:ip_t+imax-1) = 0.0_wp               ! no penetration
!      w(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,2);
!      IF ( scal_on ) THEN
!         DO is=1,inb_scal
!            s(ip_t:ip_t+imax-1,is) = BcsScalJmax%ref(1:imax,k,is)
!         ENDDO
!      ENDIF
!      ip_t = ip_t + nxy
!   ENDDO

    return
end subroutine RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_1
