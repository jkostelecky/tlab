#include "dns_const.h"

program STATE
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS
    use TLAB_PROCS
    use THERMO_VARS
    use THERMO_THERMAL
    use THERMO_ANELASTIC
    use THERMO_CALORIC

    implicit none

    real(wp) p(1), ps(1), t(1), qs(1), qv(1), qt(1), ql(1), r(1), e(1), h(1), z1(2), dummy(1), dqldqt(1), ep(1), theta(1), theta_e(1), Td(1)
    real(wp) heat1(1), heat2(1), cp1(1), cp2(1), alpha(1), as(1), bs(1)
    real(wp) r1(1), h1(1), s(3)
    integer(wi) iopt

! ###################################################################
    call TLAB_START()

    imixture = MIXT_TYPE_AIRWATER
    nondimensional = .false.
    call THERMO_INITIALIZE()
    ep = 1.0_wp
    dsmooth = 1.0_wp
    scaleheight = 1.0_wp

    write (*, *) 'Case p-t (1) or d-e (2) or p-h (3)?'
    read (*, *) iopt

    if (iopt == 1) then
        write (*, *) 'temperature (K) ?'
        read (*, *) t
        t = t/TREF !t = (t+273.15)/TREF
        write (*, *) 'pressure (bar) ?'
        read (*, *) p

    else if (iopt == 2) then
        write (*, *) 'density ?'
        read (*, *) r
        write (*, *) 'energy ?'
        read (*, *) e

    else if (iopt == 3) then
        write (*, *) 'enthalpy (kJ/kg)?'
        read (*, *) h
        write (*, *) 'pressure (bar) ?'
        read (*, *) p

    end if

    write (*, *) 'water specific humidity (g/kg) ?'
    read (*, *) qt
    qt = qt*0.001_wp

! ###################################################################
    if (iopt == 1) then
        call THERMO_POLYNOMIAL_PSAT(1, t, ps)
        qs = 1.0_wp/(p/ps - 1.0_wp)*rd_ov_rv
        qs = qs/(1.0_wp + qs)
        if (qt(1) > qs(1)) then
            qv = qs*(1 - qt)
            ql = qt - qv
        else
            qv = qt
            ql = 1.0_wp
        end if
        z1(1) = qt(1)
        z1(2) = ql(1)
        call THERMO_CALORIC_ENTHALPY(1, z1, t, h)
        call THERMO_CALORIC_ENERGY(1, z1, t, e)
        call THERMO_THERMAL_DENSITY(1, z1, p, t, r)

        s(1) = h(1); s(2:3) = z1(1:2)
        call THERMO_ANELASTIC_THETA_L(1, 1, 1, s, ep, p, theta)
        call THERMO_ANELASTIC_THETA_E(1, 1, 1, s, ep, p, theta_e)
        call THERMO_ANELASTIC_DEWPOINT(1, 1, 1, s, ep, p, r, Td, dummy)

    else if (iopt == 2) then
        z1(1) = qt(1)
        call THERMO_CALORIC_TEMPERATURE(1, z1, e, r, T, dqldqt)
        ql = z1(2)
        qv = qt - ql
        qs = qv ! initial condition for next routine
        call THERMO_THERMAL_PRESSURE(1, z1, r, t, p)
        call THERMO_POLYNOMIAL_PSAT(1, t, ps)
        qs = 1.0_wp/(p/ps - 1.0_wp)*rd_ov_rv
        qs = qs/(1.0_wp + qs)
        call THERMO_CALORIC_ENTHALPY(1, z1, t, h)

        s(1) = h(1); s(2:3) = z1(1:2)
        call THERMO_ANELASTIC_THETA_L(1, 1, 1, s, ep, p, theta)
        call THERMO_ANELASTIC_THETA_E(1, 1, 1, s, ep, p, theta_e)
        call THERMO_ANELASTIC_DEWPOINT(1, 1, 1, s, ep, p, r, Td, dummy)

    else if (iopt == 3) then
        h = h/TREF/1.007
        z1(1) = qt(1)
        call THERMO_ANELASTIC_PH(1, 1, 1, z1, h, ep, p)
        s(1) = h(1); s(2:3) = z1(1:2)
        call THERMO_ANELASTIC_TEMPERATURE(1, 1, 1, s, ep, T)
        ! CALL THERMO_AIRWATER_PH_RE(1, z1, p, h, T)
        ql(1) = z1(2)
        qv = qt - ql

        call THERMO_POLYNOMIAL_PSAT(1, T, ps)
        qs = 1.0_wp/(p/ps - 1.0_wp)*rd_ov_rv
        qs = qs/(1.0_wp + qs)
        call THERMO_THERMAL_DENSITY(1, z1, p, T, r)
        call THERMO_CALORIC_ENERGY(1, z1, T, e)
        call THERMO_ANELASTIC_THETA_L(1, 1, 1, s, ep, p, theta)
        call THERMO_ANELASTIC_THETA_E(1, 1, 1, s, ep, p, theta_e)
        call THERMO_ANELASTIC_DEWPOINT(1, 1, 1, s, ep, p, r, Td, dummy)

! check
        call THERMO_ANELASTIC_DENSITY(1, 1, 1, s, ep, p, r1)
!     r2 = p/(T*(1- qt +qv/rd_ov_rv ) )
        call THERMO_CALORIC_ENTHALPY(1, z1, T, h1)

    end if

    write (*, 1000) 'Saturation specific humidity ......:', qs
    write (*, 1000) 'Vapor specific humidity ...........:', qv
    write (*, 1000) 'Liquid specific humidity ..........:', ql
    write (*, 1000) 'Density ...........................:', r
    write (*, 1000) 'Pressure (bar) ....................:', p
    write (*, 1000) 'Saturation pressure (bar) .........:', ps
    write (*, 1000) 'Temperature (K) ...................:', t*TREF !- 273.15
    write (*, 1000) 'Dewpoint temperature (K) ..........:', Td*TREF
    write (*, 1000) 'Specific heat capacity ............:', Cd + qt*Cdv + ql*Cvl
    write (*, 1000) 'Specific energy ...................:', e
    write (*, 1000) 'Specific enthalpy .................:', h
    write (*, 1000) 'Reference latent heat (kJ/kg) .....:', -THERMO_AI(6, 1, 3)*1.007*TREF
    WRITE(*,1000) 'Latent heat (kJ/kg) ...............:', (-Cl-t*Lvl ) *1.007 *TREF
    write (*, 1000) 'Liquid-water potential T (K) ......:', theta*TREF
    write (*, 1000) 'Equivalent potential T (K) ........:', theta_e*TREF
    if (iopt == 3) then
        write (*, 1000) 'Density ...........................:', r1
        write (*, 1000) 'Specific enthalpy .................:', h1
    end if

! ###################################################################
    write (*, *) ' '
    write (*, *) 'Calculate reversal linear coefficients (1-yes/0-no) ?'
    read (*, *) iopt

    if (iopt == 1 .and. ql(1) > 1.0_wp) then
        heat1 = -Lvl -Cvl*t
        heat2 = heat1*(1.0_wp + qv/(1.0_wp - qt)) - Cdv*t

        cp1 = (1.0_wp - qt)*Cd + qv*THERMO_AI(1, 1, 1) + ql*Cl
        dummy = (heat1**2)*qv/((t**2)*cp1*GRATIO*Rv)
        cp2 = cp1*(1.0_wp + dummy*(1.0_wp + qv/(1.0_wp - qt)/rd_ov_rv))

        alpha = 1.0_wp + heat1*qv/((1.0_wp - qt)*GRATIO*Rd*t)

        as = -alpha/cp2/t
        bs = heat2*as + 1.0_wp/(1.0_wp - qt)
        write (*, *) 'Enthalpy coefficient ..........:', as
        write (*, *) 'Water fraction coefficient ....:', bs

    else if (iopt == 1 .and. ql(1) == 1.0_wp) then
        cp1 = Cd + qt*Cdv

        as = -1.0_wp/cp1/t
        bs = Cdv/cp1 - Rdv/(Rd + qt*Rdv)
        write (*, *) 'Enthalpy coefficient ..........:', as
        write (*, *) 'Water fraction coefficient ....:', bs
    end if

    stop

1000 format(A, G_FORMAT_R)

end program STATE
