#include "dns_const.h"
#include "dns_error.h"

module Chemistry
    use TLAB_CONSTANTS, only: wp, wi, pi_wp, efile, MAX_PARS
    use TLAB_TYPES, only: term_dt
    use TLAB_VARS, only: inb_scal
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
    implicit none
    private

    type(term_dt) :: chemistryProps              ! Chemistry parameters

    public :: chemistryProps
    public :: Chemistry_Initialize
    public :: Chemistry_Source

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_QUADRATIC = 1
    integer, parameter :: TYPE_QUADRATIC3 = 2
    integer, parameter :: TYPE_LAYEREDRELAXATION = 3
    integer, parameter :: TYPE_OZONE = 4

    real(wp), allocatable :: relaxation_strength(:)

contains
    !########################################################################
    !########################################################################
    subroutine Chemistry_Initialize(inifile)
        use TLAB_TYPES, only: profiles_dt
        use TLAB_VARS, only: damkohler, sbg, g
        use PROFILES
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile
        character(len=512) sRes
        integer(wi) idummy, is, j
        type(profiles_dt) prof_loc

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        call TLAB_WRITE_ASCII(bakfile, '#')
        call TLAB_WRITE_ASCII(bakfile, '#[Chemistry]')
        call TLAB_WRITE_ASCII(bakfile, '#Type=<none/quadratic/layeredrelaxation/ozone>')
        call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

        call SCANINICHAR(bakfile, inifile, 'Chemistry', 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; chemistryProps%type = TYPE_NONE
        elseif (trim(adjustl(sRes)) == 'quadratic') then; chemistryProps%type = TYPE_QUADRATIC; 
        elseif (trim(adjustl(sRes)) == 'quadratic3') then; chemistryProps%type = TYPE_QUADRATIC3; 
        elseif (trim(adjustl(sRes)) == 'layeredrelaxation') then; chemistryProps%type = TYPE_LAYEREDRELAXATION; 
        elseif (trim(adjustl(sRes)) == 'ozone') then; chemistryProps%type = TYPE_OZONE; 
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Error in Chemistry.Type.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

        if (chemistryProps%type /= EQNS_NONE) then
            chemistryProps%parameters(:) = 0.0_wp
            call SCANINICHAR(bakfile, inifile, 'Chemistry', 'Parameters', '1.0', sRes)
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, chemistryProps%parameters)

        end if

        chemistryProps%active = .false.
        do is = 1, inb_scal
            if (abs(damkohler(is)) > 0.0_wp) chemistryProps%active(is) = .true.
        end do

        select case (chemistryProps%type)
        case (TYPE_LAYEREDRELAXATION)
            prof_loc%type = PROFILE_TANH
            prof_loc%ymean = sbg(is)%ymean
            prof_loc%thick = -chemistryProps%parameters(3)*0.5_wp
            prof_loc%mean = 0.5_wp
            prof_loc%delta = 1.0_wp
            prof_loc%lslope = 0.0_wp
            prof_loc%uslope = 0.0_wp

            allocate (relaxation_strength(g(2)%size))
            do j = 1, g(2)%size
                relaxation_strength(j) = PROFILES_CALCULATE(prof_loc, g(2)%nodes(j) - chemistryProps%parameters(2))
            end do

        end select

        return
    end subroutine Chemistry_Initialize

    !########################################################################
    !########################################################################
    subroutine Chemistry_Source(locProps, nx, ny, nz, is, s, source)
        use TLAB_VARS, only: damkohler

        type(term_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, is
        real(wp), intent(in) :: s(nx, ny, nz, inb_scal)
        real(wp), intent(out) :: source(nx, ny, nz)

        ! -----------------------------------------------------------------------
        integer(wi) j
        real(wp) dummy, dummy2

        !########################################################################
        select case (locProps%type)

        case (TYPE_LAYEREDRELAXATION)
            do j = 1, ny
                source(:, j, :) = -damkohler(is)/locProps%parameters(1)*relaxation_strength(j)*s(:, j, :, is)
            end do

        case (TYPE_QUADRATIC)
            dummy = damkohler(is)*locProps%parameters(is)
            source = dummy*s(:, :, :, 2)*s(:, :, :, 3)

        case (TYPE_QUADRATIC3)
            dummy = damkohler(is)*locProps%parameters(is)

            if (is >= 1 .and. is <= 3) then
                source = dummy*s(:, :, :, 2)*s(:, :, :, 3)
            else if (is >= 4 .and. is <= 6) then
                source = dummy*s(:, :, :, 4)*s(:, :, :, 5)
            else if (is >= 7 .and. is <= 9) then
                source = dummy*s(:, :, :, 7)*s(:, :, :, 8)
            end if

        case (TYPE_OZONE)
            dummy = damkohler(is)
            if (is == 4) dummy = -dummy

            source = -locProps%parameters(1)/(1.0_wp + locProps%parameters(2)*s(:, :, :, 1))
            source = exp(source)

            if (is == 4) then
                dummy2 = 1.0_wp + locProps%parameters(3)
                source = dummy*(dummy2*s(:, :, :, 4) - source*s(:, :, :, 2)*s(:, :, :, 3))
            else
                source = dummy*(s(:, :, :, 4) - source*s(:, :, :, 2)*s(:, :, :, 3))
            end if

        end select

        return
    end subroutine Chemistry_Source

end module Chemistry
