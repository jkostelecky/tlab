!# Grid generation tool. Origin is set always to (0,0,0)
#include "dns_error.h"

#define C_FILE_LOC "INIGRID"

program INIGRID
    use TLAB_TYPES, only: grid_dt, wp
    use TLAB_CONSTANTS, only: gfile, ifile, lfile, efile
    use TLAB_PROCS
    use GRID_LOCAL
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_pro
#endif
    implicit none

    character*32 sfile, bakfile
    type(grid_dt) :: g(3)
    real(wp), allocatable :: wrk1d(:, :)
    integer(wi) idir, iseg, isize_wrk1d, n, nmax, iloc
    real(wp) scale_old, scale_new, ds
    character(len=16), parameter :: block(3) = ['IniGridOx', 'IniGridOy', 'IniGridOz']

    ! #######################################################################
    ! Initialize
    ! #######################################################################
    sfile = TRIM(ADJUSTL(gfile))//'.sts'
    bakfile = TRIM(ADJUSTL(ifile))//'.bak'

    g(1)%name = 'x'
    g(2)%name = 'y'
    g(3)%name = 'z'

    call TLAB_START()

    do idir = 1, 3
        call GRID_READBLOCK(bakfile, ifile, block(idir), g_build(idir), g(idir)%periodic)

        g(idir)%size = g_build(idir)%size(1)                    ! Calculate total number of points
        do iseg = 2, g_build(idir)%nseg
            g(idir)%size = g(idir)%size + g_build(idir)%SIZE(iseg) - 1
        end do
        if (g_build(idir)%mirrored) g(idir)%size = 2*g(idir)%size - 2

        allocate (g(idir)%nodes(g(idir)%size))                   ! memory space for the grid nodes

    end do

    isize_wrk1d = MAX(g(1)%size, MAX(g(2)%size, g(3)%size))
    allocate (wrk1d(isize_wrk1d, 8))

    ! #######################################################################
    ! Construct grid
    ! #######################################################################
    do idir = 1, 3

        iloc = 1
        if (g_build(idir)%mirrored) iloc = g(idir)%size/2   ! mirrored case; first point in array is imax/2

        g(idir)%nodes(iloc) = 0.0_wp                        ! set origin at zero

        do iseg = 1, g_build(idir)%nseg                     ! Loop over the segments that build the grid
            nmax = g_build(idir)%size(iseg)                 ! for readability below
            ! create uniform reference grid s starting at zero
            if (nmax > 1) then
                ds = (g_build(idir)%end(iseg) - g(idir)%nodes(iloc))/real(nmax - 1, wp)
                g(idir)%nodes(iloc:) = [(real(n - 1, wp), n=1, nmax)]*ds + g(idir)%nodes(iloc)
                ! do n = 1, nmax - 1
                !     g(idir)%nodes(iloc + n) = g(idir)%nodes(iloc + n - 1) + ds
                ! end do

                select case (g_build(idir)%opts(1, iseg))

                case (GTYPE_UNIFORM)
                    ! already done

                case (GTYPE_TANH)
                    call BLD_TANH(idir, iseg, g(idir)%nodes(iloc:), nmax, wrk1d)

                case (GTYPE_EXP)
                    call BLD_EXP(idir, iseg, g(idir)%nodes(iloc:), nmax, wrk1d)

                case DEFAULT
                    call BLD_THEREST(idir, iseg, g(idir)%nodes(iloc:), nmax)

                end select

                iloc = iloc + nmax - 1

            end if
        end do

        if (g_build(idir)%mirrored) call GRID_MIRROR(g(idir)%size, g(idir)%nodes)

        if (g(idir)%size > 1) then
            g(idir)%scale = g(idir)%nodes(g(idir)%size) - g(idir)%nodes(1)
        else
            g(idir)%scale = 1.0_wp
        end if

        if (g_build(idir)%fixed_scale > 0.0_wp) then                ! rescale on exact fixed value
            scale_new = g_build(idir)%fixed_scale
            scale_old = g(idir)%scale
            g(idir)%nodes = (g(idir)%nodes/scale_old)*scale_new     ! rescale nodes
            g(idir)%nodes(g(idir)%size) = scale_new                 ! avoid rounding error
            g(idir)%scale = g(idir)%nodes(g(idir)%size)             ! update scale
        end if

        if (g(idir)%periodic) g(idir)%size = g(idir)%size - 1

    end do

    ! #######################################################################
    ! Statistics
    ! #######################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
        open (20, file=sfile)

        do idir = 1, 3
            write (20, 3000) '['//TRIM(ADJUSTL(g(idir)%name))//'-direction]'

            if (g(idir)%size > 1) then
                wrk1d(2, 1) = g(idir)%nodes(2) - g(idir)%nodes(1)
                do n = 3, g(idir)%size
                    wrk1d(n, 1) = g(idir)%nodes(n) - g(idir)%nodes(n - 1)
                    wrk1d(n, 2) = wrk1d(n, 1)/wrk1d(n - 1, 1)
                end do

                write (20, 2000) 'number of points .......: ', g(idir)%size
                write (20, 1000) 'origin .................: ', g(idir)%nodes(1)
                write (20, 1000) 'end point ..............: ', g(idir)%nodes(g(idir)%size)
                write (20, 1000) 'scale ..................: ', g(idir)%scale
                write (20, 1000) 'minimum step ...........: ', MINVAL(wrk1d(2:g(idir)%size, 1))
                write (20, 1000) 'maximum step ...........: ', MAXVAL(wrk1d(2:g(idir)%size, 1))
                write (20, 1000) 'minimum stretching .....: ', MINVAL(wrk1d(3:g(idir)%size, 2))
                write (20, 1000) 'maximum stretching .....: ', MAXVAL(wrk1d(3:g(idir)%size, 2))

            else
                write (20, '(a7)') '2D case'

            end if

        end do

        close (20)

        ! #######################################################################
        ! Writing data
        ! #######################################################################
        call TLAB_WRITE_ASCII(lfile, 'Writing grid.')
  call IO_WRITE_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, g(1)%nodes, g(2)%nodes, g(3)%nodes)

#ifdef USE_MPI
    end if
#endif

    call TLAB_STOP(0)

1000 format(a25, e12.5)
2000 format(a25, i5)
3000 format(a13)

contains
    ! #######################################################################
    ! #######################################################################
    subroutine GRID_READBLOCK(bakfile, inifile, block, var, periodic)

        character(LEN=*), intent(in) :: bakfile, inifile, block
        type(grid_build_dt), intent(out) :: var
        logical, intent(OUT) :: periodic

! -------------------------------------------------------------------
        integer(wi) idummy
        character(LEN=512) sRes, str

! #######################################################################
        call TLAB_WRITE_ASCII(bakfile, '['//block//']')
        call TLAB_WRITE_ASCII(bakfile, 'segments=<number of segments>')
        call TLAB_WRITE_ASCII(bakfile, 'periodic=<yes/no>')
        call TLAB_WRITE_ASCII(bakfile, 'mirrored=<yes/no>')
        call TLAB_WRITE_ASCII(bakfile, 'fixed_scale=<value>')

        call SCANINIINT(bakfile, inifile, block, 'segments', '1', var%nseg)

        periodic = .false.
        call SCANINICHAR(bakfile, inifile, block, 'periodic', 'no', sRes)
        if (TRIM(ADJUSTL(sRes)) == 'yes') periodic = .true.

        var%mirrored = .false.
        call SCANINICHAR(bakfile, inifile, block, 'mirrored', 'no', sRes)
        if (TRIM(ADJUSTL(sRes)) == 'yes') var%mirrored = .true.

        call SCANINIREAL(bakfile, inifile, block, 'fixed_scale', '-1.0', var%fixed_scale)

        if (periodic .and. var%mirrored) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Periodicity with mirroring is not supported.')
            call TLAB_STOP(DNS_ERROR_GRID_SCALE)
        end if

! -------------------------------------------------------------------
! Loop over the segments
! -------------------------------------------------------------------
        do iseg = 1, var%nseg
            write (str, *) iseg

            call TLAB_WRITE_ASCII(bakfile, 'Segment number '//TRIM(ADJUSTL(str)))
            call TLAB_WRITE_ASCII(bakfile, 'scales_'//TRIM(ADJUSTL(str))//'=<physical end of the segment>')
            call TLAB_WRITE_ASCII(bakfile, 'points_'//TRIM(ADJUSTL(str))//'=<points in the segment>')
            call TLAB_WRITE_ASCII(bakfile, 'opts_'//TRIM(ADJUSTL(str))//'=<option>')
            call TLAB_WRITE_ASCII(bakfile, 'vals_'//TRIM(ADJUSTL(str))//'=<values>')

            call SCANINIINT(bakfile, inifile, block, 'points_'//TRIM(ADJUSTL(str)), '1', var%size(iseg))
            call SCANINIREAL(bakfile, inifile, block, 'scales_'//TRIM(ADJUSTL(str)), '-1.0', var%end(iseg))

            var%opts(:, iseg) = 0
            call SCANINICHAR(bakfile, inifile, block, 'opts_'//TRIM(ADJUSTL(str)), '1', sRes)
            if (TRIM(ADJUSTL(sRes)) == 'uniform') then; var%opts(1, iseg) = GTYPE_UNIFORM
            else if (TRIM(ADJUSTL(sRes)) == 'tanh') then; var%opts(1, iseg) = GTYPE_TANH
            else if (TRIM(ADJUSTL(sRes)) == 'exp') then; var%opts(1, iseg) = GTYPE_EXP
            else
                idummy = MAX_PARAMES
                call LIST_INTEGER(sRes, idummy, var%opts(1, iseg))
            end if

            var%vals(:, iseg) = 0
            call SCANINICHAR(bakfile, inifile, block, 'vals_'//TRIM(ADJUSTL(str)), '1.0', sRes)
            idummy = MAX_PARAMES
            call LIST_REAL(sRes, idummy, var%vals(1, iseg))

        end do

        return
    end subroutine GRID_READBLOCK

    ! #######################################################################
    ! #######################################################################
    subroutine GRID_MIRROR(imax, x)
        implicit none

        integer(wi), intent(IN) :: imax
        real(wp), intent(INOUT) :: x(imax)

        ! -----------------------------------------------------------------------
        integer(wi) i
        real(wp) offset

        ! #######################################################################
        ! Offset for even number of points
        offset = (x(imax/2 + 1) - x(imax/2))/2.0_wp
        do i = imax/2, imax
            x(i) = x(i) - offset
        end do

        ! Mirroring
        do i = 1, imax/2 - 1
            x(i) = -x(imax + 1 - i)
        end do

        ! Global translation to set origin at zero
        offset = x(1)
        x = x - offset

        return
    end subroutine GRID_MIRROR

end program INIGRID
