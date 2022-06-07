#include "types.h"
#include "dns_error.h"
#include "dns_const_mpi.h"

#define SIZEOFBYTE 1

! The offset should be converted into logical and,
! when .TRUE., read the first integer in file as offset

SUBROUTINE IO_WRITE_SUBARRAY4(iflag_mode, fname, varname, DATA, sizes, work)

  USE TLAB_TYPES,     ONLY : subarray_dt
  USE TLAB_CONSTANTS, ONLY : lfile
  USE TLAB_VARS,    ONLY : io_aux
  USE TLAB_PROCS
#ifdef USE_MPI
  USE MPI
#endif
  IMPLICIT NONE

  TINTEGER,                                   INTENT(IN)    :: iflag_mode
  CHARACTER*(*),                              INTENT(IN)    :: fname
  TINTEGER,                                   INTENT(IN)    :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
  CHARACTER*32, DIMENSION(sizes(5)),          INTENT(IN)    :: varname
  TREAL,        DIMENSION(sizes(1),sizes(5)), INTENT(IN)    :: DATA
  REAL(4),      DIMENSION(sizes(1)),          INTENT(INOUT) :: work

  ! -----------------------------------------------------------------------
  TINTEGER iv, isize
  CHARACTER*64 name

#ifdef USE_MPI
  TINTEGER :: mpio_status(MPI_STATUS_SIZE), mpio_fh, ims_err
#else
  TINTEGER :: ioffset_local
#endif

  ! #######################################################################
#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

  isize = ( sizes(3) -sizes(2) ) /sizes(4) +1

#ifdef USE_MPI
  IF ( io_aux(iflag_mode)%active ) THEN
#endif

    DO iv = 1,sizes(5)
      name = TRIM(ADJUSTL(fname))
      IF ( varname(iv) .NE. '' ) name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(varname(iv)))

      CALL TLAB_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(name))//'...')

      work(1:isize) = SNGL(DATA(sizes(2):sizes(3):sizes(4),iv))

#ifdef USE_MPI
      CALL MPI_File_open(io_aux(iflag_mode)%communicator, TRIM(ADJUSTL(name)), IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,mpio_fh, ims_err)
      DO_MPI_ERROR_CHECK
      CALL MPI_File_set_view(mpio_fh, io_aux(iflag_mode)%offset, MPI_REAL4, io_aux(iflag_mode)%subarray, 'native', MPI_INFO_NULL, ims_err)
      DO_MPI_ERROR_CHECK
      CALL MPI_File_write_all(mpio_fh, work, isize, MPI_REAL4, mpio_status, ims_err)
      DO_MPI_ERROR_CHECK
      CALL MPI_File_close(mpio_fh, ims_err)
      DO_MPI_ERROR_CHECK

#else
#include "dns_open_file.h"
      ioffset_local = io_aux(iflag_mode)%offset + 1
      WRITE(LOC_UNIT_ID,POS=ioffset_local) work(1:isize)
      CLOSE(LOC_UNIT_ID)

#endif

    ENDDO

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE IO_WRITE_SUBARRAY4

!########################################################################
!########################################################################
SUBROUTINE IO_READ_SUBARRAY8(iflag_mode, fname, varname, DATA, sizes, work)

  USE TLAB_TYPES,     ONLY : subarray_dt
  USE TLAB_CONSTANTS, ONLY : lfile, efile
#ifdef TRACE_ON
  USE TLAB_CONSTANTS, ONLY : tfile
#endif
  USE TLAB_VARS,      ONLY : io_aux
  USE TLAB_PROCS
#ifdef USE_MPI
  USE MPI
  USE TLAB_MPI_VARS,  ONLY : ims_err
#endif

  IMPLICIT NONE

  TINTEGER,                                   INTENT(IN)    :: iflag_mode
  CHARACTER*(*),                              INTENT(IN)    :: fname
  TINTEGER,                                   INTENT(IN)    :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
  CHARACTER*32, DIMENSION(sizes(5)),          INTENT(IN)    :: varname
  TREAL,        DIMENSION(sizes(1),sizes(5)), INTENT(OUT)   :: DATA
  TREAL,        DIMENSION(sizes(1)),          INTENT(INOUT) :: work

  ! -----------------------------------------------------------------------
  TINTEGER iv, isize
  CHARACTER*64 name

#ifdef USE_MPI
  INTEGER :: mpio_status(MPI_STATUS_SIZE), mpio_fh
#else
  TINTEGER :: ioffset_local
#endif

  ! #######################################################################
#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile,'ENTERING IO_READ_SUBARRAY8')
#endif
  isize = ( sizes(3) -sizes(2) ) /sizes(4) +1

  DO iv = 1,sizes(5)
     name = TRIM(ADJUSTL(fname))
     IF ( varname(iv) .NE. '' ) name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(varname(iv)))
     CALL TLAB_WRITE_ASCII(lfile, 'Reading field '//TRIM(ADJUSTL(name))//'...')
  ENDDO

#ifdef USE_MPI
  IF ( io_aux(iflag_mode)%active ) THEN
#endif

    DO iv = 1,sizes(5)
      name = TRIM(ADJUSTL(fname))
      IF ( varname(iv) .NE. '' ) name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(varname(iv)))
#ifdef USE_MPI
      CALL MPI_File_open(io_aux(iflag_mode)%communicator, TRIM(ADJUSTL(name)), MPI_MODE_RDONLY,MPI_INFO_NULL,mpio_fh, ims_err)
      DO_MPI_ERROR_CHECK
      CALL MPI_File_set_view(mpio_fh, io_aux(iflag_mode)%offset, MPI_REAL8, io_aux(iflag_mode)%subarray, 'native', MPI_INFO_NULL, ims_err)
      DO_MPI_ERROR_CHECK
      CALL MPI_File_read_all(mpio_fh, work, isize, MPI_REAL8, mpio_status, ims_err)
      DO_MPI_ERROR_CHECK
      CALL MPI_File_close(mpio_fh, ims_err)
      DO_MPI_ERROR_CHECK
#else
#include "dns_open_file.h"
      ioffset_local = io_aux(iflag_mode)%offset + 1
      READ(LOC_UNIT_ID,POS=ioffset_local) work(1:isize)
      CLOSE(LOC_UNIT_ID)
#endif
      DATA(sizes(2):sizes(3):sizes(4),iv) = work(1:isize)
   ENDDO

#ifdef USE_MPI
  ENDIF
#endif
#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile,'LEAVING IO_READ_SUBARRAY8')
#endif

  RETURN
END SUBROUTINE IO_READ_SUBARRAY8

!########################################################################
!########################################################################
SUBROUTINE IO_WRITE_SUBARRAY8(iflag_mode, fname, varname, DATA, sizes, work)

  USE TLAB_TYPES,     ONLY : subarray_dt
  USE TLAB_CONSTANTS, ONLY : lfile
  USE TLAB_VARS,    ONLY : io_aux
  USE TLAB_PROCS
#ifdef USE_MPI
  USE MPI
#endif

  IMPLICIT NONE

  TINTEGER,                                   INTENT(IN)    :: iflag_mode
  CHARACTER*(*),                              INTENT(IN)    :: fname
  TINTEGER,                                   INTENT(IN)    :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
  CHARACTER*32, DIMENSION(sizes(5)),          INTENT(IN)    :: varname
  TREAL,        DIMENSION(sizes(1),sizes(5)), INTENT(IN)    :: DATA
  TREAL,        DIMENSION(sizes(1)),          INTENT(INOUT) :: work

  ! -----------------------------------------------------------------------
  TINTEGER iv, isize
  CHARACTER*64 name

#ifdef USE_MPI
  TINTEGER :: mpio_status(MPI_STATUS_SIZE), mpio_fh, ims_err
#else
  TINTEGER :: ioffset_local
#endif

  ! #######################################################################
#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

  isize = ( sizes(3) -sizes(2) ) /sizes(4) +1

#ifdef USE_MPI
  IF ( io_aux(iflag_mode)%active ) THEN
#endif

    DO iv = 1,sizes(5)
      name = TRIM(ADJUSTL(fname))
      IF ( varname(iv) .NE. '' ) name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(varname(iv)))

      CALL TLAB_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(name))//'...')

      work(1:isize) = DATA(sizes(2):sizes(3):sizes(4),iv)

#ifdef USE_MPI
      CALL MPI_File_open(io_aux(iflag_mode)%communicator, TRIM(ADJUSTL(name)), IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,mpio_fh, ims_err)
      DO_MPI_ERROR_CHECK
      CALL MPI_File_set_view(mpio_fh, io_aux(iflag_mode)%offset, MPI_REAL8, io_aux(iflag_mode)%subarray, 'native', MPI_INFO_NULL, ims_err)
      DO_MPI_ERROR_CHECK
      CALL MPI_File_write_all(mpio_fh, work, isize, MPI_REAL8, mpio_status, ims_err)
      DO_MPI_ERROR_CHECK
      CALL MPI_File_close(mpio_fh, ims_err)
      DO_MPI_ERROR_CHECK

#else
#include "dns_open_file.h"
      ioffset_local = io_aux(iflag_mode)%offset + 1
      WRITE(LOC_UNIT_ID,POS=ioffset_local) work(1:isize)
      CLOSE(LOC_UNIT_ID)

#endif

    ENDDO

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE IO_WRITE_SUBARRAY8
