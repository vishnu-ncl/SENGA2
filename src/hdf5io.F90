module hdf5io

! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:25:50

#ifdef hdf5
use hdf5
INCLUDE 'mpif.h'
CHARACTER (LEN=60) :: h5_filename(2)
#endif
contains

!     HDF5 DMPI FILES ORGANISATION

!     THE  DMPI FILES CONSISTS OF 7 DATASETS:
!     - '/CONSTANTS' (2 DIM ARRAY)  WHICH CONTAINS THE VALUE OF
!       NXNODE, NYNODE, NZNODE, NSPEC FOR EACH MPI RANK
!     -'/DRUN': (4 DIM ARRAY), THE FIRST DIMENSION IS THE MPI RANK,
!        THE THREE OTHERS ARE NXNODE,NYNODE,NZNODE.
!     -'/URUN': IDEM
!     -'/VRUN': IDEM
!     -'/WRUN': IDEM
!     -'/ERUN': IDEM
!     -'/YRUN': (5 DIM ARRAY), THE FIRST DIMENSION IS THE MPI RANK,
!        THE FOUR OTHERS ARE NXNODE,NYNODE,NZNODE,NSPEC.

!     THE SUBROUTINES TO MANAGE THE DMPI FILES ARE:
!       - CREATE_H5DUMP_FILES
!       - READ_H5DUMP_FILES
!       - WRITE_H5_DUMPFILE


SUBROUTINE hdf5_init()
!     INITIALISE THE HDF5 LIBRARY.
!     IT IS MANDATORY TO CALL IT BEFORE CALLING ANY HDF5 SUBROUTINES.
#ifdef hdf5
integer(kind=4) :: ERR
CALL h5open_f(ERR)
IF(ERR/=0)THEN
  WRITE(*,*)"ERROR WHEN INITIALIZING HDF5 LIBRARY"
  CALL mpi_abort(mpi_comm_world,-1,ERR)
END IF
#endif
END SUBROUTINE

SUBROUTINE hdf5_close()
!       CLOSE THE HDF5 LIBRARY, TO CALL BEFORE THE END OF THE PROGRAM
!       NO CALL TO HDF5 SUBROUTINES CAN BE MADE AFTER.
#ifdef hdf5
integer(kind=4) :: ERR
CALL h5close_f(ERR)
IF(ERR/=0)THEN
  WRITE(*,*)"ERROR WHEN CLOSING HDF5 LIBRARY"
  CALL mpi_abort(mpi_comm_world,-1,ERR)
END IF
#endif
END SUBROUTINE

SUBROUTINE create_h5dump_files()
#ifdef hdf5
use com_senga
integer(kind=4) :: ERR, i
LOGICAL :: exists
INTEGER(:: hid_t) space_id, dset_id,plist_id
INTEGER(:: hsize_t) dims(5)
INTEGER(:: hid_t) dumpfile_id(2)
CHARACTER (LEN=6) :: pnproc
CHARACTER (LEN=5) :: dset_name

WRITE(pnproc,'(I6.6)')iproc
CALL h5pcreate_f(h5p_file_access_f, plist_id, ERR)
CALL h5pset_fapl_mpio_f(plist_id, mpi_comm_world, mpi_info_null, ERR)
DO i=1,2
  INQUIRE(FILE=h5_filename(i),EXIST=exists)
  IF(.NOT.exists)THEN
    CALL h5fcreate_f(h5_filename(i), h5f_acc_trunc_f,  &
        dumpfile_id(i),ERR, access_prp=plist_id)
    IF(ERR/=0)THEN
      WRITE(*,*)"ERROR WHEN CALLING H5FCREATE_F INSIDE CREATE_H5DUMP_FILES"
      CALL mpi_abort(mpi_comm_world,-1,ERR)
    END IF

    dims(1)=nproc
    dims(2)=4
    CALL h5screate_simple_f(2, dims, space_id, ERR)
    CALL h5dcreate_f(dumpfile_id(i),"/CONSTANTS",  &
        h5t_native_integer,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    CALL h5sclose_f(space_id,ERR)

!IN ORDER TO SIMPLIFY THE HDF5 CODE, EACH RANK WILL USE
!AN ARRAY OF THE SAME SIZE, SO WE USE THE BIGGEST VALUE OF
!NXNODE, NYNODE AND NZNODE
    dims(1)=nproc
    dims(2)=nxsize
    dims(3)=nysize
    dims(4)=nzsize
    dims(5)=nspec
    CALL mpi_allreduce(mpi_in_place, dims, 5, mpi_integer,  &
        mpi_max, mpi_comm_world, ERR)

    CALL h5screate_simple_f(4, dims, space_id, ERR)
    dset_name = "/DRUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    dset_name = "/URUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    dset_name = "/VRUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    dset_name = "/WRUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    dset_name = "/ERUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    CALL h5sclose_f(space_id,ERR)

    CALL h5screate_simple_f(5,dims,space_id,ERR)
    dset_name = "/YRUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    CALL h5sclose_f(space_id,ERR)
    CALL h5fclose_f(dumpfile_id(i),ERR)
  END IF
END DO
CALL h5pclose_f(plist_id, ERR)
#endif
END SUBROUTINE

! UMOD DATA I/O START
SUBROUTINE create_h5data_files()
#ifdef hdf5
use com_senga
integer(kind=4) :: ERR, i
LOGICAL :: exists
INTEGER(:: hid_t) space_id, dset_id,plist_id
INTEGER(:: hsize_t) dims(5)
INTEGER(:: hid_t) dumpfile_id(2)
CHARACTER (LEN=6) :: pnproc
CHARACTER (LEN=5) :: dset_name

WRITE(pnproc,'(I6.6)')iproc
CALL h5pcreate_f(h5p_file_access_f, plist_id, ERR)
CALL h5pset_fapl_mpio_f(plist_id, mpi_comm_world, mpi_info_null, ERR)
DO i=1,2
  INQUIRE(FILE=h5_filename(i),EXIST=exists)
  IF(.NOT.exists)THEN
    CALL h5fcreate_f(h5_filename(i), h5f_acc_trunc_f,  &
        dumpfile_id(i),ERR, access_prp=plist_id)
    IF(ERR/=0)THEN
      WRITE(*,*)"ERROR WHEN CALLING H5FCREATE_F INSIDE CREATE_H5DATA_FILE"
      CALL mpi_abort(mpi_comm_world,-1,ERR)
    END IF

    dims(1)=nproc
    dims(2)=4
    CALL h5screate_simple_f(2, dims, space_id, ERR)
    CALL h5dcreate_f(dumpfile_id(i),"/CONSTANTS",  &
        h5t_native_integer,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    CALL h5sclose_f(space_id,ERR)

!IN ORDER TO SIMPLIFY THE HDF5 CODE, EACH RANK WILL USE
!AN ARRAY OF THE SAME SIZE, SO WE USE THE BIGGEST VALUE OF
!NXNODE, NYNODE AND NZNODE
    dims(1)=nproc
    dims(2)=nxsize
    dims(3)=nysize
    dims(4)=nzsize
    dims(5)=nspec
    CALL mpi_allreduce(mpi_in_place, dims, 5, mpi_integer,  &
        mpi_max, mpi_comm_world, ERR)

    CALL h5screate_simple_f(4, dims, space_id, ERR)
    dset_name = "/DRUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    dset_name = "/URUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    dset_name = "/VRUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    dset_name = "/WRUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    dset_name = "/ERUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    CALL h5sclose_f(space_id,ERR)

    CALL h5screate_simple_f(5,dims,space_id,ERR)
    dset_name = "/YRUN"
    CALL h5dcreate_f(dumpfile_id(i),dset_name,  &
        h5t_native_double,space_id,dset_id,ERR)
    CALL h5dclose_f(dset_id,ERR)
    CALL h5sclose_f(space_id,ERR)
    CALL h5fclose_f(dumpfile_id(i),ERR)
  END IF
END DO
CALL h5pclose_f(plist_id, ERR)
#endif
END SUBROUTINE
! UMOD DATA I/O END

SUBROUTINE read_h5dump_files()
#ifdef hdf5
use com_senga
integer(kind=4) :: ERR, consts(4)
INTEGER(:: hid_t) dset_id,file_id,plist_id,space_id,memspace_id
INTEGER(:: hsize_t) count(5), dims1(1), dims3(3), dims4(4)
INTEGER(:: hsize_t) offset(5)
CHARACTER (LEN=6) :: pnproc
CHARACTER (LEN=5) :: dset_name

CALL h5pcreate_f(h5p_file_access_f, plist_id, ERR)
CALL h5pset_fapl_mpio_f(plist_id, mpi_comm_world, mpi_info_null, ERR)

WRITE(pnproc,'(I6.6)')iproc
CALL h5fopen_f(h5_filename(idflag+1),h5f_acc_rdonly_f,file_id,  &
    ERR, access_prp=plist_id)
IF(ERR/=0)THEN
  WRITE(*,*)"COULDN'T OPEN FILE:",h5_filename(idflag+1)
  CALL mpi_abort(mpi_comm_world,-1,ERR)
END IF
CALL h5pclose_f(plist_id, ERR)

CALL h5pcreate_f(h5p_dataset_xfer_f, plist_id, ERR)
CALL h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, ERR)

offset(1)=iproc
offset(2)=0
count(1)=1
count(2)=4

CALL h5screate_simple_f(2, count, memspace_id, ERR)
dims1 = (/4/)
CALL h5dopen_f(file_id, "/CONSTANTS", dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dread_f(dset_id, h5t_native_integer, consts, dims1,ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5sclose_f(space_id, ERR)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(memspace_id, ERR)

IF(nxsize /= consts(1))WRITE(6,*)'DUMP INPUT SIZE ERROR: X'
IF(nysize /= consts(2))WRITE(6,*)'DUMP INPUT SIZE ERROR: Y'
IF(nzsize /= consts(3))WRITE(6,*)'DUMP INPUT SIZE ERROR: Z'
IF(nspec /= consts(4))WRITE(6,*) 'DUMP INPUT SIZE ERROR: SPECIES'

offset(1)=iproc
offset(2)=0
offset(3)=0
offset(4)=0
count(1)=1
count(2)=nxsize
count(3)=nysize
count(4)=nzsize

CALL h5screate_simple_f(4, count, memspace_id, ERR)

dims3 = (/nxsize,nysize,nzsize/)

dset_name = "/DRUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dread_f(dset_id, h5t_native_double, drun, dims3, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

dset_name = "/URUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dread_f(dset_id, h5t_native_double, urun, dims3, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

dset_name = "/VRUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dread_f(dset_id, h5t_native_double, vrun, dims3, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

dset_name = "/WRUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dread_f(dset_id, h5t_native_double, wrun, dims3, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

dset_name = "/ERUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dread_f(dset_id, h5t_native_double, erun, dims3, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

offset(5)=0
count(5)=nspec
CALL h5sclose_f(memspace_id, ERR)
CALL h5screate_simple_f(5, count, memspace_id, ERR)
dims4 = (/nxsize,nysize,nzsize,nspec/)

dset_name = "/YRUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dread_f(dset_id, h5t_native_double, yrun, dims4, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

CALL h5pclose_f(plist_id, ERR)
CALL h5sclose_f(memspace_id, ERR)
CALL h5fclose_f(file_id,ERR)
#endif
END SUBROUTINE

SUBROUTINE write_h5_dumpfile()
#ifdef hdf5
use com_senga
integer(kind=4) :: ERR, consts(4)
INTEGER(:: hid_t) dset_id,file_id,plist_id,space_id,memspace_id
INTEGER(:: hsize_t) count(5), dims1(1), dims3(3), dims4(4)
INTEGER(:: hsize_t) offset(5)
CHARACTER (LEN=6) :: pnproc
CHARACTER (LEN=5) :: dset_name

CALL h5pcreate_f(h5p_file_access_f, plist_id, ERR)
CALL h5pset_fapl_mpio_f(plist_id, mpi_comm_world, mpi_info_null, ERR)

WRITE(pnproc,'(I6.6)')iproc
CALL h5fopen_f(h5_filename(idflag+1),h5f_acc_rdwr_f,file_id,  &
    ERR,access_prp=plist_id)
IF(ERR/=0)THEN
  WRITE(*,*)"COULDN'T OPEN FILE:",h5_filename(idflag+1)
  CALL mpi_abort(mpi_comm_world,-1,ERR)
END IF
CALL h5pclose_f(plist_id, ERR)

CALL h5pcreate_f(h5p_dataset_xfer_f, plist_id, ERR)
CALL h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, ERR)

!WRITE NXNODES, NYNODES, NZNODES, NSPEC VALUES
dims1 = (/4/)
offset(1)=iproc
offset(2)=0
count(1)=1
count(2)=4
consts(1) = nxsize
consts(2) = nysize
consts(3) = nzsize
consts(4) = nspec

CALL h5screate_simple_f(2, count, memspace_id, ERR)
CALL h5dopen_f(file_id, "/CONSTANTS", dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dwrite_f(dset_id, h5t_native_integer, consts,dims1,ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

!WRITE DRUN, URUN, VRUN, WRUN, ERUN, YRUN VALUES
offset(1)=iproc
offset(2)=0
offset(3)=0
offset(4)=0
count(1)=1
count(2)=nxsize
count(3)=nysize
count(4)=nzsize

CALL h5screate_simple_f(4, count, memspace_id, ERR)

dims3 = (/nxsize,nysize,nzsize/)

write(*, '(a)') "Using the arrays not allocated by OPS, &
                        Please implement the function in OPS first, hdf5io.F90: ID=412"
            STOP

dset_name = "/DRUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dwrite_f(dset_id, h5t_native_double, drun, dims3, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

dset_name = "/URUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dwrite_f(dset_id, h5t_native_double, urun, dims3, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

dset_name = "/VRUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dwrite_f(dset_id, h5t_native_double, vrun, dims3, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

dset_name = "/WRUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dwrite_f(dset_id, h5t_native_double, wrun, dims3, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

dset_name = "/ERUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dwrite_f(dset_id, h5t_native_double, erun, dims3, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

offset(5)=0
count(5)=nspec
CALL h5sclose_f(memspace_id, ERR)
CALL h5screate_simple_f(5, count, memspace_id, ERR)
dims4 = (/nxsize,nysize,nzsize,nspec/)

dset_name = "/YRUN"
CALL h5dopen_f(file_id, dset_name, dset_id, ERR)
CALL h5dget_space_f(dset_id, space_id, ERR)
CALL h5sselect_hyperslab_f(space_id, h5s_select_set_f, offset, count, ERR)
CALL h5dwrite_f(dset_id, h5t_native_double, yrun, dims4, ERR,  &
    file_space_id=space_id, mem_space_id=memspace_id, xfer_prp = plist_id)
CALL check_err(ERR, __LINE__)
CALL h5dclose_f(dset_id, ERR)
CALL h5sclose_f(space_id, ERR)

CALL h5pclose_f(plist_id, ERR)
CALL h5sclose_f(memspace_id, ERR)
CALL h5fclose_f(file_id,ERR)
#endif
END SUBROUTINE

SUBROUTINE check_err(ERR, line)
integer(kind=4) :: ERR, line
#ifdef hdf5
IF (ERR /=0)THEN
  WRITE(*,*)"ERROR ON LINE", line
  CALL mpi_abort(mpi_comm_world, -1, ERR)
END IF
#endif
END SUBROUTINE
END module hdf5io
