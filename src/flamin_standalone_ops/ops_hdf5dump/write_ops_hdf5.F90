subroutine write_dat_to_hdf5(dataset_name, data, with_halos)

    use hdf5
    use h5lt
    use com_senga
    implicit none

    character(len=*), intent(in) :: dataset_name
    real(kind=8), dimension(nxglbl,nyglbl,nzglbl), intent(in) :: data
    integer, intent(in) :: with_halos

    character(len=*), parameter :: fname = "flame_init.h5"

!   HDF5 variables
    integer :: hdferr
    integer(hid_t):: g_id     ! grid group identifer
    integer(hid_t):: d_id     ! data set identifier
    integer(hid_t):: a_id     ! attribute identifier
    integer(hid_t):: type_id
    integer(hid_t):: dspace   ! Dataspace identifier

    integer(hsize_t) xdim(1)
    integer(hsize_t) ydim(1)
    integer(hsize_t) zdim(1)
    logical :: group_exists, dataset_exists

    integer(hsize_t) :: dims1(1), str_dims(2)
    integer(hsize_t) :: str_size
    integer(size_t)  :: size_temp

    integer, parameter :: ndim = 3
    integer(hsize_t) gdim(ndim)
    integer(hid_t)  :: file_id

!   Attributes
    integer, dimension(3) :: base = [0,0,0]
    integer :: d_m(3)
    integer :: d_p(3)
    integer, dimension(3) :: size = [nxglbl,nyglbl,nzglbl]

    integer, dimension(1) :: block_index = [0]
    integer, dimension(1) :: data_ndim = [1]

    character(len=11) :: block_str = "SENGA_GRID"
    character(len=8) :: ops_type_str = "ops_dat"
    character(len=13) :: data_type = "real(kind=8)"

    logical :: file_exists

    if(with_halos .eq. 1) then
        d_m = [-nhalox,-nhaloy,-nhaloz]
        d_p = [nhalox,nhaloy,nhaloz]
    else
        d_m = [0,0,0]  
        d_p = [0,0,0]
    end if

    gdim(1)=nxglbl
    gdim(2)=nyglbl
    gdim(3)=nzglbl

    xdim=nxglbl
    ydim=nyglbl
    zdim=nzglbl

!   ---- initialize fortran interface ----
    call h5open_f(hdferr)

!   ---- create a new file using default properties ----
    inquire(file=fname, exist=file_exists)
    if (.not. file_exists) then
        call h5fcreate_f(fname, h5f_acc_trunc_f, file_id, hdferr)
    else
        call h5fopen_f(fname, h5f_acc_rdwr_f, file_id, hdferr)
    end if

!   =========================================================================
!   Write grid data
!   =========================================================================
!   ---- create group at the root using absolute name containing the x/y/z coordinates datasets ----
    call h5lexists_f(file_id, "SENGA_GRID", group_exists, hdferr)
    if(.not. group_exists) then
        call h5gcreate_f(file_id, "/SENGA_GRID", g_id, hdferr)
    else
        call h5gopen_f(file_id, "SENGA_GRID", g_id, hdferr)
    end if

!   =========================================================================
!   Start writing dataset and its attributes
!   =========================================================================
    call h5lexists_f(g_id, trim(dataset_name), dataset_exists, hdferr)

    if(.not. dataset_exists) then

        call h5screate_simple_f(ndim, gdim, dspace, hdferr)
        call h5dcreate_f(g_id, trim(dataset_name), H5T_IEEE_F64LE, dspace, d_id, hdferr)
        call h5dwrite_f(d_id, H5T_NATIVE_DOUBLE, data, gdim, hdferr)
        call h5sclose_f(dspace, hdferr)
!       ---- Write string attribute "ops_type" ----
        call h5ltset_attribute_string_f(g_id, trim(dataset_name), "ops_type",ops_type_str, hdferr)
!       ---- Write string attribute "block" ----
        call h5ltset_attribute_string_f(g_id, trim(dataset_name), "block",block_str, hdferr)
!       ---- Write integer attribute "block_index" ----
        size_temp = 1
        call h5ltset_attribute_int_f(g_id, trim(dataset_name), "block_index",block_index, size_temp, hdferr)
!       ---- Write integer attribute "dim" ----
        size_temp = 1
        call h5ltset_attribute_int_f(g_id, trim(dataset_name), "dim",data_ndim, size_temp, hdferr)
!       ---- Write integer attribute "size" ----
        size_temp = ndim
        call h5ltset_attribute_int_f(g_id, trim(dataset_name), "size",size, size_temp, hdferr)
!       ---- Write integer attribute "d_m" ----
        size_temp = ndim
        call h5ltset_attribute_int_f(g_id, trim(dataset_name), "d_m",d_m, size_temp, hdferr)
!       ---- Write integer attribute "d_p" ----
        size_temp = ndim
        call h5ltset_attribute_int_f(g_id, trim(dataset_name), "d_p",d_p, size_temp, hdferr)
!       ---- Write integer attribute "base" ----
        size_temp = ndim
        call h5ltset_attribute_int_f(g_id, trim(dataset_name), "base",base, size_temp, hdferr)
!       ---- Write string attribute "data_type" ----
        call h5ltset_attribute_string_f(g_id, trim(dataset_name), "type",data_type, hdferr)
!       ---- close the dataset ----
        call h5dclose_f(d_id, hdferr)
    else
        print *, "Dataset already exists"
    end if

!   ---- close group ----
    call h5gclose_f(g_id, hdferr)
!   ---- close the file ----
    call h5fclose_f(file_id, hdferr)
!   ---- close fortran interface ----
    call h5close_f(hdferr)

end subroutine write_dat_to_hdf5
