subroutine read_h5()
   use decl_var
   use hdf5
   implicit none
   integer :: i,j,k,nt,l
   character(len=4)::tstamp

   character(len=15) :: dataset
   character(len=35) :: fname

   do nt = fstart,fstop,fstep
     write(tstamp,'(I4.4)')nt
     fname="../../test_dir/drun_timestep"//tstamp//".h5"
     dataset="SENGA_GRID/DRUN"

     !initialize fortran interface
     call h5open_f(hdferr)

     !Open file and dataset using the default properties.
     call h5fopen_f(trim(fname), h5f_acc_rdonly_f, file_id, hdferr)
     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, drun, gdim, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)
     call h5fclose_f(file_id , hdferr)

     fname="../../test_dir/urun_timestep"//tstamp//".h5"
     dataset="SENGA_GRID/URUN"

     !initialize fortran interface
     call h5open_f(hdferr)

     !Open file and dataset using the default properties.
     call h5fopen_f(trim(fname), h5f_acc_rdonly_f, file_id, hdferr)
     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, urun, gdim, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)
     call h5fclose_f(file_id , hdferr)

     fname="../../test_dir/vrun_timestep"//tstamp//".h5"
     dataset="SENGA_GRID/VRUN"

     !initialize fortran interface
     call h5open_f(hdferr)

     !Open file and dataset using the default properties.
     call h5fopen_f(trim(fname), h5f_acc_rdonly_f, file_id, hdferr)
     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, vrun, gdim, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)
     call h5fclose_f(file_id , hdferr)

     fname="../../test_dir/wrun_timestep"//tstamp//".h5"
     dataset="SENGA_GRID/WRUN"

     !initialize fortran interface
     call h5open_f(hdferr)

     !Open file and dataset using the default properties.
     call h5fopen_f(trim(fname), h5f_acc_rdonly_f, file_id, hdferr)
     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, wrun, gdim, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)
     call h5fclose_f(file_id , hdferr)

     fname="../../test_dir/prun_timestep"//tstamp//".h5"
     dataset="SENGA_GRID/PRN2"

     !initialize fortran interface
     call h5open_f(hdferr)

     !Open file and dataset using the default properties.
     call h5fopen_f(trim(fname), h5f_acc_rdonly_f, file_id, hdferr)
     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, prun, gdim, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)
     call h5fclose_f(file_id , hdferr)
     !print *, "here"

     fname="../../test_dir/trun_timestep"//tstamp//".h5"
     dataset="SENGA_GRID/TRN2"

     !initialize fortran interface
     call h5open_f(hdferr)

     !Open file and dataset using the default properties.
     call h5fopen_f(trim(fname), h5f_acc_rdonly_f, file_id, hdferr)
     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, trun, gdim, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)
     call h5fclose_f(file_id , hdferr)

     fname="../../test_dir/yrun_timestep"//tstamp//".h5"
     dataset="SENGA_GRID/YRUN"

     !initialize fortran interface
     call h5open_f(hdferr)

     !Open file and dataset using the default properties.
     call h5fopen_f(trim(fname), h5f_acc_rdonly_f, file_id, hdferr)
     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, yrun, gdim+1, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)
     call h5fclose_f(file_id , hdferr)

     do i=1,nx
       do j=1,ny
         do k=1,nz
           urun(i,j,k)=urun(i,j,k)/drun(i,j,k)
           vrun(i,j,k)=vrun(i,j,k)/drun(i,j,k)
           wrun(i,j,k)=wrun(i,j,k)/drun(i,j,k)
           do l=1,nspec
             yrun(l,i,j,k) = yrun(l,i,j,k)/drun(i,j,k)
           end do
         end do
       end do
     end do

     call hdf5_write(nt)
   end do
end subroutine read_h5
