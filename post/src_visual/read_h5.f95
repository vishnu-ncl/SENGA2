subroutine read_h5()
   use decl_var
   use hdf5
   implicit none
   integer :: i,j,k,nt,ispec
   character(len=8)::tstamp
   character(len=6)::spec

   character(len=15) :: dataset
   character(len=17) :: dataset2
   character(len=35) :: fname

   do nt = fstart,fstop,fstep
     write(tstamp,'(I8.8)')nt
     fname="../../output/timestep"//tstamp//".h5"
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

     dataset="SENGA_GRID/URUN"
     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, urun, gdim, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)

     dataset="SENGA_GRID/VRUN"

     !Open dataset using the default properties.
     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, vrun, gdim, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)

     dataset="SENGA_GRID/WRUN"
     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, wrun, gdim, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)

     dataset="SENGA_GRID/PRN2"

     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, prun, gdim, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)

     dataset="SENGA_GRID/TRN2"

     call h5dopen_f(file_id, dataset, dset_id, hdferr)

     !Read the data using the default properties.
     call h5dread_f(dset_id, h5t_native_double, trun, gdim, hdferr)

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)

     !Open file and dataset using the default properties.
     do i=1,nspec
       write(spec,"(A4,I2.2)") "YRUN",i
       dataset2="SENGA_GRID/"//spec
       call h5dopen_f(file_id, dataset2, dset_id, hdferr)

       !Read the data using the default properties.
       call h5dread_f(dset_id, h5t_native_double, yrun(:,:,:,i), gdim, hdferr)
     end do

     !Close and release resources.
     call h5dclose_f(dset_id , hdferr)
     call h5fclose_f(file_id , hdferr)

     do i=1,nx
       do j=1,ny
         do k=1,nz
           urun(i,j,k)=urun(i,j,k)/drun(i,j,k)
           vrun(i,j,k)=vrun(i,j,k)/drun(i,j,k)
           wrun(i,j,k)=wrun(i,j,k)/drun(i,j,k)
           do ispec=1,nspec
             yrun(i,j,k,ispec)=yrun(i,j,k,ispec)/drun(i,j,k)
           end do
         end do
       end do
     end do

     call hdf5_write(nt)
   end do
end subroutine read_h5
