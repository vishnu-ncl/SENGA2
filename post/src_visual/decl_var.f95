module decl_var
   use hdf5
   implicit none

   integer,parameter :: nvars=6
   real(kind=8),allocatable, dimension(:,:,:) :: urun,vrun,wrun
   real(kind=8),allocatable, dimension(:,:,:) :: prun,trun,drun
   real(kind=8),allocatable, dimension(:,:,:,:) :: yrun
   real(kind=8) :: xlen,ylen,zlen,deltax,deltay,deltaz

   integer :: nx,ny,nz,nspec,nxproc,nyproc,nzproc
   integer :: fstart = 1
   integer :: fstop = 1
   integer :: fstep = 1
 
   integer :: hdferr
   integer(hid_t)  :: file_id, dset_id
   !-mod-----------------
   integer, parameter :: ndim = 3
   integer, parameter :: nxdim = 1, nydim = 1, nzdim = 1
   integer, parameter :: switch = 3
   integer(hsize_t) gdim(ndim)

   type tag
     character(len=12) wrtdst
   end type
   
   type varspace
     real, dimension(:,:,:), allocatable :: var
   end type

   type (tag) :: wrtvar(nvars)
   type (varspace) :: work(nvars)
end module decl_var
