subroutine read_cont()
   use decl_var
   use chmech
   implicit none
   character(len=20) :: fcont
   integer :: ncont,i,j,k

   fcont = "../../input/cont.dat"
   ncont = 11

   open(unit=ncont,file=fcont,status='old',action='read',form='formatted')
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)xlen,ylen,zlen
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)nx,ny,nz
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)nxproc,nyproc,nzproc
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)
   read(ncont,*)nspec
   close(ncont)

   allocate(drun(nx,ny,nz))
   allocate(urun(nx,ny,nz))
   allocate(vrun(nx,ny,nz))
   allocate(wrun(nx,ny,nz))
   allocate(prun(nx,ny,nz))
   allocate(trun(nx,ny,nz))
   allocate(yrun(nx,ny,nz,nspec))
   allocate(wrtspc(nspec))
   allocate(wrtrrt(nspec))
end subroutine read_cont
