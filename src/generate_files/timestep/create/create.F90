! Program to create HDF5 file suitable to restart OPS version
! change following as per required before running this
!   > nxglbl,nyglbl,nzglbl 
!   > nspcmx
!   > ntdump
!   > number of YRUN_1D array declaration and wherever it is used
!   > fndmpi - file which contain data from F90 version of senga2 to be converted to HDF5

PROGRAM create

    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    implicit none

!   Declare paramaters
    integer(kind=4), parameter :: nxglbl=504, nyglbl=252, nzglbl=252
    integer(kind=4), parameter :: nspcmx=9
    integer(kind=4) :: nxdmax,nydmax,nzdmax,ndspec,nsize
    integer(kind=4) :: iproc,ic,jc,kc,itime,idflag,ispec,index1d
    integer(kind=4), parameter :: ncdmpi=1
    integer(kind=4), parameter :: ncdmpo=2
    integer(kind=4), parameter :: ntdump=1000
    integer(kind=4) :: dtime

    real(kind=8) :: etime,tstep,errold,errldr    

!   SENGA dats present in file
    real(kind=8) :: drun(nxglbl,nyglbl,nzglbl),urun(nxglbl,nyglbl,nzglbl),  &
    vrun(nxglbl,nyglbl,nzglbl),wrun(nxglbl,nyglbl,nzglbl),  &
    erun(nxglbl,nyglbl,nzglbl), yrun(nxglbl,nyglbl,nzglbl,nspcmx)

    real(kind=8) :: drun_1D(nxglbl*nyglbl*nzglbl), urun_1D(nxglbl*nyglbl*nzglbl), &
                    vrun_1D(nxglbl*nyglbl*nzglbl), wrun_1D(nxglbl*nyglbl*nzglbl), &
                    erun_1D(nxglbl*nyglbl*nzglbl), &
                    yrun1_1D(nxglbl*nyglbl*nzglbl), yrun2_1D(nxglbl*nyglbl*nzglbl), &
                    yrun3_1D(nxglbl*nyglbl*nzglbl), yrun4_1D(nxglbl*nyglbl*nzglbl), &
                    yrun5_1D(nxglbl*nyglbl*nzglbl), yrun6_1D(nxglbl*nyglbl*nzglbl), &
                    yrun7_1D(nxglbl*nyglbl*nzglbl), yrun8_1D(nxglbl*nyglbl*nzglbl), &
                    yrun9_1D(nxglbl*nyglbl*nzglbl)

    character(len=60) :: fndmpi, fndmpo, fname
    character(len=8)  :: citime
    character(len=11) :: pndmpi
    character(len=4)  :: pnxdat
    parameter(pndmpi = 'output/dmpi',pnxdat = '.dat')
    character(len=3) :: pnxhdf
    parameter(pnxhdf = '.h5')
    character(len=6) :: pnproc
    character(len=1) :: pnflag

!   Declare ops_block and ops_dats
    TYPE(ops_block) :: senga_grid
    TYPE(ops_dat) :: d_drun, d_urun, d_vrun, d_wrun, d_erun
    TYPE(ops_dat) :: d_yrun(nspcmx)
    TYPE(ops_dat) :: d_temp

!   Declare stencil
    TYPE(ops_stencil) :: s3d_000

    integer(kind=4) :: a3d_000(3) = [0,0,0]
    integer(kind=4) :: d_size(3) = [nxglbl,nyglbl,nzglbl]
    integer(kind=4) :: d_base(3) = [1,1,1]
    integer(kind=4) :: d_p(3) = [0,0,0]
    integer(kind=4) :: d_m(3) = [0,0,0]

    real(kind=8), dimension(:), allocatable :: temp_real_null

    character(len=20) :: buf

    integer(kind=4) :: rangexyz(6)

!   Init OPS environment
    call ops_init(2)

    call ops_decl_block(3, senga_grid, "SENGA_GRID")

!   Create OPS dats and refer to 1D arrays so that OPS_DATS will be automatically set
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, drun_1D, d_drun, "real(kind=8)", "DRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, urun_1D, d_urun, "real(kind=8)", "URUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vrun_1D, d_vrun, "real(kind=8)", "VRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wrun_1D, d_wrun, "real(kind=8)", "WRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, erun_1D, d_erun, "real(kind=8)", "ERUN")

    ispec = 1
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun1_1D, d_yrun(ispec), "real(kind=8)", trim(buf))
    
    ispec = 2
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun2_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 3
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun3_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 4
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun4_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 5
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun5_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 6
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun6_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 7
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun7_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 8
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun8_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 9
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun9_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_temp, "real(kind=8)", "temp")

    call ops_decl_stencil( 3, 1, a3d_000, s3d_000, "0,0,0")

    call ops_partition(" ")

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_temp, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!   Read the last dump file
    fndmpi = 'h2flame_3D_timestep1770000.dat'
    OPEN(UNIT=ncdmpi,FILE=trim(fndmpi),STATUS='OLD',FORM='UNFORMATTED')
    READ(ncdmpi)nxdmax,nydmax,nzdmax,ndspec, drun,urun,vrun,wrun,erun,yrun, &
               etime,tstep,errold,errldr

    CLOSE(ncdmpi)

!   SIZE ERROR CHECK
    IF(nxdmax /= nxglbl)WRITE(6,*)'Dump input size error: x'
    IF(nydmax /= nyglbl)WRITE(6,*)'Dump input size error: y'
    IF(nzdmax /= nzglbl)WRITE(6,*)'Dump input size error: z'
    IF(ndspec /= nspcmx)WRITE(6,*)'Dump input size error: species'

    write(*,*) "nxglbl: ",nxdmax,"  nyglbl: ",nydmax,"  nzglbl:",nzdmax, "  ndspec:",ndspec
    write(*,*) "etime: ",etime,"  tstep: ",tstep,"  errold:",errold,"  errldr:",errldr

    index1d = 0
    DO kc = 1, nzglbl
        DO jc = 1, nyglbl
            DO ic = 1, nxglbl
                index1d = index1d + 1

                drun_1D(index1d) = drun(ic,jc,kc)
                urun_1D(index1d) = urun(ic,jc,kc)
                vrun_1D(index1d) = vrun(ic,jc,kc)
                wrun_1D(index1d) = wrun(ic,jc,kc)
                erun_1D(index1d) = erun(ic,jc,kc)

                yrun1_1D(index1d) = yrun(ic,jc,kc,1)
                yrun2_1D(index1d) = yrun(ic,jc,kc,2)
                yrun3_1D(index1d) = yrun(ic,jc,kc,3)
                yrun4_1D(index1d) = yrun(ic,jc,kc,4)
                yrun5_1D(index1d) = yrun(ic,jc,kc,5)
                yrun6_1D(index1d) = yrun(ic,jc,kc,6)
                yrun7_1D(index1d) = yrun(ic,jc,kc,7)
                yrun8_1D(index1d) = yrun(ic,jc,kc,8)
                yrun9_1D(index1d) = yrun(ic,jc,kc,9)

            END DO
        END DO
    END DO

!   Write etime,tstep,errold,errldr to dmpi<timestep>.dat file
    DO iproc = 0, 127
        WRITE(pnproc,'(I6.6)')iproc
        idflag = 0
        WRITE(pnflag,'(I1)')idflag
        fndmpo = pndmpi//pnproc//pnflag//pnxdat
        OPEN(UNIT=ncdmpo,FILE=fndmpo,STATUS='NEW',FORM='UNFORMATTED')
        WRITE(ncdmpo)nxdmax,nydmax,nzdmax,ndspec, &
                     etime,tstep,errold,errldr
        CLOSE(ncdmpo)

        idflag = 1
        WRITE(pnflag,'(I1)')idflag
        fndmpo = pndmpi//pnproc//pnflag//pnxdat
        OPEN(UNIT=ncdmpo,FILE=fndmpo,STATUS='NEW',FORM='UNFORMATTED')
        WRITE(ncdmpo)nxdmax,nydmax,nzdmax,ndspec, &
                     etime,tstep,errold,errldr
        CLOSE(ncdmpo)

    END DO

!   Write OPS_DATS to HDF5 file
    itime = 1770000
    dtime=INT(itime/ntdump)
    WRITE(citime,'(I8.8)') dtime
    fname = 'timestep'//citime//pnxhdf

    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))

    call ops_fetch_dat_hdf5_file(d_drun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vrun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wrun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_erun, trim(fname))

    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_yrun(ispec), trim(fname))
    END DO

!   Finalize OPS env
    call ops_exit( )

END PROGRAM create
