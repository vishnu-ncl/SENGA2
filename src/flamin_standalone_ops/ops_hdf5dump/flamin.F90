PROGRAM flamin

    use com_senga
    use, intrinsic :: ISO_C_BINDING

    implicit none

!   *************************************************************************

!   FLAMIN
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   28-DEC-2003:  CREATED
!   08-JAN-2005:  RSC INITIAL 1D LAMINAR FLAME PROFILE

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   SETS INITIAL THERMOCHEMICAL FIELD
!   1D LAMINAR FLAME PROFILE (LEFT OR RIGHT FACING)
!   SPECIAL FOR 21 STEP HYDROGEN MECHAMISM

!   *************************************************************************

!   LOCAL DATA
!   ==========
    real(kind=8) :: drun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: urun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: vrun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: wrun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: trun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: yrun(nxglbl,nyglbl,nzglbl,nspcmx)
    integer(kind=4) :: ix, ic, ispec, fnum

    character(len=20) :: buf

!   BEGIN
!   =====

!   =========================================================================

    fnum = 111
    OPEN(UNIT=fnum,FILE='h2_init.dat',FORM="FORMATTED",STATUS= "OLD")
    DO ic=1,nxglbl
        READ(fnum,*)ix,urun(ic,1,1),trun(ic,1,1),drun(ic,1,1),  &
            yrun(ic,1,1,1),yrun(ic,1,1,2),yrun(ic,1,1,3),yrun(ic,1,1,4),  &
            yrun(ic,1,1,5),yrun(ic,1,1,6),yrun(ic,1,1,7),yrun(ic,1,1,8), yrun(ic,1,1,9)
    END DO
    CLOSE(fnum)

    DO ic = 1,nxglbl
        vrun(ic,1,1) = 0.0_8
        wrun(ic,1,1) = 0.0_8
    END DO


    WRITE(buf,"(A4)") "DRUN"
    call write_dat_to_hdf5(trim(buf), drun, 0)
    WRITE(buf,"(A4)") "URUN"
    call write_dat_to_hdf5(trim(buf), urun, 0)
    WRITE(buf,"(A4)") "VRUN"
    call write_dat_to_hdf5(trim(buf), vrun, 0)
    WRITE(buf,"(A4)") "WRUN"
    call write_dat_to_hdf5(trim(buf), wrun, 0)
    WRITE(buf,"(A4)") "TRUN"
    call write_dat_to_hdf5(trim(buf), trun, 0)

    DO ispec = 1,nspcmx
        WRITE(buf,"(A4,I2.2)") "YRUN",ispec
        call write_dat_to_hdf5(trim(buf), yrun(:,:,:,ispec), 0)
    END DO

!   =========================================================================

END PROGRAM flamin
