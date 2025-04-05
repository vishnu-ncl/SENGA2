SUBROUTINE parfer

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   PARFER
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   11-MAY-2003:  CREATED
!   04-JAN-2007:  RSC REVISE PARALLEL RECEIVES
!   26-OCT-2008:  RSC/TDD BUG FIX JSTAB

!   DESCRIPTION
!   -----------
!   CARRIES OUT TRANSFER OF PARALLEL DATA

!   *************************************************************************


!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    integer(kind=4) :: rangexyz(6)

!   =========================================================================

!   Forcing data exchange
!   Making data dirty
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    call ops_par_loop(maths_kernel_set_dirty, "set_dirty", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(maths_kernel_set_dirty, "set_dirty", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(maths_kernel_set_dirty, "set_dirty", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(maths_kernel_set_dirty, "set_dirty", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(maths_kernel_set_dirty, "set_dirty", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    DO ispec = 1,nspec
        call ops_par_loop(maths_kernel_set_dirty, "set_dirty", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
    END DO

!   Trigerring halo exchange
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_halo_exchange_xdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drhs, 1, s3d_p500_to_m500_x, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_xdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_p500_to_m500_x, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_xdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_p500_to_m500_x, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_xdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_p500_to_m500_x, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_xdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_p500_to_m500_x, "real(kind=8)", OPS_READ))
    DO ispec = 1,nspec
        call ops_par_loop(maths_kernel_halo_exchange_xdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_yrhs(ispec), 1, s3d_p500_to_m500_x, "real(kind=8)", OPS_READ))
    END DO

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_halo_exchange_ydir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drhs, 1, s3d_p050_to_m050_y, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_ydir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_p050_to_m050_y, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_ydir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_p050_to_m050_y, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_ydir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_p050_to_m050_y, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_ydir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_p050_to_m050_y, "real(kind=8)", OPS_READ))
    DO ispec = 1,nspec
        call ops_par_loop(maths_kernel_halo_exchange_ydir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_yrhs(ispec), 1, s3d_p050_to_m050_y, "real(kind=8)", OPS_READ))
    END DO

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_halo_exchange_zdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drhs, 1, s3d_p005_to_m005_z, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_zdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_p005_to_m005_z, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_zdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_p005_to_m005_z, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_zdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_p005_to_m005_z, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_halo_exchange_zdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_p005_to_m005_z, "real(kind=8)", OPS_READ))
    DO ispec = 1,nspec
        call ops_par_loop(maths_kernel_halo_exchange_zdir, "halo_exchange", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_yrhs(ispec), 1, s3d_p005_to_m005_z, "real(kind=8)", OPS_READ))
    END DO

!   SPECIAL CASE OF PERIODIC TRANSFER ON SINGLE PROCESSOR
!   -----------------------------------------------------

!   X-DIRECTION
!   ONLY NEED TO CHECK ONE END
    IF(nendxl == nperi) call ops_halo_transfer(halos_grp_x)

!   Y-DIRECTION
!   ONLY NEED TO CHECK ONE END
!   NOTE EXTENDED X-LIMITS FOR Y TRANSFERS
    IF(nendyl == nperi .and. nhaloy /= 0) call ops_halo_transfer(halos_grp_y)

!   Z-DIRECTION
!   ONLY NEED TO CHECK ONE END
!   NOTE EXTENDED X- AND Y-LIMITS FOR Z TRANSFERS
    IF(nendzl == nperi .and. nhaloz /= 0) call ops_halo_transfer(halos_grp_z)

!   =========================================================================

END SUBROUTINE parfer
