SUBROUTINE rhsvel

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   RHSVEL
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   08-DEC-2002:  CREATED
!   17-APR-2013:  RSC MIXTURE AVERAGED TRANSPORT

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   COMPUTES RIGHT-HAND-SIDES FOR TIME INTEGRATION
!   OF CONTINUITY AND MOMENTUM EQUATIONS
!   EVALUATES PRESSURE WORK AND VISCOUS WORK TERMS IN ENERGY EQUATION

!   *************************************************************************


!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------

!   -------------------------------------------------------------------------


!   LOCAL DATA
!   ==========
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   CONVERT VELOCITIES
!     ------------------

!   U,V,WRHS CONTAIN RHO U,V,W: CONVERT TO U,V,W
!   U,V,W HELD IN U,V,WTMP THROUGHOUT THIS ROUTINE
!   U,V,W ARE PARALLEL
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
!    call ops_par_loop(maths_kernel_eqU, "A = B/C", senga_grid, 3, rangexyz,  &
!                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
!                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
!                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

!    call ops_par_loop(maths_kernel_eqU, "A = B/C", senga_grid, 3, rangexyz,  &
!                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
!                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
!                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

!    call ops_par_loop(maths_kernel_eqU, "A = B/C", senga_grid, 3, rangexyz,  &
!                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
!                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
!                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))
    call ops_par_loop(maths_kernel_eqU_fused, "A = B/C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

!   =========================================================================

!   COLLECT VELOCITY COMPONENTS FOR BCs
!   -----------------------------------

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
    call dfbydy(d_utmp,d_store1)
    call dfbydz(d_utmp,d_store2)

!   X-DIRECTION
    IF(fxlcnv)THEN
        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velcomp_xdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t2bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

    END IF
    IF(fxrcnv)THEN
        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velcomp_xdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t2bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
    call dfbydx(d_vtmp,d_store1)
    call dfbydz(d_vtmp,d_store2)

!   Y-DIRECTION
    IF(fylcnv)THEN
        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(boundary_kernel_velcomp_ydir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t2byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

    END IF
    IF(fyrcnv)THEN
        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velcomp_ydir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t2byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
    call dfbydx(d_wtmp,d_store1)
    call dfbydy(d_wtmp,d_store2)

!   Z-DIRECTION
    IF(fzlcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(boundary_kernel_velcomp_zdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t2bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

    END IF
    IF(fzrcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(boundary_kernel_velcomp_zdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t2bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
    call dfbydy(d_vtmp,d_store1)
    call dfbydz(d_vtmp,d_store2)

!   X-DIRECTION
    IF(fxlcnv)THEN
        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_eqA_xdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t3bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF
    IF(fxrcnv)THEN
        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_eqA_xdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t3bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
!   NC and VM changed here: ycorrection
    call dfbydx(d_wtmp,d_store1)
    call dfbydz(d_wtmp,d_store2)

!   Y-DIRECTION
    IF(fylcnv)THEN
        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(boundary_kernel_eqA_ydir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t4byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF
    IF(fyrcnv)THEN
        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_eqA_ydir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t4byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
    call dfbydx(d_utmp,d_store1)
    call dfbydy(d_utmp,d_store2)

!   Z-DIRECTION
    IF(fzlcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(boundary_kernel_eqA_zdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t3bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF
    IF(fzrcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(boundary_kernel_eqA_zdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t3bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
    call dfbydy(d_wtmp,d_store1)
    call dfbydz(d_wtmp,d_store2)

!   X-DIRECTION
    IF(fxlcnv)THEN
        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_eqA_xdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t4bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF
    IF(fxrcnv)THEN
        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_eqA_xdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t4bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
!   NC and VM changed here: ycorrection
    call dfbydx(d_utmp,d_store1)
    call dfbydz(d_utmp,d_store2)

!   Y-DIRECTION
    IF(fylcnv)THEN
        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(boundary_kernel_eqA_ydir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t3byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF
    IF(fyrcnv)THEN
        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_eqA_ydir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t3byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
    call dfbydx(d_vtmp,d_store1)
    call dfbydy(d_vtmp,d_store2)

!   Z-DIRECTION
    IF(fzlcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(boundary_kernel_eqA_zdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t4bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF
    IF(fzrcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(boundary_kernel_eqA_zdir, "COLLECT VELOCITY COMPONENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_t4bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF


!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   MOMENTUM EQUATIONS: CONVECTIVE TERMS
!   ------------------------------------
!   RHO U U
!   RHO U U IS PARALLEL
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    call ops_par_loop(maths_kernel_eqW, "A = B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ))

!   D/DX RHO U U
!   STRAIGHT INTO STORE4 FOR NOW
    call dfbydx(d_store7,d_store4)


!   U-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   (HALF) D/DX RHO U U: ALREADY IN STORE4
!                                                  STORE4 = U CONVECTIVE TERMS
!                                                         U,V,WRHS = RHO U,V,W
!    =========================================================================

!   RHO U V
!   RHO U V IS PARALLEL
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    call ops_par_loop(maths_kernel_eqW, "A = B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ))

!   D/DY RHO V U
!   D/DX RHO U V
    call dfbydy(d_store7,d_store1)
    call dfbydx(d_store7,d_store5)


!   U-EQUATION: CONVECTIVE TERMS
!   V-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   (HALF) D/DY RHO V U
!   (HALF) D/DX RHO U V: ALREADY IN STORE5
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                             STORE4,5 = U,V CONVECTIVE TERMS
!                                                        U,V,WRHS = RHO U,V,W
!   =========================================================================

!   RHO U W
!   RHO U W IS PARALLEL
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    call ops_par_loop(maths_kernel_eqW, "A = B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ))

!   D/DZ RHO W U
!   D/DX RHO U W
    call dfbydz(d_store7,d_store1)
    call dfbydx(d_store7,d_store6)

!   U-EQUATION: CONVECTIVE TERMS
!   W-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   (HALF) D/DZ RHO W U
!   (HALF) D/DX RHO U W: ALREADY IN STORE6
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                         STORE4,5,6 = U,V,W CONVECTIVE TERMS
!                                                        U,V,WRHS = RHO U,V,W
!   =========================================================================

!   RHO V V
!   RHO V V IS PARALLEL
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    call ops_par_loop(maths_kernel_eqW, "A = B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ))

!   D/DY RHO V V
    call dfbydy(d_store7,d_store1)

!   V-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   (HALF) D/DY RHO V V
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                         STORE4,5,6 = U,V,W CONVECTIVE TERMS
!                                                        U,V,WRHS = RHO U,V,W
!   =========================================================================

!   RHO V W
!   RHO V W IS PARALLEL
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    call ops_par_loop(maths_kernel_eqW, "A = B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ))

!   D/DZ RHO W V
!   D/DY RHO V W
    call dfbydz(d_store7,d_store1)
    call dfbydy(d_store7,d_store2)

!   V-EQUATION: CONVECTIVE TERMS
!   W-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   (HALF) D/DZ RHO W V
!   (HALF) D/DY RHO V W
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                         STORE4,5,6 = U,V,W CONVECTIVE TERMS
!                                                        U,V,WRHS = RHO U,V,W
!   =========================================================================

!   RHO W W
!   RHO W W IS PARALLEL
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    call ops_par_loop(maths_kernel_eqW, "A = B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ))

!   D/DZ RHO W W
    call dfbydz(d_store7,d_store1)

!   W-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   (HALF) D/DZ RHO W W
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                         STORE4,5,6 = U,V,W CONVECTIVE TERMS
!                                                        U,V,WRHS = RHO U,V,W
!   =========================================================================

!   VELOCITY NORMAL DERIVATIVES
!   ---------------------------
!   DUDX,DVDY,DWDZ
    call dfbydx(d_utmp,d_store1)
    call dfbydy(d_vtmp,d_store2)
    call dfbydz(d_wtmp,d_store3)

!                                                 STORE1,2,3 = DUDX,DVDY,DWDZ
!                                         STORE4,5,6 = U,V,W CONVECTIVE TERMS
!                                                        U,V,WRHS = RHO U,V,W
!   =========================================================================

!   COLLECT VELOCITY DERIVATIVES FOR BCs
!   ------------------------------------

!   X-DIRECTION: DUDX
    IF(fxlcnv)THEN
        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_xdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF
    IF(fxrcnv)THEN
        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_xdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF

!   Y-DIRECTION: DVDY
    IF(fylcnv)THEN
        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_ydir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF
    IF(fyrcnv)THEN
        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_ydir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF

!   Z-DIRECTION: DWDZ
    IF(fzlcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(boundary_kernel_velderiv_zdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF
    IF(fzrcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_zdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF

!   =========================================================================

!   U-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   DIV RHO U U  + U DIV RHO U
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAB, "A = B+C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_divm, 1, s3d_000, "real(kind=8)", OPS_READ))

!   RHO U DUDX
    call ops_par_loop(maths_kernel_eqW, "A = B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!   HALF DIV RHO U U + HALF RHO U DUDX + HALF U DIV RHO U
!   STORE IN URHS
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAO, "A = -half*(B+C)", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

!   V-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   DIV RHO U V + V DIV RHO U
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAB, "A = B+C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_divm, 1, s3d_000, "real(kind=8)", OPS_READ))

!   RHO V DVDY
    call ops_par_loop(maths_kernel_eqW, "A = B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!   HALF DIV RHO U V + HALF RHO V DVDY + HALF V DIV RHO U
!   STORE IN VRHS
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAO, "A = -half*(B+C)", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

!   W-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   DIV RHO U W + W DIV RHO U
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAB, "A = B+C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_divm, 1, s3d_000, "real(kind=8)", OPS_READ))

!   RHO W DWDZ
    call ops_par_loop(maths_kernel_eqW, "A = B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!   HALF DIV RHO U W + HALF RHO W DWDZ + HALF W DIV RHO U
!   STORE IN WRHS
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAO, "A = -half*(B+C)", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

!   U,V,WRHS CONTAIN U,V,W SOURCE TERMS THROUGHOUT REMAINDER OF THIS ROUTINE
!                                                STORE1,2,3 = DUDX,DVDY,DWDZ
!   ========================================================================

!   PRESSURE GRADIENTS
!   ------------------

!   DPDX,DPDY,DPDZ
    call dfbydx(d_prun,d_store4)
    call dfbydy(d_prun,d_store5)
    call dfbydz(d_prun,d_store6)

!                                              STORE1,2,3 = DUDX,DVDY,DWDZ
!                                              STORE4,5,6 = DPDX,DPDY,DPDZ
!   ======================================================================

!   COLLECT PRESSURE AND ITS GRADIENTS FOR BCs
!   ------------------------------------------

!   X-DIRECTION: P AND DPDX
    IF(fxlcnv)THEN
        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_pressure_xdir, "COLLECT PRESSURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_strpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t4bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t51bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t52bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

    END IF
    IF(fxrcnv)THEN
        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_pressure_xdir, "COLLECT PRESSURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_strpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_bcl5xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t4bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t51bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t52bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF

!   Y-DIRECTION: P AND DPDY
    IF(fylcnv)THEN
        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(boundary_kernel_pressure_ydir, "COLLECT PRESSURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_strpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t4byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t51byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t52byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

    END IF
    IF(fyrcnv)THEN
        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_pressure_ydir, "COLLECT PRESSURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_strpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_bcl5yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t4byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t51byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t52byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF

!   Z-DIRECTION: P AND DPDZ
    IF(fzlcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(boundary_kernel_pressure_zdir, "COLLECT PRESSURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_strpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t4bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t51bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t52bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

    END IF
    IF(fzrcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(boundary_kernel_pressure_zdir, "COLLECT PRESSURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_strpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_bcl5zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t4bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_t51bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t52bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))
    END IF

!   =========================================================================

!   U,V,W-EQUATIONS: PRESSURE GRADIENT TERMS
!   E-EQUATION: PRESSURE WORK TERMS
!   -------------------------------
!   DPDX,DpDY,DPDZ
!   U DPDX + V DPDY + W DPDZ
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqS, "A = A-B", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqS, "A = A-B", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqS, "A = A-B", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAL,"A = A-B*C-D*E-F*G", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                 STORE1,2,3 = DUDX,DVDY,DWDZ
!   =========================================================================

!   E-EQUATION: PRESSURE WORK TERMS
!   -------------------------------
!   P DIV U
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAK, "A = A-B*(C+D+E)", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                 STORE1,2,3 = DUDX,DVDY,DWDZ
!   =========================================================================

!   VISCOSITY
!   ---------

!   VISCOSITY IS PARALLEL
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    call ops_par_loop(maths_kernel_eqG, "A = A*var", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_gbl(prantl_ops, 1,"real(kind=8)", OPS_READ))

!                                                 STORE1,2,3 = DUDX,DVDY,DWDZ
!    -------------------------------------------------------------------------

!   MIXTURE-AVERAGED TRANSPORT
!   RSC 17-APR-2013
!   DIAGNOSTICS
!   WRITE(6,*)'RHSVEL: visc: ',FLMAVT
    IF(flmavt)THEN
        rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_difmix, 1, s3d_000, "real(kind=8)", OPS_READ))

!   DIAGNOSTICS
!   KC = 1
!   JC = 1
!   WRITE(6,'(4I5)')ITIME,IRKSTP,JC,KC
!   DO IC = ISTAB,ISTOB
!       WRITE(6,'(I5,2(1PE15.7))')IC,TRANSP(IC,JC,KC)
!   ENDDO

!   DIAGNOSTICS
!   KC = 1
!   JC = 1
!   IC = 1
!       WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
!   IC = 2
!       WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
!   IC = 500
!       WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
!   IC = 1000
!       WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
!   IC = 1001
!       WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)

    END IF

!   -------------------------------------------------------------------------

!   VISCOUS TERMS: TAUXXb,e,f
!   -------------

!   DVDY+DWDZ
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqR, "A = B+C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!   TAUXXb,e,f
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqN, "A = var1*B-var2*A", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(fthd_ops, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(tthd_ops, 1, "real(kind=8)", OPS_READ))

!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUXXb,e,f DUDX
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAD, "A = A+B*C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                 STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                         STORE6 = TAUXXb,e,f
!   =========================================================================

!   VISCOSITY GRADIENT: X COMPONENT
!   ------------------
    call dfbydx(d_transp,d_store4)

!                                                 STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                              STORE4 = DMUDX
!                                                         STORE6 = TAUXXb,e,f
!   =========================================================================

!   U-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUXX,Xb,e,f
!   U TAUXX,Xb,e,f

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqV, "A = A*B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN X: TAUXX,X TERM ZERO ON END POINTS
    IF(fxlvsn) call zeroxl(d_store6)
    IF(fxrvsn) call zeroxr(d_store6)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAA, "A = A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                 STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                              STORE4 = DMUDX
!   =========================================================================

!   VISCOUS TERMS: TAUYYb,e,f
!   -------------

!   DUDX+DWDZ
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqR, "A = B+C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!   TAUYYb,e,f
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqN, "A = var1*B-var2*A", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(fthd_ops, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(tthd_ops, 1, "real(kind=8)", OPS_READ))

!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUYYb,e,f DVDY
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAD, "A = A+B*C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                 STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                              STORE4 = DMUDX
!                                                         STORE6 = TAUYYb,e,f
!   =========================================================================

!   VISCOSITY GRADIENT: Y COMPONENT
!   ------------------
    call dfbydy(d_transp,d_store5)

!                                                 STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                          STORE4,5 = DMUDX,Y
!                                                         STORE6 = TAUYYb,e,f
!   =========================================================================

!   V-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUYY,Yb,e,f
!   V TAUYY,Yb,e,f
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqV, "A = A*B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN Y: TAUYY,Y TERM ZERO ON END POINTS
    IF(fylvsn) call zeroyl(d_store6)
    IF(fyrvsn) call zeroyr(d_store6)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAA, "A = A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                 STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                          STORE4,5 = DMUDX,Y
!   =========================================================================

!   VISCOUS TERMS: TAUZZb,e,f
!   -------------

!   DUDX+DVDY
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!   TAUZZb,e,f
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqN, "A = var1*B-var2*A", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(fthd_ops, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(tthd_ops, 1, "real(kind=8)", OPS_READ))

!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUZZb,e,f DWDZ
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAD, "A = A+B*C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                               STORE3 = DWDZ
!                                                          STORE4,5 = DMUDX,Y
!                                                         STORE1 = TAUZZb,e,f
!   =========================================================================

!   VISCOSITY GRADIENT: Z COMPONENT
!   ------------------
    call dfbydz(d_transp,d_store6)

!                                                      STORE4,5,6 = DMUDX,Y,Z
!                                                         STORE1 = TAUZZb,e,f
!   =========================================================================

!   W-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUZZ,Zb,e,f
!   W TAUZZ,Zb,e,f

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqV, "A = A*B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN Z: TAUZZ,Z TERM ZERO ON END POINTS
    IF(fzlvsn) call zerozl(d_store1)
    IF(fzrvsn) call zerozr(d_store1)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAA, "A = A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   VELOCITY CROSS-DERIVATIVES
!   --------------------------

!   DUDY
!   ----
    call dfbydy(d_utmp,d_store1)

!   COLLECT VELOCITY DERIVATIVE FOR BCs
!   -----------------------------------
!   Y-DIRECTION: DUDY
    IF(fylcnv)THEN
        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_ydir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF
    IF(fyrcnv)THEN
        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_ydir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF

!                                                               STORE1 = DUDY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   U-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   HALF RHO V DUDY

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAR, "A = A-half*B*C*D", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                               STORE1 = DUDY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   DVDX
!   ----
    call dfbydx(d_vtmp,d_store2)


!   COLLECT VELOCITY DERIVATIVES FOR BCs
!   ------------------------------------
!   X-DIRECTION: DVDX
    IF(fxlcnv)THEN
        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_xdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF
    IF(fxrcnv)THEN
        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_xdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF

!                                                               STORE1 = DUDY
!                                                               STORE2 = DVDX
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   V-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   HALF RHO U DVDX
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAR, "A = A-half*B*C*D", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                               STORE1 = DUDY
!                                                               STORE2 = DVDX
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   VELOCITY SECOND CROSS DERIVATIVES
!   ---------------------------------

!   DUDY+DVDX
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqR, "A = B+C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!   D2UDXY
    call d2fdxy(d_utmp,d_store1)

!   D2VDXY
    call d2fdxy(d_vtmp,d_store2)

!                                                             STORE1 = D2UDXY
!                                                             STORE2 = D2VDXY
!                                                        STORE7 = (DUDY+DVDX)
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   U-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUXX,Xc
!   U TAUXX,Xc
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqI, "A = var*B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(tthd_ops, 1, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN X: TAUXX,X TERMS ZERO ON END POINTS
    IF(fxlvsn) call zeroxl(d_store3)
    IF(fxrvsn) call zeroxr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqS, "A = A-B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAJ, "A = A-B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                             STORE1 = D2UDXY
!                                                             STORE2 = D2VDXY
!                                                        STORE7 = (DUDY+DVDX)
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   V-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUYY,Yc
!   V TAUYY,Yc

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqI, "A = var*B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(tthd_ops, 1, "real(kind=8)", OPS_READ))


!   BOUNDARY CONDITIONS
!   BC IN Y: TAUYY,Y TERMS ZERO ON END POINTS
    IF(fylvsn) call zeroyl(d_store3)
    IF(fyrvsn) call zeroyr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqS, "A = A-B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAJ, "A = A-B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                             STORE1 = D2UDXY
!                                                             STORE2 = D2VDXY
!                                                        STORE7 = (DUDY+DVDX)
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   D2UDY2
    call d2fdy2(d_utmp,d_store3)

!   D2UDY2+D2VDXY
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!   U-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUXY,Ya,b,c
!   U TAUXY,Ya,b,c

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAE, "A = A*B+C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN Y: TAUXY,Y TERM ZERO ON END POINTS
    IF(fylvst) call zeroyl(d_store3)
    IF(fyrvst) call zeroyr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAA, "A = A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                             STORE1 = D2UDXY
!                                                        STORE7 = (DUDY+DVDX)
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   D2VDX2
    call d2fdx2(d_vtmp,d_store3)

!   D2UDXY+D2VDX2
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!   V-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUYX,Xa,b,c
!   V TAUYX,Xa,b,c

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAE, "A = A*B+C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN X: TAUYX,X TERM ZERO ON END POINTS
    IF(fxlvst) call zeroxl(d_store3)
    IF(fxrvst) call zeroxr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAA, "A = A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                        STORE7 = (DUDY+DVDX)
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUXY(DUDY+DVDX)
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAC, "A = A+B*C*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   DUDZ
!   ----
    call dfbydz(d_utmp,d_store1)

!   COLLECT VELOCITY DERIVATIVES FOR BCs
!   ------------------------------------
!   Z-DIRECTION: DUDZ
    IF(fzlcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(boundary_kernel_velderiv_zdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF
    IF(fzrcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_zdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF

!                                                               STORE1 = DUDZ
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   U-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   HALF RHO W DUDZ

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAR, "A = A-half*B*C*D", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                               STORE1 = DUDZ
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   DWDX
!   ----
    call dfbydx(d_wtmp,d_store2)

!   COLLECT VELOCITY DERIVATIVES FOR BCs
!   ------------------------------------
!   X-DIRECTION: DWDX
    IF(fxlcnv)THEN
        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_xdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF
    IF(fxrcnv)THEN
        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_xdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF

!                                                               STORE1 = DUDZ
!                                                               STORE2 = DWDX
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   W-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   HALF RHO U DWDX

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAR, "A = A-half*B*C*D", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                               STORE1 = DUDZ
!                                                               STORE2 = DWDX
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   VELOCITY SECOND CROSS DERIVATIVES
!   ---------------------------------

!   DUDZ+DWDX
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqR, "A = B+C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!   D2UDXZ
    call d2fdxz(d_utmp,d_store1)

!   D2WDXZ
    call d2fdxz(d_wtmp,d_store2)

!                                                             STORE1 = D2UDXZ
!                                                             STORE2 = D2WDXZ
!                                                        STORE7 = (DUDZ+DWDX)
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   U-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUXX,Xd
!   U TAUXX,Xd
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqI, "A = var*B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(tthd_ops, 1, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN X: TAUXX,X TERMS ZERO ON END POINTS
    IF(fxlvsn) call zeroxl(d_store3)
    IF(fxrvsn) call zeroxr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqS, "A = A-B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAJ, "A = A-B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                             STORE1 = D2UDXZ
!                                                             STORE2 = D2WDXZ
!                                                        STORE7 = (DUDZ+DWDX)
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   W-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUZZ,Zc
!   W TAUZZ,Zc
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqI, "A = var*B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(tthd_ops, 1, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN Z: TAUZZ,Z TERMS ZERO ON END POINTS
    IF(fzlvsn) call zerozl(d_store3)
    IF(fzrvsn) call zerozr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqS, "A = A-B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAJ, "A = A-B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                             STORE1 = D2UDXZ
!                                                             STORE2 = D2WDXZ
!                                                        STORE7 = (DUDZ+DWDX)
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   D2UDZ2
    call d2fdz2(d_utmp,d_store3)

!   D2UDZ2+D2WDXZ
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!   U-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUXZ,Za,b,c
!   U TAUXZ,Za,b,c

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAE, "A = A*B+C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN Y: TAUXZ,Z TERM ZERO ON END POINTS
    IF(fzlvst) call zerozl(d_store3)
    IF(fzrvst) call zerozr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAA, "A = A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                             STORE1 = D2UDXZ
!                                                        STORE7 = (DUDZ+DWDX)
!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   D2WDX2
    call d2fdx2(d_wtmp,d_store3)

!   D2UDXZ+D2WDX2
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!   W-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUZX,Xa,b,c
!   W TAUZX,Xa,b,c

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAE, "A = A*B+C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN X: TAUZX,X TERM ZERO ON END POINTS
    IF(fxlvst) call zeroxl(d_store3)
    IF(fxrvst) call zeroxr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAA, "A = A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                        STORE7 = (DUDZ+DWDX)
!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUXZ(DUDZ+DWDX)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAC, "A = A+B*C*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   DVDZ
    call dfbydz(d_vtmp,d_store1)

!   COLLECT VELOCITY DERIVATIVES FOR BCs
!   ------------------------------------

!   Z-DIRECTION: DVDZ
    IF(fzlcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(boundary_kernel_velderiv_zdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF
    IF(fzrcnv)THEN
        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_zdir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF

!                                                               STORE1 = DVDZ
!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   V-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   HALF RHO W DVDZ
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAR, "A = A-half*B*C*D", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                               STORE1 = DVDZ
!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   DWDY
    call dfbydy(d_wtmp,d_store2)

!   COLLECT VELOCITY DERIVATIVES FOR BCs
!   ------------------------------------

!   Y-DIRECTION: DWDY
    IF(fylcnv)THEN
        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_ydir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF
    IF(fyrcnv)THEN
        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(boundary_kernel_velderiv_ydir, "COLLECT VELOCITY DERIVATIVES FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF

!                                                               STORE1 = DVDZ
!                                                               STORE2 = DWDY
!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   W-EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   HALF RHO V DWDY
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAR, "A = A-half*B*C*D", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                                STORE1 = DVDZ
!                                                               STORE2 = DWDY
!                                                          TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   VELOCITY SECOND CROSS DERIVATIVES
!   ---------------------------------

!   DVDZ+DWDY
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqR, "A = B+C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!   D2VDYZ
    call d2fdyz(d_vtmp,d_store1)

!   D2WDYZ
    call d2fdyz(d_wtmp,d_store2)

!                                                             STORE1 = D2VDYZ
!                                                             STORE2 = D2WDYZ
!                                                        STORE7 = (DVDZ+DWDY)
!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   V-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUYY,Yd
!   V TAUYY,Yd

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqI, "A = var*B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(tthd_ops, 1, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN Y: TAUYY,2 TERMS ZERO ON END POINTS
    IF(fylvsn) call zeroyl(d_store3)
    IF(fyrvsn) call zeroyr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqS, "A = A-B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAJ, "A = A-B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                             STORE1 = D2VDYZ
!                                                             STORE2 = D2WDYZ
!                                                        STORE7 = (DVDZ+DWDY)
!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   W-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUZZ,Zd
!   W TAUZZ,Zd

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqI, "A = var*B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(tthd_ops, 1, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN Z: TAUZZ,Z TERMS ZERO ON END POINTS
    IF(fzlvsn) call zerozl(d_store3)
    IF(fzrvsn) call zerozr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqS, "A = A-B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAJ, "A = A-B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                             STORE1 = D2VDYZ
!                                                             STORE2 = D2WDYZ
!                                                        STORE7 = (DVDZ+DWDY)
!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   D2VDZ2
    call d2fdz2(d_vtmp,d_store3)

!   D2VDZ2+D2WDYZ
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!   V-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUYZ,Za,b,c
!   V TAUYZ,Za,b,c

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAE, "A = A*B+C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN Z: TAUYZ,Z TERM ZERO ON END POINTS
    IF(fzlvst) call zerozl(d_store3)
    IF(fzrvst) call zerozr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAA, "A = A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                             STORE1 = D2VDYZ
!                                                        STORE7 = (DVDZ+DWDY)
!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   D2WDY2
    call d2fdy2(d_wtmp,d_store3)

!   D2VDYZ+D2WDY2
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!   W-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUZY,Ya,b,c
!   W TAUZY,Ya,b,c

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAE, "A = A*B+C*D", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN Y: TAUZY,Y TERM ZERO ON END POINTS
    IF(fylvst)call zeroyl(d_store3)
    IF(fyrvst)call zeroyr(d_store3)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqQ, "A = A+B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqAA, "A = A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                        STORE7 = (DVDZ+DWDY)
!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUYZ(DVDZ+DWDY)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAC, "A = A+B*C*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                          TRANSP = VISCOSITY
!                                                      STORE4,5,6 = DMUDX,Y,Z
!   =========================================================================

!   VELOCITY SECOND NORMAL DERIVATIVE TERMS
!    ---------------------------------------
!   D2UDX2,D2VDY2,D2WDZ2
    call d2fdx2(d_utmp,d_store1)
    call d2fdy2(d_vtmp,d_store2)
    call d2fdz2(d_wtmp,d_store3)

!   BOUNDARY CONDITIONS
!   BC IN X: TAUXX,Xa TERM ZERO ON END POINTS
    IF(fxlvsn) call zeroxl(d_store1)
    IF(fxrvsn) call zeroxr(d_store1)

!   BC IN Y: TAUYY,Ya TERM ZERO ON END POINTS
    IF(fylvsn) call zeroyl(d_store2)
    IF(fyrvsn) call zeroyr(d_store2)

!   BC IN Z: TAUZZ,Za TERM ZERO ON END POINTS
    IF(fzlvsn) call zerozl(d_store3)
    IF(fzrvsn) call zerozr(d_store3)

!   U-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUXX,Xa
!   U TAUXX,Xa
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqtau, "TAU", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(fthd_ops, 1, "real(kind=8)", OPS_READ))

!   V-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUYY,Ya
!   V TAUYY,Ya

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqtau, "TAU", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(fthd_ops, 1, "real(kind=8)", OPS_READ))

!   W-EQUATION: VISCOUS STRESS TERMS
!   E-EQUATION: VISCOUS WORK TERMS
!   ------------------------------
!   TAUZZ,Za
!   W TAUZZ,Za

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqtau, "TAU", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(fthd_ops, 1, "real(kind=8)", OPS_READ))

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   CONTINUITY EQUATION
!   -------------------
!   DIV RHO U

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqC, "A = -B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_divm, 1, s3d_000, "real(kind=8)", OPS_READ))

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

END SUBROUTINE rhsvel
