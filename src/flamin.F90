SUBROUTINE flamin

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga 

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

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   PARAMETERS
!   ==========
!   ESTIMATED FLAME LOCATION AND THICKNESS
    real(kind=8) :: clocat,cthick
    PARAMETER(clocat = 0.0025_dp, cthick = 0.0005_dp)

!   FUNCTION
!   ========
    real(kind=8) :: erfunc
    EXTERNAL erfunc

!   LOCAL DATA
!   ==========
    real(kind=8) :: trinr,u0
    real(kind=8) :: deltagx,deltagy,deltagz
    real(kind=8) :: rglocl
    real(kind=8) :: angfrx,angfry,angfrz
    integer :: jspec
    integer :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   SPECIFY INITIAL THERMOCHEMICAL FIELD HERE
!   =========================================

!   SET PRODUCT TEMPERATURE
!   -----------------------
!   REACTANT TEMPERATURE SET IN CONTROL FILE
    trinr = trin
    u0 = 34.789806_dp

!   Mixture gas constant
    rglocl = zero
    do jspec = 1,nspec
        rglocl = rglocl + rgspec(jspec)*yrin(jspec)
    end do

!   Times (constant) density
    rglocl  = drin*rglocl

!   Global indexing
!   ---------------
    deltagx = xgdlen/(real(nxglbl-1,8))
    deltagy = ygdlen/(real(nyglbl-1,8))
    deltagz = zgdlen/(real(nzglbl-1,8))

!   SET THE VELOCITY PROFILE FOR TGV
!   --------------------------------
    angfrx = 8.0_8*atan(1.0_8)/xgdlen
    angfry = 8.0_8*atan(1.0_8)/ygdlen
    angfrz = 8.0_8*atan(1.0_8)/zgdlen

    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    call ops_par_loop(flamin_kernel_set_velocity_tgv, "SET THE VELOCITY PROFILE FOR TGV",  senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(8)", OPS_WRITE), &
                    ops_arg_dat(d_vrun, 1, s3d_000, "real(8)", OPS_WRITE), &
                    ops_arg_dat(d_wrun, 1, s3d_000, "real(8)", OPS_WRITE), &
                    ops_arg_dat(d_prun, 1, s3d_000, "real(8)", OPS_WRITE), &
                    ops_arg_gbl(prin, 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(drin, 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(u0, 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(deltagx, 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(deltagy, 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(deltagz, 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(angfrx, 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(angfry, 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(angfrz, 1, "real(8)", OPS_READ), &
                    ops_arg_idx())

!   set temperature profile assuming constant density
    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    call ops_par_loop(flamin_kernel_eqA, "A = B/var", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(8)", OPS_WRITE), &
                    ops_arg_dat(d_prun, 1, s3d_000, "real(8)", OPS_READ), &
                    ops_arg_gbl(rglocl, 1, "real(8)", OPS_READ))

!   set constant density
    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    call ops_par_loop(math_kernel_eqAN, "set constant density",  senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(8)", OPS_WRITE), &
                    ops_arg_gbl(drin, 1, "real(8)", OPS_READ))

9000  FORMAT(a)

END SUBROUTINE flamin
