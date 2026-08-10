!out_arr = s3d_000,                in_arr = s3d_000
SUBROUTINE copy_kernel(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel

!A = B_multidim(ispec)
SUBROUTINE copy_kernel_mdim_to_sdim(out_arr, in_arr, ispec)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    integer(kind=4), intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE copy_kernel_mdim_to_sdim

!A_multidim(ispec) = B
SUBROUTINE copy_kernel_sdim_to_mdim(out_arr, in_arr, ispec)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    integer(kind=4), intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel_sdim_to_mdim

SUBROUTINE copy_kernel_mdim_to_sdim_int(out_arr, in_arr, ispec)

    integer(kind=4), dimension(1) :: out_arr
    integer(kind=4), dimension(1), intent(in) :: in_arr
    integer(kind=4), intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE copy_kernel_mdim_to_sdim_int

!A_multidim(ispec) = B
SUBROUTINE copy_kernel_sdim_to_mdim_int(out_arr, in_arr, ispec)

    integer(kind=4), dimension(1) :: out_arr
    integer(kind=4), dimension(1), intent(in) :: in_arr
    integer(kind=4), intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel_sdim_to_mdim_int

!out_arr = s3d_000_strid3d_yz,      in_arr = s3d_000
SUBROUTINE copy_kernel_xdir(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel_xdir

!out_arr = s3d_000_strid3d_yz,      in_arr = s3d_000_strid3d_yz
SUBROUTINE copy_kernel_xxdir(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel_xxdir

!out_arr = s3d_000_strid3d_xz,      in_arr = s3d_000
SUBROUTINE copy_kernel_ydir(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel_ydir

!out_arr = s3d_000_strid3d_xz,      in_arr = s3d_000_strid3d_xz
SUBROUTINE copy_kernel_yydir(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel_yydir

!out_arr = s3d_000_strid3d_xy,      in_arr = s3d_000
SUBROUTINE copy_kernel_zdir(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel_zdir

!out_arr = s3d_000_strid3d_xy,      in_arr = s3d_000_strid3d_xy
SUBROUTINE copy_kernel_zzdir(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel_zzdir
