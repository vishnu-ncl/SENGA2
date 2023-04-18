!out_arr = s3d_000,                in_arr = s3d_000
SUBROUTINE copy_kernel(out_arr, in_arr)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel

!out_arr = s3d_000_strid3d_yz,      in_arr = s3d_000
SUBROUTINE copy_kernel_xdir(out_arr, in_arr)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel_xdir

!out_arr = s3d_000_strid3d_xz,      in_arr = s3d_000
SUBROUTINE copy_kernel_ydir(out_arr, in_arr)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel_ydir

!out_arr = s3d_000_strid3d_xy,      in_arr = s3d_000
SUBROUTINE copy_kernel_zdir(out_arr, in_arr)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel_zdir
