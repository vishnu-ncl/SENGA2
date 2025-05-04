SUBROUTINE ops_data_init()
    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

    integer(kind=4) :: d_size(3)
    integer(kind=4) :: d_base(3) = [1,1,1] !array indexing - start from 1
    integer(kind=4) :: d_p(3) !max boundary depths for the dat in the possitive direction
    integer(kind=4) :: d_m(3) !max boundary depths for the dat in the negative direction

    real(kind=8), dimension(:), allocatable :: temp_real_null
    integer(kind=4), dimension(:), allocatable :: temp_int_null

    integer(kind=4) :: ispec,jspec,iindex,line,status,ic
    integer(kind=4) :: rangexyz(6)

    character(len=20) :: buf

    integer(kind=4) :: halo_idx
    integer(kind=4) :: iter_size(3), base_from(3), base_to(3), dir_from(3), dir_to(3)

    integer(kind=4) :: a3d_000(3) = [0,0,0]

    integer(kind=4) :: stride3d_x(3) = [1,0,0]
    integer(kind=4) :: stride3d_y(3) = [0,1,0]
    integer(kind=4) :: stride3d_z(3) = [0,0,1]

    integer(kind=4) :: stride3d_xy(3) = [1,1,0]
    integer(kind=4) :: stride3d_xz(3) = [1,0,1]
    integer(kind=4) :: stride3d_yz(3) = [0,1,1]

    character(len=10) :: pncont
    character(len=4) :: pnxdat
    parameter(pncont = 'input/cont', pnxdat = '.dat')

    character(len=60) :: fname
    character(len=3) :: pnxhdf
    parameter(pnxhdf = '.h5')
    integer(kind=4) :: dtime
    character(len=8) :: citime

!------------------------------------------------------------------------------------------------------------

    integer(kind=4) :: a3d_000_to_p400_x(15) = [0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0]
    integer(kind=4) :: a3d_000_to_m400_x(15) = [0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0]

    integer(kind=4) :: a3d_p100_to_p400_x(12) = [1,0,0, 2,0,0, 3,0,0, 4,0,0]
    integer(kind=4) :: a3d_m100_to_m400_x(12) = [-1,0,0, -2,0,0, -3,0,0, -4,0,0]

    integer(kind=4) :: a3d_p500_to_m500_x(33) = [5,0,0, 4,0,0, 3,0,0, 2,0,0, 1,0,0, 0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0, -5,0,0]

!------------------------------------------------------------------------------------------------------------
    integer(kind=4) :: a3d_000_to_p040_y(15) = [0,0,0, 0,1,0, 0,2,0, 0,3,0, 0,4,0]
    integer(kind=4) :: a3d_000_to_m040_y(15) = [0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0]

    integer(kind=4) :: a3d_p010_to_p040_y(12) = [0,1,0, 0,2,0, 0,3,0, 0,4,0]
    integer(kind=4) :: a3d_m010_to_m040_y(12) = [0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0]

    integer(kind=4) :: a3d_p050_to_m050_y(33) = [0,5,0, 0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, 0,-5,0]

!------------------------------------------------------------------------------------------------------------
    integer(kind=4) :: a3d_000_to_p004_z(15) = [0,0,0, 0,0,1, 0,0,2, 0,0,3, 0,0,4]
    integer(kind=4) :: a3d_000_to_m004_z(15) = [0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4]

    integer(kind=4) :: a3d_p001_to_p004_z(12) = [0,0,1, 0,0,2, 0,0,3, 0,0,4]
    integer(kind=4) :: a3d_m001_to_m004_z(12) = [0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4]

    integer(kind=4) :: a3d_p005_to_m005_z(33) = [0,0,5, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, 0,0,-5]

!------------------------------------------------------------------------------------------------------------

    integer(kind=4) :: a3d_p000_m440_mixed_xy(39) = [0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, -1,0,0, -1,-1,0, -2,0,0, -2,-2,0, -3,0,0, -3,-3,0, -4,0,0, -4,-4,0]
    integer(kind=4) :: a3d_p010_m430_mixed_xy(75) = [0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, -1,1,0, -1,0,0, -1,-1,0, -1,-2,0, -1,-3,0, -2,1,0, -2,0,0, -2,-1,0, -2,-2,0, -2,-3,0, -3,1,0, -3,0,0, -3,-1,0, -3,-2,0, -3,-3,0, -4,1,0, -4,0,0, -4,-1,0, -4,-2,0, -4,-3,0]
    integer(kind=4) :: a3d_p020_m420_mixed_xy(60) = [0,2,0, 0,1,0, 0,-1,0, 0,-2,0, -1,2,0, -1,1,0, -1,-1,0, -1,-2,0, -2,2,0, -2,1,0, -2,-1,0, -2,-2,0, -3,2,0, -3,1,0, -3,-1,0, -3,-2,0, -4,2,0, -4,1,0, -4,-1,0, -4,-2,0]
    integer(kind=4) :: a3d_p030_m410_mixed_xy(75) = [0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, -1,3,0, -1,2,0, -1,1,0, -1,0,0, -1,-1,0, -2,3,0, -2,2,0, -2,1,0, -2,0,0, -2,-1,0, -3,3,0, -3,2,0, -3,1,0, -3,0,0, -3,-1,0, -4,3,0, -4,2,0, -4,1,0, -4,0,0, -4,-1,0]
    integer(kind=4) :: a3d_p040_m400_mixed_xy(39) = [0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, -1,1,0, -1,0,0, -2,2,0, -2,0,0, -3,3,0, -3,0,0, -4,4,0, -4,0,0]

    integer(kind=4) :: a3d_p100_m300_mixed_xy(39) = [1,0,0, 1,-1,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, -1,1,0, -1,0,0, -2,2,0, -2,0,0, -3,3,0, -3,0,0]
    integer(kind=4) :: a3d_p100_m340_mixed_xy(75) = [1,0,0, 1,-1,0, 1,-2,0, 1,-3,0, 1,-4,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, -1,0,0, -1,-1,0, -1,-2,0, -1,-3,0, -1,-4,0, -2,0,0, -2,-1,0, -2,-2,0, -2,-3,0, -2,-4,0, -3,0,0, -3,-1,0, -3,-2,0, -3,-3,0, -3,-4,0]
    integer(kind=4) :: a3d_p110_m330_mixed_xy(39) = [1,1,0, 1,0,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, -1,0,0, -1,-1,0, -2,0,0, -2,-2,0, -3,0,0, -3,-3,0]
    integer(kind=4) :: a3d_p120_m320_mixed_xy(60) = [1,2,0, 1,1,0, 1,-1,0, 1,-2,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0, -1,2,0, -1,1,0, -1,-1,0, -1,-2,0, -2,2,0, -2,1,0, -2,-1,0, -2,-2,0, -3,2,0, -3,1,0, -3,-1,0, -3,-2,0]
    integer(kind=4) :: a3d_p140_m300_mixed_xy(75) = [1,4,0, 1,3,0, 1,2,0, 1,1,0, 1,0,0, 0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, -1,4,0, -1,3,0, -1,2,0, -1,1,0, -1,0,0, -2,4,0, -2,3,0, -2,2,0, -2,1,0, -2,0,0, -3,4,0, -3,3,0, -3,2,0, -3,1,0, -3,0,0]

    integer(kind=4) :: a3d_p200_m240_mixed_xy(60) = [2,0,0, 2,-1,0, 2,-2,0, 2,-3,0, 2,-4,0, 1,0,0, 1,-1,0, 1,-2,0, 1,-3,0, 1,-4,0, -1,0,0, -1,-1,0, -1,-2,0, -1,-3,0, -1,-4,0, -2,0,0, -2,-1,0, -2,-2,0, -2,-3,0, -2,-4,0]
    integer(kind=4) :: a3d_p210_m230_mixed_xy(60) = [2,1,0, 2,0,0, 2,-1,0, 2,-2,0, 2,-3,0, 1,1,0, 1,0,0, 1,-1,0, 1,-2,0, 1,-3,0, -1,1,0, -1,0,0, -1,-1,0, -1,-2,0, -1,-3,0, -2,1,0, -2,0,0, -2,-1,0, -2,-2,0, -2,-3,0]
    integer(kind=4) :: a3d_p220_m220_mixed_xy(24) = [2,2,0, 2,-2,0, 1,1,0, 1,-1,0, -1,1,0, -1,-1,0, -2,2,0, -2,-2,0]
    integer(kind=4) :: a3d_p230_m210_mixed_xy(60) = [2,3,0, 2,2,0, 2,1,0, 2,0,0, 2,-1,0, 1,3,0, 1,2,0, 1,1,0, 1,0,0, 1,-1,0, -1,3,0, -1,2,0, -1,1,0, -1,0,0, -1,-1,0, -2,3,0, -2,2,0, -2,1,0, -2,0,0, -2,-1,0]
    integer(kind=4) :: a3d_p240_m200_mixed_xy(60) = [2,4,0, 2,3,0, 2,2,0, 2,1,0, 2,0,0, 1,4,0, 1,3,0, 1,2,0, 1,1,0, 1,0,0, -1,4,0, -1,3,0, -1,2,0, -1,1,0, -1,0,0, -2,4,0, -2,3,0, -2,2,0, -2,1,0, -2,0,0]

    integer(kind=4) :: a3d_p300_m100_mixed_xy(39) = [3,0,0, 3,-3,0, 2,0,0, 2,-2,0, 1,0,0, 1,-1,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, -1,1,0, -1,0,0]
    integer(kind=4) :: a3d_p300_m140_mixed_xy(75) = [3,0,0, 3,-1,0, 3,-2,0, 3,-3,0, 3,-4,0, 2,0,0, 2,-1,0, 2,-2,0, 2,-3,0, 2,-4,0, 1,0,0, 1,-1,0, 1,-2,0, 1,-3,0, 1,-4,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, -1,0,0, -1,-1,0, -1,-2,0, -1,-3,0, -1,-4,0]
    integer(kind=4) :: a3d_p320_m120_mixed_xy(60) = [3,2,0, 3,1,0, 3,-1,0, 3,-2,0, 2,2,0, 2,1,0, 2,-1,0, 2,-2,0, 1,2,0, 1,1,0, 1,-1,0, 1,-2,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0, -1,2,0, -1,1,0, -1,-1,0, -1,-2,0]
    integer(kind=4) :: a3d_p330_m110_mixed_xy(39) = [3,3,0, 3,0,0, 2,2,0, 2,0,0, 1,1,0, 1,0,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, -1,0,0, -1,-1,0]
    integer(kind=4) :: a3d_p330_m330_mixed_xy(36) = [3,3,0, 3,-3,0, 2,2,0, 2,-2,0, 1,1,0, 1,-1,0, -1,1,0, -1,-1,0, -2,2,0, -2,-2,0, -3,3,0, -3,-3,0]
    integer(kind=4) :: a3d_p340_m100_mixed_xy(75) = [3,4,0, 3,3,0, 3,2,0, 3,1,0, 3,0,0, 2,4,0, 2,3,0, 2,2,0, 2,1,0, 2,0,0, 1,4,0, 1,3,0, 1,2,0, 1,1,0, 1,0,0, 0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, -1,4,0, -1,3,0, -1,2,0, -1,1,0, -1,0,0]

    integer(kind=4) :: a3d_p400_p040_mixed_xy(39) = [4,0,0, 4,-4,0, 3,0,0, 3,-3,0, 2,0,0, 2,-2,0, 1,0,0, 1,-1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0]
    integer(kind=4) :: a3d_p410_p030_mixed_xy(75) = [4,1,0, 4,0,0, 4,-1,0, 4,-2,0, 4,-3,0, 3,1,0, 3,0,0, 3,-1,0, 3,-2,0, 3,-3,0, 2,1,0, 2,0,0, 2,-1,0, 2,-2,0, 2,-3,0, 1,1,0, 1,0,0, 1,-1,0, 1,-2,0, 1,-3,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0]
    integer(kind=4) :: a3d_p420_m020_mixed_xy(60) = [4,2,0, 4,1,0, 4,-1,0, 4,-2,0, 3,2,0, 3,1,0, 3,-1,0, 3,-2,0, 2,2,0, 2,1,0, 2,-1,0, 2,-2,0, 1,2,0, 1,1,0, 1,-1,0, 1,-2,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0]
    integer(kind=4) :: a3d_p430_m010_mixed_xy(75) = [4,3,0, 4,2,0, 4,1,0, 4,0,0, 4,-1,0, 3,3,0, 3,2,0, 3,1,0, 3,0,0, 3,-1,0, 2,3,0, 2,2,0, 2,1,0, 2,0,0, 2,-1,0, 1,3,0, 1,2,0, 1,1,0, 1,0,0, 1,-1,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0]
    integer(kind=4) :: a3d_p440_p000_mixed_xy(39) = [4,4,0, 4,0,0, 3,3,0, 3,0,0, 2,2,0, 2,0,0, 1,1,0, 1,0,0, 0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0]
    integer(kind=4) :: a3d_p440_m440_mixed_xy(48) = [4,4,0, 4,-4,0, 3,3,0, 3,-3,0, 2,2,0, 2,-2,0, 1,1,0, 1,-1,0, -1,1,0, -1,-1,0, -2,2,0, -2,-2,0, -3,3,0, -3,-3,0, -4,4,0, -4,-4,0]

    integer(kind=4) :: a3d_p550_m550_mixed_xy(60) = [5,5,0, 5,-5,0, 4,4,0, 4,-4,0, 3,3,0, 3,-3,0, 2,2,0, 2,-2,0, 1,1,0, 1,-1,0, -1,1,0, -1,-1,0, -2,2,0, -2,-2,0, -3,3,0, -3,-3,0, -4,4,0, -4,-4,0, -5,5,0, -5,-5,0]

!------------------------------------------------------------------------------------------------------------

    integer(kind=4) :: a3d_p000_m404_mixed_xz(39) = [0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, -1,0,0, -1,0,-1, -2,0,0, -2,0,-2, -3,0,0, -3,0,-3, -4,0,0, -4,0,-4]
    integer(kind=4) :: a3d_p001_m403_mixed_xz(75) = [0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, -1,0,1, -1,0,0, -1,0,-1, -1,0,-2, -1,0,-3, -2,0,1, -2,0,0, -2,0,-1, -2,0,-2, -2,0,-3, -3,0,1, -3,0,0, -3,0,-1, -3,0,-2, -3,0,-3, -4,0,1, -4,0,0, -4,0,-1, -4,0,-2, -4,0,-3]
    integer(kind=4) :: a3d_p002_m402_mixed_xz(60) = [0,0,2, 0,0,1, 0,0,-1, 0,0,-2, -1,0,2, -1,0,1, -1,0,-1, -1,0,-2, -2,0,2, -2,0,1, -2,0,-1, -2,0,-2, -3,0,2, -3,0,1, -3,0,-1, -3,0,-2, -4,0,2, -4,0,1, -4,0,-1, -4,0,-2]
    integer(kind=4) :: a3d_p003_m401_mixed_xz(75) = [0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, -1,0,3, -1,0,2, -1,0,1, -1,0,0, -1,0,-1, -2,0,3, -2,0,2, -2,0,1, -2,0,0, -2,0,-1, -3,0,3, -3,0,2, -3,0,1, -3,0,0, -3,0,-1, -4,0,3, -4,0,2, -4,0,1, -4,0,0, -4,0,-1]
    integer(kind=4) :: a3d_p004_m400_mixed_xz(39) = [0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, -1,0,1, -1,0,0, -2,0,2, -2,0,0, -3,0,3, -3,0,0, -4,0,4, -4,0,0]

    integer(kind=4) :: a3d_p100_m300_mixed_xz(39) = [1,0,0, 1,0,-1, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, -1,0,1, -1,0,0, -2,0,2, -2,0,0, -3,0,3, -3,0,0]
    integer(kind=4) :: a3d_p100_m304_mixed_xz(75) = [1,0,0, 1,0,-1, 1,0,-2, 1,0,-3, 1,0,-4, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, -1,0,0, -1,0,-1, -1,0,-2, -1,0,-3, -1,0,-4, -2,0,0, -2,0,-1, -2,0,-2, -2,0,-3, -2,0,-4, -3,0,0, -3,0,-1, -3,0,-2, -3,0,-3, -3,0,-4]
    integer(kind=4) :: a3d_p101_m303_mixed_xz(39) = [1,0,1, 1,0,0, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, -1,0,0, -1,0,-1, -2,0,0, -2,0,-2, -3,0,0, -3,0,-3]
    integer(kind=4) :: a3d_p102_m302_mixed_xz(60) = [1,0,2, 1,0,1, 1,0,-1, 1,0,-2, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2, -1,0,2, -1,0,1, -1,0,-1, -1,0,-2, -2,0,2, -2,0,1, -2,0,-1, -2,0,-2, -3,0,2, -3,0,1, -3,0,-1, -3,0,-2]
    integer(kind=4) :: a3d_p104_m300_mixed_xz(75) = [1,0,4, 1,0,3, 1,0,2, 1,0,1, 1,0,0, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, -1,0,4, -1,0,3, -1,0,2, -1,0,1, -1,0,0, -2,0,4, -2,0,3, -2,0,2, -2,0,1, -2,0,0, -3,0,4, -3,0,3, -3,0,2, -3,0,1, -3,0,0]

    integer(kind=4) :: a3d_p200_m204_mixed_xz(60) = [2,0,0, 2,0,-1, 2,0,-2, 2,0,-3, 2,0,-4, 1,0,0, 1,0,-1, 1,0,-2, 1,0,-3, 1,0,-4, -1,0,0, -1,0,-1, -1,0,-2, -1,0,-3, -1,0,-4, -2,0,0, -2,0,-1, -2,0,-2, -2,0,-3, -2,0,-4]
    integer(kind=4) :: a3d_p201_m203_mixed_xz(60) = [2,0,1, 2,0,0, 2,0,-1, 2,0,-2, 2,0,-3, 1,0,1, 1,0,0, 1,0,-1, 1,0,-2, 1,0,-3, -1,0,1, -1,0,0, -1,0,-1, -1,0,-2, -1,0,-3, -2,0,1, -2,0,0, -2,0,-1, -2,0,-2, -2,0,-3]
    integer(kind=4) :: a3d_p202_m202_mixed_xz(24) = [2,0,2, 2,0,-2, 1,0,1, 1,0,-1, -1,0,1, -1,0,-1, -2,0,2, -2,0,-2]
    integer(kind=4) :: a3d_p203_m201_mixed_xz(60) = [2,0,3, 2,0,2, 2,0,1, 2,0,0, 2,0,-1, 1,0,3, 1,0,2, 1,0,1, 1,0,0, 1,0,-1, -1,0,3, -1,0,2, -1,0,1, -1,0,0, -1,0,-1, -2,0,3, -2,0,2, -2,0,1, -2,0,0, -2,0,-1]
    integer(kind=4) :: a3d_p204_m200_mixed_xz(60) = [2,0,4, 2,0,3, 2,0,2, 2,0,1, 2,0,0, 1,0,4, 1,0,3, 1,0,2, 1,0,1, 1,0,0, -1,0,4, -1,0,3, -1,0,2, -1,0,1, -1,0,0, -2,0,4, -2,0,3, -2,0,2, -2,0,1, -2,0,0]

    integer(kind=4) :: a3d_p300_m100_mixed_xz(39) = [3,0,0, 3,0,-3, 2,0,0, 2,0,-2, 1,0,0, 1,0,-1, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, -1,0,1, -1,0,0]
    integer(kind=4) :: a3d_p300_m104_mixed_xz(75) = [3,0,0, 3,0,-1, 3,0,-2, 3,0,-3, 3,0,-4, 2,0,0, 2,0,-1, 2,0,-2, 2,0,-3, 2,0,-4, 1,0,0, 1,0,-1, 1,0,-2, 1,0,-3, 1,0,-4, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, -1,0,0, -1,0,-1, -1,0,-2, -1,0,-3, -1,0,-4]
    integer(kind=4) :: a3d_p302_m102_mixed_xz(60) = [3,0,2, 3,0,1, 3,0,-1, 3,0,-2, 2,0,2, 2,0,1, 2,0,-1, 2,0,-2, 1,0,2, 1,0,1, 1,0,-1, 1,0,-2, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2, -1,0,2, -1,0,1, -1,0,-1, -1,0,-2]
    integer(kind=4) :: a3d_p303_m101_mixed_xz(39) = [3,0,3, 3,0,0, 2,0,2, 2,0,0, 1,0,1, 1,0,0, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, -1,0,0, -1,0,-1]
    integer(kind=4) :: a3d_p303_m303_mixed_xz(36) = [3,0,3, 3,0,-3, 2,0,2, 2,0,-2, 1,0,1, 1,0,-1, -1,0,1, -1,0,-1, -2,0,2, -2,0,-2, -3,0,3, -3,0,-3]
    integer(kind=4) :: a3d_p304_m100_mixed_xz(75) = [3,0,4, 3,0,3, 3,0,2, 3,0,1, 3,0,0, 2,0,4, 2,0,3, 2,0,2, 2,0,1, 2,0,0, 1,0,4, 1,0,3, 1,0,2, 1,0,1, 1,0,0, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, -1,0,4, -1,0,3, -1,0,2, -1,0,1, -1,0,0]

    integer(kind=4) :: a3d_p400_p004_mixed_xz(39) = [4,0,0, 4,0,-4, 3,0,0, 3,0,-3, 2,0,0, 2,0,-2, 1,0,0, 1,0,-1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4]
    integer(kind=4) :: a3d_p401_p003_mixed_xz(75) = [4,0,1, 4,0,0, 4,0,-1, 4,0,-2, 4,0,-3, 3,0,1, 3,0,0, 3,0,-1, 3,0,-2, 3,0,-3, 2,0,1, 2,0,0, 2,0,-1, 2,0,-2, 2,0,-3, 1,0,1, 1,0,0, 1,0,-1, 1,0,-2, 1,0,-3, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3]
    integer(kind=4) :: a3d_p402_m002_mixed_xz(60) = [4,0,2, 4,0,1, 4,0,-1, 4,0,-2, 3,0,2, 3,0,1, 3,0,-1, 3,0,-2, 2,0,2, 2,0,1, 2,0,-1, 2,0,-2, 1,0,2, 1,0,1, 1,0,-1, 1,0,-2, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2]
    integer(kind=4) :: a3d_p403_m001_mixed_xz(75) = [4,0,3, 4,0,2, 4,0,1, 4,0,0, 4,0,-1, 3,0,3, 3,0,2, 3,0,1, 3,0,0, 3,0,-1, 2,0,3, 2,0,2, 2,0,1, 2,0,0, 2,0,-1, 1,0,3, 1,0,2, 1,0,1, 1,0,0, 1,0,-1, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1]
    integer(kind=4) :: a3d_p404_p000_mixed_xz(39) = [4,0,4, 4,0,0, 3,0,3, 3,0,0, 2,0,2, 2,0,0, 1,0,1, 1,0,0, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0]
    integer(kind=4) :: a3d_p404_m404_mixed_xz(48) = [4,0,4, 4,0,-4, 3,0,3, 3,0,-3, 2,0,2, 2,0,-2, 1,0,1, 1,0,-1, -1,0,1, -1,0,-1, -2,0,2, -2,0,-2, -3,0,3, -3,0,-3, -4,0,4, -4,0,-4]

    integer(kind=4) :: a3d_p505_m505_mixed_xz(60) = [5,0,5, 5,0,-5, 4,0,4, 4,0,-4, 3,0,3, 3,0,-3, 2,0,2, 2,0,-2, 1,0,1, 1,0,-1, -1,0,1, -1,0,-1, -2,0,2, -2,0,-2, -3,0,3, -3,0,-3, -4,0,4, -4,0,-4, -5,0,5, -5,0,-5]

!------------------------------------------------------------------------------------------------------------

    integer(kind=4) :: a3d_p000_m044_mixed_yz(39) = [0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, 0,-1,0, 0,-1,-1, 0,-2,0, 0,-2,-2, 0,-3,0, 0,-3,-3, 0,-4,0, 0,-4,-4]
    integer(kind=4) :: a3d_p001_m043_mixed_yz(75) = [0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,-1,1, 0,-1,0, 0,-1,-1, 0,-1,-2, 0,-1,-3, 0,-2,1, 0,-2,0, 0,-2,-1, 0,-2,-2, 0,-2,-3, 0,-3,1, 0,-3,0, 0,-3,-1, 0,-3,-2, 0,-3,-3, 0,-4,1, 0,-4,0, 0,-4,-1, 0,-4,-2, 0,-4,-3]
    integer(kind=4) :: a3d_p002_m042_mixed_yz(60) = [0,0,2, 0,0,1, 0,0,-1, 0,0,-2, 0,-1,2, 0,-1,1, 0,-1,-1, 0,-1,-2, 0,-2,2, 0,-2,1, 0,-2,-1, 0,-2,-2, 0,-3,2, 0,-3,1, 0,-3,-1, 0,-3,-2, 0,-4,2, 0,-4,1, 0,-4,-1, 0,-4,-2]
    integer(kind=4) :: a3d_p003_m041_mixed_yz(75) = [0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,-1,3, 0,-1,2, 0,-1,1, 0,-1,0, 0,-1,-1, 0,-2,3, 0,-2,2, 0,-2,1, 0,-2,0, 0,-2,-1, 0,-3,3, 0,-3,2, 0,-3,1, 0,-3,0, 0,-3,-1, 0,-4,3, 0,-4,2, 0,-4,1, 0,-4,0, 0,-4,-1]
    integer(kind=4) :: a3d_p004_m040_mixed_yz(39) = [0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,-1,1, 0,-1,0, 0,-2,2, 0,-2,0, 0,-3,3, 0,-3,0, 0,-4,4, 0,-4,0]

    integer(kind=4) :: a3d_p010_m030_mixed_yz(39) = [0,1,0, 0,1,-1, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,-1,1, 0,-1,0, 0,-2,2, 0,-2,0, 0,-3,3, 0,-3,0]
    integer(kind=4) :: a3d_p010_m034_mixed_yz(75) = [0,1,0, 0,1,-1, 0,1,-2, 0,1,-3, 0,1,-4, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, 0,-1,0, 0,-1,-1, 0,-1,-2, 0,-1,-3, 0,-1,-4, 0,-2,0, 0,-2,-1, 0,-2,-2, 0,-2,-3, 0,-2,-4, 0,-3,0, 0,-3,-1, 0,-3,-2, 0,-3,-3, 0,-3,-4]
    integer(kind=4) :: a3d_p011_m033_mixed_yz(39) = [0,1,1, 0,1,0, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,-1,0, 0,-1,-1, 0,-2,0, 0,-2,-2, 0,-3,0, 0,-3,-3]
    integer(kind=4) :: a3d_p012_m032_mixed_yz(60) = [0,1,2, 0,1,1, 0,1,-1, 0,1,-2, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2, 0,-1,2, 0,-1,1, 0,-1,-1, 0,-1,-2, 0,-2,2, 0,-2,1, 0,-2,-1, 0,-2,-2, 0,-3,2, 0,-3,1, 0,-3,-1, 0,-3,-2]
    integer(kind=4) :: a3d_p014_m030_mixed_yz(75) = [0,1,4, 0,1,3, 0,1,2, 0,1,1, 0,1,0, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,-1,4, 0,-1,3, 0,-1,2, 0,-1,1, 0,-1,0, 0,-2,4, 0,-2,3, 0,-2,2, 0,-2,1, 0,-2,0, 0,-3,4, 0,-3,3, 0,-3,2, 0,-3,1, 0,-3,0]

    integer(kind=4) :: a3d_p020_m024_mixed_yz(60) = [0,2,0, 0,2,-1, 0,2,-2, 0,2,-3, 0,2,-4, 0,1,0, 0,1,-1, 0,1,-2, 0,1,-3, 0,1,-4, 0,-1,0, 0,-1,-1, 0,-1,-2, 0,-1,-3, 0,-1,-4, 0,-2,0, 0,-2,-1, 0,-2,-2, 0,-2,-3, 0,-2,-4]
    integer(kind=4) :: a3d_p021_m023_mixed_yz(60) = [0,2,1, 0,2,0, 0,2,-1, 0,2,-2, 0,2,-3, 0,1,1, 0,1,0, 0,1,-1, 0,1,-2, 0,1,-3, 0,-1,1, 0,-1,0, 0,-1,-1, 0,-1,-2, 0,-1,-3, 0,-2,1, 0,-2,0, 0,-2,-1, 0,-2,-2, 0,-2,-3]
    integer(kind=4) :: a3d_p022_m022_mixed_yz(24) = [0,2,2, 0,2,-2, 0,1,1, 0,1,-1, 0,-1,1, 0,-1,-1, 0,-2,2, 0,-2,-2]
    integer(kind=4) :: a3d_p023_m021_mixed_yz(60) = [0,2,3, 0,2,2, 0,2,1, 0,2,0, 0,2,-1, 0,1,3, 0,1,2, 0,1,1, 0,1,0, 0,1,-1, 0,-1,3, 0,-1,2, 0,-1,1, 0,-1,0, 0,-1,-1, 0,-2,3, 0,-2,2, 0,-2,1, 0,-2,0, 0,-2,-1]
    integer(kind=4) :: a3d_p024_m020_mixed_yz(60) = [0,2,4, 0,2,3, 0,2,2, 0,2,1, 0,2,0, 0,1,4, 0,1,3, 0,1,2, 0,1,1, 0,1,0, 0,-1,4, 0,-1,3, 0,-1,2, 0,-1,1, 0,-1,0, 0,-2,4, 0,-2,3, 0,-2,2, 0,-2,1, 0,-2,0]

    integer(kind=4) :: a3d_p030_m010_mixed_yz(39) = [0,3,0, 0,3,-3, 0,2,0, 0,2,-2, 0,1,0, 0,1,-1, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,-1,1, 0,-1,0]
    integer(kind=4) :: a3d_p030_m014_mixed_yz(75) = [0,3,0, 0,3,-1, 0,3,-2, 0,3,-3, 0,3,-4, 0,2,0, 0,2,-1, 0,2,-2, 0,2,-3, 0,2,-4, 0,1,0, 0,1,-1, 0,1,-2, 0,1,-3, 0,1,-4, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, 0,-1,0, 0,-1,-1, 0,-1,-2, 0,-1,-3, 0,-1,-4]
    integer(kind=4) :: a3d_p032_m012_mixed_yz(60) = [0,3,2, 0,3,1, 0,3,-1, 0,3,-2, 0,2,2, 0,2,1, 0,2,-1, 0,2,-2, 0,1,2, 0,1,1, 0,1,-1, 0,1,-2, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2, 0,-1,2, 0,-1,1, 0,-1,-1, 0,-1,-2]
    integer(kind=4) :: a3d_p033_m011_mixed_yz(39) = [0,3,3, 0,3,0, 0,2,2, 0,2,0, 0,1,1, 0,1,0, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,-1,0, 0,-1,-1]
    integer(kind=4) :: a3d_p033_m033_mixed_yz(36) = [0,3,3, 0,3,-3, 0,2,2, 0,2,-2, 0,1,1, 0,1,-1, 0,-1,1, 0,-1,-1, 0,-2,2, 0,-2,-2, 0,-3,3, 0,-3,-3]
    integer(kind=4) :: a3d_p034_m010_mixed_yz(75) = [0,3,4, 0,3,3, 0,3,2, 0,3,1, 0,3,0, 0,2,4, 0,2,3, 0,2,2, 0,2,1, 0,2,0, 0,1,4, 0,1,3, 0,1,2, 0,1,1, 0,1,0, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,-1,4, 0,-1,3, 0,-1,2, 0,-1,1, 0,-1,0]

    integer(kind=4) :: a3d_p040_p004_mixed_yz(39) = [0,4,0, 0,4,-4, 0,3,0, 0,3,-3, 0,2,0, 0,2,-2, 0,1,0, 0,1,-1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4]
    integer(kind=4) :: a3d_p041_p003_mixed_yz(75) = [0,4,1, 0,4,0, 0,4,-1, 0,4,-2, 0,4,-3, 0,3,1, 0,3,0, 0,3,-1, 0,3,-2, 0,3,-3, 0,2,1, 0,2,0, 0,2,-1, 0,2,-2, 0,2,-3, 0,1,1, 0,1,0, 0,1,-1, 0,1,-2, 0,1,-3, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3]
    integer(kind=4) :: a3d_p042_m002_mixed_yz(60) = [0,4,2, 0,4,1, 0,4,-1, 0,4,-2, 0,3,2, 0,3,1, 0,3,-1, 0,3,-2, 0,2,2, 0,2,1, 0,2,-1, 0,2,-2, 0,1,2, 0,1,1, 0,1,-1, 0,1,-2, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2]
    integer(kind=4) :: a3d_p043_m001_mixed_yz(75) = [0,4,3, 0,4,2, 0,4,1, 0,4,0, 0,4,-1, 0,3,3, 0,3,2, 0,3,1, 0,3,0, 0,3,-1, 0,2,3, 0,2,2, 0,2,1, 0,2,0, 0,2,-1, 0,1,3, 0,1,2, 0,1,1, 0,1,0, 0,1,-1, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1]
    integer(kind=4) :: a3d_p044_p000_mixed_yz(39) = [0,4,4, 0,4,0, 0,3,3, 0,3,0, 0,2,2, 0,2,0, 0,1,1, 0,1,0, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0]
    integer(kind=4) :: a3d_p044_m044_mixed_yz(48) = [0,4,4, 0,4,-4, 0,3,3, 0,3,-3, 0,2,2, 0,2,-2, 0,1,1, 0,1,-1, 0,-1,1, 0,-1,-1, 0,-2,2, 0,-2,-2, 0,-3,3, 0,-3,-3, 0,-4,4, 0,-4,-4]

    integer(kind=4) :: a3d_p055_m055_mixed_yz(60) = [0,5,5, 0,5,-5, 0,4,4, 0,4,-4, 0,3,3, 0,3,-3, 0,2,2, 0,2,-2, 0,1,1, 0,1,-1, 0,-1,1, 0,-1,-1, 0,-2,2, 0,-2,-2, 0,-3,3, 0,-3,-3, 0,-4,4, 0,-4,-4, 0,-5,5, 0,-5,-5]

!------------------------------------------------------------------------------------------------------------

!   *-----------------------------------------OPS Declarations-----------------------------------------------*

!   Declare OPS Block
    call ops_decl_block(3, senga_grid, "SENGA_GRID")

!   Declare OPS Dats
    d_size = [nxglbl, nyglbl, nzglbl]
    d_m    = [0,0,0]
    d_p    = [0,0,0]

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store1, "real(kind=8)", "STORE1")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store2, "real(kind=8)", "STORE2")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store3, "real(kind=8)", "STORE3")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store4, "real(kind=8)", "STORE4")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store5, "real(kind=8)", "STORE5")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store6, "real(kind=8)", "STORE6")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_divm, "real(kind=8)", "DIVM")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ucor, "real(kind=8)", "UCOR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vcor, "real(kind=8)", "VCOR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wcor, "real(kind=8)", "WCOR")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd1x, "real(kind=8)", "WD1X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd1x, "real(kind=8)", "PD1X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td1x, "real(kind=8)", "TD1X")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd1y, "real(kind=8)", "WD1Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd1y, "real(kind=8)", "PD1Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td1y, "real(kind=8)", "TD1Y")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd1z, "real(kind=8)", "WD1Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd1z, "real(kind=8)", "PD1Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td1z, "real(kind=8)", "TD1Z")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd2x, "real(kind=8)", "WD2X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd2x, "real(kind=8)", "PD2X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td2x, "real(kind=8)", "TD2X")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd2y, "real(kind=8)", "WD2Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd2y, "real(kind=8)", "PD2Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td2y, "real(kind=8)", "TD2Y")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd2z, "real(kind=8)", "WD2Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd2z, "real(kind=8)", "PD2Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td2z, "real(kind=8)", "TD2Z")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ufxl, "real(kind=8)", "UFXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vfxl, "real(kind=8)", "VFXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wfxl, "real(kind=8)", "WFXL")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_drun, "real(kind=8)", "DRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_urun, "real(kind=8)", "URUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vrun, "real(kind=8)", "VRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wrun, "real(kind=8)", "WRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_erun, "real(kind=8)", "ERUN")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_derr, "real(kind=8)", "DERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_uerr, "real(kind=8)", "UERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_verr, "real(kind=8)", "VERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_werr, "real(kind=8)", "WERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_eerr, "real(kind=8)", "EERR")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d2prun, "real(kind=8)", "PRN2")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d2trun, "real(kind=8)", "TRN2")

!---------------------------------------MULTI-DIM DAT--------------------------------------------------------

    d_size = [nxglbl, nyglbl, nzglbl]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    DO ispec = 1,nspcmx
        WRITE(buf,"(A4,I2.2)") "YRUN",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yrun(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A4,I2.2)") "YERR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yerr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A4,I2.2)") "RATE",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_rate(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A4,I2.2)") "RRTE",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_rrte(ispec), "real(kind=8)", trim(buf))
    END DO

    d_size = [ nxglbl, nyglbl, nzglbl]
    d_m    = [-nhalox,-nhaloy,-nhaloz]
    d_p    = [ nhalox, nhaloy, nhaloz]
    DO iindex = 1,nintmx
        WRITE(buf,"(A6,I2.2)") "ITNDEX",iindex
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_int_null, d_itndex(iindex), "integer(kind=4)", trim(buf))
    END DO
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_yrhs_mdim, "real(kind=8)", "YRHS-MDIM")

    DO ispec = 1,nspcmx
        WRITE(buf,"(A4,I2.2)") "YRHS",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yrhs(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "CTRANS",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ctrans(ispec), "real(kind=8)", trim(buf))
    END DO

    call ops_decl_dat(senga_grid, nctmax+1, d_size, d_base, d_m, d_p, temp_real_null, d_tcoeff, "real(kind=8)", "TCOEFF")
    call ops_decl_dat(senga_grid, nctmax, d_size, d_base, d_m, d_p, temp_real_null, d_tderiv, "real(kind=8)", "TDERIV")

    d_size = [1, nyglbl, nzglbl]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    DO ispec = 1,nspcmx
        WRITE(buf,"(A6,I2.2)") "BCLYXL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyxl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "BCLYXR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyxr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRYXL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryxl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRYXR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryxr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "DYDTXL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtxl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "DYDTXR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtxr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "RATEXL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ratexl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "RATEXR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ratexr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRHXL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhxl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRHXR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhxr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "T6BXL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6bxl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "T6BXR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6bxr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "TT6XL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6xl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "TT6XR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6xr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "YINF1",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yinf1(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STYLXL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stlyxl(ispec), "real(kind=8)", trim(buf))
    END DO

    d_size = [nxglbl, 1, nzglbl]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    DO ispec = 1,nspcmx
        WRITE(buf,"(A6,I2.2)") "BCLYYL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyyl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "BCLYYR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyyr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRYYL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryyl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRYYR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryyr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "DYDTYL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtyl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "DYDTYR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtyr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "RATEYL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_rateyl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "RATEYR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_rateyr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRHYL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhyl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRHYR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhyr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "T6BYL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6byl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "T6BYR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6byr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "TT6YL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6yl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "TT6YR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6yr(ispec), "real(kind=8)", trim(buf))
    END DO

    d_size = [nxglbl, nyglbl, 1]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    DO ispec = 1,nspcmx
        WRITE(buf,"(A6,I2.2)") "BCLYZL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyzl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "BCLYZR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyzr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRYZL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryzl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRYZR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryzr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "DYDTZL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtzl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "DYDTZR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtzr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "RATEZL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ratezl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "RATEZR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ratezr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRHZL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhzl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A6,I2.2)") "STRHZR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhzr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "T6BZL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6bzl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "T6BZR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6bzr(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "TT6ZL",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6zl(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A5,I2.2)") "TT6ZR",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6zr(ispec), "real(kind=8)", trim(buf))
    END DO

!---------------------------------------WITH HALOS-----------------------------------------------------------

    d_size = [nxglbl, nyglbl, nzglbl]
    d_m    = [-nhalox,-nhaloy,-nhaloz]
    d_p    = [nhalox,nhaloy,nhaloz]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_drhs, "real(kind=8)", "DRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_urhs, "real(kind=8)", "URHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vrhs, "real(kind=8)", "VRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wrhs, "real(kind=8)", "WRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_erhs, "real(kind=8)", "ERHS")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_utmp, "real(kind=8)", "UTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vtmp, "real(kind=8)", "VTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wtmp, "real(kind=8)", "WTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_prun, "real(kind=8)", "PRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_trun, "real(kind=8)", "TRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_transp, "real(kind=8)", "TRANSP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store7, "real(kind=8)", "STORE7")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wmomix, "real(kind=8)", "WMOMIX")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_difmix, "real(kind=8)", "DIFMIX")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tdrmix, "real(kind=8)", "TDRMIX")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_combo1, "real(kind=8)", "COMBO1")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_combo2, "real(kind=8)", "COMBO2")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_combo3, "real(kind=8)", "COMBO3")

!-----------------------------------------Boundary YZ--------------------------------------------------------

    d_size = [1,nyglbl,nzglbl]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_totyxl, "real(kind=8)", "TOTYXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1xl, "real(kind=8)", "BCL1XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1xr, "real(kind=8)", "BCL1XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2xl, "real(kind=8)", "BCL2XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2xr, "real(kind=8)", "BCL2XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3xl, "real(kind=8)", "BCL3XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3xr, "real(kind=8)", "BCL3XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4xl, "real(kind=8)", "BCL4XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4xr, "real(kind=8)", "BCL4XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5xl, "real(kind=8)", "BCL5XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5xr, "real(kind=8)", "BCL5XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltxl, "real(kind=8)", "BCLTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltxr, "real(kind=8)", "BCLTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struxl, "real(kind=8)", "STRUXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struxr, "real(kind=8)", "STRUXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvxl, "real(kind=8)", "STRVXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvxr, "real(kind=8)", "STRVXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwxl, "real(kind=8)", "STRWXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwxr, "real(kind=8)", "STRWXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpxl, "real(kind=8)", "STRPXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpxr, "real(kind=8)", "STRPXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdxl, "real(kind=8)", "STRDXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdxr, "real(kind=8)", "STRDXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtxl, "real(kind=8)", "STRTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtxr, "real(kind=8)", "STRTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strexl, "real(kind=8)", "STREXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strexr, "real(kind=8)", "STREXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgxl, "real(kind=8)", "STRGXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgxr, "real(kind=8)", "STRGXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrxl, "real(kind=8)", "STRRXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrxr, "real(kind=8)", "STRRXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtxl, "real(kind=8)", "DUDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtxr, "real(kind=8)", "DUDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtxl, "real(kind=8)", "DVDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtxr, "real(kind=8)", "DVDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtxl, "real(kind=8)", "DWDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtxr, "real(kind=8)", "DWDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtxl, "real(kind=8)", "DTDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtxr, "real(kind=8)", "DTDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtxl, "real(kind=8)", "DDDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtxr, "real(kind=8)", "DDDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouxl, "real(kind=8)", "ACOUXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouxr, "real(kind=8)", "ACOUXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2xl, "real(kind=8)", "OVA2XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2xr, "real(kind=8)", "OVA2XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1xl, "real(kind=8)", "GAM1XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1xr, "real(kind=8)", "GAM1XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmxl, "real(kind=8)", "OVGMXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmxr, "real(kind=8)", "OVGMXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtxl, "real(kind=8)", "SYDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtxr, "real(kind=8)", "SYDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpxl, "real(kind=8)", "SORPXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpxr, "real(kind=8)", "SORPXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1bxl, "real(kind=8)", "T1BXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1bxr, "real(kind=8)", "T1BXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2bxl, "real(kind=8)", "T2BXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2bxr, "real(kind=8)", "T2BXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3bxl, "real(kind=8)", "T3BXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3bxr, "real(kind=8)", "T3BXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4bxl, "real(kind=8)", "T4BXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4bxr, "real(kind=8)", "T4BXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51bxl, "real(kind=8)", "T51BXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51bxr, "real(kind=8)", "T51BXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52bxl, "real(kind=8)", "T52BXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52bxr, "real(kind=8)", "T52BXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_uinf1, "real(kind=8)", "UINF1")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vinf1, "real(kind=8)", "VINF1")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_winf1, "real(kind=8)", "WINF1")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ustead, "real(kind=8)", "USTEAD")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1xl, "real(kind=8)", "TT1XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1xr, "real(kind=8)", "TT1XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2xl, "real(kind=8)", "TT2XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2xr, "real(kind=8)", "TT2XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3xl, "real(kind=8)", "TT3XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3xr, "real(kind=8)", "TT3XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4xl, "real(kind=8)", "TT4XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4xr, "real(kind=8)", "TT4XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5xl, "real(kind=8)", "TT5XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5xr, "real(kind=8)", "TT5XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stluxl, "real(kind=8)", "STLUXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stlvxl, "real(kind=8)", "STLVXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stlwxl, "real(kind=8)", "STLWXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stltxl, "real(kind=8)", "STLTXL")

!-----------------------------------------Boundary XZ--------------------------------------------------------

    d_size = [nxglbl,1,nzglbl]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1yl, "real(kind=8)", "BCL1YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1yr, "real(kind=8)", "BCL1YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2yl, "real(kind=8)", "BCL2YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2yr, "real(kind=8)", "BCL2YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3yl, "real(kind=8)", "BCL3YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3yr, "real(kind=8)", "BCL3YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4yl, "real(kind=8)", "BCL4YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4yr, "real(kind=8)", "BCL4YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5yl, "real(kind=8)", "BCL5YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5yr, "real(kind=8)", "BCL5YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltyl, "real(kind=8)", "BCLTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltyr, "real(kind=8)", "BCLTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struyl, "real(kind=8)", "STRUYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struyr, "real(kind=8)", "STRUYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvyl, "real(kind=8)", "STRVYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvyr, "real(kind=8)", "STRVYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwyl, "real(kind=8)", "STRWYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwyr, "real(kind=8)", "STRWYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpyl, "real(kind=8)", "STRPYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpyr, "real(kind=8)", "STRPYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdyl, "real(kind=8)", "STRDYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdyr, "real(kind=8)", "STRDYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtyl, "real(kind=8)", "STRTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtyr, "real(kind=8)", "STRTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_streyl, "real(kind=8)", "STREYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_streyr, "real(kind=8)", "STREYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgyl, "real(kind=8)", "STRGYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgyr, "real(kind=8)", "STRGYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strryl, "real(kind=8)", "STRRYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strryr, "real(kind=8)", "STRRYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtyl, "real(kind=8)", "DUDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtyr, "real(kind=8)", "DUDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtyl, "real(kind=8)", "DVDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtyr, "real(kind=8)", "DVDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtyl, "real(kind=8)", "DWDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtyr, "real(kind=8)", "DWDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtyl, "real(kind=8)", "DTDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtyr, "real(kind=8)", "DTDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtyl, "real(kind=8)", "DDDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtyr, "real(kind=8)", "DDDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouyl, "real(kind=8)", "ACOUYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouyr, "real(kind=8)", "ACOUYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2yl, "real(kind=8)", "OVA2YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2yr, "real(kind=8)", "OVA2YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1yl, "real(kind=8)", "GAM1YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1yr, "real(kind=8)", "GAM1YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmyl, "real(kind=8)", "OVGMYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmyr, "real(kind=8)", "OVGMYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtyl, "real(kind=8)", "SYDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtyr, "real(kind=8)", "SYDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpyl, "real(kind=8)", "SORPYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpyr, "real(kind=8)", "SORPYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1byl, "real(kind=8)", "T1BYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1byr, "real(kind=8)", "T1BYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2byl, "real(kind=8)", "T2BYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2byr, "real(kind=8)", "T2BYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3byl, "real(kind=8)", "T3BYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3byr, "real(kind=8)", "T3BYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4byl, "real(kind=8)", "T4BYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4byr, "real(kind=8)", "T4BYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51byl, "real(kind=8)", "T51BYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51byr, "real(kind=8)", "T51BYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52byl, "real(kind=8)", "T52BYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52byr, "real(kind=8)", "T52BYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1yl, "real(kind=8)", "TT1YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1yr, "real(kind=8)", "TT1YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2yl, "real(kind=8)", "TT2YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2yr, "real(kind=8)", "TT2YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3yl, "real(kind=8)", "TT3YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3yr, "real(kind=8)", "TT3YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4yl, "real(kind=8)", "TT4YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4yr, "real(kind=8)", "TT4YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5yl, "real(kind=8)", "TT5YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5yr, "real(kind=8)", "TT5YR")

!-----------------------------------------Boundary XY--------------------------------------------------------

    d_size = [nxglbl,nyglbl,1]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1zl, "real(kind=8)", "BCL1ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1zr, "real(kind=8)", "BCL1ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2zl, "real(kind=8)", "BCL2ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2zr, "real(kind=8)", "BCL2ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3zl, "real(kind=8)", "BCL3ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3zr, "real(kind=8)", "BCL3ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4zl, "real(kind=8)", "BCL4ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4zr, "real(kind=8)", "BCL4ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5zl, "real(kind=8)", "BCL5ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5zr, "real(kind=8)", "BCL5ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltzl, "real(kind=8)", "BCLTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltzr, "real(kind=8)", "BCLTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struzl, "real(kind=8)", "STRUZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struzr, "real(kind=8)", "STRUZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvzl, "real(kind=8)", "STRVZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvzr, "real(kind=8)", "STRVZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwzl, "real(kind=8)", "STRWZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwzr, "real(kind=8)", "STRWZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpzl, "real(kind=8)", "STRPZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpzr, "real(kind=8)", "STRPZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdzl, "real(kind=8)", "STRDZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdzr, "real(kind=8)", "STRDZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtzl, "real(kind=8)", "STRTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtzr, "real(kind=8)", "STRTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strezl, "real(kind=8)", "STREZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strezr, "real(kind=8)", "STREZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgzl, "real(kind=8)", "STRGZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgzr, "real(kind=8)", "STRGZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrzl, "real(kind=8)", "STRRZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrzr, "real(kind=8)", "STRRZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtzl, "real(kind=8)", "DUDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtzr, "real(kind=8)", "DUDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtzl, "real(kind=8)", "DVDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtzr, "real(kind=8)", "DVDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtzl, "real(kind=8)", "DWDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtzr, "real(kind=8)", "DWDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtzl, "real(kind=8)", "DTDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtzr, "real(kind=8)", "DTDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtzl, "real(kind=8)", "DDDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtzr, "real(kind=8)", "DDDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouzl, "real(kind=8)", "ACOUZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouzr, "real(kind=8)", "ACOUZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2zl, "real(kind=8)", "OVA2ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2zr, "real(kind=8)", "OVA2ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1zl, "real(kind=8)", "GAM1ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1zr, "real(kind=8)", "GAM1ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmzl, "real(kind=8)", "OVGMZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmzr, "real(kind=8)", "OVGMZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtzl, "real(kind=8)", "SYDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtzr, "real(kind=8)", "SYDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpzl, "real(kind=8)", "SORPZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpzr, "real(kind=8)", "SORPZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1bzl, "real(kind=8)", "T1BZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1bzr, "real(kind=8)", "T1BZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2bzl, "real(kind=8)", "T2BZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2bzr, "real(kind=8)", "T2BZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3bzl, "real(kind=8)", "T3BZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3bzr, "real(kind=8)", "T3BZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4bzl, "real(kind=8)", "T4BZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4bzr, "real(kind=8)", "T4BZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51bzl, "real(kind=8)", "T51BZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51bzr, "real(kind=8)", "T51BZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52bzl, "real(kind=8)", "T52BZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52bzr, "real(kind=8)", "T52BZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1zl, "real(kind=8)", "TT1ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1zr, "real(kind=8)", "TT1ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2zl, "real(kind=8)", "TT2ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2zr, "real(kind=8)", "TT2ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3zl, "real(kind=8)", "TT3ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3zr, "real(kind=8)", "TT3ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4zl, "real(kind=8)", "TT4ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4zr, "real(kind=8)", "TT4ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5zl, "real(kind=8)", "TT5ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5zr, "real(kind=8)", "TT5ZR")

!------------------------------------Only X-direction--------------------------------------------------------

    d_size = [nxglbl,1,1]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_crin, "real(kind=8)", "CRIN")

!------------------------------------Inflow Filter YZ--------------------------------------------------------

    d_size = [1,nyglbl,nzglbl]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ufilt, "real(kind=8)", "UFILT")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vfilt, "real(kind=8)", "VFILT")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wfilt, "real(kind=8)", "WFILT")
    DO ispec = 1,nspcmx
        WRITE(buf,"(A5,I2.2)") "YFILT",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yfilt(ispec), "real(kind=8)", trim(buf))
    END DO


    d_size = [1,nyglbl,nzglbl]
    d_m    = [0,0,-nfz]
    d_p    = [0,0,nfz]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ufold, "real(kind=8)", "UFOLD")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vfold, "real(kind=8)", "VFOLD")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wfold, "real(kind=8)", "WFOLD")
    DO ispec = 1,nspcmx
        WRITE(buf,"(A5,I2.2)") "YFOLD",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yfold(ispec), "real(kind=8)", trim(buf))
    END DO


    d_size = [1,nyglbl,nzglbl]
    d_m    = [0,-nfy,-nfz]
    d_p    = [0,nfy,nfz]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_urand, "real(kind=8)", "URAND")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vrand, "real(kind=8)", "VRAND")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wrand, "real(kind=8)", "WRAND")
    DO ispec = 1,nspcmx
        WRITE(buf,"(A5,I2.2)") "YRAND",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yrand(ispec), "real(kind=8)", trim(buf))
    END DO

!------------------------------------Check for COLD/WARM start-----------------------------------------------

    fncont = pncont//pnxdat
!   Open the input file and find if its cold start/restart
    OPEN(UNIT=nccont,FILE=fncont,STATUS='OLD',FORM='FORMATTED')
    DO line = 1,16
        READ(nccont,*)
    END DO
    READ(nccont,*)tstep,ntime1,ntime,nstpsw
    READ(nccont,*)
    READ(nccont,*)
    READ(nccont,*)ntdump,ntrept,ntstat,ndofmt
    READ(nccont,*)
    READ(nccont,*)
!   COLD START SWITCH (0=COLD START, 1=RESTART), DUMP INPUT FORMAT
    READ(nccont,*)ncdmpi,ndifmt
    DO line = 1,6
        READ(nccont,*)
    END DO
    READ(nccont,*)inflam
    DO line = 1,5
        READ(nccont,*)
    END DO
    READ(nccont,*)nspreq
    DO ispec = 1, nspreq
        READ(nccont,*)jspec,yrin(ispec)
    END DO
    DO line = 1,5
        READ(nccont,*)
    END DO
    READ(nccont,*)ngbcxl,(nxlprm(ic),ic=1,nbcpri), (rxlprm(ic),ic=1,nbcprr)

    CLOSE(nccont)

!   ====================
!   WARM START
!   ====================
    IF(ncdmpi == 1) THEN

        dtime=INT((ntime1-1)/ntdump)
        WRITE(citime,'(I8.8)') dtime

        fname = 'output/timestep'//citime//pnxhdf

!       -----------------------------------------------------------------------
!       THIS BLOCK MAY BE MODIFIED AS REQUIRED
!       TO BLEND INITIAL VELOCITY AND SCALAR FIELDS
!       WITH PREVIOUSLY DUMPED DATA
!       -----------------------------------------------------------------------

!       RESTART FROM FULL DUMP FILES
!       ----------------------------
!       READ THE DATA FROM DUMP INPUT FILE
!       NOTE THAT URUN,VRUN,WRUN,ERUN AND YRUN ARE ALL IN CONSERVATIVE FORM
        IF (ops_is_root() == 1) THEN
            WRITE(*,*) "Warm Start: start step -> ", ntime1
            WRITE(*,*) "Reading from dumped file: ", trim(fname)
        END IF

        call ops_decl_dat_hdf5(d_drun_dump, senga_grid, 1, "real(kind=8)", "DRUN", trim(fname), status)
        call ops_decl_dat_hdf5(d_urun_dump, senga_grid, 1, "real(kind=8)", "URUN", trim(fname), status)
        call ops_decl_dat_hdf5(d_vrun_dump, senga_grid, 1, "real(kind=8)", "VRUN", trim(fname), status)
        call ops_decl_dat_hdf5(d_wrun_dump, senga_grid, 1, "real(kind=8)", "WRUN", trim(fname), status)
        call ops_decl_dat_hdf5(d_erun_dump, senga_grid, 1, "real(kind=8)", "ERUN", trim(fname), status)

        DO ispec = 1,nspcmx
            WRITE(buf,"(A4,I2.2)") "YRUN",ispec
             call ops_decl_dat_hdf5(d_yrun_dump(ispec), senga_grid, 1, "real(kind=8)", trim(buf), trim(fname), status)
        END DO

    END IF

    IF(nxlprm(1) == 4) THEN
        IF(ncdmpi == 1) THEN

            fname = 'output/inflow'//pnxhdf

            call ops_decl_dat_hdf5(d_uinf2, senga_grid, 1, "real(kind=8)", "UINF2", trim(fname), status)
            call ops_decl_dat_hdf5(d_vinf2, senga_grid, 1, "real(kind=8)", "VINF2", trim(fname), status)
            call ops_decl_dat_hdf5(d_winf2, senga_grid, 1, "real(kind=8)", "WINF2", trim(fname), status)

            DO ispec = 1,nspcmx
                WRITE(buf,"(A5,I2.2)") "YINF2",ispec
                call ops_decl_dat_hdf5(d_yinf2(ispec), senga_grid, 1, "real(kind=8)", trim(buf), trim(fname), status)
            END DO

        ELSE
            d_size = [1,nyglbl,nzglbl]
            d_m    = [0,0,0]
            d_p    = [0,0,0]

            call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_uinf2, "real(kind=8)", "UINF2")
            call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vinf2, "real(kind=8)", "VINF2")
            call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_winf2, "real(kind=8)", "WINF2")
            DO ispec = 1,nspcmx
                WRITE(buf,"(A5,I2.2)") "YINF2",ispec
                call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yinf2(ispec), "real(kind=8)", trim(buf))
            END DO

        END IF
    END IF

    IF(inflam == 1) THEN

        fname = 'input/flame_init'//pnxhdf

        call ops_decl_dat_hdf5(d_drun_dump, senga_grid, 1, "real(kind=8)", "DRUN", trim(fname), status)
        call ops_decl_dat_hdf5(d_urun_dump, senga_grid, 1, "real(kind=8)", "URUN", trim(fname), status)
        call ops_decl_dat_hdf5(d_vrun_dump, senga_grid, 1, "real(kind=8)", "VRUN", trim(fname), status)
        call ops_decl_dat_hdf5(d_wrun_dump, senga_grid, 1, "real(kind=8)", "WRUN", trim(fname), status)
        call ops_decl_dat_hdf5(d_trun_dump, senga_grid, 1, "real(kind=8)", "TRUN", trim(fname), status)

        DO ispec = 1,nspcmx
            WRITE(buf,"(A4,I2.2)") "YRUN",ispec
             call ops_decl_dat_hdf5(d_yrun_dump(ispec), senga_grid, 1, "real(kind=8)", trim(buf), trim(fname), status)
        END DO

    END IF

!------------------------------------OPS Reduction Handles---------------------------------------------------

    call ops_decl_reduction_handle(8, h_erdtot, "real(kind=8)", "erdtot")
    call ops_decl_reduction_handle(8, h_erutot, "real(kind=8)", "erutot")
    call ops_decl_reduction_handle(8, h_ervtot, "real(kind=8)", "ervtot")
    call ops_decl_reduction_handle(8, h_erwtot, "real(kind=8)", "erwtot")
    call ops_decl_reduction_handle(8, h_eretot, "real(kind=8)", "eretot")
    call ops_decl_reduction_handle(8, h_erytot, "real(kind=8)", "erytot")
    call ops_decl_reduction_handle(8, h_tket, "real(kind=8)", "tket")
    call ops_decl_reduction_handle(8, h_tkes, "real(kind=8)", "tkes")
    call ops_decl_reduction_handle(8, h_enstro, "real(kind=8)", "enstro")
    call ops_decl_reduction_handle(8, h_ubart, "real(kind=8)", "ubart")
    call ops_decl_reduction_handle(8, h_vbart, "real(kind=8)", "vbart")
    call ops_decl_reduction_handle(8, h_wbart, "real(kind=8)", "wbart")
    call ops_decl_reduction_handle(8, h_uvart, "real(kind=8)", "uvart")
    call ops_decl_reduction_handle(8, h_vvart, "real(kind=8)", "vvart")
    call ops_decl_reduction_handle(8, h_wvart, "real(kind=8)", "wvart")
    call ops_decl_reduction_handle(8, h_umean, "real(kind=8)", "umean")
    call ops_decl_reduction_handle(8, h_denom, "real(kind=8)", "denom")

!------------------------------------OPS Stencil-------------------------------------------------------------

    call ops_decl_stencil( 3, 1, a3d_000, s3d_000, "0,0,0")

    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_x, s3d_000_strid3d_x, "stride 3D X dir")
    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_y, s3d_000_strid3d_y, "stride 3D Y dir")
    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_z, s3d_000_strid3d_z, "stride 3D Z dir")

    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_xy, s3d_000_strid3d_xy, "stride 3D XY dir")
    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_xz, s3d_000_strid3d_xz, "stride 3D XZ dir")
    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_yz, s3d_000_strid3d_yz, "stride 3D YZ dir")

!------------------------------------------------------------------------------------------------------------

    call ops_decl_stencil( 3, 5, a3d_000_to_p400_x, s3d_000_to_p400_x, "0,0,0 to 4,0,0")
    call ops_decl_stencil( 3, 5, a3d_000_to_m400_x, s3d_000_to_m400_x, "0,0,0 to -4,0,0")

    call ops_decl_stencil( 3, 4, a3d_p100_to_p400_x, s3d_p100_to_p400_x, "1,0,0 to 4,0,0")
    call ops_decl_stencil( 3, 4, a3d_m100_to_m400_x, s3d_m100_to_m400_x, "-1,0,0 to -4,0,0")

    call ops_decl_stencil( 3, 11, a3d_p500_to_m500_x, s3d_p500_to_m500_x, "5,0,0 to -5,0,0")

!------------------------------------------------------------------------------------------------------------

    call ops_decl_stencil( 3, 5, a3d_000_to_p040_y, s3d_000_to_p040_y, "0,0,0 to 0,4,0")
    call ops_decl_stencil( 3, 5, a3d_000_to_m040_y, s3d_000_to_m040_y, "0,0,0 to  0,-4,0")

    call ops_decl_stencil( 3, 4, a3d_p010_to_p040_y, s3d_p010_to_p040_y, "0,1,0 to 0,4,0")
    call ops_decl_stencil( 3, 4, a3d_m010_to_m040_y, s3d_m010_to_m040_y, "0,-1,0 to 0,-4,0")

    call ops_decl_stencil( 3, 11, a3d_p050_to_m050_y, s3d_p050_to_m050_y, "0,5,0 to  0,-5,0")

!------------------------------------------------------------------------------------------------------------

    call ops_decl_stencil( 3, 5, a3d_000_to_p004_z, s3d_000_to_p004_z, "0,0,0 to 0,0,4")
    call ops_decl_stencil( 3, 5, a3d_000_to_m004_z, s3d_000_to_m004_z, "0,0,0 to  0,0,-4")

    call ops_decl_stencil( 3, 4, a3d_p001_to_p004_z, s3d_p001_to_p004_z, "0,0,1 to 0,0,4")
    call ops_decl_stencil( 3, 4, a3d_m001_to_m004_z, s3d_m001_to_m004_z, "0,0,-1 to 0,0,-4")

    call ops_decl_stencil( 3, 11, a3d_p005_to_m005_z, s3d_p005_to_m005_z, "0,0,5 to  0,0,-5")

!------------------------------------------------------------------------------------------------------------

    call ops_decl_stencil( 3, 13, a3d_p000_m440_mixed_xy, s3d_p000_m440_mixed_xy, "0,0,0 to -4,-4,0")
    call ops_decl_stencil( 3, 25, a3d_p010_m430_mixed_xy, s3d_p010_m430_mixed_xy, "0,1,0 to -4,-3,0")
    call ops_decl_stencil( 3, 20, a3d_p020_m420_mixed_xy, s3d_p020_m420_mixed_xy, "0,2,0 to -4,-2,0")
    call ops_decl_stencil( 3, 25, a3d_p030_m410_mixed_xy, s3d_p030_m410_mixed_xy, "0,3,0 to -4,-1,0")
    call ops_decl_stencil( 3, 13, a3d_p040_m400_mixed_xy, s3d_p040_m400_mixed_xy, "0,4,0 to -4,0,0")

    call ops_decl_stencil( 3, 13, a3d_p100_m300_mixed_xy, s3d_p100_m300_mixed_xy, "1,0,0 to -3,0,0")
    call ops_decl_stencil( 3, 25, a3d_p100_m340_mixed_xy, s3d_p100_m340_mixed_xy, "1,0,0 to -3,-4,0")
    call ops_decl_stencil( 3, 13, a3d_p110_m330_mixed_xy, s3d_p110_m330_mixed_xy, "1,1,0 to -3,-3,0")
    call ops_decl_stencil( 3, 20, a3d_p120_m320_mixed_xy, s3d_p120_m320_mixed_xy, "1,2,0 to -3,-2,0")
    call ops_decl_stencil( 3, 25, a3d_p140_m300_mixed_xy, s3d_p140_m300_mixed_xy, "1,4,0 to -3,0,0")

    call ops_decl_stencil( 3, 20, a3d_p200_m240_mixed_xy, s3d_p200_m240_mixed_xy, "2,0,0 to -2,-4,0")
    call ops_decl_stencil( 3, 20, a3d_p210_m230_mixed_xy, s3d_p210_m230_mixed_xy, "2,1,0 to -2,-3,0")
    call ops_decl_stencil( 3, 8,  a3d_p220_m220_mixed_xy, s3d_p220_m220_mixed_xy, "2,2,0 to -2,-2,0")
    call ops_decl_stencil( 3, 20, a3d_p230_m210_mixed_xy, s3d_p230_m210_mixed_xy, "2,3,0 to -2,-1,0")
    call ops_decl_stencil( 3, 20, a3d_p240_m200_mixed_xy, s3d_p240_m200_mixed_xy, "2,4,0 to -2,0,0")

    call ops_decl_stencil( 3, 13, a3d_p300_m100_mixed_xy, s3d_p300_m100_mixed_xy, "3,0,0 to -1,0,0")
    call ops_decl_stencil( 3, 25, a3d_p300_m140_mixed_xy, s3d_p300_m140_mixed_xy, "3,0,0 to -1,-4,0")
    call ops_decl_stencil( 3, 20, a3d_p320_m120_mixed_xy, s3d_p320_m120_mixed_xy, "3,2,0 to -1,-2,0")
    call ops_decl_stencil( 3, 13, a3d_p330_m110_mixed_xy, s3d_p330_m110_mixed_xy, "3,3,0 to -1,-1,0")
    call ops_decl_stencil( 3, 12, a3d_p330_m330_mixed_xy, s3d_p330_m330_mixed_xy, "3,3,0 to -3,-3,0")
    call ops_decl_stencil( 3, 25, a3d_p340_m100_mixed_xy, s3d_p340_m100_mixed_xy, "3,4,0 to -1,0,0")

    call ops_decl_stencil( 3, 13, a3d_p400_p040_mixed_xy, s3d_p400_p040_mixed_xy, "4,0,0 to 0,-4,0")
    call ops_decl_stencil( 3, 25, a3d_p410_p030_mixed_xy, s3d_p410_p030_mixed_xy, "4,1,0 to 0,-3,0")
    call ops_decl_stencil( 3, 20, a3d_p420_m020_mixed_xy, s3d_p420_m020_mixed_xy, "4,2,0 to 0,-2,0")
    call ops_decl_stencil( 3, 25, a3d_p430_m010_mixed_xy, s3d_p430_m010_mixed_xy, "4,3,0 to 0,-1,0")
    call ops_decl_stencil( 3, 13, a3d_p440_p000_mixed_xy, s3d_p440_p000_mixed_xy, "4,4,0 to 0,0,0")
    call ops_decl_stencil( 3, 16, a3d_p440_m440_mixed_xy, s3d_p440_m440_mixed_xy, "4,4,0 to -4,-4,0")

    call ops_decl_stencil( 3, 20, a3d_p550_m550_mixed_xy, s3d_p550_m550_mixed_xy, "5,5,0 to -5,-5,0")

!------------------------------------------------------------------------------------------------------------

    call ops_decl_stencil( 3, 13, a3d_p000_m404_mixed_xz, s3d_p000_m404_mixed_xz, "0,0,0 to -4,0,-4")
    call ops_decl_stencil( 3, 25, a3d_p001_m403_mixed_xz, s3d_p001_m403_mixed_xz, "0,0,1 to -4,0,-3")
    call ops_decl_stencil( 3, 20, a3d_p002_m402_mixed_xz, s3d_p002_m402_mixed_xz, "0,0,2 to -4,0,-2")
    call ops_decl_stencil( 3, 25, a3d_p003_m401_mixed_xz, s3d_p003_m401_mixed_xz, "0,0,3 to -4,0,-1")
    call ops_decl_stencil( 3, 13, a3d_p004_m400_mixed_xz, s3d_p004_m400_mixed_xz, "0,0,4 to -4,0,0")

    call ops_decl_stencil( 3, 13, a3d_p100_m300_mixed_xz, s3d_p100_m300_mixed_xz, "1,0,0 to -3,0,0")
    call ops_decl_stencil( 3, 25, a3d_p100_m304_mixed_xz, s3d_p100_m304_mixed_xz, "1,0,0 to -3,0,-4")
    call ops_decl_stencil( 3, 13, a3d_p101_m303_mixed_xz, s3d_p101_m303_mixed_xz, "1,0,1 to -3,0,-3")
    call ops_decl_stencil( 3, 20, a3d_p102_m302_mixed_xz, s3d_p102_m302_mixed_xz, "1,0,2 to -3,0,-2")
    call ops_decl_stencil( 3, 25, a3d_p104_m300_mixed_xz, s3d_p104_m300_mixed_xz, "1,0,4 to -3,0,0")

    call ops_decl_stencil( 3, 20, a3d_p200_m204_mixed_xz, s3d_p200_m204_mixed_xz, "2,0,0 to -2,0,-4")
    call ops_decl_stencil( 3, 20, a3d_p201_m203_mixed_xz, s3d_p201_m203_mixed_xz, "2,0,1 to -2,0,-3")
    call ops_decl_stencil( 3, 8,  a3d_p202_m202_mixed_xz, s3d_p202_m202_mixed_xz, "2,0,2 to -2,0,-2")
    call ops_decl_stencil( 3, 20, a3d_p203_m201_mixed_xz, s3d_p203_m201_mixed_xz, "2,0,3 to -2,0,-1")
    call ops_decl_stencil( 3, 20, a3d_p204_m200_mixed_xz, s3d_p204_m200_mixed_xz, "2,0,4 to -2,0,0")

    call ops_decl_stencil( 3, 13, a3d_p300_m100_mixed_xz, s3d_p300_m100_mixed_xz, "3,0,0 to -1,0,0")
    call ops_decl_stencil( 3, 25, a3d_p300_m104_mixed_xz, s3d_p300_m104_mixed_xz, "3,0,0 to -1,0,-4")
    call ops_decl_stencil( 3, 20, a3d_p302_m102_mixed_xz, s3d_p302_m102_mixed_xz, "3,0,2 to -1,0,-2")
    call ops_decl_stencil( 3, 13, a3d_p303_m101_mixed_xz, s3d_p303_m101_mixed_xz, "3,0,3 to -1,0,-1")
    call ops_decl_stencil( 3, 12, a3d_p303_m303_mixed_xz, s3d_p303_m303_mixed_xz, "3,0,3 to -3,0,-3")
    call ops_decl_stencil( 3, 25, a3d_p304_m100_mixed_xz, s3d_p304_m100_mixed_xz, "3,0,4 to -1,0,0")

    call ops_decl_stencil( 3, 13, a3d_p400_p004_mixed_xz, s3d_p400_p004_mixed_xz, "4,0,0 to 0,0,-4")
    call ops_decl_stencil( 3, 25, a3d_p401_p003_mixed_xz, s3d_p401_p003_mixed_xz, "4,0,1 to 0,0,-3")
    call ops_decl_stencil( 3, 20, a3d_p402_m002_mixed_xz, s3d_p402_m002_mixed_xz, "4,0,2 to 0,0,-2")
    call ops_decl_stencil( 3, 25, a3d_p403_m001_mixed_xz, s3d_p403_m001_mixed_xz, "4,0,3 to 0,0,-1")
    call ops_decl_stencil( 3, 13, a3d_p404_p000_mixed_xz, s3d_p404_p000_mixed_xz, "4,0,4 to 0,0,0")
    call ops_decl_stencil( 3, 16, a3d_p404_m404_mixed_xz, s3d_p404_m404_mixed_xz, "4,0,4 to -4,0,-4")

    call ops_decl_stencil( 3, 20, a3d_p505_m505_mixed_xz, s3d_p505_m505_mixed_xz, "5,0,5 to -5,0,-5")

!------------------------------------------------------------------------------------------------------------

    call ops_decl_stencil( 3, 13, a3d_p000_m044_mixed_yz, s3d_p000_m044_mixed_yz, "0,0,0 to 0,-4,-4")
    call ops_decl_stencil( 3, 25, a3d_p001_m043_mixed_yz, s3d_p001_m043_mixed_yz, "0,0,1 to 0,-4,-3")
    call ops_decl_stencil( 3, 20, a3d_p002_m042_mixed_yz, s3d_p002_m042_mixed_yz, "0,0,2 to 0,-4,-2")
    call ops_decl_stencil( 3, 25, a3d_p003_m041_mixed_yz, s3d_p003_m041_mixed_yz, "0,0,3 to 0,-4,-1")
    call ops_decl_stencil( 3, 13, a3d_p004_m040_mixed_yz, s3d_p004_m040_mixed_yz, "0,0,4 to 0,-4,0")

    call ops_decl_stencil( 3, 13, a3d_p010_m030_mixed_yz, s3d_p010_m030_mixed_yz, "0,1,0 to 0,-3,0")
    call ops_decl_stencil( 3, 25, a3d_p010_m034_mixed_yz, s3d_p010_m034_mixed_yz, "0,1,0 to 0,-3,-4")
    call ops_decl_stencil( 3, 13, a3d_p011_m033_mixed_yz, s3d_p011_m033_mixed_yz, "0,1,1 to 0,-3,-3")
    call ops_decl_stencil( 3, 20, a3d_p012_m032_mixed_yz, s3d_p012_m032_mixed_yz, "0,1,2 to 0,-3,-2")
    call ops_decl_stencil( 3, 25, a3d_p014_m030_mixed_yz, s3d_p014_m030_mixed_yz, "0,1,4 to 0,-3,0")

    call ops_decl_stencil( 3, 20, a3d_p020_m024_mixed_yz, s3d_p020_m024_mixed_yz, "0,2,0 to 0,-2,-4")
    call ops_decl_stencil( 3, 20, a3d_p021_m023_mixed_yz, s3d_p021_m023_mixed_yz, "0,2,1 to 0,-2,-3")
    call ops_decl_stencil( 3, 8,  a3d_p022_m022_mixed_yz, s3d_p022_m022_mixed_yz, "0,2,2 to 0,-2,-2")
    call ops_decl_stencil( 3, 20, a3d_p023_m021_mixed_yz, s3d_p023_m021_mixed_yz, "0,2,3 to 0,-2,-1")
    call ops_decl_stencil( 3, 20, a3d_p024_m020_mixed_yz, s3d_p024_m020_mixed_yz, "0,2,4 to 0,-2,0")

    call ops_decl_stencil( 3, 13, a3d_p030_m010_mixed_yz, s3d_p030_m010_mixed_yz, "0,3,0 to 0,-1,0")
    call ops_decl_stencil( 3, 25, a3d_p030_m014_mixed_yz, s3d_p030_m014_mixed_yz, "0,3,0 to 0,-1,-4")
    call ops_decl_stencil( 3, 20, a3d_p032_m012_mixed_yz, s3d_p032_m012_mixed_yz, "0,3,2 to 0,-1,-2")
    call ops_decl_stencil( 3, 13, a3d_p033_m011_mixed_yz, s3d_p033_m011_mixed_yz, "0,3,3 to 0,-1,-1")
    call ops_decl_stencil( 3, 12, a3d_p033_m033_mixed_yz, s3d_p033_m033_mixed_yz, "0,3,3 to 0,-3,-3")
    call ops_decl_stencil( 3, 25, a3d_p034_m010_mixed_yz, s3d_p034_m010_mixed_yz, "0,3,4 to 0,-1,0")

    call ops_decl_stencil( 3, 13, a3d_p040_p004_mixed_yz, s3d_p040_p004_mixed_yz, "0,4,0 to 0,0,-4")
    call ops_decl_stencil( 3, 25, a3d_p041_p003_mixed_yz, s3d_p041_p003_mixed_yz, "0,4,1 to 0,0,-3")
    call ops_decl_stencil( 3, 20, a3d_p042_m002_mixed_yz, s3d_p042_m002_mixed_yz, "0,4,2 to 0,0,-2")
    call ops_decl_stencil( 3, 25, a3d_p043_m001_mixed_yz, s3d_p043_m001_mixed_yz, "0,4,3 to 0,0,-1")
    call ops_decl_stencil( 3, 13, a3d_p044_p000_mixed_yz, s3d_p044_p000_mixed_yz, "0,4,4 to 0,0,0")
    call ops_decl_stencil( 3, 16, a3d_p044_m044_mixed_yz, s3d_p044_m044_mixed_yz, "0,4,4 to 0,-4,-4")

    call ops_decl_stencil( 3, 20, a3d_p055_m055_mixed_yz, s3d_p055_m055_mixed_yz, "0,5,5 to 0,-5,-5")

!------------------------------------------------------------------------------------------------------------

!---------------------OPS DECL HALO for periodic transfer on single processor--------------------------------

    dir_from = [1,2,3]
    dir_to   = [1,2,3]

!   X-DIRECTION : RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
    iter_size = [nhalox,nyglbl,nzglbl]
    base_from = [1,1,1]
    base_to   = [nxglbl+1,1,1]

    halo_idx = 0

    halo_idx = halo_idx+1
    call ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    DO ispec = 1,nspcmx
        halo_idx = halo_idx+1
        call ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))
    END DO

!   X-DIRECTION : LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
    iter_size = [nhalox,nyglbl,nzglbl]
    base_from = [nxglbl-nhalox+1,1,1]
    base_to   = [1-nhalox,1,1]

    halo_idx = halo_idx+1
    call ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    DO ispec = 1,nspcmx
        halo_idx = halo_idx+1
        call ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))
    END DO

    call ops_decl_halo_group(halo_idx, halos_x, halos_grp_x)

!------------------------------------------------------------------------------------------------------------

!   Y-DIRECTION : RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
    iter_size = [nxglbl+2*nhalox,nhaloy,nzglbl]
    base_from = [1-nhalox,1,1]
    base_to   = [1-nhalox,nyglbl+1,1]

    halo_idx = 0

    halo_idx = halo_idx+1
    call ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    DO ispec = 1,nspcmx
        halo_idx = halo_idx+1
        call ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))
    END DO

!   Y-DIRECTION : LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
    iter_size = [nxglbl+2*nhalox,nhaloy,nzglbl]
    base_from = [1-nhalox,nyglbl-nhaloy+1,1]
    base_to   = [1-nhalox,1-nhaloy,1]

    halo_idx = halo_idx+1
    call ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    DO ispec = 1,nspcmx
        halo_idx = halo_idx+1
        call ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))
    END DO

    call ops_decl_halo_group(halo_idx, halos_y, halos_grp_y)

!------------------------------------------------------------------------------------------------------------

!   Z-DIRECTION : RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
    iter_size = [nxglbl+2*nhalox,nyglbl+2*nhaloy,nhaloz]
    base_from = [1-nhalox,1-nhaloy,1]
    base_to   = [1-nhalox,1-nhaloy,nzglbl+1]

    halo_idx = 0

    halo_idx = halo_idx+1
    call ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    DO ispec = 1,nspcmx
        halo_idx = halo_idx+1
        call ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))
    END DO

!   Z-DIRECTION : LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
    iter_size = [nxglbl+2*nhalox,nyglbl+2*nhaloy,nhaloz]
    base_from = [1-nhalox,1-nhaloy,nzglbl-nhaloz+1]
    base_to   = [1-nhalox,1-nhaloy,1-nhaloz]

    halo_idx = halo_idx+1
    call ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    halo_idx = halo_idx+1
    call ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    DO ispec = 1,nspcmx
        halo_idx = halo_idx+1
        call ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))
    END DO

    call ops_decl_halo_group(halo_idx, halos_z, halos_grp_z)

!------------------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------------------
    call ops_partition(" ")
!------------------------------------------------------------------------------------------------------------

!---------------------------------First touch - OPS DATS without halos---------------------------------------
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_divm, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd1x, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd1x, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td1x, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd1y, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd1y, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td1y, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd1z, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd1z, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td1z, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd2x, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd2x, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td2x, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd2y, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd2y, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td2y, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd2z, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd2z, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td2z, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_ufxl, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vfxl, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wfxl, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!-------------------------------First touch - OPS DATS without halos-----------------------------------------
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wmomix, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_difmix, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!------------------------------------First touch - OPS DATS Boundary YZ--------------------------------------
    rangexyz = [1,1,1,nyglbl,1,nzglbl]

    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strexl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t1bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t2bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t3bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t4bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t51bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t52bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]

    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strexr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strrxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t1bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t2bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t3bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t4bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t51bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_t52bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

!--------------------------First touch - OPS DATS Boundary XZ------------------------------------------------
    rangexyz = [1,nxglbl,1,1,1,nzglbl]

    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_streyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))


    rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]

    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_streyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strryr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

!------------------------------First touch - OPS DATS Boundary XY--------------------------------------------
    rangexyz = [1,nxglbl,1,nyglbl,1,1]

    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strezl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))


    rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]

    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strezr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strrzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

!--------------------------------First touch - MULTI-DIM DAT-------------------------------------------------
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_rrte(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))

    END DO

    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    DO iindex = 1,nintmx
        call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
    END DO

    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_ctrans(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
    END DO

    rangexyz = [1,1,1,nyglbl,1,nzglbl]
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_ratexl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    END DO

    rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_ratexr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    END DO

    rangexyz = [1,nxglbl,1,1,1,nzglbl]
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_rateyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    END DO

    rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_rateyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    END DO

    rangexyz = [1,nxglbl,1,nyglbl,1,1]
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_ratezl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    END DO

    rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_ratezr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    END DO

!------------------------------------------------------------------------------------------------------------

#ifdef OPS_WITH_CUDAFOR

    ncofmx_opsconstant = ncofmx
    ntinmx_opsconstant = ntinmx
    nspcmx_opsconstant = nspcmx
    nssmax_opsconstant = nssmax
    nstpmx_opsconstant = nstpmx
    ndcfmx_opsconstant = ndcfmx
    nvcfmx_opsconstant = nvcfmx
    nccfmx_opsconstant = nccfmx
    nrkmax_opsconstant = nrkmax
    ncbcsz_opsconstant = ncbcsz
    nbcprr_opsconstant = nbcprr
    nspimx_opsconstant = nspimx
    ntbase_opsconstant = ntbase
    nintmx_opsconstant = nintmx
    nctmax_opsconstant = nctmax
    nctmm1_opsconstant = nctmm1
    nrsmax_opsconstant = nrsmax
    nbcpri_opsconstant = nbcpri
    ncfrmx_opsconstant = ncfrmx

#endif

    call ops_decl_const("ncofmx", 1, "integer(kind=4)", ncofmx)
    call ops_decl_const("ntinmx", 1, "integer(kind=4)", ntinmx)
    call ops_decl_const("nspcmx", 1, "integer(kind=4)", nspcmx)
    call ops_decl_const("nssmax", 1, "integer(kind=4)", nssmax)
    call ops_decl_const("nstpmx", 1, "integer(kind=4)", nstpmx)
    call ops_decl_const("ndcfmx", 1, "integer(kind=4)", ndcfmx)
    call ops_decl_const("nvcfmx", 1, "integer(kind=4)", nvcfmx)
    call ops_decl_const("nccfmx", 1, "integer(kind=4)", nccfmx)
    call ops_decl_const("nrkmax", 1, "integer(kind=4)", nrkmax)
    call ops_decl_const("ncbcsz", 1, "integer(kind=4)", ncbcsz)
    call ops_decl_const("nbcprr", 1, "integer(kind=4)", nbcprr)
    call ops_decl_const("nspimx", 1, "integer(kind=4)", nspimx)
    call ops_decl_const("ntbase", 1, "integer(kind=4)", ntbase)
    call ops_decl_const("nintmx", 1, "integer(kind=4)", nintmx)
    call ops_decl_const("nctmax", 1, "integer(kind=4)", nctmax)
    call ops_decl_const("nctmm1", 1, "integer(kind=4)", nctmm1)
    call ops_decl_const("nrsmax", 1, "integer(kind=4)", nrsmax)
    call ops_decl_const("nbcpri", 1, "integer(kind=4)", nbcpri)
    call ops_decl_const("ncfrmx", 1, "integer(kind=4)", ncfrmx)

END SUBROUTINE ops_data_init
