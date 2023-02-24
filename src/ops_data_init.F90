SUBROUTINE ops_data_init()
    use OPS_Fortran_Reference

    use data_types
    use com_senga
    use com_ops_senga

    INTEGER :: d_size(3)
    INTEGER d_base(3) /1,1,1/ !array indexing - start from 1
    INTEGER :: d_p(3) !max boundary depths for the dat in the possitive direction
    INTEGER :: d_m(3) !max boundary depths for the dat in the negative direction

    real(8), dimension(:), allocatable :: temp_real_null
    integer, dimension(:), allocatable :: temp_int_null

    INTEGER :: ispec
    INTEGER :: rangexyz(6)

    INTEGER :: halo_idx
    INTEGER :: iter_size(3), base_from(3), base_to(3), dir_from(3), dir_to(3)

    INTEGER a3d_000(3)    /0,0,0/

    integer stride3d_x(3) /1,0,0/
    integer stride3d_y(3) /0,1,0/
    integer stride3d_z(3) /0,0,1/
    
    integer stride3d_xy(3) /1,1,0/
    integer stride3d_xz(3) /1,0,1/
    integer stride3d_yz(3) /0,1,1/

!------------------------------------------------------------------------------------------------------------

    INTEGER a3d_000_to_p400_x(15) /0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0/
    INTEGER a3d_000_to_m400_x(15) /0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0/

    INTEGER a3d_000_to_p500_x(18) /0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0/
    INTEGER a3d_000_to_m500_x(18) /0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0, -5,0,0/

    INTEGER a3d_p100_to_p400_x(12) /1,0,0, 2,0,0, 3,0,0, 4,0,0/
    INTEGER a3d_m100_to_m400_x(12) /-1,0,0, -2,0,0, -3,0,0, -4,0,0/

    INTEGER a3d_p200_to_m200_x_no000(12) /2,0,0, 1,0,0, -1,0,0, -2,0,0/
    INTEGER a3d_p300_to_m300_x_no000(18) /3,0,0, 2,0,0, 1,0,0, -1,0,0, -2,0,0, -3,0,0/
    INTEGER a3d_p400_to_m400_x_no000(24) /4,0,0, 3,0,0, 2,0,0, 1,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0/
    INTEGER a3d_p500_to_m500_x_no000(30) /5,0,0, 4,0,0, 3,0,0, 2,0,0, 1,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0, -5,0,0/

    INTEGER a3d_p200_to_m200_x(15) /2,0,0, 1,0,0, 0,0,0, -1,0,0, -2,0,0/
    INTEGER a3d_p300_to_m300_x(21) /3,0,0, 2,0,0, 1,0,0, 0,0,0, -1,0,0, -2,0,0, -3,0,0/
    INTEGER a3d_p400_to_m400_x(27) /4,0,0, 3,0,0, 2,0,0, 1,0,0, 0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0/
    INTEGER a3d_p500_to_m500_x(33) /5,0,0, 4,0,0, 3,0,0, 2,0,0, 1,0,0, 0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0, -5,0,0/

    INTEGER a3d_p300_to_m100_x(15) /3,0,0, 2,0,0, 1,0,0, 0,0,0, -1,0,0/
    INTEGER a3d_p100_to_m300_x(15) /1,0,0, 0,0,0, -1,0,0, -2,0,0, -3,0,0/

    INTEGER a3d_p400_to_m100_x(18) /4,0,0, 3,0,0, 2,0,0, 1,0,0, 0,0,0, -1,0,0/
    INTEGER a3d_p100_to_m400_x(18) /1,0,0, 0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0/
!------------------------------------------------------------------------------------------------------------
    INTEGER a3d_000_to_p040_y(15) /0,0,0, 0,1,0, 0,2,0, 0,3,0, 0,4,0/
    INTEGER a3d_000_to_m040_y(15) /0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0/

    INTEGER a3d_000_to_p050_y(18) /0,0,0, 0,1,0, 0,2,0, 0,3,0, 0,4,0, 0,5,0/
    INTEGER a3d_000_to_m050_y(18) /0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, 0,-5,0/

    INTEGER a3d_p010_to_p040_y(12) /0,1,0, 0,2,0, 0,3,0, 0,4,0/
    INTEGER a3d_m010_to_m040_y(12) /0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0/

    INTEGER a3d_p020_to_m020_y_no000(12) /0,2,0, 0,1,0, 0,-1,0, 0,-2,0/
    INTEGER a3d_p030_to_m030_y_no000(18) /0,3,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0, 0,-3,0/
    INTEGER a3d_p040_to_m040_y_no000(24) /0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0/
    INTEGER a3d_p050_to_m050_y_no000(30) /0,5,0, 0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, 0,-5,0/

    INTEGER a3d_p020_to_m020_y(15) /0,2,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0/
    INTEGER a3d_p030_to_m030_y(21) /0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0/
    INTEGER a3d_p040_to_m040_y(27) /0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0/
    INTEGER a3d_p050_to_m050_y(33) /0,5,0, 0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, 0,-5,0/

    INTEGER a3d_p030_to_m010_y(15) /0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0/
    INTEGER a3d_p010_to_m030_y(15) /0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0/

    INTEGER a3d_p040_to_m010_y(18) /0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0/
    INTEGER a3d_p010_to_m040_y(18) /0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0/
!------------------------------------------------------------------------------------------------------------
    INTEGER a3d_000_to_p004_z(15) /0,0,0, 0,0,1, 0,0,2, 0,0,3, 0,0,4/
    INTEGER a3d_000_to_m004_z(15) /0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4/

    INTEGER a3d_000_to_p005_z(18) /0,0,0, 0,0,1, 0,0,2, 0,0,3, 0,0,4, 0,0,5/
    INTEGER a3d_000_to_m005_z(18) /0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, 0,0,-5/

    INTEGER a3d_p001_to_p004_z(12) /0,0,1, 0,0,2, 0,0,3, 0,0,4/
    INTEGER a3d_m001_to_m004_z(12) /0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4/

    INTEGER a3d_p002_to_m002_z_no000(12) /0,0,2, 0,0,1, 0,0,-1, 0,0,-2/
    INTEGER a3d_p003_to_m003_z_no000(18) /0,0,3, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2, 0,0,-3/
    INTEGER a3d_p004_to_m004_z_no000(24) /0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4/
    INTEGER a3d_p005_to_m005_z_no000(30) /0,0,5, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, 0,0,-5/

    INTEGER a3d_p002_to_m002_z(15) /0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2/
    INTEGER a3d_p003_to_m003_z(21) /0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3/
    INTEGER a3d_p004_to_m004_z(27) /0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4/
    INTEGER a3d_p005_to_m005_z(33) /0,0,5, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, 0,0,-5/

    INTEGER a3d_p003_to_m001_z(15) /0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1/
    INTEGER a3d_p001_to_m003_z(15) /0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3/

    INTEGER a3d_p004_to_m001_z(18) /0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1/
    INTEGER a3d_p001_to_m004_z(18) /0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4/

!------------------------------------------------------------------------------------------------------------

    INTEGER a3d_p000_m440_mixed_xy(39) /0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, -1,0,0, -1,-1,0, -2,0,0, -2,-2,0, -3,0,0, -3,-3,0, -4,0,0, -4,-4,0/
    INTEGER a3d_p010_m430_mixed_xy(75) /0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, -1,1,0, -1,0,0, -1,-1,0, -1,-2,0, -1,-3,0, -2,1,0, -2,0,0, -2,-1,0, -2,-2,0, -2,-3,0, -3,1,0, -3,0,0, -3,-1,0, -3,-2,0, -3,-3,0, -4,1,0, -4,0,0, -4,-1,0, -4,-2,0, -4,-3,0/
    INTEGER a3d_p020_m420_mixed_xy(60) /0,2,0, 0,1,0, 0,-1,0, 0,-2,0, -1,2,0, -1,1,0, -1,-1,0, -1,-2,0, -2,2,0, -2,1,0, -2,-1,0, -2,-2,0, -3,2,0, -3,1,0, -3,-1,0, -3,-2,0, -4,2,0, -4,1,0, -4,-1,0, -4,-2,0/
    INTEGER a3d_p030_m410_mixed_xy(75) /0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, -1,3,0, -1,2,0, -1,1,0, -1,0,0, -1,-1,0, -2,3,0, -2,2,0, -2,1,0, -2,0,0, -2,-1,0, -3,3,0, -3,2,0, -3,1,0, -3,0,0, -3,-1,0, -4,3,0, -4,2,0, -4,1,0, -4,0,0, -4,-1,0/
    INTEGER a3d_p040_m400_mixed_xy(39) /0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, -1,1,0, -1,0,0, -2,2,0, -2,0,0, -3,3,0, -3,0,0, -4,4,0, -4,0,0/

    INTEGER a3d_p100_m300_mixed_xy(39) /1,0,0, 1,-1,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, -1,1,0, -1,0,0, -2,2,0, -2,0,0, -3,3,0, -3,0,0/
    INTEGER a3d_p100_m340_mixed_xy(75) /1,0,0, 1,-1,0, 1,-2,0, 1,-3,0, 1,-4,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, -1,0,0, -1,-1,0, -1,-2,0, -1,-3,0, -1,-4,0, -2,0,0, -2,-1,0, -2,-2,0, -2,-3,0, -2,-4,0, -3,0,0, -3,-1,0, -3,-2,0, -3,-3,0, -3,-4,0/
    INTEGER a3d_p110_m330_mixed_xy(39) /1,1,0, 1,0,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, -1,0,0, -1,-1,0, -2,0,0, -2,-2,0, -3,0,0, -3,-3,0/
    INTEGER a3d_p120_m320_mixed_xy(60) /1,2,0, 1,1,0, 1,-1,0, 1,-2,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0, -1,2,0, -1,1,0, -1,-1,0, -1,-2,0, -2,2,0, -2,1,0, -2,-1,0, -2,-2,0, -3,2,0, -3,1,0, -3,-1,0, -3,-2,0/
    INTEGER a3d_p140_m300_mixed_xy(75) /1,4,0, 1,3,0, 1,2,0, 1,1,0, 1,0,0, 0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, -1,4,0, -1,3,0, -1,2,0, -1,1,0, -1,0,0, -2,4,0, -2,3,0, -2,2,0, -2,1,0, -2,0,0, -3,4,0, -3,3,0, -3,2,0, -3,1,0, -3,0,0/

    INTEGER a3d_p200_m240_mixed_xy(60) /2,0,0, 2,-1,0, 2,-2,0, 2,-3,0, 2,-4,0, 1,0,0, 1,-1,0, 1,-2,0, 1,-3,0, 1,-4,0, -1,0,0, -1,-1,0, -1,-2,0, -1,-3,0, -1,-4,0, -2,0,0, -2,-1,0, -2,-2,0, -2,-3,0, -2,-4,0/
    INTEGER a3d_p210_m230_mixed_xy(60) /2,1,0, 2,0,0, 2,-1,0, 2,-2,0, 2,-3,0, 1,1,0, 1,0,0, 1,-1,0, 1,-2,0, 1,-3,0, -1,1,0, -1,0,0, -1,-1,0, -1,-2,0, -1,-3,0, -2,1,0, -2,0,0, -2,-1,0, -2,-2,0, -2,-3,0/
    INTEGER a3d_p220_m220_mixed_xy(24) /2,2,0, 2,-2,0, 1,1,0, 1,-1,0, -1,1,0, -1,-1,0, -2,2,0, -2,-2,0/
    INTEGER a3d_p230_m210_mixed_xy(60) /2,3,0, 2,2,0, 2,1,0, 2,0,0, 2,-1,0, 1,3,0, 1,2,0, 1,1,0, 1,0,0, 1,-1,0, -1,3,0, -1,2,0, -1,1,0, -1,0,0, -1,-1,0, -2,3,0, -2,2,0, -2,1,0, -2,0,0, -2,-1,0/
    INTEGER a3d_p240_m200_mixed_xy(60) /2,4,0, 2,3,0, 2,2,0, 2,1,0, 2,0,0, 1,4,0, 1,3,0, 1,2,0, 1,1,0, 1,0,0, -1,4,0, -1,3,0, -1,2,0, -1,1,0, -1,0,0, -2,4,0, -2,3,0, -2,2,0, -2,1,0, -2,0,0/

    INTEGER a3d_p300_m100_mixed_xy(39) /3,0,0, 3,-3,0, 2,0,0, 2,-2,0, 1,0,0, 1,-1,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, -1,1,0, -1,0,0/
    INTEGER a3d_p300_m140_mixed_xy(75) /3,0,0, 3,-1,0, 3,-2,0, 3,-3,0, 3,-4,0, 2,0,0, 2,-1,0, 2,-2,0, 2,-3,0, 2,-4,0, 1,0,0, 1,-1,0, 1,-2,0, 1,-3,0, 1,-4,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, -1,0,0, -1,-1,0, -1,-2,0, -1,-3,0, -1,-4,0/
    INTEGER a3d_p320_m120_mixed_xy(60) /3,2,0, 3,1,0, 3,-1,0, 3,-2,0, 2,2,0, 2,1,0, 2,-1,0, 2,-2,0, 1,2,0, 1,1,0, 1,-1,0, 1,-2,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0, -1,2,0, -1,1,0, -1,-1,0, -1,-2,0/
    INTEGER a3d_p330_m110_mixed_xy(39) /3,3,0, 3,0,0, 2,2,0, 2,0,0, 1,1,0, 1,0,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, -1,0,0, -1,-1,0/
    INTEGER a3d_p330_m330_mixed_xy(36) /3,3,0, 3,-3,0, 2,2,0, 2,-2,0, 1,1,0, 1,-1,0, -1,1,0, -1,-1,0, -2,2,0, -2,-2,0, -3,3,0, -3,-3,0/
    INTEGER a3d_p340_m100_mixed_xy(75) /3,4,0, 3,3,0, 3,2,0, 3,1,0, 3,0,0, 2,4,0, 2,3,0, 2,2,0, 2,1,0, 2,0,0, 1,4,0, 1,3,0, 1,2,0, 1,1,0, 1,0,0, 0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, -1,4,0, -1,3,0, -1,2,0, -1,1,0, -1,0,0/

    INTEGER a3d_p400_p040_mixed_xy(39) /4,0,0, 4,-4,0, 3,0,0, 3,-3,0, 2,0,0, 2,-2,0, 1,0,0, 1,-1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0/
    INTEGER a3d_p410_p030_mixed_xy(75) /4,1,0, 4,0,0, 4,-1,0, 4,-2,0, 4,-3,0, 3,1,0, 3,0,0, 3,-1,0, 3,-2,0, 3,-3,0, 2,1,0, 2,0,0, 2,-1,0, 2,-2,0, 2,-3,0, 1,1,0, 1,0,0, 1,-1,0, 1,-2,0, 1,-3,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0/
    INTEGER a3d_p420_m020_mixed_xy(60) /4,2,0, 4,1,0, 4,-1,0, 4,-2,0, 3,2,0, 3,1,0, 3,-1,0, 3,-2,0, 2,2,0, 2,1,0, 2,-1,0, 2,-2,0, 1,2,0, 1,1,0, 1,-1,0, 1,-2,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0/
    INTEGER a3d_p430_m010_mixed_xy(75) /4,3,0, 4,2,0, 4,1,0, 4,0,0, 4,-1,0, 3,3,0, 3,2,0, 3,1,0, 3,0,0, 3,-1,0, 2,3,0, 2,2,0, 2,1,0, 2,0,0, 2,-1,0, 1,3,0, 1,2,0, 1,1,0, 1,0,0, 1,-1,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0/
    INTEGER a3d_p440_p000_mixed_xy(39) /4,4,0, 4,0,0, 3,3,0, 3,0,0, 2,2,0, 2,0,0, 1,1,0, 1,0,0, 0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0/
    INTEGER a3d_p440_m440_mixed_xy(48) /4,4,0, 4,-4,0, 3,3,0, 3,-3,0, 2,2,0, 2,-2,0, 1,1,0, 1,-1,0, -1,1,0, -1,-1,0, -2,2,0, -2,-2,0, -3,3,0, -3,-3,0, -4,4,0, -4,-4,0/

    INTEGER a3d_p550_m550_mixed_xy(60) /5,5,0, 5,-5,0, 4,4,0, 4,-4,0, 3,3,0, 3,-3,0, 2,2,0, 2,-2,0, 1,1,0, 1,-1,0, -1,1,0, -1,-1,0, -2,2,0, -2,-2,0, -3,3,0, -3,-3,0, -4,4,0, -4,-4,0, -5,5,0, -5,-5,0/

!------------------------------------------------------------------------------------------------------------

    INTEGER a3d_p000_m404_mixed_xz(39) /0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, -1,0,0, -1,0,-1, -2,0,0, -2,0,-2, -3,0,0, -3,0,-3, -4,0,0, -4,0,-4/
    INTEGER a3d_p001_m403_mixed_xz(75) /0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, -1,0,1, -1,0,0, -1,0,-1, -1,0,-2, -1,0,-3, -2,0,1, -2,0,0, -2,0,-1, -2,0,-2, -2,0,-3, -3,0,1, -3,0,0, -3,0,-1, -3,0,-2, -3,0,-3, -4,0,1, -4,0,0, -4,0,-1, -4,0,-2, -4,0,-3/
    INTEGER a3d_p002_m402_mixed_xz(60) /0,0,2, 0,0,1, 0,0,-1, 0,0,-2, -1,0,2, -1,0,1, -1,0,-1, -1,0,-2, -2,0,2, -2,0,1, -2,0,-1, -2,0,-2, -3,0,2, -3,0,1, -3,0,-1, -3,0,-2, -4,0,2, -4,0,1, -4,0,-1, -4,0,-2/
    INTEGER a3d_p003_m401_mixed_xz(75) /0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, -1,0,3, -1,0,2, -1,0,1, -1,0,0, -1,0,-1, -2,0,3, -2,0,2, -2,0,1, -2,0,0, -2,0,-1, -3,0,3, -3,0,2, -3,0,1, -3,0,0, -3,0,-1, -4,0,3, -4,0,2, -4,0,1, -4,0,0, -4,0,-1/
    INTEGER a3d_p004_m400_mixed_xz(39) /0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, -1,0,1, -1,0,0, -2,0,2, -2,0,0, -3,0,3, -3,0,0, -4,0,4, -4,0,0/

    INTEGER a3d_p100_m300_mixed_xz(39) /1,0,0, 1,0,-1, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, -1,0,1, -1,0,0, -2,0,2, -2,0,0, -3,0,3, -3,0,0/
    INTEGER a3d_p100_m304_mixed_xz(75) /1,0,0, 1,0,-1, 1,0,-2, 1,0,-3, 1,0,-4, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, -1,0,0, -1,0,-1, -1,0,-2, -1,0,-3, -1,0,-4, -2,0,0, -2,0,-1, -2,0,-2, -2,0,-3, -2,0,-4, -3,0,0, -3,0,-1, -3,0,-2, -3,0,-3, -3,0,-4/
    INTEGER a3d_p101_m303_mixed_xz(39) /1,0,1, 1,0,0, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, -1,0,0, -1,0,-1, -2,0,0, -2,0,-2, -3,0,0, -3,0,-3/
    INTEGER a3d_p102_m302_mixed_xz(60) /1,0,2, 1,0,1, 1,0,-1, 1,0,-2, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2, -1,0,2, -1,0,1, -1,0,-1, -1,0,-2, -2,0,2, -2,0,1, -2,0,-1, -2,0,-2, -3,0,2, -3,0,1, -3,0,-1, -3,0,-2/
    INTEGER a3d_p104_m300_mixed_xz(75) /1,0,4, 1,0,3, 1,0,2, 1,0,1, 1,0,0, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, -1,0,4, -1,0,3, -1,0,2, -1,0,1, -1,0,0, -2,0,4, -2,0,3, -2,0,2, -2,0,1, -2,0,0, -3,0,4, -3,0,3, -3,0,2, -3,0,1, -3,0,0/

    INTEGER a3d_p200_m204_mixed_xz(60) /2,0,0, 2,0,-1, 2,0,-2, 2,0,-3, 2,0,-4, 1,0,0, 1,0,-1, 1,0,-2, 1,0,-3, 1,0,-4, -1,0,0, -1,0,-1, -1,0,-2, -1,0,-3, -1,0,-4, -2,0,0, -2,0,-1, -2,0,-2, -2,0,-3, -2,0,-4/
    INTEGER a3d_p201_m203_mixed_xz(60) /2,0,1, 2,0,0, 2,0,-1, 2,0,-2, 2,0,-3, 1,0,1, 1,0,0, 1,0,-1, 1,0,-2, 1,0,-3, -1,0,1, -1,0,0, -1,0,-1, -1,0,-2, -1,0,-3, -2,0,1, -2,0,0, -2,0,-1, -2,0,-2, -2,0,-3/
    INTEGER a3d_p202_m202_mixed_xz(24) /2,0,2, 2,0,-2, 1,0,1, 1,0,-1, -1,0,1, -1,0,-1, -2,0,2, -2,0,-2/
    INTEGER a3d_p203_m201_mixed_xz(60) /2,0,3, 2,0,2, 2,0,1, 2,0,0, 2,0,-1, 1,0,3, 1,0,2, 1,0,1, 1,0,0, 1,0,-1, -1,0,3, -1,0,2, -1,0,1, -1,0,0, -1,0,-1, -2,0,3, -2,0,2, -2,0,1, -2,0,0, -2,0,-1/
    INTEGER a3d_p204_m200_mixed_xz(60) /2,0,4, 2,0,3, 2,0,2, 2,0,1, 2,0,0, 1,0,4, 1,0,3, 1,0,2, 1,0,1, 1,0,0, -1,0,4, -1,0,3, -1,0,2, -1,0,1, -1,0,0, -2,0,4, -2,0,3, -2,0,2, -2,0,1, -2,0,0/

    INTEGER a3d_p300_m100_mixed_xz(39) /3,0,0, 3,0,-3, 2,0,0, 2,0,-2, 1,0,0, 1,0,-1, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, -1,0,1, -1,0,0/
    INTEGER a3d_p300_m104_mixed_xz(75) /3,0,0, 3,0,-1, 3,0,-2, 3,0,-3, 3,0,-4, 2,0,0, 2,0,-1, 2,0,-2, 2,0,-3, 2,0,-4, 1,0,0, 1,0,-1, 1,0,-2, 1,0,-3, 1,0,-4, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, -1,0,0, -1,0,-1, -1,0,-2, -1,0,-3, -1,0,-4/
    INTEGER a3d_p302_m102_mixed_xz(60) /3,0,2, 3,0,1, 3,0,-1, 3,0,-2, 2,0,2, 2,0,1, 2,0,-1, 2,0,-2, 1,0,2, 1,0,1, 1,0,-1, 1,0,-2, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2, -1,0,2, -1,0,1, -1,0,-1, -1,0,-2/
    INTEGER a3d_p303_m101_mixed_xz(39) /3,0,3, 3,0,0, 2,0,2, 2,0,0, 1,0,1, 1,0,0, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, -1,0,0, -1,0,-1/
    INTEGER a3d_p303_m303_mixed_xz(36) /3,0,3, 3,0,-3, 2,0,2, 2,0,-2, 1,0,1, 1,0,-1, -1,0,1, -1,0,-1, -2,0,2, -2,0,-2, -3,0,3, -3,0,-3/
    INTEGER a3d_p304_m100_mixed_xz(75) /3,0,4, 3,0,3, 3,0,2, 3,0,1, 3,0,0, 2,0,4, 2,0,3, 2,0,2, 2,0,1, 2,0,0, 1,0,4, 1,0,3, 1,0,2, 1,0,1, 1,0,0, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, -1,0,4, -1,0,3, -1,0,2, -1,0,1, -1,0,0/

    INTEGER a3d_p400_p004_mixed_xz(39) /4,0,0, 4,0,-4, 3,0,0, 3,0,-3, 2,0,0, 2,0,-2, 1,0,0, 1,0,-1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4/
    INTEGER a3d_p401_p003_mixed_xz(75) /4,0,1, 4,0,0, 4,0,-1, 4,0,-2, 4,0,-3, 3,0,1, 3,0,0, 3,0,-1, 3,0,-2, 3,0,-3, 2,0,1, 2,0,0, 2,0,-1, 2,0,-2, 2,0,-3, 1,0,1, 1,0,0, 1,0,-1, 1,0,-2, 1,0,-3, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3/
    INTEGER a3d_p402_m002_mixed_xz(60) /4,0,2, 4,0,1, 4,0,-1, 4,0,-2, 3,0,2, 3,0,1, 3,0,-1, 3,0,-2, 2,0,2, 2,0,1, 2,0,-1, 2,0,-2, 1,0,2, 1,0,1, 1,0,-1, 1,0,-2, 0,0,2, 0,0,1, 0,0,-1, 0,0,-2/
    INTEGER a3d_p403_m001_mixed_xz(75) /4,0,3, 4,0,2, 4,0,1, 4,0,0, 4,0,-1, 3,0,3, 3,0,2, 3,0,1, 3,0,0, 3,0,-1, 2,0,3, 2,0,2, 2,0,1, 2,0,0, 2,0,-1, 1,0,3, 1,0,2, 1,0,1, 1,0,0, 1,0,-1, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1/
    INTEGER a3d_p404_p000_mixed_xz(39) /4,0,4, 4,0,0, 3,0,3, 3,0,0, 2,0,2, 2,0,0, 1,0,1, 1,0,0, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0/
    INTEGER a3d_p404_m404_mixed_xz(48) /4,0,4, 4,0,-4, 3,0,3, 3,0,-3, 2,0,2, 2,0,-2, 1,0,1, 1,0,-1, -1,0,1, -1,0,-1, -2,0,2, -2,0,-2, -3,0,3, -3,0,-3, -4,0,4, -4,0,-4/

    INTEGER a3d_p505_m505_mixed_xz(60) /5,0,5, 5,0,-5, 4,0,4, 4,0,-4, 3,0,3, 3,0,-3, 2,0,2, 2,0,-2, 1,0,1, 1,0,-1, -1,0,1, -1,0,-1, -2,0,2, -2,0,-2, -3,0,3, -3,0,-3, -4,0,4, -4,0,-4, -5,0,5, -5,0,-5/    

!------------------------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------------------------

!   *-----------------------------------------OPS Declarations-----------------------------------------------*

!   Declare OPS Block
    call ops_decl_block(3, senga_grid, "SENGA_GRID")

!   Declare OPS Dats
    d_size = (/nxglbl, nyglbl, nzglbl/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store1, "real(8)", "STORE1")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store2, "real(8)", "STORE2")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store3, "real(8)", "STORE3")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store4, "real(8)", "STORE4")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store5, "real(8)", "STORE5")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store6, "real(8)", "STORE6")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_divm, "real(8)", "DIVM")
    
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ucor, "real(8)", "UCOR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vcor, "real(8)", "VCOR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wcor, "real(8)", "WCOR")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd1x, "real(8)", "WD1X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd1x, "real(8)", "PD1X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td1x, "real(8)", "TD1X")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd1y, "real(8)", "WD1Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd1y, "real(8)", "PD1Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td1y, "real(8)", "TD1Y")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd1z, "real(8)", "WD1Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd1z, "real(8)", "PD1Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td1z, "real(8)", "TD1Z")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd2x, "real(8)", "WD2X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd2x, "real(8)", "PD2X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td2x, "real(8)", "TD2X")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd2y, "real(8)", "WD2Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd2y, "real(8)", "PD2Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td2y, "real(8)", "TD2Y")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd2z, "real(8)", "WD2Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd2z, "real(8)", "PD2Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td2z, "real(8)", "TD2Z")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ufxl, "real(8)", "UFXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vfxl, "real(8)", "VFXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wfxl, "real(8)", "WFXL")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_drun, "real(8)", "DRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_urun, "real(8)", "URUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vrun, "real(8)", "VRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wrun, "real(8)", "WRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_erun, "real(8)", "ERUN")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_derr, "real(8)", "DERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_uerr, "real(8)", "UERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_verr, "real(8)", "VERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_werr, "real(8)", "WERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_eerr, "real(8)", "EERR")

!---------------------------------------MULTI-DIM DAT--------------------------------------------------------

    d_size = (/nxglbl, nyglbl, nzglbl/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_yrun, "real(8)", "YRUN")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_yerr, "real(8)", "YERR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_rate, "real(8)", "RATE")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_rrte, "real(8)", "RRTE")

    d_size = (/nxglbl, nyglbl, nzglbl/)
    d_m = (/-nhalox,-nhaloy,-nhaloz/)
    d_p = (/nhalox,nhaloy,nhaloz/)
    call ops_decl_dat(senga_grid, nintmx, d_size, d_base, d_m, d_p, temp_int_null, d_itndex, "integer", "ITNDEX")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_yrhs, "real(8)", "YRHS")

    d_size = (/1, nyglbl, nzglbl/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_bclyxl, "real(8)", "BCLYXL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_bclyxr, "real(8)", "BCLYXR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_stryxl, "real(8)", "STRYXL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_stryxr, "real(8)", "STRYXR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_dydtxl, "real(8)", "DYDTXL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_dydtxr, "real(8)", "DYDTXR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_ratexl, "real(8)", "RATEXL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_ratexr, "real(8)", "RATEXR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_strhxl, "real(8)", "STRHXL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_strhxr, "real(8)", "STRHXR")

    d_size = (/nxglbl, 1, nzglbl/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_bclyyl, "real(8)", "BCLYYL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_bclyyr, "real(8)", "BCLYYR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_stryyl, "real(8)", "STRYYL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_stryyr, "real(8)", "STRYYR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_dydtyl, "real(8)", "DYDTYL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_dydtyr, "real(8)", "DYDTYR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_rateyl, "real(8)", "RATEYL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_rateyr, "real(8)", "RATEYR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_strhyl, "real(8)", "STRHYL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_strhyr, "real(8)", "STRHYR")

    d_size = (/nxglbl, nyglbl, 1/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_bclyzl, "real(8)", "BCLYZL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_bclyzr, "real(8)", "BCLYZR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_stryzl, "real(8)", "STRYZL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_stryzr, "real(8)", "STRYZR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_dydtzl, "real(8)", "DYDTZL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_dydtzr, "real(8)", "DYDTZR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_ratezl, "real(8)", "RATEZL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_ratezr, "real(8)", "RATEZR")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_strhzl, "real(8)", "STRHZL")
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_strhzr, "real(8)", "STRHZR")

!---------------------------------------WITH HALOS-----------------------------------------------------------

    d_size = (/nxglbl, nyglbl, nzglbl/)
    d_m = (/-nhalox,-nhaloy,-nhaloz/)
    d_p = (/nhalox,nhaloy,nhaloz/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_drhs, "real(8)", "DRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_urhs, "real(8)", "URHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vrhs, "real(8)", "VRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wrhs, "real(8)", "WRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_erhs, "real(8)", "ERHS")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_utmp, "real(8)", "UTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vtmp, "real(8)", "VTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wtmp, "real(8)", "WTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_prun, "real(8)", "PRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_trun, "real(8)", "TRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_transp, "real(8)", "TRANSP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store7, "real(8)", "STORE7")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wmomix, "real(8)", "WMOMIX")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_difmix, "real(8)", "DIFMIX")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tdrmix, "real(8)", "TDRMIX")

!-----------------------------------------Boundary YZ--------------------------------------------------------

    d_size = (/1,nyglbl,nzglbl/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1xl, "real(8)", "BCL1XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1xr, "real(8)", "BCL1XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2xl, "real(8)", "BCL2XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2xr, "real(8)", "BCL2XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3xl, "real(8)", "BCL3XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3xr, "real(8)", "BCL3XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4xl, "real(8)", "BCL4XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4xr, "real(8)", "BCL4XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5xl, "real(8)", "BCL5XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5xr, "real(8)", "BCL5XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltxl, "real(8)", "BCLTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltxr, "real(8)", "BCLTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struxl, "real(8)", "STRUXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struxr, "real(8)", "STRUXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvxl, "real(8)", "STRVXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvxr, "real(8)", "STRVXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwxl, "real(8)", "STRWXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwxr, "real(8)", "STRWXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpxl, "real(8)", "STRPXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpxr, "real(8)", "STRPXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdxl, "real(8)", "STRDXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdxr, "real(8)", "STRDXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtxl, "real(8)", "STRTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtxr, "real(8)", "STRTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strexl, "real(8)", "STREXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strexr, "real(8)", "STREXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgxl, "real(8)", "STRGXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgxr, "real(8)", "STRGXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrxl, "real(8)", "STRRXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrxr, "real(8)", "STRRXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtxl, "real(8)", "DUDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtxr, "real(8)", "DUDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtxl, "real(8)", "DVDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtxr, "real(8)", "DVDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtxl, "real(8)", "DWDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtxr, "real(8)", "DWDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtxl, "real(8)", "DTDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtxr, "real(8)", "DTDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtxl, "real(8)", "DDDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtxr, "real(8)", "DDDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouxl, "real(8)", "ACOUXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouxr, "real(8)", "ACOUXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2xl, "real(8)", "OVA2XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2xr, "real(8)", "OVA2XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1xl, "real(8)", "GAM1XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1xr, "real(8)", "GAM1XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmxl, "real(8)", "OVGMXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmxr, "real(8)", "OVGMXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtxl, "real(8)", "SYDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtxr, "real(8)", "SYDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpxl, "real(8)", "SORPXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpxr, "real(8)", "SORPXR")

!-----------------------------------------Boundary XZ--------------------------------------------------------

    d_size = (/nxglbl,1,nzglbl/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1yl, "real(8)", "BCL1YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1yr, "real(8)", "BCL1YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2yl, "real(8)", "BCL2YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2yr, "real(8)", "BCL2YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3yl, "real(8)", "BCL3YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3yr, "real(8)", "BCL3YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4yl, "real(8)", "BCL4YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4yr, "real(8)", "BCL4YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5yl, "real(8)", "BCL5YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5yr, "real(8)", "BCL5YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltyl, "real(8)", "BCLTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltyr, "real(8)", "BCLTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struyl, "real(8)", "STRUYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struyr, "real(8)", "STRUYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvyl, "real(8)", "STRVYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvyr, "real(8)", "STRVYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwyl, "real(8)", "STRWYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwyr, "real(8)", "STRWYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpyl, "real(8)", "STRPYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpyr, "real(8)", "STRPYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdyl, "real(8)", "STRDYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdyr, "real(8)", "STRDYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtyl, "real(8)", "STRTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtyr, "real(8)", "STRTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_streyl, "real(8)", "STREYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_streyr, "real(8)", "STREYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgyl, "real(8)", "STRGYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgyr, "real(8)", "STRGYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strryl, "real(8)", "STRRYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strryr, "real(8)", "STRRYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtyl, "real(8)", "DUDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtyr, "real(8)", "DUDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtyl, "real(8)", "DVDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtyr, "real(8)", "DVDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtyl, "real(8)", "DWDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtyr, "real(8)", "DWDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtyl, "real(8)", "DTDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtyr, "real(8)", "DTDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtyl, "real(8)", "DDDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtyr, "real(8)", "DDDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouyl, "real(8)", "ACOUYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouyr, "real(8)", "ACOUYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2yl, "real(8)", "OVA2YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2yr, "real(8)", "OVA2YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1yl, "real(8)", "GAM1YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1yr, "real(8)", "GAM1YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmyl, "real(8)", "OVGMYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmyr, "real(8)", "OVGMYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtyl, "real(8)", "SYDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtyr, "real(8)", "SYDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpyl, "real(8)", "SORPYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpyr, "real(8)", "SORPYR")

!-----------------------------------------Boundary XY--------------------------------------------------------

    d_size = (/nxglbl,nyglbl,1/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1zl, "real(8)", "BCL1ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1zr, "real(8)", "BCL1ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2zl, "real(8)", "BCL2ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2zr, "real(8)", "BCL2ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3zl, "real(8)", "BCL3ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3zr, "real(8)", "BCL3ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4zl, "real(8)", "BCL4ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4zr, "real(8)", "BCL4ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5zl, "real(8)", "BCL5ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5zr, "real(8)", "BCL5ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltzl, "real(8)", "BCLTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltzr, "real(8)", "BCLTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struzl, "real(8)", "STRUZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struzr, "real(8)", "STRUZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvzl, "real(8)", "STRVZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvzr, "real(8)", "STRVZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwzl, "real(8)", "STRWZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwzr, "real(8)", "STRWZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpzl, "real(8)", "STRPZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpzr, "real(8)", "STRPZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdzl, "real(8)", "STRDZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdzr, "real(8)", "STRDZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtzl, "real(8)", "STRTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtzr, "real(8)", "STRTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strezl, "real(8)", "STREZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strezr, "real(8)", "STREZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgzl, "real(8)", "STRGZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgzr, "real(8)", "STRGZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrzl, "real(8)", "STRRZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrzr, "real(8)", "STRRZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtzl, "real(8)", "DUDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtzr, "real(8)", "DUDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtzl, "real(8)", "DVDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtzr, "real(8)", "DVDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtzl, "real(8)", "DWDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtzr, "real(8)", "DWDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtzl, "real(8)", "DTDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtzr, "real(8)", "DTDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtzl, "real(8)", "DDDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtzr, "real(8)", "DDDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouzl, "real(8)", "ACOUZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouzr, "real(8)", "ACOUZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2zl, "real(8)", "OVA2ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2zr, "real(8)", "OVA2ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1zl, "real(8)", "GAM1ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1zr, "real(8)", "GAM1ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmzl, "real(8)", "OVGMZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmzr, "real(8)", "OVGMZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtzl, "real(8)", "SYDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtzr, "real(8)", "SYDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpzl, "real(8)", "SORPZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpzr, "real(8)", "SORPZR")

!------------------------------------Only X-direction--------------------------------------------------------

    d_size = (/nxglbl,1,1/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_crin, "real(8)", "CRIN")

!------------------------------------X=Y=Z=1-----------------------------------------------------------------

    d_size = (/1,1,1/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_prefer, "real(8)", "PREFER")

!------------------------------------OPS Reduction Handles---------------------------------------------------

    call ops_decl_reduction_handle(8, h_erdtot, "real(8)", "erdtot")
    call ops_decl_reduction_handle(8, h_erutot, "real(8)", "erutot")
    call ops_decl_reduction_handle(8, h_ervtot, "real(8)", "ervtot")
    call ops_decl_reduction_handle(8, h_erwtot, "real(8)", "erwtot")
    call ops_decl_reduction_handle(8, h_eretot, "real(8)", "eretot")
    call ops_decl_reduction_handle(8, h_erytot, "real(8)", "erytot")
    call ops_decl_reduction_handle(8, h_tket, "real(8)", "tket")
    call ops_decl_reduction_handle(8, h_ubart, "real(8)", "ubart")
    call ops_decl_reduction_handle(8, h_vbart, "real(8)", "vbart")
    call ops_decl_reduction_handle(8, h_wbart, "real(8)", "wbart")
    call ops_decl_reduction_handle(8, h_uvart, "real(8)", "uvart")
    call ops_decl_reduction_handle(8, h_vvart, "real(8)", "vvart")
    call ops_decl_reduction_handle(8, h_wvart, "real(8)", "wvart")
    call ops_decl_reduction_handle(8, h_prefer, "real(8)", "prefer")

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
    
    call ops_decl_stencil( 3, 6, a3d_000_to_p500_x, s3d_000_to_p500_x, "0,0,0 to 5,0,0")
    call ops_decl_stencil( 3, 6, a3d_000_to_m500_x, s3d_000_to_m500_x, "0,0,0 to -5,0,0")

    call ops_decl_stencil( 3, 4, a3d_p100_to_p400_x, s3d_p100_to_p400_x, "1,0,0 to 4,0,0")
    call ops_decl_stencil( 3, 4, a3d_m100_to_m400_x, s3d_m100_to_m400_x, "-1,0,0 to -4,0,0")

    call ops_decl_stencil( 3,  4, a3d_p200_to_m200_x_no000, s3d_p200_to_m200_x_no000, "2,0,0 to -2,0,0")
    call ops_decl_stencil( 3,  6, a3d_p300_to_m300_x_no000, s3d_p300_to_m300_x_no000, "3,0,0 to -3,0,0")
    call ops_decl_stencil( 3,  8, a3d_p400_to_m400_x_no000, s3d_p400_to_m400_x_no000, "4,0,0 to -4,0,0")
    call ops_decl_stencil( 3, 10, a3d_p500_to_m500_x_no000, s3d_p500_to_m500_x_no000, "5,0,0 to -5,0,0")
    
    call ops_decl_stencil( 3,  5, a3d_p200_to_m200_x, s3d_p200_to_m200_x, "2,0,0 to -2,0,0")
    call ops_decl_stencil( 3,  7, a3d_p300_to_m300_x, s3d_p300_to_m300_x, "3,0,0 to -3,0,0")
    call ops_decl_stencil( 3,  9, a3d_p400_to_m400_x, s3d_p400_to_m400_x, "4,0,0 to -4,0,0")
    call ops_decl_stencil( 3, 11, a3d_p500_to_m500_x, s3d_p500_to_m500_x, "5,0,0 to -5,0,0")

    call ops_decl_stencil( 3,  5, a3d_p300_to_m100_x, s3d_p300_to_m100_x, "3,0,0 to -1,0,0")
    call ops_decl_stencil( 3,  5, a3d_p100_to_m300_x, s3d_p100_to_m300_x, "1,0,0 to -3,0,0")

    call ops_decl_stencil( 3,  6, a3d_p400_to_m100_x, s3d_p400_to_m100_x, "4,0,0 to -1,0,0")
    call ops_decl_stencil( 3,  6, a3d_p100_to_m400_x, s3d_p100_to_m400_x, "1,0,0 to -4,0,0")

!------------------------------------------------------------------------------------------------------------

    call ops_decl_stencil( 3, 5, a3d_000_to_p040_y, s3d_000_to_p040_y, "0,0,0 to 0,4,0")
    call ops_decl_stencil( 3, 5, a3d_000_to_m040_y, s3d_000_to_m040_y, "0,0,0 to  0,-4,0")

    call ops_decl_stencil( 3, 6, a3d_000_to_p050_y, s3d_000_to_p050_y, "0,0,0 to 0,5,0")
    call ops_decl_stencil( 3, 6, a3d_000_to_m050_y, s3d_000_to_m050_y, "0,0,0 to  0,-5,0")

    call ops_decl_stencil( 3, 4, a3d_p010_to_p040_y, s3d_p010_to_p040_y, "0,1,0 to 0,4,0")
    call ops_decl_stencil( 3, 4, a3d_m010_to_m040_y, s3d_m010_to_m040_y, "0,-1,0 to 0,-4,0")

    call ops_decl_stencil( 3,  4, a3d_p020_to_m020_y_no000, s3d_p020_to_m020_y_no000, "0,2,0 to  0,-2,0")
    call ops_decl_stencil( 3,  6, a3d_p030_to_m030_y_no000, s3d_p030_to_m030_y_no000, "0,3,0 to  0,-3,0")
    call ops_decl_stencil( 3,  8, a3d_p040_to_m040_y_no000, s3d_p040_to_m040_y_no000, "0,4,0 to  0,-4,0")
    call ops_decl_stencil( 3, 10, a3d_p050_to_m050_y_no000, s3d_p050_to_m050_y_no000, "0,5,0 to  0,-5,0")

    call ops_decl_stencil( 3,  5, a3d_p020_to_m020_y, s3d_p020_to_m020_y, "0,2,0 to  0,-2,0")
    call ops_decl_stencil( 3,  7, a3d_p030_to_m030_y, s3d_p030_to_m030_y, "0,3,0 to  0,-3,0")
    call ops_decl_stencil( 3,  9, a3d_p040_to_m040_y, s3d_p040_to_m040_y, "0,4,0 to  0,-4,0")
    call ops_decl_stencil( 3, 11, a3d_p050_to_m050_y, s3d_p050_to_m050_y, "0,5,0 to  0,-5,0")

    call ops_decl_stencil( 3,  5, a3d_p030_to_m010_y, s3d_p030_to_m010_y, "0,3,0 to  0,-1,0")
    call ops_decl_stencil( 3,  5, a3d_p010_to_m030_y, s3d_p010_to_m030_y, "0,1,0 to  0,-3,0")

    call ops_decl_stencil( 3,  6, a3d_p040_to_m010_y, s3d_p040_to_m010_y, "0,4,0 to  0,-1,0")
    call ops_decl_stencil( 3,  6, a3d_p010_to_m040_y, s3d_p010_to_m040_y, "0,1,0 to  0,-4,0")

!------------------------------------------------------------------------------------------------------------

    call ops_decl_stencil( 3, 5, a3d_000_to_p004_z, s3d_000_to_p004_z, "0,0,0 to 0,0,4")
    call ops_decl_stencil( 3, 5, a3d_000_to_m004_z, s3d_000_to_m004_z, "0,0,0 to  0,0,-4")

    call ops_decl_stencil( 3, 6, a3d_000_to_p005_z, s3d_000_to_p005_z, "0,0,0 to 0,0,5")
    call ops_decl_stencil( 3, 6, a3d_000_to_m005_z, s3d_000_to_m005_z, "0,0,0 to  0,0,-5")

    call ops_decl_stencil( 3, 4, a3d_p001_to_p004_z, s3d_p001_to_p004_z, "0,0,1 to 0,0,4")
    call ops_decl_stencil( 3, 4, a3d_m001_to_m004_z, s3d_m001_to_m004_z, "0,0,-1 to 0,0,-4")

    call ops_decl_stencil( 3,  4, a3d_p002_to_m002_z_no000, s3d_p002_to_m002_z_no000, "0,0,2 to  0,0,-2")
    call ops_decl_stencil( 3,  6, a3d_p003_to_m003_z_no000, s3d_p003_to_m003_z_no000, "0,0,3 to  0,0,-3")
    call ops_decl_stencil( 3,  8, a3d_p004_to_m004_z_no000, s3d_p004_to_m004_z_no000, "0,0,4 to  0,0,-4")
    call ops_decl_stencil( 3, 10, a3d_p005_to_m005_z_no000, s3d_p005_to_m005_z_no000, "0,0,5 to  0,0,-5")

    call ops_decl_stencil( 3,  5, a3d_p002_to_m002_z, s3d_p002_to_m002_z, "0,0,2 to  0,0,-2")
    call ops_decl_stencil( 3,  7, a3d_p003_to_m003_z, s3d_p003_to_m003_z, "0,0,3 to  0,0,-3")
    call ops_decl_stencil( 3,  9, a3d_p004_to_m004_z, s3d_p004_to_m004_z, "0,0,4 to  0,0,-4")
    call ops_decl_stencil( 3, 11, a3d_p005_to_m005_z, s3d_p005_to_m005_z, "0,0,5 to  0,0,-5")

    call ops_decl_stencil( 3,  5, a3d_p003_to_m001_z, s3d_p003_to_m001_z, "0,0,3 to  0,0,-1")
    call ops_decl_stencil( 3,  5, a3d_p001_to_m003_z, s3d_p001_to_m003_z, "0,0,1 to  0,0,-3")

    call ops_decl_stencil( 3,  6, a3d_p004_to_m001_z, s3d_p004_to_m001_z, "0,0,4 to  0,0,-1")
    call ops_decl_stencil( 3,  6, a3d_p001_to_m004_z, s3d_p001_to_m004_z, "0,0,1 to  0,0,-4")

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
    call ops_decl_stencil( 3, 12, a3d_p330_m330_mixed_xy, s3d_p330_m330_mixed_xy, "3,3,0 to -3,-3,0")
    call ops_decl_stencil( 3, 13, a3d_p330_m110_mixed_xy, s3d_p330_m110_mixed_xy, "3,3,0 to -1,-1,0")
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



!------------------------------------------------------------------------------------------------------------

!---------------------OPS DECL HALO for periodic transfer on single processor--------------------------------

    dir_from = (/1,2,3/)
    dir_to = (/1,2,3/)

!   X-DIRECTION : RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
    iter_size = (/nhalox,nyglbl,nzglbl/)
    base_from = (/1,1,1/)
    base_to = (/nxglbl+1,1,1/)

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

    halo_idx = halo_idx+1
    call ops_decl_halo(d_yrhs, d_yrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

!   X-DIRECTION : LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
    iter_size = (/nhalox,nyglbl,nzglbl/)
    base_from = (/nxglbl-nhalox+1,1,1/)
    base_to = (/1-nhalox,1,1/)

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

    halo_idx = halo_idx+1
    call ops_decl_halo(d_yrhs, d_yrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))


    call ops_decl_halo_group(12, halos_x, halos_grp_x)

!------------------------------------------------------------------------------------------------------------

!   Y-DIRECTION : RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
    iter_size = (/nxglbl+2*nhalox,nhaloy,nzglbl/)
    base_from = (/1-nhalox,1,1/)
    base_to = (/1-nhalox,nyglbl+1,1/)

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

    halo_idx = halo_idx+1
    call ops_decl_halo(d_yrhs, d_yrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

!   Y-DIRECTION : LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
    iter_size = (/nxglbl+2*nhalox,nhaloy,nzglbl/)
    base_from = (/1-nhalox,nyglbl-nhaloy+1,1/)
    base_to = (/1-nhalox,1-nhaloy,1/)

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

    halo_idx = halo_idx+1
    call ops_decl_halo(d_yrhs, d_yrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))


    call ops_decl_halo_group(12, halos_y, halos_grp_y)

!------------------------------------------------------------------------------------------------------------

!   Z-DIRECTION : RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
    iter_size = (/nxglbl+2*nhalox,nyglbl+2*nhaloy,nhaloz/)
    base_from = (/1-nhalox,1-nhaloy,1/)
    base_to = (/1-nhalox,1-nhaloy,nzglbl+1/)

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

    halo_idx = halo_idx+1
    call ops_decl_halo(d_yrhs, d_yrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

!   Z-DIRECTION : LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
    iter_size = (/nxglbl+2*nhalox,nyglbl+2*nhaloy,nhaloz/)
    base_from = (/1-nhalox,1-nhaloy,nzglbl-nhaloz+1/)
    base_to = (/1-nhalox,1-nhaloy,1-nhaloz/)

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

    halo_idx = halo_idx+1
    call ops_decl_halo(d_yrhs, d_yrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))


    call ops_decl_halo_group(12, halos_z, halos_grp_z)

!------------------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------------------
    call ops_partition(" ")
!------------------------------------------------------------------------------------------------------------

!---------------------------------First touch - OPS DATS without halos---------------------------------------
    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_divm, 1, s3d_000, "real(8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_ucor, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vcor, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wcor, 1, s3d_000, "real(8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd1x, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd1x, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td1x, 1, s3d_000, "real(8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd1y, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd1y, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td1y, 1, s3d_000, "real(8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd1z, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd1z, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td1z, 1, s3d_000, "real(8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd2x, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd2x, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td2x, 1, s3d_000, "real(8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd2y, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd2y, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td2y, 1, s3d_000, "real(8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wd2z, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_pd2z, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_td2z, 1, s3d_000, "real(8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_ufxl, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vfxl, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wfxl, 1, s3d_000, "real(8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_drun, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_urun, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vrun, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wrun, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_erun, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_derr, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_uerr, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_verr, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_werr, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_eerr, 1, s3d_000, "real(8)", OPS_WRITE))

!-------------------------------First touch - OPS DATS without halos-----------------------------------------
    rangexyz = (/1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz/)

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_prun, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_transp, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wmomix, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_difmix, 1, s3d_000, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tdrmix, 1, s3d_000, "real(8)", OPS_WRITE))

!------------------------------------First touch - OPS DATS Boundary YZ--------------------------------------
    rangexyz = (/1,1,1,nyglbl,1,nzglbl/)

    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strexl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))

    rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)

    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5xr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strexr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strrxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_WRITE))

!--------------------------First touch - OPS DATS Boundary XZ------------------------------------------------
    rangexyz = (/1,nxglbl,1,1,1,nzglbl/)

    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_streyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))


    rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,nzglbl/)

    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5yr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_streyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strryr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_WRITE))

!------------------------------First touch - OPS DATS Boundary XY--------------------------------------------
    rangexyz = (/1,nxglbl,1,nyglbl,1,1/)

    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strezl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))


    rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)

    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcl5zr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_bcltzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strpzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strezr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strgzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strrzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dudtzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dvdtzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dwdtzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dtdtzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_dddtzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ovgmzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sydtzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_WRITE))

!--------------------------------First touch - MULTI-DIM DAT-------------------------------------------------
    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_MD, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_yrun, 9, s3d_000, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_yerr, 9, s3d_000, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_rate, 9, s3d_000, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_rrte, 9, s3d_000, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))

    END DO

    rangexyz=(/1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz/)
    DO ispec = 1,nintmx
        call ops_par_loop(set_zero_kernel_MD_int, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
    END DO
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_MD, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
    END DO

    rangexyz = (/1,1,1,nyglbl,1,nzglbl/)
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_MD_xdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyxl, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_xdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryxl, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_xdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtxl, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_xdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_ratexl, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_xdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhxl, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
    END DO

    rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_MD_xdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyxr, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_xdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryxr, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_xdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtxr, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_xdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_ratexr, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_xdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhxr, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
    END DO

    rangexyz = (/1,nxglbl,1,1,1,nzglbl/)
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_MD_ydir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyyl, 9, s3d_000_strid3d_xz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_ydir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryyl, 9, s3d_000_strid3d_xz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_ydir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtyl, 9, s3d_000_strid3d_xz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_ydir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_rateyl, 9, s3d_000_strid3d_xz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_ydir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhyl, 9, s3d_000_strid3d_xz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
    END DO

    rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,nzglbl/)
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_MD_ydir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyyr, 9, s3d_000_strid3d_xz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_ydir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryyr, 9, s3d_000_strid3d_xz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_ydir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtyr, 9, s3d_000_strid3d_xz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_ydir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_rateyr, 9, s3d_000_strid3d_xz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_ydir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhyr, 9, s3d_000_strid3d_xz, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
    END DO

    rangexyz = (/1,nxglbl,1,nyglbl,1,1/)
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_MD_zdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyzl, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_zdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryzl, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_zdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtzl, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_zdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_ratezl, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_zdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhzl, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
    END DO

    rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
    DO ispec = 1,nspcmx
        call ops_par_loop(set_zero_kernel_MD_zdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_bclyzr, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_zdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_stryzr, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_zdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_dydtzr, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_zdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_ratezr, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
        call ops_par_loop(set_zero_kernel_MD_zdir, "set_zero_multidim", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_strhzr, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                          ops_arg_gbl(ispec, 1, "integer", OPS_READ))
    END DO

!------------------------------------------------------------------------------------------------------------

END SUBROUTINE ops_data_init
