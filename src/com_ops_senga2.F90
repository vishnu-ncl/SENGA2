MODULE com_ops_senga

    use OPS_Fortran_Reference
    
    implicit none

    !Senga2 OPS vars

    TYPE(ops_block) :: senga_grid

    TYPE(ops_dat) :: d_store1, d_store2, d_store3, d_store4, d_store5, d_store6
    
    TYPE(ops_dat) :: d_wd1x, d_pd1x, d_td1x
    TYPE(ops_dat) :: d_wd1y, d_pd1y, d_td1y
    TYPE(ops_dat) :: d_wd1z, d_pd1z, d_td1z

    TYPE(ops_dat) :: d_wd2x, d_pd2x, d_td2x
    TYPE(ops_dat) :: d_wd2y, d_pd2y, d_td2y
    TYPE(ops_dat) :: d_wd2z, d_pd2z, d_td2z

    TYPE(ops_dat) :: d_urhs, d_drhs, d_erhs, d_wrhs, d_vrhs
    TYPE(ops_dat) :: d_utmp, d_vtmp, d_wtmp, d_trun, d_transp, d_store7

    TYPE(ops_stencil) :: s3d_000
    TYPE(ops_stencil) :: s3d_000_to_p400_x, s3d_000_to_m400_x
    TYPE(ops_stencil) :: s3d_000_to_p500_x, s3d_000_to_m500_x
    TYPE(ops_stencil) :: s3d_p200_to_m200_x, s3d_p300_to_m300_x, s3d_p400_to_m400_x, s3d_p500_to_m500_x
    TYPE(ops_stencil) :: s3d_p300_to_m100_x, s3d_p100_to_m300_x
    TYPE(ops_stencil) :: s3d_p400_to_m100_x, s3d_p100_to_m400_x

    TYPE(ops_stencil) :: s3d_000_to_p040_y, s3d_000_to_m040_y
    TYPE(ops_stencil) :: s3d_000_to_p050_y, s3d_000_to_m050_y
    TYPE(ops_stencil) :: s3d_p020_to_m020_y, s3d_p030_to_m030_y, s3d_p040_to_m040_y, s3d_p050_to_m050_y
    TYPE(ops_stencil) :: s3d_p030_to_m010_y, s3d_p010_to_m030_y
    TYPE(ops_stencil) :: s3d_p040_to_m010_y, s3d_p010_to_m040_y

    TYPE(ops_stencil) :: s3d_000_to_p004_z, s3d_000_to_m004_z
    TYPE(ops_stencil) :: s3d_000_to_p005_z, s3d_000_to_m005_z
    TYPE(ops_stencil) :: s3d_p002_to_m002_z, s3d_p003_to_m003_z, s3d_p004_to_m004_z, s3d_p005_to_m005_z
    TYPE(ops_stencil) :: s3d_p003_to_m001_z, s3d_p001_to_m003_z
    TYPE(ops_stencil) :: s3d_p004_to_m001_z, s3d_p001_to_m004_z

END MODULE com_ops_senga
