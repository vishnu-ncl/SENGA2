MODULE com_ops_senga

    use OPS_Fortran_Reference
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    implicit none

!------------------------------------------------------------------------------------------------------------
!   Senga2 OPS vars

!------------------------------------------------------------------------------------------------------------
!   OPS block
    TYPE(ops_block) :: senga_grid

!------------------------------------------------------------------------------------------------------------
!   OPS dats
    TYPE(ops_dat) :: d_store1, d_store2, d_store3, d_store4, d_store5, d_store6, d_divm

    TYPE(ops_dat) :: d_drun, d_urun, d_vrun, d_wrun, d_erun
    TYPE(ops_dat) :: d_drhs, d_urhs, d_vrhs, d_wrhs, d_erhs


    TYPE(ops_dat) :: d_utmp, d_vtmp, d_wtmp, d_prun, d_trun, d_transp, d_store7

END MODULE com_ops_senga
