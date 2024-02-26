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

!   =========================================================================

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
