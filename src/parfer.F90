SUBROUTINE parfer

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
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
    IF(nendxl == nperi) THEN

        call ops_halo_transfer(halos_grp_x)

    END IF


!   Y-DIRECTION
!   ONLY NEED TO CHECK ONE END
!   NOTE EXTENDED X-LIMITS FOR Y TRANSFERS
    IF(nendyl == nperi)THEN

        IF(nhaloy /= 0) call ops_halo_transfer(halos_grp_y)

    END IF

!   Z-DIRECTION
!   ONLY NEED TO CHECK ONE END
!   NOTE EXTENDED X- AND Y-LIMITS FOR Z TRANSFERS
    IF(nendzl == nperi) THEN

        IF(nhaloz /= 0) call ops_halo_transfer(halos_grp_z)

    END IF

!   =========================================================================

END SUBROUTINE parfer
