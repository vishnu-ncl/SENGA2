SUBROUTINE ardums(array,nxphys,nyphys,nzphys,nspec,ispec)
 
!   *************************************************************************

!   ARDUMS
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   11-SEP-2005:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   DIAGNOSTIC ROUTINE
!   PRINTS THE CONTENTS OF THE SPECIFIED SPECIES ARRAY

!   *************************************************************************
    use data_types

!   ARGUMENTS
!   =========

    real(8), intent(IN OUT)       :: array(nxphys,nyphys,nzphys,nspec)
    integer(4), intent(IN OUT)             :: nxphys
    integer(4), intent(IN)                 :: nyphys
    integer(4), intent(IN)                 :: nzphys
    integer(4), intent(IN OUT)             :: nspec
    integer(4), intent(IN OUT)             :: ispec

!   LOCAL DATA
!   ==========
    integer(4) :: ic,jc,kc

!   BEGIN
!   =====

    WRITE(6,*)

    DO kc = nzphys,1,-1
        WRITE(6,'(" Z-PLANE ",I10)')kc
        DO jc = nyphys,1,-1
            WRITE(6,'(I5,21(1PE10.3))') jc,(array(ic,jc,kc,ispec),ic = 1,nxphys)
        END DO
    END DO

!   ================================================================================

END SUBROUTINE ardums
