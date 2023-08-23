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

!   ARGUMENTS
!   =========

    integer(kind=4), intent(IN OUT)             :: nxphys
    integer(kind=4), intent(IN)                 :: nyphys
    integer(kind=4), intent(IN)                 :: nzphys
    integer(kind=4), intent(IN OUT)             :: nspec
    integer(kind=4), intent(IN OUT)             :: ispec
    real(kind=8), intent(IN OUT)       :: array(nxphys,nyphys,nzphys,nspec)

!   LOCAL DATA
!   ==========
    integer(kind=4) :: ic,jc,kc

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
