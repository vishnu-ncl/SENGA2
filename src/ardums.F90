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

    real(kind=8), intent(IN OUT)       :: array(nxphys,nyphys,nzphys,nspec)
    integer, intent(IN OUT)             :: nxphys
    integer, intent(IN)                 :: nyphys
    integer, intent(IN)                 :: nzphys
    integer, intent(IN OUT)             :: nspec
    integer, intent(IN OUT)             :: ispec

!   LOCAL DATA
!   ==========
    integer :: ic,jc,kc

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
