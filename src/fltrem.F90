SUBROUTINE fltrem
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-05  Time: 22:58:53

!     *************************************************************************

!     FLTREM
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     30-AUG-2009:  CREATED
!     08-AUG-2012:  RSC EVALUATE ALL SPECIES

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     CARRIES OUT SPATIAL FILTERING

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc,ispec


!     BEGIN
!     =====

!     =========================================================================

!     DENSITY
CALL filtrx(drhs)
CALL filtry(drhs)
CALL filtrz(drhs)

!     U-VELOCITY
CALL filtrx(urhs)
CALL filtry(urhs)
CALL filtrz(urhs)

!     V-VELOCITY
CALL filtrx(vrhs)
CALL filtry(vrhs)
CALL filtrz(vrhs)

!     W-VELOCITY
CALL filtrx(wrhs)
CALL filtry(wrhs)
CALL filtrz(wrhs)

!     INTERNAL ENERGY
CALL filtrx(erhs)
CALL filtry(erhs)
CALL filtrz(erhs)

!     SPECIES MASS FRACTION
!     RSC 08-AUG-2012 EVALUATE ALL SPECIES
!      DO ISPEC = 1,NSPM1
DO ispec = 1,nspec
  
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        store7(ic,jc,kc) = yrhs(ic,jc,kc,ispec)
        
      END DO
    END DO
  END DO
  
  CALL filtrx(store7)
  CALL filtry(store7)
  CALL filtrz(store7)
  
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        yrhs(ic,jc,kc,ispec) = store7(ic,jc,kc)
        
      END DO
    END DO
  END DO
  
END DO

!     =========================================================================


RETURN
END SUBROUTINE fltrem
