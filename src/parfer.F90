SUBROUTINE parfer
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-04  Time: 20:55:55

!     *************************************************************************

!     PARFER
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     11-MAY-2003:  CREATED
!     04-JAN-2007:  RSC REVISE PARALLEL RECEIVES
!     26-OCT-2008:  RSC/TDD BUG FIX JSTAB

!     DESCRIPTION
!     -----------
!     CARRIES OUT TRANSFER OF PARALLEL DATA

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
INTEGER :: ispec,ncount


!     BEGIN
!     =====

!     =========================================================================

!     X-WISE PARALLEL TRANSFER
!     ------------------------
IF(proddx)THEN
  
!       ODD-NUMBERED PROCESSOR - RL SR RR SL
!       -----------------------------------------------------------------------
  
  IF(prgoxl)THEN
    
!         RECEIVE FROM THE LEFT
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(drhs,parray,istalo,istolo,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(urhs,parray,istalo,istolo,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(vrhs,parray,istalo,istolo,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(wrhs,parray,istalo,istolo,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(erhs,parray,istalo,istolo,jstal,jstol,kstal,kstol)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,ixprom,itgxrl)
      CALL inhals(yrhs,ispec,parray,istalo,istolo, jstal,jstol,kstal,kstol)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoxr)THEN
    
!         SEND TO THE RIGHT
    
    ncount = nhalox*nynode*nznode
    
    CALL exhalo(drhs,parray,istari,istori,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(urhs,parray,istari,istori,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(vrhs,parray,istari,istori,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(wrhs,parray,istari,istori,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(erhs,parray,istari,istori,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istari,istori, jstal,jstol,kstal,kstol)
      CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(drhs,parray,istaro,istoro,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(urhs,parray,istaro,istoro,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(vrhs,parray,istaro,istoro,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(wrhs,parray,istaro,istoro,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(erhs,parray,istaro,istoro,jstal,jstol,kstal,kstol)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,ixprop,itgxrr)
      CALL inhals(yrhs,ispec,parray,istaro,istoro, jstal,jstol,kstal,kstol)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoxl)THEN
    
!         SEND TO THE LEFT
    
    ncount = nhalox*nynode*nznode
    
    CALL exhalo(drhs,parray,istali,istoli,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(urhs,parray,istali,istoli,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(vrhs,parray,istali,istoli,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(wrhs,parray,istali,istoli,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(erhs,parray,istali,istoli,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istali,istoli, jstal,jstol,kstal,kstol)
      CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
ELSE
  
!       EVEN-NUMBERED PROCESSOR - SR RL SL RR
!       -----------------------------------------------------------------------
  
  IF(prgoxr)THEN
    
!         SEND TO THE RIGHT
    
    ncount = nhalox*nynode*nznode
    
    CALL exhalo(drhs,parray,istari,istori,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(urhs,parray,istari,istori,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(vrhs,parray,istari,istori,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(wrhs,parray,istari,istori,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(erhs,parray,istari,istori,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istari,istori, jstal,jstol,kstal,kstol)
      CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoxl)THEN
    
!         RECEIVE FROM THE LEFT
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(drhs,parray,istalo,istolo,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(urhs,parray,istalo,istolo,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(vrhs,parray,istalo,istolo,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(wrhs,parray,istalo,istolo,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(erhs,parray,istalo,istolo,jstal,jstol,kstal,kstol)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,ixprom,itgxrl)
      CALL inhals(yrhs,ispec,parray,istalo,istolo, jstal,jstol,kstal,kstol)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         SEND TO THE LEFT
    
    ncount = nhalox*nynode*nznode
    
    CALL exhalo(drhs,parray,istali,istoli,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(urhs,parray,istali,istoli,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(vrhs,parray,istali,istoli,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(wrhs,parray,istali,istoli,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(erhs,parray,istali,istoli,jstal,jstol,kstal,kstol)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istali,istoli, jstal,jstol,kstal,kstol)
      CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoxr)THEN
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(drhs,parray,istaro,istoro,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(urhs,parray,istaro,istoro,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(vrhs,parray,istaro,istoro,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(wrhs,parray,istaro,istoro,jstal,jstol,kstal,kstol)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(erhs,parray,istaro,istoro,jstal,jstol,kstal,kstol)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,ixprop,itgxrr)
      CALL inhals(yrhs,ispec,parray,istaro,istoro, jstal,jstol,kstal,kstol)
      
    END DO
    
  END IF
  
END IF

!     =========================================================================

!     Y-WISE PARALLEL TRANSFER
!     ------------------------
!     NOTE EXTENDED X-LIMITS FOR Y TRANSFERS

IF(proddy)THEN
  
!       ODD-NUMBERED PROCESSOR - RL SR RR SL
!       -----------------------------------------------------------------------
  
  IF(prgoyl)THEN
    
!         RECEIVE FROM THE LEFT
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(drhs,parray,istab,istob,jstalo,jstolo,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(urhs,parray,istab,istob,jstalo,jstolo,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(vrhs,parray,istab,istob,jstalo,jstolo,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(wrhs,parray,istab,istob,jstalo,jstolo,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(erhs,parray,istab,istob,jstalo,jstolo,kstal,kstol)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,iyprom,itgyrl)
      CALL inhals(yrhs,ispec,parray,istab,istob, jstalo,jstolo,kstal,kstol)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoyr)THEN
    
!         SEND TO THE RIGHT
    
    ncount = nhaloy*nxnbig*nznode
    
    CALL exhalo(drhs,parray,istab,istob,jstari,jstori,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(urhs,parray,istab,istob,jstari,jstori,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(vrhs,parray,istab,istob,jstari,jstori,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(wrhs,parray,istab,istob,jstari,jstori,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(erhs,parray,istab,istob,jstari,jstori,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istab,istob, jstari,jstori,kstal,kstol)
      CALL p_send(parray,nparay,ncount,iyprop,itgysr)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(drhs,parray,istab,istob,jstaro,jstoro,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(urhs,parray,istab,istob,jstaro,jstoro,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(vrhs,parray,istab,istob,jstaro,jstoro,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(wrhs,parray,istab,istob,jstaro,jstoro,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(erhs,parray,istab,istob,jstaro,jstoro,kstal,kstol)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,iyprop,itgyrr)
      CALL inhals(yrhs,ispec,parray,istab,istob, jstaro,jstoro,kstal,kstol)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoyl)THEN
    
!         SEND TO THE LEFT
    
    ncount = nhaloy*nxnbig*nznode
    
    CALL exhalo(drhs,parray,istab,istob,jstali,jstoli,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(urhs,parray,istab,istob,jstali,jstoli,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(vrhs,parray,istab,istob,jstali,jstoli,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(wrhs,parray,istab,istob,jstali,jstoli,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(erhs,parray,istab,istob,jstali,jstoli,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istab,istob, jstali,jstoli,kstal,kstol)
      CALL p_send(parray,nparay,ncount,iyprom,itgysl)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
ELSE
  
!       EVEN-NUMBERED PROCESSOR - SR RL SL RR
!       -----------------------------------------------------------------------
  
  IF(prgoyr)THEN
    
!         SEND TO THE RIGHT
    
    ncount = nhaloy*nxnbig*nznode
    
    CALL exhalo(drhs,parray,istab,istob,jstari,jstori,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(urhs,parray,istab,istob,jstari,jstori,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(vrhs,parray,istab,istob,jstari,jstori,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(wrhs,parray,istab,istob,jstari,jstori,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(erhs,parray,istab,istob,jstari,jstori,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istab,istob, jstari,jstori,kstal,kstol)
      CALL p_send(parray,nparay,ncount,iyprop,itgysr)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoyl)THEN
    
!         RECEIVE FROM THE LEFT
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(drhs,parray,istab,istob,jstalo,jstolo,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(urhs,parray,istab,istob,jstalo,jstolo,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(vrhs,parray,istab,istob,jstalo,jstolo,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(wrhs,parray,istab,istob,jstalo,jstolo,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(erhs,parray,istab,istob,jstalo,jstolo,kstal,kstol)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,iyprom,itgyrl)
      CALL inhals(yrhs,ispec,parray,istab,istob, jstalo,jstolo,kstal,kstol)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         SEND TO THE LEFT
    
    ncount = nhaloy*nxnbig*nznode
    
    CALL exhalo(drhs,parray,istab,istob,jstali,jstoli,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(urhs,parray,istab,istob,jstali,jstoli,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(vrhs,parray,istab,istob,jstali,jstoli,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(wrhs,parray,istab,istob,jstali,jstoli,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(erhs,parray,istab,istob,jstali,jstoli,kstal,kstol)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istab,istob, jstali,jstoli,kstal,kstol)
      CALL p_send(parray,nparay,ncount,iyprom,itgysl)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoyr)THEN
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(drhs,parray,istab,istob,jstaro,jstoro,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(urhs,parray,istab,istob,jstaro,jstoro,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(vrhs,parray,istab,istob,jstaro,jstoro,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(wrhs,parray,istab,istob,jstaro,jstoro,kstal,kstol)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(erhs,parray,istab,istob,jstaro,jstoro,kstal,kstol)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,iyprop,itgyrr)
      CALL inhals(yrhs,ispec,parray,istab,istob, jstaro,jstoro,kstal,kstol)
      
    END DO
    
  END IF
  
END IF

!     =========================================================================

!     Z-WISE PARALLEL TRANSFER
!     ------------------------
!     NOTE EXTENDED X- AND Y-LIMITS FOR Y AND Z TRANSFERS

IF(proddz)THEN
  
!       ODD-NUMBERED PROCESSOR - RL SR RR SL
!       -----------------------------------------------------------------------
  
  IF(prgozl)THEN
    
!         RECEIVE FROM THE LEFT
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(drhs,parray,istab,istob,jstab,jstob,kstalo,kstolo)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(urhs,parray,istab,istob,jstab,jstob,kstalo,kstolo)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(vrhs,parray,istab,istob,jstab,jstob,kstalo,kstolo)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(wrhs,parray,istab,istob,jstab,jstob,kstalo,kstolo)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(erhs,parray,istab,istob,jstab,jstob,kstalo,kstolo)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,izprom,itgzrl)
      CALL inhals(yrhs,ispec,parray,istab,istob, jstab,jstob,kstalo,kstolo)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgozr)THEN
    
!         SEND TO THE RIGHT
    
    ncount = nhaloz*nxnbig*nynbig
    
    CALL exhalo(drhs,parray,istab,istob,jstab,jstob,kstari,kstori)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(urhs,parray,istab,istob,jstab,jstob,kstari,kstori)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(vrhs,parray,istab,istob,jstab,jstob,kstari,kstori)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(wrhs,parray,istab,istob,jstab,jstob,kstari,kstori)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(erhs,parray,istab,istob,jstab,jstob,kstari,kstori)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istab,istob, jstab,jstob,kstari,kstori)
      CALL p_send(parray,nparay,ncount,izprop,itgzsr)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(drhs,parray,istab,istob,jstab,jstob,kstaro,kstoro)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(urhs,parray,istab,istob,jstab,jstob,kstaro,kstoro)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(vrhs,parray,istab,istob,jstab,jstob,kstaro,kstoro)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(wrhs,parray,istab,istob,jstab,jstob,kstaro,kstoro)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(erhs,parray,istab,istob,jstab,jstob,kstaro,kstoro)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,izprop,itgzrr)
      CALL inhals(yrhs,ispec,parray,istab,istob, jstab,jstob,kstaro,kstoro)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgozl)THEN
    
!         SEND TO THE LEFT
    
    ncount = nhaloz*nxnbig*nynbig
    
    CALL exhalo(drhs,parray,istab,istob,jstab,jstob,kstali,kstoli)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(urhs,parray,istab,istob,jstab,jstob,kstali,kstoli)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(vrhs,parray,istab,istob,jstab,jstob,kstali,kstoli)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(wrhs,parray,istab,istob,jstab,jstob,kstali,kstoli)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(erhs,parray,istab,istob,jstab,jstob,kstali,kstoli)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istab,istob, jstab,jstob,kstali,kstoli)
      CALL p_send(parray,nparay,ncount,izprom,itgzsl)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
ELSE
  
!       EVEN-NUMBERED PROCESSOR - SR RL SL RR
!       -----------------------------------------------------------------------
  
  IF(prgozr)THEN
    
!         SEND TO THE RIGHT
    
    ncount = nhaloz*nxnbig*nynbig
    
    CALL exhalo(drhs,parray,istab,istob,jstab,jstob,kstari,kstori)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(urhs,parray,istab,istob,jstab,jstob,kstari,kstori)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(vrhs,parray,istab,istob,jstab,jstob,kstari,kstori)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(wrhs,parray,istab,istob,jstab,jstob,kstari,kstori)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(erhs,parray,istab,istob,jstab,jstob,kstari,kstori)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istab,istob, jstab,jstob,kstari,kstori)
      CALL p_send(parray,nparay,ncount,izprop,itgzsr)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgozl)THEN
    
!         RECEIVE FROM THE LEFT
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(drhs,parray,istab,istob,jstab,jstob,kstalo,kstolo)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(urhs,parray,istab,istob,jstab,jstob,kstalo,kstolo)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(vrhs,parray,istab,istob,jstab,jstob,kstalo,kstolo)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(wrhs,parray,istab,istob,jstab,jstob,kstalo,kstolo)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(erhs,parray,istab,istob,jstab,jstob,kstalo,kstolo)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,izprom,itgzrl)
      CALL inhals(yrhs,ispec,parray,istab,istob, jstab,jstob,kstalo,kstolo)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         SEND TO THE LEFT
    
    ncount = nhaloz*nxnbig*nynbig
    
    CALL exhalo(drhs,parray,istab,istob,jstab,jstob,kstali,kstoli)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(urhs,parray,istab,istob,jstab,jstob,kstali,kstoli)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
!         RSC/TDD BUG FIX JSTAB
    CALL exhalo(vrhs,parray,istab,istob,jstab,jstob,kstali,kstoli)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(wrhs,parray,istab,istob,jstab,jstob,kstali,kstoli)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(erhs,parray,istab,istob,jstab,jstob,kstali,kstoli)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,istab,istob, jstab,jstob,kstali,kstoli)
      CALL p_send(parray,nparay,ncount,izprom,itgzsl)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgozr)THEN
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(drhs,parray,istab,istob,jstab,jstob,kstaro,kstoro)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(urhs,parray,istab,istob,jstab,jstob,kstaro,kstoro)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(vrhs,parray,istab,istob,jstab,jstob,kstaro,kstoro)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(wrhs,parray,istab,istob,jstab,jstob,kstaro,kstoro)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(erhs,parray,istab,istob,jstab,jstob,kstaro,kstoro)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,izprop,itgzrr)
      CALL inhals(yrhs,ispec,parray,istab,istob, jstab,jstob,kstaro,kstoro)
      
    END DO
    
  END IF
  
END IF

!     =========================================================================

!     SPECIAL CASE OF PERIODIC TRANSFER ON SINGLE PROCESSOR
!     -----------------------------------------------------

!     X-DIRECTION
!     ONLY NEED TO CHECK ONE END
IF(nendxl == nperi)THEN
  
  CALL cxhalo(drhs,jstal,jstol,kstal,kstol)
  CALL cxhalo(urhs,jstal,jstol,kstal,kstol)
  CALL cxhalo(vrhs,jstal,jstol,kstal,kstol)
  CALL cxhalo(wrhs,jstal,jstol,kstal,kstol)
  CALL cxhalo(erhs,jstal,jstol,kstal,kstol)
  CALL cxhals(yrhs,jstal,jstol,kstal,kstol)
  
END IF

!     Y-DIRECTION
!     ONLY NEED TO CHECK ONE END
!     NOTE EXTENDED X-LIMITS FOR Y TRANSFERS
IF(nendyl == nperi)THEN
  
  CALL cyhalo(drhs,istab,istob,kstal,kstol)
  CALL cyhalo(urhs,istab,istob,kstal,kstol)
  CALL cyhalo(vrhs,istab,istob,kstal,kstol)
  CALL cyhalo(wrhs,istab,istob,kstal,kstol)
  CALL cyhalo(erhs,istab,istob,kstal,kstol)
  CALL cyhals(yrhs,istab,istob,kstal,kstol)
  
END IF

!     Z-DIRECTION
!     ONLY NEED TO CHECK ONE END
!     NOTE EXTENDED X- AND Y-LIMITS FOR Z TRANSFERS
IF(nendzl == nperi)THEN
  
  CALL czhalo(drhs,istab,istob,jstab,jstob)
  CALL czhalo(urhs,istab,istob,jstab,jstob)
  CALL czhalo(vrhs,istab,istob,jstab,jstob)
  CALL czhalo(wrhs,istab,istob,jstab,jstob)
  CALL czhalo(erhs,istab,istob,jstab,jstob)
  CALL czhals(yrhs,istab,istob,jstab,jstob)
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE parfer
