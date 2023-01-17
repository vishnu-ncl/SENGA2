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
    CALL inhalo(drhs,parray,1-nhalox,0,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(urhs,parray,1-nhalox,0,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(vrhs,parray,1-nhalox,0,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(wrhs,parray,1-nhalox,0,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(erhs,parray,1-nhalox,0,1,nysize,1,nzsize)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,ixprom,itgxrl)
      CALL inhals(yrhs,ispec,parray,1-nhalox,0, 1,nysize,1,nzsize)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoxr)THEN
    
!         SEND TO THE RIGHT
    
    ncount = nhalox*nynode*nznode
    
    CALL exhalo(drhs,parray,nxsize+1-nhalox,nxsize,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(urhs,parray,nxsize+1-nhalox,nxsize,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(vrhs,parray,nxsize+1-nhalox,nxsize,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(wrhs,parray,nxsize+1-nhalox,nxsize,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(erhs,parray,nxsize+1-nhalox,nxsize,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,nxsize+1-nhalox,nxsize, 1,nysize,1,nzsize)
      CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(drhs,parray,nxsize+1,nxsize+nhalox,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(urhs,parray,nxsize+1,nxsize+nhalox,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(vrhs,parray,nxsize+1,nxsize+nhalox,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(wrhs,parray,nxsize+1,nxsize+nhalox,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(erhs,parray,nxsize+1,nxsize+nhalox,1,nysize,1,nzsize)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,ixprop,itgxrr)
      CALL inhals(yrhs,ispec,parray,nxsize+1,nxsize+nhalox, 1,nysize,1,nzsize)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoxl)THEN
    
!         SEND TO THE LEFT
    
    ncount = nhalox*nynode*nznode
    
    CALL exhalo(drhs,parray,1,nhalox,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(urhs,parray,1,nhalox,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(vrhs,parray,1,nhalox,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(wrhs,parray,1,nhalox,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(erhs,parray,1,nhalox,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,1,nhalox, 1,nysize,1,nzsize)
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
    
    CALL exhalo(drhs,parray,nxsize+1-nhalox,nxsize,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(urhs,parray,nxsize+1-nhalox,nxsize,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(vrhs,parray,nxsize+1-nhalox,nxsize,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(wrhs,parray,nxsize+1-nhalox,nxsize,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    CALL exhalo(erhs,parray,nxsize+1-nhalox,nxsize,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,nxsize+1-nhalox,nxsize, 1,nysize,1,nzsize)
      CALL p_send(parray,nparay,ncount,ixprop,itgxsr)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoxl)THEN
    
!         RECEIVE FROM THE LEFT
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(drhs,parray,1-nhalox,0,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(urhs,parray,1-nhalox,0,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(vrhs,parray,1-nhalox,0,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(wrhs,parray,1-nhalox,0,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprom,itgxrl)
    CALL inhalo(erhs,parray,1-nhalox,0,1,nysize,1,nzsize)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,ixprom,itgxrl)
      CALL inhals(yrhs,ispec,parray,1-nhalox,0, 1,nysize,1,nzsize)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         SEND TO THE LEFT
    
    ncount = nhalox*nynode*nznode
    
    CALL exhalo(drhs,parray,1,nhalox,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(urhs,parray,1,nhalox,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(vrhs,parray,1,nhalox,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(wrhs,parray,1,nhalox,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    CALL exhalo(erhs,parray,1,nhalox,1,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,1,nhalox, 1,nysize,1,nzsize)
      CALL p_send(parray,nparay,ncount,ixprom,itgxsl)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoxr)THEN
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(drhs,parray,nxsize+1,nxsize+nhalox,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(urhs,parray,nxsize+1,nxsize+nhalox,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(vrhs,parray,nxsize+1,nxsize+nhalox,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(wrhs,parray,nxsize+1,nxsize+nhalox,1,nysize,1,nzsize)
    
    CALL p_recv(parray,nparay,ixprop,itgxrr)
    CALL inhalo(erhs,parray,nxsize+1,nxsize+nhalox,1,nysize,1,nzsize)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,ixprop,itgxrr)
      CALL inhals(yrhs,ispec,parray,nxsize+1,nxsize+nhalox, 1,nysize,1,nzsize)
      
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
    CALL inhalo(drhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,0,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(urhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,0,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,0,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,0,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(erhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,0,1,nzsize)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,iyprom,itgyrl)
      CALL inhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1-nhaloy,0,1,nzsize)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoyr)THEN
    
!         SEND TO THE RIGHT
    
    ncount = nhaloy*nxnbig*nznode
    
    CALL exhalo(drhs,parray,1-nhalox,nxsize+nhalox,nysize+1-nhaloy,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(urhs,parray,1-nhalox,nxsize+nhalox,nysize+1-nhaloy,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(vrhs,parray,1-nhalox,nxsize+nhalox,nysize+1-nhaloy,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(wrhs,parray,1-nhalox,nxsize+nhalox,nysize+1-nhaloy,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(erhs,parray,1-nhalox,nxsize+nhalox,nysize+1-nhaloy,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, nysize+1-nhaloy,nysize,1,nzsize)
      CALL p_send(parray,nparay,ncount,iyprop,itgysr)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(drhs,parray,1-nhalox,nxsize+nhalox,nysize+1,nysize+nhaloy,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(urhs,parray,1-nhalox,nxsize+nhalox,nysize+1,nysize+nhaloy,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(vrhs,parray,1-nhalox,nxsize+nhalox,nysize+1,nysize+nhaloy,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(wrhs,parray,1-nhalox,nxsize+nhalox,nysize+1,nysize+nhaloy,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(erhs,parray,1-nhalox,nxsize+nhalox,nysize+1,nysize+nhaloy,1,nzsize)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,iyprop,itgyrr)
      CALL inhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, nysize+1,nysize+nhaloy,1,nzsize)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoyl)THEN
    
!         SEND TO THE LEFT
    
    ncount = nhaloy*nxnbig*nznode
    
    CALL exhalo(drhs,parray,1-nhalox,nxsize+nhalox,1,nhaloy,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(urhs,parray,1-nhalox,nxsize+nhalox,1,nhaloy,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1,nhaloy,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1,nhaloy,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(erhs,parray,1-nhalox,nxsize+nhalox,1,nhaloy,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1,nhaloy,1,nzsize)
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
    
    CALL exhalo(drhs,parray,1-nhalox,nxsize+nhalox,nysize+1-nhaloy,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(urhs,parray,1-nhalox,nxsize+nhalox,nysize+1-nhaloy,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(vrhs,parray,1-nhalox,nxsize+nhalox,nysize+1-nhaloy,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(wrhs,parray,1-nhalox,nxsize+nhalox,nysize+1-nhaloy,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    CALL exhalo(erhs,parray,1-nhalox,nxsize+nhalox,nysize+1-nhaloy,nysize,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprop,itgysr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, nysize+1-nhaloy,nysize,1,nzsize)
      CALL p_send(parray,nparay,ncount,iyprop,itgysr)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoyl)THEN
    
!         RECEIVE FROM THE LEFT
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(drhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,0,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(urhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,0,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,0,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,0,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprom,itgyrl)
    CALL inhalo(erhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,0,1,nzsize)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,iyprom,itgyrl)
      CALL inhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1-nhaloy,0,1,nzsize)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         SEND TO THE LEFT
    
    ncount = nhaloy*nxnbig*nznode
    
    CALL exhalo(drhs,parray,1-nhalox,nxsize+nhalox,1,nhaloy,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(urhs,parray,1-nhalox,nxsize+nhalox,1,nhaloy,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1,nhaloy,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1,nhaloy,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    CALL exhalo(erhs,parray,1-nhalox,nxsize+nhalox,1,nhaloy,1,nzsize)
    CALL p_send(parray,nparay,ncount,iyprom,itgysl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1,nhaloy,1,nzsize)
      CALL p_send(parray,nparay,ncount,iyprom,itgysl)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgoyr)THEN
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(drhs,parray,1-nhalox,nxsize+nhalox,nysize+1,nysize+nhaloy,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(urhs,parray,1-nhalox,nxsize+nhalox,nysize+1,nysize+nhaloy,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(vrhs,parray,1-nhalox,nxsize+nhalox,nysize+1,nysize+nhaloy,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(wrhs,parray,1-nhalox,nxsize+nhalox,nysize+1,nysize+nhaloy,1,nzsize)
    
    CALL p_recv(parray,nparay,iyprop,itgyrr)
    CALL inhalo(erhs,parray,1-nhalox,nxsize+nhalox,nysize+1,nysize+nhaloy,1,nzsize)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,iyprop,itgyrr)
      CALL inhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, nysize+1,nysize+nhaloy,1,nzsize)
      
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
    CALL inhalo(drhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1-nhaloz,0)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(urhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1-nhaloz,0)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1-nhaloz,0)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1-nhaloz,0)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(erhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1-nhaloz,0)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,izprom,itgzrl)
      CALL inhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1-nhaloy,nysize+nhaloy,1-nhaloz,0)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgozr)THEN
    
!         SEND TO THE RIGHT
    
    ncount = nhaloz*nxnbig*nynbig
    
    CALL exhalo(drhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(urhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(erhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
      CALL p_send(parray,nparay,ncount,izprop,itgzsr)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(drhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(urhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(erhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,izprop,itgzrr)
      CALL inhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgozl)THEN
    
!         SEND TO THE LEFT
    
    ncount = nhaloz*nxnbig*nynbig
    
    CALL exhalo(drhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1,nhaloz)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(urhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1,nhaloz)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1,nhaloz)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1,nhaloz)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(erhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1,nhaloz)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1-nhaloy,nysize+nhaloy,1,nhaloz)
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
    
    CALL exhalo(drhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(urhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    CALL exhalo(erhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
    CALL p_send(parray,nparay,ncount,izprop,itgzsr)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1-nhaloy,nysize+nhaloy,nzsize+1-nhaloz,nzsize)
      CALL p_send(parray,nparay,ncount,izprop,itgzsr)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgozl)THEN
    
!         RECEIVE FROM THE LEFT
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(drhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1-nhaloz,0)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(urhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1-nhaloz,0)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1-nhaloz,0)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1-nhaloz,0)
    
    CALL p_recv(parray,nparay,izprom,itgzrl)
    CALL inhalo(erhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1-nhaloz,0)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,izprom,itgzrl)
      CALL inhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1-nhaloy,nysize+nhaloy,1-nhaloz,0)
      
    END DO
    
!         ---------------------------------------------------------------------
    
!         SEND TO THE LEFT
    
    ncount = nhaloz*nxnbig*nynbig
    
    CALL exhalo(drhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1,nhaloz)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(urhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1,nhaloz)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
!         RSC/TDD BUG FIX JSTAB
    CALL exhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1,nhaloz)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1,nhaloz)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    CALL exhalo(erhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,1,nhaloz)
    CALL p_send(parray,nparay,ncount,izprom,itgzsl)
    
    DO ispec = 1,nspec
      
      CALL exhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1-nhaloy,nysize+nhaloy,1,nhaloz)
      CALL p_send(parray,nparay,ncount,izprom,itgzsl)
      
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
  IF(prgozr)THEN
    
!         RECEIVE FROM THE RIGHT
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(drhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(urhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(vrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(wrhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
    
    CALL p_recv(parray,nparay,izprop,itgzrr)
    CALL inhalo(erhs,parray,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
    
    DO ispec = 1,nspec
      
      CALL p_recv(parray,nparay,izprop,itgzrr)
      CALL inhals(yrhs,ispec,parray,1-nhalox,nxsize+nhalox, 1-nhaloy,nysize+nhaloy,nzsize+1,nzsize+nhaloz)
      
    END DO
    
  END IF
  
END IF

!     =========================================================================

!     SPECIAL CASE OF PERIODIC TRANSFER ON SINGLE PROCESSOR
!     -----------------------------------------------------

!     X-DIRECTION
!     ONLY NEED TO CHECK ONE END
IF(nendxl == nperi)THEN
  
  CALL cxhalo(drhs,1,nysize,1,nzsize)
  CALL cxhalo(urhs,1,nysize,1,nzsize)
  CALL cxhalo(vrhs,1,nysize,1,nzsize)
  CALL cxhalo(wrhs,1,nysize,1,nzsize)
  CALL cxhalo(erhs,1,nysize,1,nzsize)
  CALL cxhals(yrhs,1,nysize,1,nzsize)
  
END IF

!     Y-DIRECTION
!     ONLY NEED TO CHECK ONE END
!     NOTE EXTENDED X-LIMITS FOR Y TRANSFERS
IF(nendyl == nperi)THEN
  
  CALL cyhalo(drhs,1-nhalox,nxsize+nhalox,1,nzsize)
  CALL cyhalo(urhs,1-nhalox,nxsize+nhalox,1,nzsize)
  CALL cyhalo(vrhs,1-nhalox,nxsize+nhalox,1,nzsize)
  CALL cyhalo(wrhs,1-nhalox,nxsize+nhalox,1,nzsize)
  CALL cyhalo(erhs,1-nhalox,nxsize+nhalox,1,nzsize)
  CALL cyhals(yrhs,1-nhalox,nxsize+nhalox,1,nzsize)
  
END IF

!     Z-DIRECTION
!     ONLY NEED TO CHECK ONE END
!     NOTE EXTENDED X- AND Y-LIMITS FOR Z TRANSFERS
IF(nendzl == nperi)THEN
  
  CALL czhalo(drhs,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy)
  CALL czhalo(urhs,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy)
  CALL czhalo(vrhs,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy)
  CALL czhalo(wrhs,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy)
  CALL czhalo(erhs,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy)
  CALL czhals(yrhs,1-nhalox,nxsize+nhalox,1-nhaloy,nysize+nhaloy)
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE parfer
