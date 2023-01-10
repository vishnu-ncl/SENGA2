SUBROUTINE fftixl(npencl,ilproc,nlproc,nlsize,ngsize,nprocs)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:25:37

!     *************************************************************************

!     FFTIXL
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPERTMENT

!     CHANGE RECORD
!     -------------
!     07-APR-2006:  CREATED
!     04-JAN-2007:  RSC REVISE PARALLEL RECEIVES

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     FORWARD FFT FOR PARALLEL INVERSE DFT FOR INLET VELOCITY FIELD
!     X-LEFT

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     PARAMETER
!     =========

INTEGER, INTENT(IN)                      :: npencl
INTEGER, INTENT(IN OUT)                  :: ilproc
INTEGER, INTENT(IN)                      :: nlproc
INTEGER, INTENT(IN)                      :: nlsize
INTEGER, INTENT(IN)                      :: ngsize
INTEGER, INTENT(IN)                      :: nprocs(0:nlproc)

INTEGER, PARAMETER :: iforw = 1


!     ARGUMENTS
!     =========




!     LOCAL DATA
!     ==========
real(kind=dp) :: pcount
INTEGER :: nstotl(0:nprmax)
INTEGER :: nsize2,iofset
INTEGER :: nlprm1,nlsiz2,ngsiz2,ncount
INTEGER :: irproc,irtag
INTEGER :: ipencl,idirec,ix,icount,icproc


!     BEGIN
!     =====

!     =========================================================================

nlprm1 = nlproc-1
nlsiz2 = nlsize*2
ngsiz2 = ngsize*2

!     =========================================================================

!     CHECK FOR LEFTMOST PROCESSOR
IF(ilproc == 0) THEN
  
!       =======================================================================
  
!       LEFTMOST PROCESSOR OF THE ROW
!       -----------------------------
  
!       GATHER ALL PHYSICAL-SPACE DATA AND DO THE FFT
  
!       PHYSICAL-SPACE DATA ON LOCAL PROCESSOR
  DO ipencl = 1, npencl
    DO idirec = 1,3
      DO ix = 1,nlsiz2
        
        fftrow(ix,idirec,ipencl) = ftpart(ix,idirec,ipencl)
        
      END DO
    END DO
  END DO
  nstotl(0) = nlsiz2
  
!       RECEIVE PHYSICAL-SPACE DATA FROM REMOTE PROCESSORS
!       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
  iofset = 0
  DO icproc = 1,nlprm1
    
    irproc = nprocs(icproc)
    irtag = irproc*nproc+iproc
    
    CALL p_recv(pcount,1,irproc,irtag)
    CALL p_recv(parray,nparay,irproc,irtag)
    
    ncount = INT(pcount)
    nsize2 = ncount/npencl/3
    iofset = iofset + nstotl(icproc-1)
    
    icount = 0
    DO ipencl = 1,npencl
      DO idirec = 1,3
        DO ix = 1,nsize2
          
          icount = icount + 1
          fftrow(ix+iofset,idirec,ipencl) = parray(icount)
          
        END DO
      END DO
    END DO
    nstotl(icproc) = nsize2
    
  END DO
  
!       -----------------------------------------------------------------------
  
!       DO THE FOURIER TRANSFORM
  CALL fftgin(ngsize,iforw)
  
  DO ipencl = 1,npencl
    DO idirec = 1,3
      
      DO ix = 1,ngsiz2
        fftinx(ix) = fftrow(ix,idirec,ipencl)
      END DO
      
      CALL fftgen(fftinx,ngsize)
      
      DO ix = 1,ngsiz2
        fftrow(ix,idirec,ipencl) = fftinx(ix)
      END DO
      
    END DO
  END DO
  
!       -----------------------------------------------------------------------
  
!       COMPRESS FT BY REMOVING IMAGINARY PART OF ZEROTH FOURIER COEFFICIENT
  DO ipencl = 1, npencl
    DO idirec = 1,3
      
      DO ix = 2,ngsize
        fftrow(ix,idirec,ipencl) = fftrow(ix+1,idirec,ipencl)
      END DO
      
    END DO
  END DO
  
!       -----------------------------------------------------------------------
  
!       FOURIER-SPACE DATA ON LOCAL PROCESSOR
  DO ipencl = 1, npencl
    DO idirec = 1,3
      
      DO ix = 1,nlsize
        ftpart(ix,idirec,ipencl) = fftrow(ix,idirec,ipencl)
      END DO
      
    END DO
  END DO
  
!       SEND FOURIER-SPACE DATA BACK TO REMOTE PROCESSORS
!       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
  iofset = 0
  DO icproc = 1,nlprm1
    
    iofset = iofset + nstotl(icproc-1)/2
    nsize2 = nstotl(icproc)/2
    
    icount = 0
    DO ipencl = 1,npencl
      DO idirec = 1,3
        DO ix = 1,nsize2
          
          icount = icount + 1
          parray(icount) = fftrow(ix+iofset,idirec,ipencl)
          
        END DO
      END DO
    END DO
    
    ncount = icount
    pcount = REAL(ncount,kind=dp)
    irproc = nprocs(icproc)
    irtag = iproc*nproc+irproc
    CALL p_send(pcount,1,1,irproc,irtag)
    CALL p_send(parray,nparay,ncount,irproc,irtag)
    
  END DO
  
!       =======================================================================
  
ELSE
  
!       =======================================================================
  
!       NOT THE LEFTMOST PROCESSOR
!       --------------------------
!       IDENTIFY THE LEFTMOST PROCESSOR
  irproc = nprocs(0)
  
!       SEND THE PHYSICAL-SPACE DATA TO THE LEFTMOST PROCRESSOR
!       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
  icount = 0
  DO ipencl = 1,npencl
    DO idirec = 1,3
      DO ix = 1,nlsiz2
        
        icount = icount + 1
        parray(icount) = ftpart(ix,idirec,ipencl)
        
      END DO
    END DO
  END DO
  
  ncount = icount
  pcount = REAL(icount,kind=dp)
  irtag = iproc*nproc+irproc
  CALL p_send(pcount,1,1,irproc,irtag)
  CALL p_send(parray,nparay,ncount,irproc,irtag)
  
!       RECEIVE THE FOURIER-SPACE DATA FROM THE LEFTMOST PROCESSOR
!       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
  irtag = irproc*nproc+iproc
  CALL p_recv(pcount,1,irproc,irtag)
  CALL p_recv(parray,nparay,irproc,irtag)
  
  icount = 0
  DO ipencl = 1,npencl
    DO idirec = 1,3
      DO ix = 1,nlsize
        
        icount = icount + 1
        ftpart(ix,idirec,ipencl) = parray(icount)
        
      END DO
    END DO
  END DO
  
!       =======================================================================
  
END IF
!     LEFTMOST PROCESSOR

!     =========================================================================


RETURN
END SUBROUTINE fftixl
