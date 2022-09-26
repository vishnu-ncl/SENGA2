SUBROUTINE fftsym(ndosym,npencl,ilproc,nlproc,nlsize,ngsize,  &
        nprocs)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:25:40

!     *************************************************************************

!     FFTSYM
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPERTMENT

!     CHANGE RECORD
!     -------------
!     25-MAY-2003:  CREATED
!     04-JAN-2007:  RSC REVISE PARALLEL RECEIVES

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     SENGA INITIAL TURBULENCE GENERATOR
!     PROCESSES THE DATA AND IMPOSES SYMMETRY FOR PARALLEL INVERSE FFT

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     PARAMETERS
!     ==========

INTEGER, INTENT(IN OUT)                  :: ndosym
INTEGER, INTENT(IN)                      :: npencl
INTEGER, INTENT(IN OUT)                  :: ilproc
INTEGER, INTENT(IN)                      :: nlproc
INTEGER, INTENT(IN)                      :: nlsize
INTEGER, INTENT(IN)                      :: ngsize
INTEGER, INTENT(IN)                      :: nprocs(0:nlproc)

INTEGER, PARAMETER :: invrs=-1


!     ARGUMENTS
!     =========




!     LOCAL DATA
!     ==========
real(kind=dp) :: pcount
INTEGER :: nstotl(0:nprmax)
INTEGER :: nsize2,iofset
INTEGER :: nlprm1,nlsiz2,ngsiz2,ncount
INTEGER :: nsymmt,nsymod,nsymax,nodbal
INTEGER :: irproc,irtag
INTEGER :: ipencl,idirec,ix,iix,iiy,icount,icproc


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
  
!       GATHER ALL FOURIER-SPACE DATA AND DO THE FFT
  
!       FOURIER-SPACE DATA ON LOCAL PROCESSOR
  DO ipencl = 1, npencl
    DO idirec = 1,3
      DO ix = 1,nlsiz2
        
        fftrow(ix,idirec,ipencl) = ftpart(ix,idirec,ipencl)
        
      END DO
    END DO
  END DO
  nstotl(0) = nlsiz2
  
!       RECEIVE FOURIER-SPACE DATA FROM REMOTE PROCESSORS
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
  
!       SYMMETRY CONDITION
  IF(ndosym == 1)THEN
    
    nsymmt = ngsize
    nsymod = MOD(nsymmt,2)
    nsymax = nsymmt/2 + nsymod
    nodbal = nsymax+1
    nsymmt = 2*(nsymmt+2)
    
    DO ipencl = 1,npencl
      DO idirec = 1,3
        
!             IGNORE THE ZEROTH MODE (IX=1)
        
        DO ix = 2,nsymax
          
          iix = 2*ix
          iiy = nsymmt-iix
          
          fftrow(iiy-1,idirec,ipencl) =  fftrow(iix-1,idirec,ipencl)
          fftrow(iiy,idirec,ipencl)   = -fftrow(iix,idirec,ipencl)
          
        END DO
        
!             ODDBALL MODE (IF ANY)
        IF(nsymod == 0)THEN
          iiy = 2*nodbal
          fftrow(iiy-1,idirec,ipencl) = zero
          fftrow(iiy,idirec,ipencl) = zero
        END IF
        
      END DO
    END DO
    
  END IF
!       SYMMETRY CONDITION
  
!       -----------------------------------------------------------------------
  
!       DO THE FOURIER TRANSFORM
  CALL fftgin(ngsize,invrs)
  
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
  
!       PHYSICAL-SPACE DATA ON LOCAL PROCESSOR
  DO ipencl = 1, npencl
    DO idirec = 1,3
      DO ix = 1,nlsiz2
        
        ftpart(ix,idirec,ipencl) = fftrow(ix,idirec,ipencl)
        
      END DO
    END DO
  END DO
  
!       SEND PHYSICAL-SPACE DATA BACK TO REMOTE PROCESSORS
!       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
  iofset = 0
  DO icproc = 1,nlprm1
    
    iofset = iofset + nstotl(icproc-1)
    nsize2 = nstotl(icproc)
    
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
  
!       SEND THE FOURIER-SPACE DATA TO THE LEFTMOST PROCESSOR
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
  pcount = REAL(ncount,kind=dp)
  irtag = iproc*nproc+irproc
  CALL p_send(pcount,1,1,irproc,irtag)
  CALL p_send(parray,nparay,ncount,irproc,irtag)
  
!       RECEIVE THE PHYSICAL-SPACE DATA FROM THE LEFTMOST PROCESSOR
!       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
  irtag = irproc*nproc+iproc
  CALL p_recv(pcount,1,irproc,irtag)
  CALL p_recv(parray,nparay,irproc,irtag)
  
  icount = 0
  DO ipencl = 1,npencl
    DO idirec = 1,3
      DO ix = 1,nlsiz2
        
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
END SUBROUTINE fftsym
