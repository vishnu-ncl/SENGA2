SUBROUTINE fftf3d(carrre,carrim,nxphys,nyphys,nzphys,  &
        nx,ny,nz,iforw)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:25:18

!     *************************************************************************

!     FFTF3D
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     29-DEC-1998:  CREATED

!     DESCRIPTION
!     -----------
!     CARRIES OUT AN FFT IN 3D

!     REFERENCES
!     ----------
!     1) NUMERICAL RECIPES pp451-453

!     *************************************************************************


!     PARAMETERS
!     ==========

INTEGER, INTENT(IN OUT)                  :: nxphys
INTEGER, INTENT(IN OUT)                  :: nyphys
INTEGER, INTENT(IN OUT)                  :: nzphys
INTEGER, INTENT(IN)                      :: nx
INTEGER, INTENT(IN)                      :: ny
INTEGER, INTENT(IN)                      :: nz
INTEGER, INTENT(IN OUT)                  :: iforw
real(kind=8),INTENT(IN OUT)         :: carrre(nxphys,nyphys,nzphys)
real(kind=8),INTENT(IN OUT)         :: carrim(nxphys,nyphys,nzphys)

INTEGER, PARAMETER :: nftmax=1024


!     ARGUMENTS
!     =========





!     LOCAL DATA
!     ==========
real(kind=8):: tarray(nftmax)
INTEGER :: ix,jx,kx


!     BEGIN
!     =====

!     FT ON I-INDEX
!     -------------
CALL fftgin(nx,iforw)

DO kx = 1,nz
  DO jx = 1,ny
    
    DO ix = 1,nx
      tarray(2*ix-1) = carrre(ix,jx,kx)
      tarray(2*ix)   = carrim(ix,jx,kx)
    END DO
    
    CALL fftgen(tarray,nx)
    
    DO ix = 1,nx
      carrre(ix,jx,kx) = tarray(2*ix-1)
      carrim(ix,jx,kx) = tarray(2*ix)
    END DO
    
  END DO
END DO


!     FT ON J-INDEX
!     -------------
CALL fftgin(ny,iforw)

DO kx = 1,nz
  DO ix = 1,nx
    
    DO jx = 1,ny
      tarray(2*jx-1) = carrre(ix,jx,kx)
      tarray(2*jx)   = carrim(ix,jx,kx)
    END DO
    
    CALL fftgen(tarray,ny)
    
    DO jx = 1,ny
      carrre(ix,jx,kx) = tarray(2*jx-1)
      carrim(ix,jx,kx) = tarray(2*jx)
    END DO
    
  END DO
END DO


!     FT ON K-INDEX
!     -------------
CALL fftgin(nz,iforw)

DO jx = 1,ny
  DO ix = 1,nx
    
    DO kx = 1,nz
      tarray(2*kx-1) = carrre(ix,jx,kx)
      tarray(2*kx)   = carrim(ix,jx,kx)
    END DO
    
    CALL fftgen(tarray,nz)
    
    DO kx = 1,nz
      carrre(ix,jx,kx) = tarray(2*kx-1)
      carrim(ix,jx,kx) = tarray(2*kx)
    END DO
    
  END DO
END DO


RETURN
END SUBROUTINE fftf3d
