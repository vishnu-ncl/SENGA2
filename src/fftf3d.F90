SUBROUTINE fftf3d(carrre,carrim,nxphys,nyphys,nzphys,  &
        nx,ny,nz,iforw)

!   *************************************************************************

!   FFTF3D
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   29-DEC-1998:  CREATED

!   DESCRIPTION
!   -----------
!   CARRIES OUT AN FFT IN 3D

!   REFERENCES
!   ----------
!   1) NUMERICAL RECIPES pp451-453

!   *************************************************************************

!   PARAMETERS
!   ==========
    integer(kind=4), parameter :: nftmax=1024

    integer(kind=4), intent(in out)                  :: nxphys
    integer(kind=4), intent(in out)                  :: nyphys
    integer(kind=4), intent(in out)                  :: nzphys
    integer(kind=4), intent(in)                      :: nx
    integer(kind=4), intent(in)                      :: ny
    integer(kind=4), intent(in)                      :: nz
    integer(kind=4), intent(in out)                  :: iforw
    real(kind=8), intent(in out)         :: carrre(nxphys,nyphys,nzphys)
    real(kind=8), intent(in out)         :: carrim(nxphys,nyphys,nzphys)


!   ARGUMENTS
!   =========

!   LOCAL DATA
!   ==========
    real(kind=8) :: tarray(nftmax)
    integer(kind=4) :: ix,jx,kx

!   BEGIN
!   =====

!   FT ON I-INDEX
!   -------------
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

!   FT ON J-INDEX
!   -------------
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

!   FT ON K-INDEX
!   -------------
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

!   =========================================================================

END SUBROUTINE fftf3d
