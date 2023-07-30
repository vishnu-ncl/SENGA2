SUBROUTINE buftxl

! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:15:32

!     *************************************************************************

!     BUFTXL
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     02-APR-2006:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     INLET TURBULENT VELOCITY FIELD
!     FORWARD FFT WITH FOLDING IN X

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
integer(kind=4) :: jpen(npenmx),kpen(npenmx)
integer(kind=4) :: ic,jc,kc,ipencl
integer(kind=4) :: iim,iic
integer(kind=4) :: npencl


!     BEGIN
!     =====

!     =========================================================================

!     CARRY OUT A FORWARD FOURIER TRANSFORM
!     =====================================

!     =========================================================================

!     X-DIRECTION
!     -----------

!     PENCIL COUNTER
ipencl = 0

!     LOOP OVER X-LINES
DO kc = 1,nzsize
  DO jc = 1,nysize

!         PENCIL INDEXING
    ipencl = ipencl + 1
    jpen(ipencl) = jc
    kpen(ipencl) = kc

!         ASSEMBLE FFT DATA
    DO ic = 1,nxsize

      iic = 2*ic
      iim = iic-1

        write(*, '(a)') "Using the arrays not allocated by OPS, &
                        Please implement the function in OPS first, buftxl.F90: ID=76"
            STOP

      ftpart(iim,1,ipencl) = urun(ic,jc,kc)
      ftpart(iic,1,ipencl) = utmp(ic,jc,kc)
      ftpart(iim,2,ipencl) = vrun(ic,jc,kc)
      ftpart(iic,2,ipencl) = vtmp(ic,jc,kc)
      ftpart(iim,3,ipencl) = wrun(ic,jc,kc)
      ftpart(iic,3,ipencl) = wtmp(ic,jc,kc)

    END DO

!         CHECK FOR MAX NO OF PENCILS
    IF(ipencl == npenmx)THEN

!           DO THE FFT
      npencl = ipencl
      CALL fftixl(npencl,ixproc,nxproc,nxsize,nxglbl,nprocx)

!           RESTORE TRANSFORMED DATA
      DO ipencl = 1,npencl

        DO ic = 1,nxsize

            write(*, '(a)') "Using the arrays not allocated by OPS, &
                        Please implement the function in OPS first, buftxl.F90: ID=101"
            STOP

          ufxl(ic,jpen(ipencl),kpen(ipencl)) = ftpart(ic,1,ipencl)
          vfxl(ic,jpen(ipencl),kpen(ipencl)) = ftpart(ic,2,ipencl)
          wfxl(ic,jpen(ipencl),kpen(ipencl)) = ftpart(ic,3,ipencl)

        END DO

      END DO
      ipencl = 0

    END IF

  END DO
END DO
!     END OF LOOP OVER X-LINES

!     CHECK FOR ANY REMAINING STORED PENCILS
IF(ipencl /= 0)THEN

!       DO THE FFT
  npencl = ipencl
  CALL fftixl(npencl,ixproc,nxproc,nxsize,nxglbl,nprocx)

!       RESTORE TRANSFORMED DATA
  DO ipencl = 1,npencl

    DO ic = 1,nxsize

        write(*, '(a)') "Using the arrays not allocated by OPS, &
                        Please implement the function in OPS first, buftxl.F90: ID=132"
        STOP

      ufxl(ic,jpen(ipencl),kpen(ipencl)) = ftpart(ic,1,ipencl)
      vfxl(ic,jpen(ipencl),kpen(ipencl)) = ftpart(ic,2,ipencl)
      wfxl(ic,jpen(ipencl),kpen(ipencl)) = ftpart(ic,3,ipencl)

    END DO

  END DO

END IF

!     =========================================================================


RETURN
END SUBROUTINE buftxl
