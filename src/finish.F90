SUBROUTINE finish
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:16:50

!     *************************************************************************

!     FINISH
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     09-MAR-2003:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     TERMINATES THE PROGRAM

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
#ifdef HDF5
use hdf5io
#endif

use com_senga
!     -------------------------------------------------------------------------

!     BEGIN
!     =====

!     =========================================================================

!     FINAL REPORT
!     ------------
IF(iproc == 0)THEN
  
  OPEN(UNIT=ncrept,FILE=fnrept,STATUS='OLD',FORM='FORMATTED')
  
!       GO TO EOF
  1000    CONTINUE
  READ(ncrept,'(A)',END=1010)
  GO TO 1000
!       END OF LOOP 1000
  
  1010    BACKSPACE(ncrept)
  WRITE(ncrept,*)
  WRITE(ncrept,*)'Normal termination of program'
  WRITE(ncrept,*)'-----------------------------'
  CLOSE(ncrept)
  
END IF

!     =========================================================================
#ifdef HDF5
CALL hdf5_close
#endif
!     TERMINATE PARALLEL MESSAGE PASSING
CALL p_stop

!     =========================================================================

RETURN
END SUBROUTINE finish
