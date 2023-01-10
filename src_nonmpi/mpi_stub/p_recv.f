      SUBROUTINE P_RECV(PLOCAL,NPARAY,IDPROC,ITAG)

C     *************************************************************************
C
C     P_RECV
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     13-APR-1997:  CREATED
C     11-MAY-2003:  RSC MODIFIED FOR SENGA2
C     03-JAN-2007:  RSC WORKAROUND FOR LAM 7.1.1 ON SUSE
C     04-JAN-2007:  RSC REVISE PARALLEL RECEIVES
C     11-AUG-2009:  RSC STUB VERSION
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     PARALLEL RECEIVE ROUTINE
C     STUB FOR NON-PARALLEL SYSTEMS
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'mpif.h'
C     -------------------------------------------------------------------------


C     ARGUMENTS
C     =========
      INTEGER NPARAY,IDPROC,ITAG
      DOUBLE PRECISION PLOCAL(NPARAY)


C     BEGIN
C     =====

C     =========================================================================

C     NO ACTION REQUIRED

C     =========================================================================


      RETURN
      END
