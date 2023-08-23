SUBROUTINE chkarr(array,  &
        nchkxl,nchkxr,nchkyl,nchkyr,nchkzl,nchkzr,  &
        istac,istoc,jstac,jstoc,kstac,kstoc)

!   *************************************************************************

!   CHKARR
!   ======

!   AUTHOR
!   ------
!   R.S.CANT - CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   04-JUL-2004:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   DIAGNOSTIC ROUTINE
!   CHECKS UNIFORMITY OF THE SPECIFIED ARRAY

!   *************************************************************************

!   ARGUMENTS
!   =========

    integer(kind=4), intent(IN OUT)     :: nchkxl
    integer(kind=4), intent(IN OUT)     :: nchkxr
    integer(kind=4), intent(IN OUT)     :: nchkyl
    integer(kind=4), intent(IN OUT)     :: nchkyr
    integer(kind=4), intent(IN OUT)     :: nchkzl
    integer(kind=4), intent(IN OUT)     :: nchkzr
    integer(kind=4), intent(IN)         :: istac
    integer(kind=4), intent(IN)         :: istoc
    integer(kind=4), intent(IN)         :: jstac
    integer(kind=4), intent(IN)         :: jstoc
    integer(kind=4), intent(IN)         :: kstac
    integer(kind=4), intent(IN)         :: kstoc
    real(kind=8), intent(IN)   :: array(nchkxl:nchkxr,nchkyl:nchkyr,nchkzl:nchkzr)

!   LOCAL DATA
!   ==========
    real(kind=8) :: arrmin,arrmax
    integer(kind=4) :: ic,jc,kc
    integer(kind=4) :: icmin,jcmin,kcmin
    integer(kind=4) :: icmax,jcmax,kcmax

!   BEGIN
!   =====

    WRITE(6,*)nchkxl,nchkxr,nchkyl,nchkyr,nchkzl,nchkzr
    WRITE(6,*)istac,istoc,jstac,jstoc,kstac,kstoc

!   =========================================================================

    arrmin = array(istac,jstac,kstac)
    icmin = istac
    jcmin = jstac
    kcmin = kstac
    arrmax = array(istac,jstac,kstac)
    icmax = istac
    jcmax = jstac
    kcmax = kstac

    DO kc = kstac,kstoc
        DO jc = jstac,jstoc
            DO ic = istac,istoc

                IF(array(ic,jc,kc) < arrmin) THEN
                    arrmin = array(ic,jc,kc)
                    icmin = ic
                    jcmin = jc
                    kcmin = kc
                END IF
                IF(array(ic,jc,kc) > arrmax) THEN
                    arrmax = array(ic,jc,kc)
                    icmax = ic
                    jcmax = jc
                    kcmax = kc
                END IF

            END DO
        END DO
    END DO

    WRITE(6,'(2(1PE12.4,3I5))')arrmin,icmin,jcmin,kcmin, arrmax,icmax,jcmax,kcmax

!   =========================================================================

END SUBROUTINE chkarr
