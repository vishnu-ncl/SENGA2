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
    use data_types

!   ARGUMENTS
!   =========

    real(kind=dp), intent(IN)   :: array(nchkxl:nchkxr,nchkyl:nchkyr,nchkzl:nchkzr)
    integer, intent(IN OUT)     :: nchkxl
    integer, intent(IN OUT)     :: nchkxr
    integer, intent(IN OUT)     :: nchkyl
    integer, intent(IN OUT)     :: nchkyr
    integer, intent(IN OUT)     :: nchkzl
    integer, intent(IN OUT)     :: nchkzr
    integer, intent(IN)         :: istac
    integer, intent(IN)         :: istoc
    integer, intent(IN)         :: jstac
    integer, intent(IN)         :: jstoc
    integer, intent(IN)         :: kstac
    integer, intent(IN)         :: kstoc

!   LOCAL DATA
!   ==========
    real(kind=dp) :: arrmin,arrmax
    integer :: ic,jc,kc
    integer :: icmin,jcmin,kcmin
    integer :: icmax,jcmax,kcmax

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
