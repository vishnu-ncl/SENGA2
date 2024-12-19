SUBROUTINE bcutxl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2024-11-08  Time: 05:06:40

!     *************************************************************************

!     BCUTXL
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     30-DEC-2003:  CREATED
!     04-JAN-2007:  RSC REVISE PARALLEL RECEIVES

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
!     AND THEIR TIME DERIVATIVES

!     X-DIRECTION LEFT-HAND END

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
!KA   FIX INFLOW BUG, BTIME IS DEFINED IN COM_SENGA2.H
!KA      DOUBLE PRECISION BTIME
DOUBLE PRECISION :: fornow,argmnt,argval,realkx
DOUBLE PRECISION :: cosval,sinval,costht,sintht
DOUBLE PRECISION :: pcount
INTEGER :: ic,jc,kc
INTEGER :: iic,iim,kx,kxbase
INTEGER :: icproc,ncount,irproc,irtag
!     VM: SYNTHETIC DIGITAL FILTERING METHOD
DOUBLE PRECISION :: vfac,dvfdt,coflow

!     FY - FOR NON-REFLECTING INLOW (INFLOW OPTION 4)
DOUBLE PRECISION :: deltagy,ycoord
INTEGER :: igofsty,iy
DOUBLE PRECISION :: lambda
PARAMETER(lambda = 3.74D-3) ! SOUND WAVELENGTH
DOUBLE PRECISION :: pulrat
PARAMETER(pulrat = 0.1D0) ! RATIO OF PULSE WIDTH TO SOUND WAVELENGTH
DOUBLE PRECISION :: ptly
PARAMETER(ptly = pulrat*lambda) ! BASE WIDTH OF PULSE
DOUBLE PRECISION :: widthp
PARAMETER(widthp = 0.5D0) ! ADDITIONAL WIDTH PARAMETER - DEFAULT 0.5
DOUBLE PRECISION :: slope
PARAMETER(slope = 2.0D4) ! LARGER VALUE RESULTS IN SHARPER SLOPE
DOUBLE PRECISION :: hyptan ! FOR HYPERBOLIC TANGENT PROFILE

!     BEGIN
!     =====

!     =========================================================================

!KA   THIS WAS MOVED TO BOUNDT & BOUNTT TO FIX INFLOW SCANNING LOCATION
!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)
!KA      BTIME = ETIME + RKTIM(IRKSTP)

!     =========================================================================

!     CONSTANT U-VELOCITY
!     PARAMETER I1=1, R1=U-VELOCITY
IF(nxlprm(1) == 1)THEN
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      struxl(jc,kc) = rxlprm(1)
      strvxl(jc,kc) = zero
      strwxl(jc,kc) = zero
      
      dudtxl(jc,kc) = zero
      dvdtxl(jc,kc) = zero
      dwdtxl(jc,kc) = zero
      
!           WRITE ETIME, PRESSURE
      IF (irkstp == nrkstp) THEN
        OPEN(UNIT=16,FILE="vars.dat",ACCESS='APPEND')
        WRITE(16,'(3E20.9E3)')etime,(strpxl(jc,1)-prin),strvel
        CLOSE(16)
      END IF
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     SINUSOIDAL U-VELOCITY
!     PARAMETER I1=2, R1=AMPLITUDE, R2=PERIOD

!     FY - EDITED FOR USE WITH NON-REFLECTING INFLOW

IF(nxlprm(1) == 2)THEN
  
  deltagy = ygdlen/REAL(nyglbl-1)
  
  igofsty = 0
  DO icproc = 0, iyproc-1
    igofsty = igofsty + npmapy(icproc)
  END DO
  
  fornow = two*pi/rxlprc(8)
  argmnt = fornow*(etime+rxlprc(9)*rxlprc(8))
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
!             GLOBAL Y COORD
      iy = igofsty + jc
      
!             NXLPRM(4)=0 - OPTION FOR STANDARD FLAT SINUSOIDAL INFLOW VELOCITY IN TIME ON WHOLE XL FACE
      IF (nxlprm(4) == 0) THEN
        
        struxl(jc,kc) = rxlprm(1)+rxlprc(7)*SIN(argmnt)
        strvxl(jc,kc) = zero
        strwxl(jc,kc) = zero
        
        dudtxl(jc,kc) = fornow*rxlprc(7)*COS(argmnt)
        dvdtxl(jc,kc) = zero
        dwdtxl(jc,kc) = zero
        
!               WRITE ETIME, PRESSURE
        IF (irkstp == nrkstp) THEN
          OPEN(UNIT=16,FILE="vars.dat",ACCESS='APPEND')
          WRITE(16,'(3E20.9E3)')etime,(strpxl(jc,1)-prin),strvel
          CLOSE(16)
        END IF
        
!             NXLPRM(4)=4 - OPTION FOR SINUSOIDAL VELOCITY IN TIME FOR PART OF XL FACE - ACTS AS A POINT SOURCE
      ELSE IF (nxlprm(4) == 4) THEN
        
!               DEFINE HYPERBOLIC TANGENT PROFILE
        hyptan=(0.5D0*TANH(slope*(DBLE(iy-1)*deltagy-  &
            (ygdlen/2.0D0-widthp*ptly)))+0.5D0)*  &
            (-0.5D0*TANH(slope*(DBLE(iy-1)*deltagy-  &
            (ygdlen/2.0D0+widthp*ptly)))+0.5D0)
        
!               SET VELOCITY ON FACE BASED ON PROFILE
        struxl(jc,kc) = rxlprm(1)+hyptan*rxlprc(7)*SIN(argmnt)
        strvxl(jc,kc) = zero
        strwxl(jc,kc) = zero
        
        dudtxl(jc,kc) = hyptan*fornow*rxlprc(7)*COS(argmnt)
        dvdtxl(jc,kc) = zero
        dwdtxl(jc,kc) = zero
        
!               WRITE ETIME, PRESSURE
        IF (irkstp == nrkstp) THEN
          OPEN(UNIT=16,FILE="vars.dat",ACCESS='APPEND')
          WRITE(16,'(3E20.9E3)')etime,(strpxl(jc,1)-prin),strvel
          CLOSE(16)
        END IF
        
        
      END IF
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     GENERATING TURBULENT FIELD USING SYNTHETIC DIGITAL FILTERING
!     METHOD
!     VM: NXLPRM(1)=4 IMPLIES THAT THE VELOCITY SYTHETIC DIGITAL FILTERING
!     IS ON
IF(nxlprm(1) == 4)THEN
  
!        VFAC=MIN(1.0D0,0.1+ETIME/0.0001)
  vfac=one
  
  coflow = rxlprm(4)
  
  IF (vfac < 1.0D0) THEN
    dvfdt=10000.0D0
  ELSE
    dvfdt=0.0D0
  END IF
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      struxl(jc,kc) = (rxlprm(1)*ustead(jc,kc)+ uinf2(jc,kc))*vfac+coflow
      strvxl(jc,kc) = vinf2(jc,kc)*vfac
      strwxl(jc,kc) = winf2(jc,kc)*vfac
      
      
      dudtxl(jc,kc) = (uinf2(jc,kc)-uinf1(jc,kc))/tstep  &
          *vfac+(rxlprm(1)*ustead(jc,kc)+uinf2(jc,kc)) *dvfdt
      dvdtxl(jc,kc) = (vinf2(jc,kc)-vinf1(jc,kc))/tstep  &
          *vfac+vinf2(jc,kc)*dvfdt
      dwdtxl(jc,kc) = (winf2(jc,kc)-winf1(jc,kc))/tstep  &
          *vfac+winf2(jc,kc)*dvfdt
      
      
    END DO
  END DO
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE bcutxl
