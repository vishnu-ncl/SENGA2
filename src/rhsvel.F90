SUBROUTINE rhsvel
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-02  Time: 15:08:24

!     *************************************************************************

!     RHSVEL
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     08-DEC-2002:  CREATED
!     17-APR-2013:  RSC MIXTURE AVERAGED TRANSPORT

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     COMPUTES RIGHT-HAND-SIDES FOR TIME INTEGRATION
!     OF CONTINUITY AND MOMENTUM EQUATIONS
!     EVALUATES PRESSURE WORK AND VISCOUS WORK TERMS IN ENERGY EQUATION

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
use com_ops_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
REAL(KIND=dp) :: fornow,prefer
INTEGER :: ic,jc,kc


!     BEGIN
!     =====

!     =========================================================================

!     CONVERT VELOCITIES
!     ------------------

!     U,V,WRHS CONTAIN RHO U,V,W: CONVERT TO U,V,W
!     U,V,W HELD IN U,V,WTMP THROUGHOUT THIS ROUTINE
!     U,V,W ARE PARALLEL
DO kc = kstab,kstob
  DO jc = jstab,jstob
    DO ic = istab,istob
      
      utmp(ic,jc,kc) = urhs(ic,jc,kc)/drhs(ic,jc,kc)
      vtmp(ic,jc,kc) = vrhs(ic,jc,kc)/drhs(ic,jc,kc)
      wtmp(ic,jc,kc) = wrhs(ic,jc,kc)/drhs(ic,jc,kc)
      
    END DO
  END DO
END DO

!     =========================================================================

!     COLLECT VELOCITY COMPONENTS FOR BCs
!     -----------------------------------

!     X-DIRECTION
IF(fxlcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      struxl(jc,kc) = utmp(istal,jc,kc)
      strvxl(jc,kc) = vtmp(istal,jc,kc)
      strwxl(jc,kc) = wtmp(istal,jc,kc)
      
    END DO
  END DO
END IF
IF(fxrcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      struxr(jc,kc) = utmp(istol,jc,kc)
      strvxr(jc,kc) = vtmp(istol,jc,kc)
      strwxr(jc,kc) = wtmp(istol,jc,kc)
      
    END DO
  END DO
END IF

!     Y-DIRECTION
IF(fylcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      struyl(ic,kc) = utmp(ic,jstal,kc)
      strvyl(ic,kc) = vtmp(ic,jstal,kc)
      strwyl(ic,kc) = wtmp(ic,jstal,kc)
      
    END DO
  END DO
END IF
IF(fyrcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      struyr(ic,kc) = utmp(ic,jstol,kc)
      strvyr(ic,kc) = vtmp(ic,jstol,kc)
      strwyr(ic,kc) = wtmp(ic,jstol,kc)
      
    END DO
  END DO
END IF

!     Z-DIRECTION
IF(fzlcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      struzl(ic,jc) = utmp(ic,jc,kstal)
      strvzl(ic,jc) = vtmp(ic,jc,kstal)
      strwzl(ic,jc) = wtmp(ic,jc,kstal)
      
    END DO
  END DO
END IF
IF(fzrcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      struzr(ic,jc) = utmp(ic,jc,kstol)
      strvzr(ic,jc) = vtmp(ic,jc,kstol)
      strwzr(ic,jc) = wtmp(ic,jc,kstol)
      
    END DO
  END DO
END IF

!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================

!     MOMENTUM EQUATIONS: CONVECTIVE TERMS
!     ------------------------------------
!     RHO U U
!     RHO U U IS PARALLEL
DO kc = kstab,kstob
  DO jc = jstab,jstob
    DO ic = istab,istob
      
      store7(ic,jc,kc) = urhs(ic,jc,kc)*utmp(ic,jc,kc)
      
    END DO
  END DO
END DO

!     D/DX RHO U U
!     STRAIGHT INTO STORE4 FOR NOW
CALL dfbydx(d_store7,d_store4)


!     U-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     (HALF) D/DX RHO U U: ALREADY IN STORE4
!                                                   STORE4 = U CONVECTIVE TERMS
!                                                          U,V,WRHS = RHO U,V,W
!     =========================================================================

!     RHO U V
!     RHO U V IS PARALLEL
DO kc = kstab,kstob
  DO jc = jstab,jstob
    DO ic = istab,istob
      
      store7(ic,jc,kc) = urhs(ic,jc,kc)*vtmp(ic,jc,kc)
      
    END DO
  END DO
END DO

!     D/DY RHO V U
!     D/DX RHO U V
CALL dfbydy(d_store7,d_store1)
CALL dfbydx(d_store7,d_store5)


!     U-EQUATION: CONVECTIVE TERMS
!     V-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     (HALF) D/DY RHO V U
!     (HALF) D/DX RHO U V: ALREADY IN STORE5
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store4(ic,jc,kc) = store4(ic,jc,kc) + store1(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                               STORE4,5 = U,V CONVECTIVE TERMS
!                                                          U,V,WRHS = RHO U,V,W
!     =========================================================================

!     RHO U W
!     RHO U W IS PARALLEL
DO kc = kstab,kstob
  DO jc = jstab,jstob
    DO ic = istab,istob
      
      store7(ic,jc,kc) = urhs(ic,jc,kc)*wtmp(ic,jc,kc)
      
    END DO
  END DO
END DO

!     D/DZ RHO W U
!     D/DX RHO U W
CALL dfbydz(d_store7,d_store1)
CALL dfbydx(d_store7,d_store6)


!     U-EQUATION: CONVECTIVE TERMS
!     W-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     (HALF) D/DZ RHO W U
!     (HALF) D/DX RHO U W: ALREADY IN STORE6
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store4(ic,jc,kc) = store4(ic,jc,kc) + store1(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                           STORE4,5,6 = U,V,W CONVECTIVE TERMS
!                                                          U,V,WRHS = RHO U,V,W
!     =========================================================================

!     RHO V V
!     RHO V V IS PARALLEL
DO kc = kstab,kstob
  DO jc = jstab,jstob
    DO ic = istab,istob
      
      store7(ic,jc,kc) = vrhs(ic,jc,kc)*vtmp(ic,jc,kc)
      
    END DO
  END DO
END DO

!     D/DY RHO V V
CALL dfbydy(d_store7,d_store1)


!     V-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     (HALF) D/DY RHO V V
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store5(ic,jc,kc) = store5(ic,jc,kc) + store1(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                           STORE4,5,6 = U,V,W CONVECTIVE TERMS
!                                                          U,V,WRHS = RHO U,V,W
!     =========================================================================

!     RHO V W
!     RHO V W IS PARALLEL
DO kc = kstab,kstob
  DO jc = jstab,jstob
    DO ic = istab,istob
      
      store7(ic,jc,kc) = vrhs(ic,jc,kc)*wtmp(ic,jc,kc)
      
    END DO
  END DO
END DO

!     D/DZ RHO W V
!     D/DY RHO V W
CALL dfbydz(d_store7,d_store1)
CALL dfbydy(d_store7,d_store2)


!     V-EQUATION: CONVECTIVE TERMS
!     W-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     (HALF) D/DZ RHO W V
!     (HALF) D/DY RHO V W
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store5(ic,jc,kc) = store5(ic,jc,kc) + store1(ic,jc,kc)
      store6(ic,jc,kc) = store6(ic,jc,kc) + store2(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                           STORE4,5,6 = U,V,W CONVECTIVE TERMS
!                                                          U,V,WRHS = RHO U,V,W
!     =========================================================================

!     RHO W W
!     RHO W W IS PARALLEL
DO kc = kstab,kstob
  DO jc = jstab,jstob
    DO ic = istab,istob
      
      store7(ic,jc,kc) = wrhs(ic,jc,kc)*wtmp(ic,jc,kc)
      
    END DO
  END DO
END DO

!     D/DZ RHO W W
CALL dfbydz(d_store7,d_store1)


!     W-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     (HALF) D/DZ RHO W W
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store6(ic,jc,kc) = store6(ic,jc,kc) + store1(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                           STORE4,5,6 = U,V,W CONVECTIVE TERMS
!                                                          U,V,WRHS = RHO U,V,W
!     =========================================================================

!     VELOCITY NORMAL DERIVATIVES
!     ---------------------------
!     DUDX,DVDY,DWDZ
CALL dfbydx(d_utmp,d_store1)
CALL dfbydy(d_vtmp,d_store2)
CALL dfbydz(d_wtmp,d_store3)
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!                                           STORE4,5,6 = U,V,W CONVECTIVE TERMS
!                                                          U,V,WRHS = RHO U,V,W
!     =========================================================================

!     COLLECT VELOCITY DERIVATIVES FOR BCs
!     ------------------------------------

!     X-DIRECTION: DUDX
IF(fxlcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      bcl1xl(jc,kc) = store1(istal,jc,kc)
      
    END DO
  END DO
END IF
IF(fxrcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      bcl1xr(jc,kc) = store1(istol,jc,kc)
      
    END DO
  END DO
END IF

!     Y-DIRECTION: DVDY
IF(fylcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      bcl1yl(ic,kc) = store2(ic,jstal,kc)
      
    END DO
  END DO
END IF
IF(fyrcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      bcl1yr(ic,kc) = store2(ic,jstol,kc)
      
    END DO
  END DO
END IF

!     Z-DIRECTION: DWDZ
IF(fzlcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      bcl1zl(ic,jc) = store3(ic,jc,kstal)
      
    END DO
  END DO
END IF
IF(fzrcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      bcl1zr(ic,jc) = store3(ic,jc,kstol)
      
    END DO
  END DO
END IF

!     =========================================================================

!     U-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     DIV RHO U U  + U DIV RHO U
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store7(ic,jc,kc) = store4(ic,jc,kc) + utmp(ic,jc,kc)*divm(ic,jc,kc)
      
    END DO
  END DO
END DO

!     RHO U DUDX
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store4(ic,jc,kc) = urhs(ic,jc,kc)*store1(ic,jc,kc)
      
    END DO
  END DO
END DO

!     HALF DIV RHO U U + HALF RHO U DUDX + HALF U DIV RHO U
!     STORE IN URHS
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      urhs(ic,jc,kc) = -half*(store4(ic,jc,kc) + store7(ic,jc,kc))
      
    END DO
  END DO
END DO


!     V-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     DIV RHO U V + V DIV RHO U
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store7(ic,jc,kc) = store5(ic,jc,kc) + vtmp(ic,jc,kc)*divm(ic,jc,kc)
      
    END DO
  END DO
END DO

!     RHO V DVDY
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store5(ic,jc,kc) = vrhs(ic,jc,kc)*store2(ic,jc,kc)
      
    END DO
  END DO
END DO

!     HALF DIV RHO U V + HALF RHO V DVDY + HALF V DIV RHO U
!     STORE IN VRHS
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      vrhs(ic,jc,kc) = -half*(store5(ic,jc,kc) + store7(ic,jc,kc))
      
    END DO
  END DO
END DO


!     W-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     DIV RHO U W + W DIV RHO U
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store7(ic,jc,kc) = store6(ic,jc,kc) + wtmp(ic,jc,kc)*divm(ic,jc,kc)
      
    END DO
  END DO
END DO

!     RHO W DWDZ
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store6(ic,jc,kc) = wrhs(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO

!     HALF DIV RHO U W + HALF RHO W DWDZ + HALF W DIV RHO U
!     STORE IN WRHS
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      wrhs(ic,jc,kc) = -half*(store6(ic,jc,kc) + store7(ic,jc,kc))
      
    END DO
  END DO
END DO
!     U,V,WRHS CONTAIN U,V,W SOURCE TERMS THROUGHOUT REMAINDER OF THIS ROUTINE
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!     =========================================================================

!     PRESSURE GRADIENTS
!     ------------------

prefer = prun(ipref,jpref,kpref)
DO kc = kstab,kstob
  DO jc = jstab,jstob
    DO ic = istab,istob
      
      store7(ic,jc,kc) = prun(ic,jc,kc) - prefer
      
    END DO
  END DO
END DO

!     8DX,8DY,8DZ
CALL dfbydx(d_store7,d_store4)
CALL dfbydy(d_store7,d_store5)
CALL dfbydz(d_store7,d_store6)
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                   STORE4,5,6 = 8DX,8DY,8DZ
!     =========================================================================

!     COLLECT PRESSURE AND ITS GRADIENTS FOR BCs
!     ------------------------------------------

!     X-DIRECTION: P AND 8DX
IF(fxlcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      strpxl(jc,kc) = prun(istal,jc,kc)
      bcl5xl(jc,kc) = store4(istal,jc,kc)
      
    END DO
  END DO
END IF
IF(fxrcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      strpxr(jc,kc) = prun(istol,jc,kc)
      bcl5xr(jc,kc) = store4(istol,jc,kc)
      
    END DO
  END DO
END IF

!     Y-DIRECTION: P AND 8DY
IF(fylcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      strpyl(ic,kc) = prun(ic,jstal,kc)
      bcl5yl(ic,kc) = store5(ic,jstal,kc)
      
    END DO
  END DO
END IF
IF(fyrcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      strpyr(ic,kc) = prun(ic,jstol,kc)
      bcl5yr(ic,kc) = store5(ic,jstol,kc)
      
    END DO
  END DO
END IF

!     Z-DIRECTION: P AND 8DZ
IF(fzlcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      strpzl(ic,jc) = prun(ic,jc,kstal)
      bcl5zl(ic,jc) = store6(ic,jc,kstal)
      
    END DO
  END DO
END IF
IF(fzrcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      strpzr(ic,jc) = prun(ic,jc,kstol)
      bcl5zr(ic,jc) = store6(ic,jc,kstol)
      
    END DO
  END DO
END IF

!     =========================================================================

!     U,V,W-EQUATIONS: PRESSURE GRADIENT TERMS
!     E-EQUATION: PRESSURE WORK TERMS
!     -------------------------------
!     8DX,8DY,8DZ
!     U 8DX + V 8DY + W 8DZ

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      urhs(ic,jc,kc) = urhs(ic,jc,kc) - store4(ic,jc,kc)
      vrhs(ic,jc,kc) = vrhs(ic,jc,kc) - store5(ic,jc,kc)
      wrhs(ic,jc,kc) = wrhs(ic,jc,kc) - store6(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) - utmp(ic,jc,kc)*store4(ic,jc,kc)  &
          - vtmp(ic,jc,kc)*store5(ic,jc,kc) - wtmp(ic,jc,kc)*store6(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!     =========================================================================

!     E-EQUATION: PRESSURE WORK TERMS
!     -------------------------------
!     P DIV U

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      erhs(ic,jc,kc) = erhs(ic,jc,kc) - prun(ic,jc,kc)*(store1(ic,jc,kc)  &
          + store2(ic,jc,kc) + store3(ic,jc,kc))
      
    END DO
  END DO
END DO
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!     =========================================================================

!     VISCOSITY
!     ---------

!     VISCOSITY IS PARALLEL
DO kc = kstab,kstob
  DO jc = jstab,jstob
    DO ic = istab,istob
      
      transp(ic,jc,kc) = transp(ic,jc,kc)*prantl
      
    END DO
  END DO
END DO
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!     -------------------------------------------------------------------------

!     MIXTURE-AVERAGED TRANSPORT
!     RSC 17-APR-2013
!C     DIAGNOSTICS
!      WRITE(6,*)'RHSVEL: visc: ',FLMAVT
IF(flmavt)THEN
  
  
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        transp(ic,jc,kc) = difmix(ic,jc,kc)
        
      END DO
    END DO
  END DO
  
!C       DIAGNOSTICS
!        KC = 1
!        JC = 1
!        WRITE(6,'(4I5)')ITIME,IRKSTP,JC,KC
!        DO IC = ISTAB,ISTOB
!          WRITE(6,'(I5,2(1PE15.7))')IC,TRANSP(IC,JC,KC)
!        ENDDO
  
!C       DIAGNOSTICS
!        KC = 1
!        JC = 1
!        IC = 1
!        WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
!        IC = 2
!        WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
!        IC = 500
!        WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
!        IC = 1000
!        WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
!        IC = 1001
!        WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
  
END IF

!     -------------------------------------------------------------------------

!     VISCOUS TERMS: TAUXXb,e,f
!     -------------

!     DVDY+DWDZ
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store6(ic,jc,kc) = store2(ic,jc,kc)+store3(ic,jc,kc)
      
    END DO
  END DO
END DO


!     TAUXXb,e,f
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store6(ic,jc,kc) = fthd*store1(ic,jc,kc) - tthd*store6(ic,jc,kc)
      
    END DO
  END DO
END DO


!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUXXb,e,f DUDX
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      erhs(ic,jc,kc) = erhs(ic,jc,kc)  &
          + transp(ic,jc,kc)*store6(ic,jc,kc)*store1(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                           STORE6 = TAUXXb,e,f
!     =========================================================================

!     VISCOSITY GRADIENT: X COMPONENT
!     ------------------
CALL dfbydx(d_transp,d_store4)
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                                STORE4 = DMUDX
!                                                           STORE6 = TAUXXb,e,f
!     =========================================================================

!     U-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUXX,Xb,e,f
!     U TAUXX,Xb,e,f

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store6(ic,jc,kc) = store6(ic,jc,kc)*store4(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN X: TAUXX,X TERM ZERO ON END POINTS
IF(fxlvsn)CALL zeroxl(d_store6)
IF(fxrvsn)CALL zeroxr(d_store6)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      urhs(ic,jc,kc) = urhs(ic,jc,kc) + store6(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + utmp(ic,jc,kc)*store6(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                                STORE4 = DMUDX
!     =========================================================================

!     VISCOUS TERMS: TAUYYb,e,f
!     -------------

!     DUDX+DWDZ
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store6(ic,jc,kc) = store1(ic,jc,kc)+store3(ic,jc,kc)
      
    END DO
  END DO
END DO


!     TAUYYb,e,f
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store6(ic,jc,kc) = fthd*store2(ic,jc,kc) - tthd*store6(ic,jc,kc)
      
    END DO
  END DO
END DO


!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUYYb,e,f DVDY
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      erhs(ic,jc,kc) = erhs(ic,jc,kc)  &
          + transp(ic,jc,kc)*store6(ic,jc,kc)*store2(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                                STORE4 = DMUDX
!                                                           STORE6 = TAUYYb,e,f
!     =========================================================================

!     VISCOSITY GRADIENT: Y COMPONENT
!     ------------------
CALL dfbydy(d_transp,d_store5)
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                            STORE4,5 = DMUDX,Y
!                                                           STORE6 = TAUYYb,e,f
!     =========================================================================

!     V-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUYY,Yb,e,f
!     V TAUYY,Yb,e,f

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store6(ic,jc,kc) = store6(ic,jc,kc)*store5(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN Y: TAUYY,Y TERM ZERO ON END POINTS
IF(fylvsn)CALL zeroyl(d_store6)
IF(fyrvsn)CALL zeroyr(d_store6)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      vrhs(ic,jc,kc) = vrhs(ic,jc,kc) + store6(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + vtmp(ic,jc,kc)*store6(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
!                                                            STORE4,5 = DMUDX,Y
!     =========================================================================

!     VISCOUS TERMS: TAUZZb,e,f
!     -------------

!     DUDX+DVDY
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store1(ic,jc,kc) = store1(ic,jc,kc)+store2(ic,jc,kc)
      
    END DO
  END DO
END DO


!     TAUZZb,e,f
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store1(ic,jc,kc) = fthd*store3(ic,jc,kc) - tthd*store1(ic,jc,kc)
      
    END DO
  END DO
END DO


!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUZZb,e,f DWDZ
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      erhs(ic,jc,kc) = erhs(ic,jc,kc)  &
          + transp(ic,jc,kc)*store1(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                                 STORE3 = DWDZ
!                                                            STORE4,5 = DMUDX,Y
!                                                           STORE1 = TAUZZb,e,f
!     =========================================================================

!     VISCOSITY GRADIENT: Z COMPONENT
!     ------------------
CALL dfbydz(d_transp,d_store6)

!                                                        STORE4,5,6 = DMUDX,Y,Z
!                                                           STORE1 = TAUZZb,e,f
!     =========================================================================

!     W-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUZZ,Zb,e,f
!     W TAUZZ,Zb,e,f

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store1(ic,jc,kc) = store1(ic,jc,kc)*store6(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN Z: TAUZZ,Z TERM ZERO ON END POINTS
IF(fzlvsn)CALL zerozl(d_store1)
IF(fzrvsn)CALL zerozr(d_store1)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      wrhs(ic,jc,kc) = wrhs(ic,jc,kc) + store1(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + wtmp(ic,jc,kc)*store1(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     VELOCITY CROSS-DERIVATIVES
!     --------------------------

!     DUDY
!     ----
CALL dfbydy(d_utmp,d_store1)


!     COLLECT VELOCITY DERIVATIVE FOR BCs
!     -----------------------------------
!     Y-DIRECTION: DUDY
IF(fylcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      bcl3yl(ic,kc) = store1(ic,jstal,kc)
      
    END DO
  END DO
END IF
IF(fyrcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      bcl3yr(ic,kc) = store1(ic,jstol,kc)
      
    END DO
  END DO
END IF
!                                                                 STORE1 = DUDY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     U-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     HALF RHO V DUDY

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      urhs(ic,jc,kc) = urhs(ic,jc,kc)  &
          - half*drhs(ic,jc,kc)*vtmp(ic,jc,kc)*store1(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                                 STORE1 = DUDY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     DVDX
!     ----
CALL dfbydx(d_vtmp,d_store2)


!     COLLECT VELOCITY DERIVATIVES FOR BCs
!     ------------------------------------
!     X-DIRECTION: DVDX
IF(fxlcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      bcl3xl(jc,kc) = store2(istal,jc,kc)
      
    END DO
  END DO
END IF
IF(fxrcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      bcl3xr(jc,kc) = store2(istol,jc,kc)
      
    END DO
  END DO
END IF
!                                                                 STORE1 = DUDY
!                                                                 STORE2 = DVDX
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     V-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     HALF RHO U DVDX

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      vrhs(ic,jc,kc) = vrhs(ic,jc,kc)  &
          - half*drhs(ic,jc,kc)*utmp(ic,jc,kc)*store2(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                                 STORE1 = DUDY
!                                                                 STORE2 = DVDX
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     VELOCITY SECOND CROSS DERIVATIVES
!     ---------------------------------

!     DUDY+DVDX
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store7(ic,jc,kc) = store1(ic,jc,kc) + store2(ic,jc,kc)
      
    END DO
  END DO
END DO

!     D2UDXY
CALL d2fdxy(d_utmp,d_store1)

!     D2VDXY
CALL d2fdxy(d_vtmp,d_store2)
!                                                               STORE1 = D2UDXY
!                                                               STORE2 = D2VDXY
!                                                          STORE7 = (DUDY+DVDX)
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     U-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUXX,Xc
!     U TAUXX,Xc

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = tthd*transp(ic,jc,kc)*store2(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN X: TAUXX,X TERMS ZERO ON END POINTS
IF(fxlvsn)CALL zeroxl(d_store3)
IF(fxrvsn)CALL zeroxr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      urhs(ic,jc,kc) = urhs(ic,jc,kc) - store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) - utmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                               STORE1 = D2UDXY
!                                                               STORE2 = D2VDXY
!                                                          STORE7 = (DUDY+DVDX)
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     V-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUYY,Yc
!     V TAUYY,Yc

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = tthd*transp(ic,jc,kc)*store1(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN Y: TAUYY,Y TERMS ZERO ON END POINTS
IF(fylvsn)CALL zeroyl(d_store3)
IF(fyrvsn)CALL zeroyr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      vrhs(ic,jc,kc) = vrhs(ic,jc,kc) - store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) - vtmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                               STORE1 = D2UDXY
!                                                               STORE2 = D2VDXY
!                                                          STORE7 = (DUDY+DVDX)
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     D2UDY2
CALL d2fdy2(d_utmp,d_store3)

!     D2UDY2+D2VDXY
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = store3(ic,jc,kc) + store2(ic,jc,kc)
      
    END DO
  END DO
END DO


!     U-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUXY,Ya,b,c
!     U TAUXY,Ya,b,c

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = transp(ic,jc,kc)*store3(ic,jc,kc)  &
          + store7(ic,jc,kc)*store5(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN Y: TAUXY,Y TERM ZERO ON END POINTS
IF(fylvst)CALL zeroyl(d_store3)
IF(fyrvst)CALL zeroyr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      urhs(ic,jc,kc) = urhs(ic,jc,kc) + store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + utmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                               STORE1 = D2UDXY
!                                                          STORE7 = (DUDY+DVDX)
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     D2VDX2
CALL d2fdx2(d_vtmp,d_store3)

!     D2UDXY+D2VDX2
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = store3(ic,jc,kc) + store1(ic,jc,kc)
      
    END DO
  END DO
END DO


!     V-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUYX,Xa,b,c
!     V TAUYX,Xa,b,c

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = transp(ic,jc,kc)*store3(ic,jc,kc)  &
          + store7(ic,jc,kc)*store4(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN X: TAUYX,X TERM ZERO ON END POINTS
IF(fxlvst)CALL zeroxl(d_store3)
IF(fxrvst)CALL zeroxr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      vrhs(ic,jc,kc) = vrhs(ic,jc,kc) + store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + vtmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                          STORE7 = (DUDY+DVDX)
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUXY(DUDY+DVDX)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      erhs(ic,jc,kc) = erhs(ic,jc,kc)  &
          + transp(ic,jc,kc)*store7(ic,jc,kc)*store7(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     DUDZ
!     ----
CALL dfbydz(d_utmp,d_store1)


!     COLLECT VELOCITY DERIVATIVES FOR BCs
!     ------------------------------------
!     Z-DIRECTION: DUDZ
IF(fzlcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      bcl3zl(ic,jc) = store1(ic,jc,kstal)
      
    END DO
  END DO
END IF
IF(fzrcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      bcl3zr(ic,jc) = store1(ic,jc,kstol)
      
    END DO
  END DO
END IF
!                                                                 STORE1 = DUDZ
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     U-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     HALF RHO W DUDZ

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      urhs(ic,jc,kc) = urhs(ic,jc,kc)  &
          - half*drhs(ic,jc,kc)*wtmp(ic,jc,kc)*store1(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                                 STORE1 = DUDZ
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     DWDX
!     ----
CALL dfbydx(d_wtmp,d_store2)


!     COLLECT VELOCITY DERIVATIVES FOR BCs
!     ------------------------------------
!     X-DIRECTION: DWDX
IF(fxlcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      bcl4xl(jc,kc) = store2(istal,jc,kc)
      
    END DO
  END DO
END IF
IF(fxrcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      bcl4xr(jc,kc) = store2(istol,jc,kc)
      
    END DO
  END DO
END IF
!                                                                 STORE1 = DUDZ
!                                                                 STORE2 = DWDX
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     W-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     HALF RHO U DWDX

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      wrhs(ic,jc,kc) = wrhs(ic,jc,kc)  &
          - half*drhs(ic,jc,kc)*utmp(ic,jc,kc)*store2(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                                 STORE1 = DUDZ
!                                                                 STORE2 = DWDX
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     VELOCITY SECOND CROSS DERIVATIVES
!     ---------------------------------

!     DUDZ+DWDX
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store7(ic,jc,kc) = store1(ic,jc,kc) + store2(ic,jc,kc)
      
    END DO
  END DO
END DO

!     D2UDXZ
CALL d2fdxz(d_utmp,d_store1)

!     D2WDXZ
CALL d2fdxz(d_wtmp,d_store2)
!                                                               STORE1 = D2UDXZ
!                                                               STORE2 = D2WDXZ
!                                                          STORE7 = (DUDZ+DWDX)
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     U-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUXX,Xd
!     U TAUXX,Xd

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = tthd*transp(ic,jc,kc)*store2(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN X: TAUXX,X TERMS ZERO ON END POINTS
IF(fxlvsn)CALL zeroxl(d_store3)
IF(fxrvsn)CALL zeroxr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      urhs(ic,jc,kc) = urhs(ic,jc,kc) - store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) - utmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                               STORE1 = D2UDXZ
!                                                               STORE2 = D2WDXZ
!                                                          STORE7 = (DUDZ+DWDX)
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     W-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUZZ,Zc
!     W TAUZZ,Zc

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = tthd*transp(ic,jc,kc)*store1(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN Z: TAUZZ,Z TERMS ZERO ON END POINTS
IF(fzlvsn)CALL zerozl(d_store3)
IF(fzrvsn)CALL zerozr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      wrhs(ic,jc,kc) = wrhs(ic,jc,kc) - store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) - wtmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                               STORE1 = D2UDXZ
!                                                               STORE2 = D2WDXZ
!                                                          STORE7 = (DUDZ+DWDX)
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     D2UDZ2
CALL d2fdz2(d_utmp,d_store3)

!     D2UDZ2+D2WDXZ
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = store3(ic,jc,kc) + store2(ic,jc,kc)
      
    END DO
  END DO
END DO


!     U-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUXZ,Za,b,c
!     U TAUXZ,Za,b,c

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = transp(ic,jc,kc)*store3(ic,jc,kc)  &
          + store7(ic,jc,kc)*store6(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN Y: TAUXZ,Z TERM ZERO ON END POINTS
IF(fzlvst)CALL zerozl(d_store3)
IF(fzrvst)CALL zerozr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      urhs(ic,jc,kc) = urhs(ic,jc,kc) + store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + utmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                               STORE1 = D2UDXZ
!                                                          STORE7 = (DUDZ+DWDX)
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     D2WDX2
CALL d2fdx2(d_wtmp,d_store3)

!     D2UDXZ+D2WDX2
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = store3(ic,jc,kc) + store1(ic,jc,kc)
      
    END DO
  END DO
END DO


!     W-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUZX,Xa,b,c
!     W TAUZX,Xa,b,c

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = transp(ic,jc,kc)*store3(ic,jc,kc)  &
          + store7(ic,jc,kc)*store4(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN X: TAUZX,X TERM ZERO ON END POINTS
IF(fxlvst)CALL zeroxl(d_store3)
IF(fxrvst)CALL zeroxr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      wrhs(ic,jc,kc) = wrhs(ic,jc,kc) + store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + wtmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                          STORE7 = (DUDZ+DWDX)
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUXZ(DUDZ+DWDX)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      erhs(ic,jc,kc) = erhs(ic,jc,kc)  &
          + transp(ic,jc,kc)*store7(ic,jc,kc)*store7(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     DVDZ
CALL dfbydz(d_vtmp,d_store1)


!     COLLECT VELOCITY DERIVATIVES FOR BCs
!     ------------------------------------

!     Z-DIRECTION: DVDZ
IF(fzlcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      bcl4zl(ic,jc) = store1(ic,jc,kstal)
      
    END DO
  END DO
END IF
IF(fzrcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      bcl4zr(ic,jc) = store1(ic,jc,kstol)
      
    END DO
  END DO
END IF
!                                                                 STORE1 = DVDZ
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     V-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     HALF RHO W DVDZ

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      vrhs(ic,jc,kc) = vrhs(ic,jc,kc)  &
          - half*drhs(ic,jc,kc)*wtmp(ic,jc,kc)*store1(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                                 STORE1 = DVDZ
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     DWDY
CALL dfbydy(d_wtmp,d_store2)


!     COLLECT VELOCITY DERIVATIVES FOR BCs
!     ------------------------------------

!     Y-DIRECTION: DWDY
IF(fylcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      bcl4yl(ic,kc) = store2(ic,jstal,kc)
      
    END DO
  END DO
END IF
IF(fyrcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      bcl4yr(ic,kc) = store2(ic,jstol,kc)
      
    END DO
  END DO
END IF
!                                                                 STORE1 = DVDZ
!                                                                 STORE2 = DWDY
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     W-EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     HALF RHO V DWDY

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      wrhs(ic,jc,kc) = wrhs(ic,jc,kc)  &
          - half*drhs(ic,jc,kc)*vtmp(ic,jc,kc)*store2(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                                 STORE1 = DVDZ
!                                                                 STORE2 = DWDY
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     VELOCITY SECOND CROSS DERIVATIVES
!     ---------------------------------

!     DVDZ+DWDY
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store7(ic,jc,kc) = store1(ic,jc,kc) + store2(ic,jc,kc)
      
    END DO
  END DO
END DO

!     D2VDYZ
CALL d2fdyz(d_vtmp,d_store1)

!     D2WDYZ
CALL d2fdyz(d_wtmp,d_store2)

!                                                               STORE1 = D2VDYZ
!                                                               STORE2 = D2WDYZ
!                                                          STORE7 = (DVDZ+DWDY)
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     V-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUYY,Yd
!     V TAUYY,Yd

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = tthd*transp(ic,jc,kc)*store2(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN Y: TAUYY,2 TERMS ZERO ON END POINTS
IF(fylvsn)CALL zeroyl(d_store3)
IF(fyrvsn)CALL zeroyr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      vrhs(ic,jc,kc) = vrhs(ic,jc,kc) - store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) - vtmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                               STORE1 = D2VDYZ
!                                                               STORE2 = D2WDYZ
!                                                          STORE7 = (DVDZ+DWDY)
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     W-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUZZ,Zd
!     W TAUZZ,Zd

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = tthd*transp(ic,jc,kc)*store1(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN Z: TAUZZ,Z TERMS ZERO ON END POINTS
IF(fzlvsn)CALL zerozl(d_store3)
IF(fzrvsn)CALL zerozr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      wrhs(ic,jc,kc) = wrhs(ic,jc,kc) - store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) - wtmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                               STORE1 = D2VDYZ
!                                                               STORE2 = D2WDYZ
!                                                          STORE7 = (DVDZ+DWDY)
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     D2VDZ2
CALL d2fdz2(d_vtmp,d_store3)

!     D2VDZ2+D2WDYZ
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = store3(ic,jc,kc) + store2(ic,jc,kc)
      
    END DO
  END DO
END DO


!     V-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUYZ,Za,b,c
!     V TAUYZ,Za,b,c

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = transp(ic,jc,kc)*store3(ic,jc,kc)  &
          + store7(ic,jc,kc)*store6(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN Z: TAUYZ,Z TERM ZERO ON END POINTS
IF(fzlvst)CALL zerozl(d_store3)
IF(fzrvst)CALL zerozr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      vrhs(ic,jc,kc) = vrhs(ic,jc,kc) + store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + vtmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                               STORE1 = D2VDYZ
!                                                          STORE7 = (DVDZ+DWDY)
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     D2WDY2
CALL d2fdy2(d_wtmp,d_store3)

!     D2VDYZ+D2WDY2
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = store3(ic,jc,kc) + store1(ic,jc,kc)
      
    END DO
  END DO
END DO


!     W-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUZY,Ya,b,c
!     W TAUZY,Ya,b,c

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store3(ic,jc,kc) = transp(ic,jc,kc)*store3(ic,jc,kc)  &
          + store7(ic,jc,kc)*store5(ic,jc,kc)
      
    END DO
  END DO
END DO

!     BOUNDARY CONDITIONS
!     BC IN Y: TAUZY,Y TERM ZERO ON END POINTS
IF(fylvst)CALL zeroyl(d_store3)
IF(fyrvst)CALL zeroyr(d_store3)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      wrhs(ic,jc,kc) = wrhs(ic,jc,kc) + store3(ic,jc,kc)
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + wtmp(ic,jc,kc)*store3(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                          STORE7 = (DVDZ+DWDY)
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUYZ(DVDZ+DWDY)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      erhs(ic,jc,kc) = erhs(ic,jc,kc)  &
          + transp(ic,jc,kc)*store7(ic,jc,kc)*store7(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                            TRANSP = VISCOSITY
!                                                        STORE4,5,6 = DMUDX,Y,Z
!     =========================================================================

!     VELOCITY SECOND NORMAL DERIVATIVE TERMS
!     ---------------------------------------
!     D2UDX2,D2VDY2,D2WDZ2
CALL d2fdx2(d_utmp,d_store1)
CALL d2fdy2(d_vtmp,d_store2)
CALL d2fdz2(d_wtmp,d_store3)

!     BOUNDARY CONDITIONS
!     BC IN X: TAUXX,Xa TERM ZERO ON END POINTS
IF(fxlvsn)CALL zeroxl(d_store1)
IF(fxrvsn)CALL zeroxr(d_store1)
!     BC IN Y: TAUYY,Ya TERM ZERO ON END POINTS
IF(fylvsn)CALL zeroyl(d_store2)
IF(fyrvsn)CALL zeroyr(d_store2)
!     BC IN Z: TAUZZ,Za TERM ZERO ON END POINTS
IF(fzlvsn)CALL zerozl(d_store3)
IF(fzrvsn)CALL zerozr(d_store3)


!     U-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUXX,Xa
!     U TAUXX,Xa

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fornow = transp(ic,jc,kc)*fthd*store1(ic,jc,kc)
      urhs(ic,jc,kc) = urhs(ic,jc,kc) + fornow
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*utmp(ic,jc,kc)
      
    END DO
  END DO
END DO

!     V-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUYY,Ya
!     V TAUYY,Ya

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fornow = transp(ic,jc,kc)*fthd*store2(ic,jc,kc)
      vrhs(ic,jc,kc) = vrhs(ic,jc,kc) + fornow
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*vtmp(ic,jc,kc)
      
    END DO
  END DO
END DO


!     W-EQUATION: VISCOUS STRESS TERMS
!     E-EQUATION: VISCOUS WORK TERMS
!     ------------------------------
!     TAUZZ,Za
!     W TAUZZ,Za

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fornow = transp(ic,jc,kc)*fthd*store3(ic,jc,kc)
      wrhs(ic,jc,kc) = wrhs(ic,jc,kc) + fornow
      erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*wtmp(ic,jc,kc)
      
    END DO
  END DO
END DO

!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================

!     CONTINUITY EQUATION
!     -------------------
!     DIV RHO U

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      drhs(ic,jc,kc) = -divm(ic,jc,kc)
      
    END DO
  END DO
END DO

!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================


RETURN
END SUBROUTINE rhsvel
