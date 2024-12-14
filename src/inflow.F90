SUBROUTINE inflow
 
    use OPS_Fortran_Reference
    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!     *************************************************************************

!     INFLOW
!     ======

!     AUTHOR
!     ------
!     M.KLEIN  --   UNIVERSITÄT DER BUNDESWEHR MÜNCHEN

!     CHANGE RECORD
!     -------------
!     27-MAR-2015:  CREATED
!     10-APR-2019:  CHANGED FOR SENGA2 (M.A.)

!     DESCRIPTION
!     -----------

!     *************************************************************************
!      IMPLICIT NONE


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
INCLUDE 'mpif.h'
INCLUDE 'com_senga2.h'
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
REAL :: ran1        ! Random number Generator

DOUBLE PRECISION :: norm,ay(-nfy:nfy),az(-nfz:nfz)  ! filter coefficients

LOGICAL :: periy,periz!,TRADIT
PARAMETER (periy=.true.,periz=.false.)!,TRADIT=.TRUE.)


INTEGER :: ic,jc,kc,jc2,kc2,yp,zp
INTEGER :: jfiltstart,kfiltstart
INTEGER :: tspc(2),nspc,ispec

!     VM: CALCULATING USTEAD
DOUBLE PRECISION :: rad
DOUBLE PRECISION :: zdist,ydist
INTEGER :: kg,jg
INTEGER :: yoffset,zoffset,ymidpnt,zmidpnt
!     VM: DEBUG
INTEGER :: snap,icproc
CHARACTER (LEN=40) :: fname
CHARACTER (LEN=4) :: psnap

    real(kind=8) :: urms,yrms,phi,delx,umean,lenx,leny,lenz,tstp,theta
    real(kind=8) :: sum_umean, sum_denom
    real(kind=8) :: yref(nspec), yref_ispec
    integer(kind=4) :: rangexyz(6)

!   LNX, LNY, LNZ: LENGTH SCALES IN TERMS OF DELX
!   THETA: INDUCED TIME SCALE
!   PHI: REQUIRED WEIGHTING PARAMETER
!   NFY,NFZ: FILTER SIZE: NF=2*LN
    delx = xgdlen/REAL(nxglbl-1,kind=8)
    umean = rxlprm(1)
    urms = rxlprm(2)
    yrms = 0.050_8
    lenx = REAL(lnx,kind=8)
    leny = REAL(lny,kind=8)
    lenz = REAL(lnz,kind=8)
    tstp = tstep
    theta = lenx*delx / umean
    phi = 1.0_8-tstp / theta

    zp = INT(iproc/(nxproc*nyproc))
    yp = INT((iproc-zp*nyproc*nxproc)/nxproc)

    jfiltstart = yp*nysize
    kfiltstart = zp*nzsize

    tspc = [1,2]
    nspc = 2
    DO ispec = 1,nspec
        yref(ispec) = 0.0_8
    END DO
    yref(1) = 1.0_8
    yref(2) = 0.233_8

!     ------------------------------------------------------------------

!C    URMS=MAX(5.0D0,3.0D0+FLOAT(ITIME)/10000.D0*1.0D0)
!     CHECK RAMP UP OF VELOCITY
!      IF (NXLPRM(2).EQ.1)THEN
!         TO BE DONE
!      END IF
pi = four*ATAN(1.0D0)
!     CALCULATE FILTER COEFFICIENTS
norm=0.0
DO jc=-nfy,nfy
  ay(jc)=EXP(-pi*jc*jc/(2*lny*lny))
  norm=norm+ay(jc)**2
END DO
norm=SQRT(norm)
DO jc=-nfy,nfy
  ay(jc)=ay(jc)/norm
END DO

norm=0.0
DO kc=-nfz,nfz
  az(kc)=EXP(-pi*kc*kc/(2*lnz*lnz))
  norm=norm+az(kc)**2
END DO
norm=SQRT(norm)
DO kc=-nfz,nfz
  az(kc)=az(kc)/norm
END DO

intran = itime

!   INITIALIZE RANDOM ARRAYS
    rangexyz = [1,nxglbl,1-nfy,nyglbl+nfy,1-nfz,nzglbl+nfz]
    call ops_par_loop(, "", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(intran, 1, "integer(kind=4)", OPS_READ))


DO jc=1-nfy,nyglbl+nfy
  DO kc=1-nfz,nzglbl+nfz
    urand(jc,kc)=ran1(intran)
    vrand(jc,kc)=ran1(intran)
    wrand(jc,kc)=ran1(intran)
  END DO
END DO
IF(nxlprm(2)==1)THEN
  DO jc=1-nfy,nyglbl+nfy
    DO kc=1-nfz,nzglbl+nfz
      yrand(jc,kc,:)=0.0D0
      DO ispc = 1,nspc
        ispec = tspc(ispc)
        yrand(jc,kc,ispec)=ran1(intran)
      END DO
    END DO
  END DO
END IF

IF (periy) THEN   ! overwrite in a periodic manner
  DO jc=-nfy+1,0
    DO kc=1-nfz,nzglbl+nfz
      urand(jc,kc)=urand(jc+nyglbl,kc)
      vrand(jc,kc)=vrand(jc+nyglbl,kc)
      wrand(jc,kc)=wrand(jc+nyglbl,kc)
    END DO
  END DO
  DO jc=nyglbl+1,nyglbl+nfy
    DO kc=1-nfz,nzglbl+nfz
      urand(jc,kc)=urand(jc-nyglbl,kc)
      vrand(jc,kc)=vrand(jc-nyglbl,kc)
      wrand(jc,kc)=wrand(jc-nyglbl,kc)
    END DO
  END DO
  IF(nxlprm(2)==1)THEN
    DO jc=-nfy+1,0
      DO kc=1-nfz,nzglbl+nfz
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yrand(jc,kc,ispec)=yrand(jc+nyglbl,kc,ispec)
        END DO
      END DO
    END DO
    DO jc=nyglbl+1,nyglbl+nfy
      DO kc=1-nfz,nzglbl+nfz
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yrand(jc,kc,ispec)=yrand(jc-nyglbl,kc,ispec)
        END DO
      END DO
    END DO
  END IF
END IF

IF (periz) THEN   ! overwrite in a periodic manner
  DO kc=-nfz+1,0
    DO jc=1-nfy,nyglbl+nfy
      urand(jc,kc)=urand(jc,kc+nzglbl)
      vrand(jc,kc)=vrand(jc,kc+nzglbl)
      wrand(jc,kc)=wrand(jc,kc+nzglbl)
    END DO
  END DO
  DO kc=nzglbl+1,nzglbl+nfz
    DO jc=1-nfy,nyglbl+nfy
      urand(jc,kc)=urand(jc,kc-nzglbl)
      vrand(jc,kc)=vrand(jc,kc-nzglbl)
      wrand(jc,kc)=wrand(jc,kc-nzglbl)
    END DO
  END DO
  IF(nxlprm(2)==1)THEN
    DO kc=-nfz+1,0
      DO jc=1-nfy,nyglbl+nfy
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yrand(jc,kc,ispec)=yrand(jc,kc+nzglbl,ispec)
        END DO
      END DO
    END DO
    DO kc=nzglbl+1,nzglbl+nfz
      DO jc=1-nfy,nyglbl+nfy
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yrand(jc,kc,ispec)=yrand(jc,kc-nzglbl,ispec)
        END DO
      END DO
    END DO
  END IF
END IF

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ufilt, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vfilt, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wfilt, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    DO ispec = 1,nspec
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yfilt(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    END DO

    rangexyz = [1,nxglbl,1,nyglbl,1-nfz,nzglbl+nfz]
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ufold, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vfold, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wfold, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    DO ispec = 1,nspec
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yfold(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    END DO

!     What is TRADIT?
DO jc=jstal,jstol
  DO kc=kstal-nfz,kstol+nfz
    DO jc2=-nfy,nfy
      
      ufold(jc,kc)=ufold(jc,kc)+  &
          urand(jc+jfiltstart+jc2,kc+kfiltstart)*ay(jc2)
      vfold(jc,kc)=vfold(jc,kc)+  &
          vrand(jc+jfiltstart+jc2,kc+kfiltstart)*ay(jc2)
      wfold(jc,kc)=wfold(jc,kc)+  &
          wrand(jc+jfiltstart+jc2,kc+kfiltstart)*ay(jc2)
    END DO
  END DO
END DO
IF(nxlprm(2)==1)THEN
  DO jc=jstal,jstol
    DO kc=kstal-nfz,kstol+nfz
      DO jc2=-nfy,nfy
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yfold(jc,kc,ispec)=yfold(jc,kc,ispec)+  &
              yrand(jc+jfiltstart+jc2,kc+kfiltstart,ispec)*ay(jc2)
        END DO
      END DO
    END DO
  END DO
END IF

DO jc=jstal,jstol
  DO kc=kstal,kstol
    DO kc2=-nfz,nfz
      ufilt(jc,kc)=ufilt(jc,kc)+ ufold(jc,kc+kc2)*az(kc2)
      vfilt(jc,kc)=vfilt(jc,kc)+ vfold(jc,kc+kc2)*az(kc2)
      wfilt(jc,kc)=wfilt(jc,kc)+ wfold(jc,kc+kc2)*az(kc2)
    END DO
  END DO
END DO
IF(nxlprm(2)==1)THEN
  DO jc=jstal,jstol
    DO kc=kstal,kstol
      DO kc2=-nfz,nfz
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yfilt(jc,kc,ispec)=yfilt(jc,kc,ispec)+  &
              yfold(jc,kc+kc2,ispec)*az(kc2)
        END DO
      END DO
    END DO
  END DO
END IF

!   Turbulence intensity
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(inflow_kernel_eqH, "inflow_kernel_eqH", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ufilt, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_vfilt, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_wfilt, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(urms, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(phi, 1, "real(kind=8)", OPS_READ))

    IF(nxlprm(2)==1)THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        DO ispc = 1,nspc
            ispec = tspc(ispc)
            yref_ispec = yref(ispec)
            call ops_par_loop(inflow_kernel_eqI, "inflow_kernel_eqI", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_yfilt(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(yrms, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(yref_ispec, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(phi, 1, "real(kind=8)", OPS_READ))
        END DO
    END IF

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(copy_kernel_xxdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_uinf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_uinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
    call ops_par_loop(copy_kernel_xxdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vinf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_vinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
    call ops_par_loop(copy_kernel_xxdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_winf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_winf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

    IF(nxlprm(2)==1)THEN
        DO ispc = 1,nspc
            ispec = tspc(ispc)
            call ops_par_loop(copy_kernel_xxdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_yinf1(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_yinf2(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
        END DO 
    END IF

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(inflow_kernel_eqA, "inflow_kernel_eqA", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_uinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_vinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_winf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_ufilt, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vfilt, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wfilt, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_uinf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vinf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_winf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(phi, 1, "real(kind=8)", OPS_READ))

    IF(nxlprm(2) == 1) THEN
        DO ispc = 1,nspc
            ispec = tspc(ispc)
            rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(inflow_kernel_eqB, "inflow_kernel_eqB", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_yinf2(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_yfilt(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_yinf1(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(phi, 1, "real(kind=8)", OPS_READ))
        END DO
    END IF

!   VM: CALCULATING USTEAD
    ymidpnt = nyglbl/2
    zmidpnt = 1

    IF (nyglbl == 1) THEN
        deltay = 0.0_8
    ELSE
        deltay = ygdlen/REAL(nyglbl-1,kind=8)
    END IF

    IF (nzglbl == 1) THEN
        deltaz = 0.0_8
    ELSE
        deltaz = zgdlen/REAL(nzglbl-1,kind=8)
    END IF

!    kg = zoffset+kc
!    jg = yoffset+jc
!    ydist = deltay*DBLE((jg-ymidpnt))
!    zdist = deltaz*DBLE((kg-zmidpnt))
!    rad=SQRT(ydist**2+zdist**2)

!    ustead(jc,kc)=0.5_8*(1.0_8-TANH((rad-exlprm(3)*ygdlen)/(0.001_8*ygdlen)))
!    ustead(jc,kc)=1.0_8   -> current assignment

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(inflow_kernel_eqC, "inflow_kernel_eqC", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ustead, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(inflow_kernel_eqD, "inflow_kernel_eqD", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_uinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_vinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_winf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_ustead, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

    IF(nxlprm(2)==1)THEN
        DO ispc = 1,nspc
            ispec = tspc(ispc)
            rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(inflow_kernel_eqE, "inflow_kernel_eqE", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_yinf2(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_ustead, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
        END DO
    END IF

!   Calculate the average speed in the entire domain, which is then
!   subtracted from u_new
    sum_umean = 0.0_8
    sum_denom = 0.0_8
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(inflow_kernel_eqF, "inflow_kernel_eqF", senga_grid, 3, rangexyz,  &
                  &  ops_arg_dat(d_ustead, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                  &  ops_arg_reduce(h_denom, 1, "real(kind=8)", OPS_INC))
    call ops_reduction_result(h_denom, sum_denom)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(inflow_kernel_eqF, "inflow_kernel_eqF", senga_grid, 3, rangexyz,  &
                  &  ops_arg_dat(d_uinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                  &  ops_arg_reduce(h_umean, 1, "real(kind=8)", OPS_INC))
    call ops_reduction_result(h_umean, sum_umean)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(inflow_kernel_eqG, "inflow_kernel_eqG", senga_grid, 3, rangexyz,  & 
                    ops_arg_dat(d_uinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                    ops_arg_gbl(sum_umean, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(sum_denom, 1, "real(kind=8)", OPS_READ))


    sum_umean = 0.0_8
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(inflow_kernel_eqF, "inflow_kernel_eqF", senga_grid, 3, rangexyz,  &
                  &  ops_arg_dat(d_vinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                  &  ops_arg_reduce(h_umean, 1, "real(kind=8)", OPS_INC))
    call ops_reduction_result(h_umean, sum_umean)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(inflow_kernel_eqG, "inflow_kernel_eqG", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                    ops_arg_gbl(sum_umean, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(sum_denom, 1, "real(kind=8)", OPS_READ))


    sum_umean = 0.0_8
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(inflow_kernel_eqF, "inflow_kernel_eqF", senga_grid, 3, rangexyz,  &
                  &  ops_arg_dat(d_winf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                  &  ops_arg_reduce(h_umean, 1, "real(kind=8)", OPS_INC))
    call ops_reduction_result(h_umean, sum_umean)

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(inflow_kernel_eqG, "inflow_kernel_eqG", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_winf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                    ops_arg_gbl(sum_umean, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(sum_denom, 1, "real(kind=8)", OPS_READ))    

!   -------------------------------------------------------------------------

END SUBROUTINE inflow

!     ==================================================================

FUNCTION ran1(idum)

INTEGER, INTENT(OUT)                     :: idum

REAL :: ran1,
INTEGER, PARAMETER :: ia=16807
INTEGER, PARAMETER :: im=2147483647
REAL, PARAMETER :: am=1./im
INTEGER, PARAMETER :: iq=127773
INTEGER, PARAMETER :: ir=2836
INTEGER, PARAMETER :: ntab=32
INTEGER, PARAMETER :: ndiv=1+(im-1)/ntab
REAL, PARAMETER :: eps=1.2E-7
REAL, PARAMETER :: rnmx=1.-eps
INTEGER :: j,k,iv(ntab),iy
SAVE iv,iy
DATA iv /ntab*0/, iy /0/

IF (idum <= 0.OR.iy == 0) THEN
  idum=MAX(-idum,1)
  DO  j=ntab+8,1,-1
    k=idum/iq
    idum=ia*(idum-k*iq)-ir*k
    IF (idum < 0) idum=idum+im
    IF (j <= ntab) iv(j)=idum
  END DO
  iy=iv(1)
END IF
k=idum/iq
idum=ia*(idum-k*iq)-ir*k
IF (idum < 0) idum=idum+im
j=1+iy/ndiv
iy=iv(j)
iv(j)=idum
ran1=MIN(am*iy,rnmx)
ran1=(ran1*2.0-1.0)/0.577
RETURN
END FUNCTION ran1
!     ==================================================================

