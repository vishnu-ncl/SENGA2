SUBROUTINE output
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-04  Time: 21:21:36

!     *************************************************************************

!     OUTPUT
!     ======

!     AUTHOR
!     ------
!     R.S.CANT - CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     28-DEC-2003:  RSC MODIFIED FOR SENGA2
!     11-JUL-2009:  RSC ADD A DUMP FORMAT SWITCH; REPORT THE DUMP
!     29-AUG-2009:  RSC UPDATE NUMBER OF PROCESSORS

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     PROCESSES THE RESULTS

!     FILES: REPORT FILE       - UNIT NCREPT - NAME FNREPT  - FORMATTED
!            DUMP OUTPUT FILES - UNIT NCDMPO - NAMES FNDMPO - UNFORMATTED
!            STATISTICS FILE   - UNIT NCSTAT - NAME FNSTAT  - FORMATTED

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
#ifdef HDF5
use hdf5io
#endif
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL PARAMETERS
!     ================
!     DIAGNOSTICS
CHARACTER (LEN=4) :: pnxres
PARAMETER(pnxres = '.res')
CHARACTER (LEN=2) :: pnxryf
PARAMETER(pnxryf = 'yf')
INTEGER :: ncdiag
PARAMETER(ncdiag = 11)
INTEGER :: iddump

!     LOCAL DATA
!     ==========
!     DIAGNOSTICS
real(kind=8) :: deltag,fornow
real(kind=8) :: ttemp(nxsize,nysize,nzsize)
real(kind=8) :: ptemp(nxsize,nysize,nzsize)
real(kind=8) :: ytemp(nxsize,nysize,nzsize,nspec)
real(kind=8) :: utgv(nxsize,nysize,nzsize)
real(kind=8) :: vtgv(nxsize,nysize,nzsize)
real(kind=8) :: wtgv(nxsize,nysize,nzsize)
real(kind=8) :: vort1(nxsize,nysize,nzsize)
real(kind=8) :: vort2(nxsize,nysize,nzsize)
real(kind=8) :: vort3(nxsize,nysize,nzsize)
real(kind=8) :: enstro(nxsize,nysize,nzsize)
real(kind=8) :: tkel(nxsize,nysize,nzsize)

real(kind=8) :: tkes,tkeg,enstrs,enstrg

INTEGER :: ispec
INTEGER :: ic,jc,kc
INTEGER :: igofst,ix
CHARACTER (LEN=21) :: fndiag
!     RSC UPDATE NUMBER OF PROCESSORS
CHARACTER (LEN=6) :: pnproc
CHARACTER (LEN=11) :: strqty
CHARACTER (LEN=2) :: strspc
LOGICAL :: bcflag

CHARACTER (LEN=60) :: fname
CHARACTER (LEN=4) :: proc
CHARACTER (LEN=5) :: ipdump

!     BEGIN
!     =====

!     =========================================================================

!     REPORT OUTPUT
!     =============
!ops IF(MOD(itime,ntrept) == 0)THEN
  
!       REPORT ON PROCESSOR NO.1 ONLY
!       ------
!ops  IF(iproc == 0)THEN
    
!ops    OPEN(UNIT=ncrept,FILE=fnrept,STATUS='OLD',FORM='FORMATTED')
    
!         GO TO EOF
!ops    1000      CONTINUE
!ops    READ(ncrept,9000,END=1010)
!ops    GO TO 1000
!ops    1010      BACKSPACE(ncrept)
    
!ops    WRITE(ncrept,9100)itime
!ops    WRITE(ncrept,9110)etime,tstep
!ops    CLOSE(ncrept)
    
!ops  END IF
  
!       =======================================================================
  
!       DIAGNOSTICS
!ops  jc = MAX(nyglbl/2,1)
!ops  kc = MAX(nzglbl/2,1)
  
!ops  WRITE(pnproc,'(I6.6)')iproc
  
!       GLOBAL INDEXING
!ops  deltag = xgdlen/(REAL(nxglbl-1))
  
!ops  igofst = 0
!ops  DO ic = 0, ixproc-1
!ops    igofst = igofst + npmapx(ic)
!ops  END DO
  
!        STRQTY = 'output/pres'
!        FNDIAG = STRQTY//PNPROC//PNXRES
!        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C       GO TO EOF
!8000    CONTINUE
!          READ(NCDIAG,9000,END=8001)
!          GOTO 8000
!8001    BACKSPACE(NCDIAG)
!        DO IC = ISTAL,ISTOL
!          IX = IGOFST + IC
!          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,PRUN(IC,JC,KC)
!        ENDDO
!        WRITE(NCDIAG,*)
!        CLOSE(NCDIAG)
  
!        STRQTY = 'output/uvel'
!        FNDIAG = STRQTY//PNPROC//PNXRES
!        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C       GO TO EOF
!8010    CONTINUE
!          READ(NCDIAG,9000,END=8011)
!          GOTO 8010
!8011    BACKSPACE(NCDIAG)
!        DO IC = ISTAL,ISTOL
!          IX = IGOFST + IC
!          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,
!     +                      URUN(IC,JC,KC)/DRUN(IC,JC,KC)
!        ENDDO
!        WRITE(NCDIAG,*)
!        CLOSE(NCDIAG)
  
!        STRQTY = 'output/temp'
!        FNDIAG = STRQTY//PNPROC//PNXRES
!        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C       GO TO EOF
!8020    CONTINUE
!          READ(NCDIAG,9000,END=8021)
!          GOTO 8020
!8021    BACKSPACE(NCDIAG)
!        DO IC = ISTAL,ISTOL
!          IX = IGOFST + IC
!          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,TRUN(IC,JC,KC)
!        ENDDO
!        WRITE(NCDIAG,*)
!        CLOSE(NCDIAG)
  
!        STRQTY = 'output/dens'
!        FNDIAG = STRQTY//PNPROC//PNXRES
!        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C       GO TO EOF
!8030    CONTINUE
!          READ(NCDIAG,9000,END=8031)
!          GOTO 8030
!8031    BACKSPACE(NCDIAG)
!        DO IC = ISTAL,ISTOL
!          IX = IGOFST + IC
!          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,DRUN(IC,JC,KC)
!        ENDDO
!        WRITE(NCDIAG,*)
!        CLOSE(NCDIAG)
  
!        STRQTY = 'output/ener'
!        FNDIAG = STRQTY//PNPROC//PNXRES
!        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C       GO TO EOF
!8040    CONTINUE
!          READ(NCDIAG,9000,END=8041)
!          GOTO 8040
!8041    BACKSPACE(NCDIAG)
!        DO IC = ISTAL,ISTOL
!          IX = IGOFST + IC
!          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,
!     +                      ERUN(IC,JC,KC)/DRUN(IC,JC,KC)
!        ENDDO
!        WRITE(NCDIAG,*)
!        CLOSE(NCDIAG)
  
!        DO ISPEC = 1, NSPEC
  
!          WRITE(STRSPC,'(I2.2)')ISPEC
!          STRQTY = 'output/'//PNXRYF//STRSPC
!          FNDIAG = STRQTY//PNPROC//PNXRES
!          OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C         GO TO EOF
!8050      CONTINUE
!            READ(NCDIAG,9000,END=8051)
!            GOTO 8050
!8051      BACKSPACE(NCDIAG)
!          DO IC = ISTAL,ISTOL
!            IX = IGOFST + IC
!              FORNOW = YRUN(IC,JC,KC,ISPEC)/DRUN(IC,JC,KC)
!              FORNOW = MAX(FORNOW,1.0D-30)
!            WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,FORNOW
!C     +                        YRUN(IC,JC,KC,ISPEC)/DRUN(IC,JC,KC)
!          ENDDO
!          WRITE(NCDIAG,*)
!          CLOSE(NCDIAG)
  
!        ENDDO
  
!C        ISPEC = NSPEC
!C        WRITE(STRSPC,'(I2.2)')ISPEC
!C        STRQTY = PNXRYF//STRSPC
!C        FNDIAG = STRQTY//PNPROC//PNXRES
!C        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!CC       GO TO EOF
!C8050    CONTINUE
!C          READ(NCDIAG,9000,END=8051)
!C          GOTO 8050
!C8051    BACKSPACE(NCDIAG)
!C        DO IC = ISTAL,ISTOL
!C          IX = IGOFST + IC
!C          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,
!C     +                      YRUN(IC,JC,KC,ISPEC)/DRUN(IC,JC,KC)
!C        ENDDO
!C        WRITE(NCDIAG,*)
!C        CLOSE(NCDIAG)
  
!ops END IF

!     =========================================================================
!     UMOD START
!     DATA OUTPUT FOR POST-PROCESSING
!     ======
!ops  IF(MOD(itime,ntdump) == 0)THEN
!ops  iddump=itime/ntdump
!ops  WRITE(ipdump,'(I5.5)') iddump
!ops  WRITE(proc,'(I4.4)') iproc
  
!ops  fname = 'output/out'//ipdump//proc//pnxres
  
  
!ops  DO kc = 1,nzsize
!ops    DO jc = 1,nysize
!ops      DO ic = 1,nxsize
!ops        DO ispec =1,nspec
!ops          ytemp(ic,jc,kc,ispec)=yrun(ic,jc,kc,ispec)/drun(ic,jc,kc)
!ops        END DO
!ops        ttemp(ic,jc,kc)=trun(ic,jc,kc)
!ops        ptemp(ic,jc,kc)=prun(ic,jc,kc)
!ops      END DO
!ops    END DO
!ops  END DO
  
!ops  OPEN(UNIT=16,FILE=trim(fname),FORM='UNFORMATTED',STATUS='NEW')
!ops  WRITE(16)drun,urun/drun,vrun/drun,wrun/drun,erun/drun,ttemp,  &
!ops      ptemp,ytemp,rrte,etime
!ops  CLOSE(16)
!ops  END IF



!     UMOD END
!     =====
!     FULL DUMP OUTPUT
!     ================
  IF(MOD(itime,ntdump) == 0)THEN
  
!       CARRY OUT A FULL DUMP
!       ---------------------
!       USE THE DUMP FILE INDICATED BY IDFLAG
!       RSC 11-JUL-2009 ADD A DUMP FORMAT SWITCH
#ifndef HDF5
  IF(ndofmt == 0)THEN
   
    write(*,*) "Timestep: ", itime, "  Writing output to: ",fndmpo(idflag+1)
 
!         UNFORMATTED DUMP OUTPUT
    OPEN(UNIT=ncdmpo,FILE=fndmpo(idflag+1),STATUS='OLD', FORM='UNFORMATTED')
    REWIND(ncdmpo)
    WRITE(ncdmpo)nxnode,nynode,nznode,nspec, drun,urun,vrun,wrun,erun,yrun,  &
        etime,tstep,errold,errldr
    CLOSE(ncdmpo)
    
  ELSE
    
!         FORMATTED DUMP OUTPUT
    OPEN(UNIT=ncdmpo,FILE=fndmpo(idflag+1),STATUS='OLD', FORM='FORMATTED')
    REWIND(ncdmpo)
    WRITE(ncdmpo,*)nxnode,nynode,nznode,nspec
    DO kc = 1,nznode
      DO jc = 1,nynode
        DO ic = 1,nxnode
          WRITE(ncdmpo,*)drun(ic,jc,kc),  &
              urun(ic,jc,kc),vrun(ic,jc,kc),wrun(ic,jc,kc), erun(ic,jc,kc),  &
              (yrun(ic,jc,kc,ispec),ispec=1,nspec)
        END DO
      END DO
    END DO
    WRITE(ncdmpo,*)etime,tstep,errold,errldr
    CLOSE(ncdmpo)
    
  END IF
  
#else
  CALL write_h5_dumpfile
#endif

!       REPORT THE DUMP
!       RSC 11-JUL-2009
IF(iproc == 0)THEN
  
  OPEN(UNIT=ncrept,FILE=fnrept,STATUS='OLD',FORM='FORMATTED')
 3000      CONTINUE
 READ(ncrept,9000,END=3010)
 GO TO 3000
 3010      BACKSPACE(ncrept)
 WRITE(ncrept,9120)fndmpo(idflag+1)
CLOSE(ncrept)
  
END IF

!       RESET THE DUMP FLAG
idflag = MOD(idflag+1,2)

END IF

!     =========================================================================

!     DUMP BC INFORMATION AS REQUIRED
!     ===============================
!ops IF(MOD(itime,ntdump) == 0)THEN
  
!ops  bcflag = (nsbcxl == nsbci2).OR.(nsbcxl == nsbci3)
!ops  bcflag = bcflag.AND.(nxlprm(1) == 3)
  
!ops  IF(bcflag)THEN
    
!         DUMP THE INLET TURBULENT VELOCITY FIELD
!ops    OPEN(UNIT=nctixl,FILE=fntixl,STATUS='OLD', FORM='UNFORMATTED')
!ops    REWIND(nctixl)
!ops    WRITE(nctixl)ufxl,vfxl,wfxl,slocxl,svelxl,bvelxl
!ops    CLOSE(nctixl)
    
!ops  END IF
  
!ops END IF
!     ==========
!     umod tgv stats 
!ops deltagx = xgdlen/(real(nxglbl-1))
!ops deltagy = ygdlen/(real(nyglbl-1))
!ops deltagz = zgdlen/(real(nzglbl-1))

!ops do kc=kstal,kstol
!ops  do jc=jstal,jstol
!ops    do ic=istal,istol
!ops      utgv(ic,jc,kc)=urun(ic,jc,kc)/drun(ic,jc,kc)
!ops      vtgv(ic,jc,kc)=vrun(ic,jc,kc)/drun(ic,jc,kc)
!ops      wtgv(ic,jc,kc)=wrun(ic,jc,kc)/drun(ic,jc,kc)
!ops    end do
!ops  end do
!ops end do


!kinetic energy
!====================

!ops do kc = kstal, kstol
!ops  do jc = jstal, jstol
!ops    do ic = istal, istol

!ops      tkel(ic,jc,kc)=(utgv(ic,jc,kc)*utgv(ic,jc,kc))+ &
!ops                    (vtgv(ic,jc,kc)*vtgv(ic,jc,kc))+ &
!ops                    (wtgv(ic,jc,kc)*wtgv(ic,jc,kc))
!ops
!ops    end do
!ops  end do
!ops end do

!ops tkes=0.0d0

!ops do kc = kstal, kstol
!ops  do jc = jstal, jstol
!ops    do ic = istal, istol

!ops     tkes=tkes+(0.5d0*drun(ic,jc,kc)*tkel(ic,jc,kc) &
!ops               *deltagx*deltagy*deltagz)

!ops    end do
!ops  end do
!ops end do
!ops call p_summ(tkes,tkeg)

!fornow = drin*xgdlen*ygdlen*zgdlen
!fornow = 1.d0/fornow
!tkeg=tkeg*fornow
!
!===================================================
!enstrophy calculation      


!ops do kc = kstal, kstol
!ops do jc = jstal, jstol 
!ops  do ic = istal, istol

!ops   vort1(ic,jc,kc)= dwtgvdy(ic,jc,kc)-dvtgvdz(ic,jc,kc)
!ops   vort2(ic,jc,kc)= dutgvdz(ic,jc,kc)-dwtgvdx(ic,jc,kc)
!ops   vort3(ic,jc,kc)= dvtgvdx(ic,jc,kc)-dutgvdy(ic,jc,kc)

!ops   enstro(ic,jc,kc)=(vort1(ic,jc,kc)*vort1(ic,jc,kc))+&
!ops              (vort2(ic,jc,kc)*vort2(ic,jc,kc))+&
!ops              (vort3(ic,jc,kc)*vort3(ic,jc,kc)) 

!ops  end do
!ops end do
!ops end do

!ops enstrs=0.0d0

!ops do kc = kstal, kstol
!ops do jc = jstal, jstol 
!ops  do ic = istal, istol

!ops   enstrs=enstrs+(0.5d0*drun(ic,jc,kc)*enstro(ic,jc,kc) &
!ops            *deltagx*deltagy*deltagz)

!ops  end do
!ops end do
!ops end do

!ops call p_summ(enstrs,enstrg)

!ops fornow = drin*xgdlen*ygdlen*zgdlen
!ops fornow = 1.d0/fornow
!enstrg=enstrg*fornow
!ops if(iproc.eq.0)print*,tkeg, enstrg*fornow,enstrs 
!ops if(iproc.eq.0) then
!ops  open(44,file='output/tgv_stat.res',access='append')
!ops    write(44,'(3e20.9)') etime,tkeg*fornow,enstrg*fornow 
!ops  close(44)
!ops endif
!umod tgv stats
!========================

!     =========================================================================

!     TIME STEP HISTORY
IF(iproc == 0)THEN
  WRITE(*,'(I7,1PE12.4,I5)')itime,tstep,inderr
END IF

!**ashutosh**!IF(iproc == 5)THEN
!**ashutosh**!  WRITE(*,'(a,I7,a,F16.8)')  &
!**ashutosh**!      'test_drhs: (step=',itime,') value: ',drhs(5,jstal,kstal)
!**ashutosh**!  WRITE(*,'(a,I7,a,F16.8)')  &
!**ashutosh**!      'test_erhs: (step=',itime,') value: ',erhs(6,jstal,kstal)
!**ashutosh**!  WRITE(*,'(a,I7,a,F16.8)')  &
!**ashutosh**!      'test_urhs: (step=',itime,') value: ',urhs(7,jstal,kstal)
!**ashutosh**!END IF

!     =========================================================================

!     STATISTICS ON THE FLY
!     =====================

!     STATISTICS MASTER SWITCH
!     ------------------------
!ops IF(ntstat >= 0)THEN
  
!     =========================================================================
  
!       OUTPUT STATISTICS
!       =================
  
!       STATISTICS ON ONE PROCESSOR ONLY
!       ----------
!ops  IF(iproc == 0)THEN
    
!ops    IF(MOD(itime,ntstat) == 0)THEN
      
!ops      OPEN(UNIT=ncstat,FILE=fnstat,STATUS='OLD',FORM='FORMATTED')
      
!           GO TO EOF
!ops      2000        CONTINUE
!ops      READ(ncstat,9200,END=2010)
!ops      GO TO 2000
!ops      2010        BACKSPACE(ncstat)
      
!ops      WRITE(ncstat,9100)itime
      
!ops      CLOSE(ncstat)
      
!ops    END IF
    
!ops  END IF
  
!       RESET STORAGE INDEX
!ops  itstat = 0
  
!     =========================================================================
  
!     STATISTICS MASTER SWITCH
!ops END IF

!     =========================================================================


!ops RETURN

 9000  FORMAT(a)
 9100  FORMAT('Time step number: ',i7)
 9110  FORMAT('Elapsed time: ',1PE12.4,';',2X,'next time step:',1PE12.4)
 9120  FORMAT('Dump completed: ',a)
 9200  FORMAT(i5,/ 5X,6(1PE12.4)/  &
    5X,6(1PE12.4)/ 5X,6(1PE12.4)/  &
   5X,4(1PE12.4)/ 5X,3(1PE12.4)/  &
    5X,3(1PE12.4)/ 5X,3(1PE12.4)/  &
    5X,3(1PE12.4)/ 5X,3(1PE12.4)/  &
    5X,2(1PE12.4))
 9300  FORMAT(2(1PE15.7))

END SUBROUTINE output
