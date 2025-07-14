SUBROUTINE indata

    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   INDATA
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   12-MAR-2003:  RSC UPDATED FOR SENGA2
!   11-JUL-2009:  RSC ADD A DUMP FORMAT SWITCH
!   29-AUG-2009:  RSC UPDATE NUMBER OF PROCESSORS
!   06-JAN-2013:  RSC MIXTURE-AVERAGED TRANSPORT
!   14-JUL-2013:  RSC RADIATION TREATMENT

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   INITIALISES ALL DATA

!   FILES: CHEMICAL DATA FILE  - UNIT NCCHEM - FORMATTED INPUT
!          TRANSPORT DATA FILE - UNIT NCDIFF - FORMATTED INPUT
!          RADIATION DATA FILE - UNIT NCRADN - FORMATTED INPUT
!          REPORT FILE         - UNIT NCREPT - FORMATTED OUTPUT
!          STATISTICS FILE     - UNIT NCSTAT - FORMATTED OUTPUT
!          DUMP OUTPUT FILES   - UNIT NCDMPO - UNFORMATTED
!          DUMP INPUT FILES    - UNIT NCDMPI - UNFORMATTED  (RESTART ONLY)

!   *************************************************************************


!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
#ifdef HDF5
    use hdf5io
#endif
    use com_espect
!   -------------------------------------------------------------------------

!   PARAMETERS
!   ==========
!   FILENAME COMPONENTS
    CHARACTER (LEN=10) :: pncont,pnchem,pndiff,pnradn
    CHARACTER (LEN=11) :: pnrept,pnstat,pndmpi,pndmpo
    CHARACTER (LEN=4) :: pnxdat,pnxres
    PARAMETER(pncont = 'input/cont', pnchem = 'input/chem',  &
              pndiff = 'input/diff', pnradn = 'input/radn',  &
              pnrept = 'output/rept', pnstat = 'output/stat',  &
              pndmpi = 'output/dmpi', pndmpo = 'output/dmpo',  &
              pnxdat = '.dat', pnxres = '.res')

!   REPORT FILE LINE LENGTH
    integer(kind=4), parameter :: lenlin = 79

!   THERMO DATA INTERVAL-MATCHING TOLERANCE
    real(kind=8), parameter :: tttol = 0.0010_8


!   LOCAL DATA
!   ==========
    real(kind=8) :: dyrin(nspcmx)
    real(kind=8) :: ctrans(nspcmx,4)
    real(kind=8) :: dtrans(nspcmx)
    real(kind=8) :: rtin,durin,dvrin,dwrin,derin
    real(kind=8) :: ttemp(5),ttold(5)
    real(kind=8) :: fornow,tmplog,tsold
    real(kind=8) :: combo1,combo2,combo3
    integer(kind=4) :: iindex,ipower,icoeff,icoef1,icoef2
    integer(kind=4) :: ic,jc,kc,ispec,istep,itint,icp
    integer(kind=4) :: jspec,ncount
    integer(kind=4) :: nxdmax,nydmax,nzdmax,ndspec
    integer(kind=4) :: rangexyz(6)
    CHARACTER (LEN=132) :: strout
    CHARACTER (LEN=79) :: clstar,cldash
    CHARACTER (LEN=10) :: streac
    CHARACTER (LEN=5) :: strstp
    CHARACTER (LEN=4) :: straro
    CHARACTER (LEN=1) :: strcof
!   RSC UPDATE NUMBER OF PROCESSORS
    CHARACTER (LEN=6) :: pnproc
    CHARACTER (LEN=1) :: pnflag
    LOGICAL :: fxdump
    LOGICAL :: flgreq

    integer(kind=4) :: dtime
    character(len=60) :: fname
    character(len=3) :: pnxhdf
    parameter(pnxhdf = '.h5')
    character(len=8) :: citime

!   BEGIN
!   =====

!   =========================================================================

!   SET UP HDF5 I/O
!   ===============
#ifdef HDF5
    call hdf5_init
#endif
!   SET UP FILE I/O
!   ===============

!   UNIT NUMBERS
    nccont = 1
    ncchem = 2
    ncdiff = 3
    ncradn = 4
    ncdmpo = 7
    ncrept = 8
    ncstat = 9

!   BUILD THE FILENAMES
    WRITE(pnproc,'(I6.6)')iproc
    fncont = pncont//pnxdat
    fnchem = pnchem//pnxdat
    fndiff = pndiff//pnxdat
    fnradn = pnradn//pnxdat
    fnrept = pnrept//pnxres
    fnstat = pnstat//pnxres
    idflag = 0
    WRITE(pnflag,'(I1)')idflag
    fndmpo(1) = pndmpi//pnproc//pnflag//pnxdat
    idflag = 1
    WRITE(pnflag,'(I1)')idflag
    fndmpo(2) = pndmpi//pnproc//pnflag//pnxdat

!   =========================================================================

!   GET THE RUN CONTROL DATA
!   ========================
    call contin

    IF(nxlprm(1) == 4) THEN
!       VM: SYNTHETIC DIGITAL FILTERING METHOD
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_uinf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_vinf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_winf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        DO ispec=1,nspec
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_yinf1(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        END DO

!       ATTENTION WHEN RESTART: THE SAME FILE EXTENSION MUST BE ENTERED HERE AS IN INFLOW!
        IF (ncdmpi == 1)THEN
!           ACTUALLY, INTRAN0.DAT OR INTRAN1.DAT SHOULD BE HERE!
            OPEN(UNIT=16,FILE='output/intran.dat',FORM='FORMATTED')
            READ(16,*)intran
            CLOSE(16)
            passer=REAL(intran,kind=8)
        ELSE
            intran=-1
        END IF

        IF(ncdmpi == 1) THEN
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
            DO ispec=1,nspec
                call ops_par_loop(copy_kernel_xxdir, "copy", senga_grid, 3, rangexyz, &
                                ops_arg_dat(d_yinf1(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                                ops_arg_dat(d_yinf2(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
            END DO
        ELSE
            rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_uinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_vinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_winf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
            DO ispec=1,nspec
                call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yinf2(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
            END DO
        END IF

    END IF

!   =========================================================================

!   WRITE THE INITIAL REPORT
!   ========================

    IF(iproc == 0)THEN

!       SET UP SEPARATOR LINES
        DO ic = 1,lenlin
            clstar(ic:ic) = '*'
            cldash(ic:ic) = '-'
        END DO

        OPEN(UNIT=ncrept,FILE=fnrept,STATUS='REPLACE',FORM='FORMATTED')

!       WRITE A HEADER
        WRITE(ncrept,9000)clstar
        WRITE(ncrept,9010)
        WRITE(ncrept,9020)
        WRITE(ncrept,9010)
        WRITE(ncrept,9030)
        WRITE(ncrept,9010)
        WRITE(ncrept,9040)
        WRITE(ncrept,9010)
        WRITE(ncrept,9000)clstar
        WRITE(ncrept,*)

        WRITE(ncrept,*)'==========='
        WRITE(ncrept,*)'REPORT FILE'
        WRITE(ncrept,*)'==========='
        WRITE(ncrept,*)

        WRITE(ncrept,9000)cldash

!       DOMAIN DECOMPOSITION DATA AND CHECKING
        WRITE(ncrept,*)
        WRITE(ncrept,*)'Domain decomposition data:'
        WRITE(ncrept,*)'-------------------------'
        WRITE(ncrept,*)

!       GLOBAL GRID SIZES
        WRITE(ncrept,*)'NXGLBL,NYGLBL,NZGLBL'
        WRITE(ncrept,9200)nxglbl,nyglbl,nzglbl
        WRITE(ncrept,*)

!       CHECK GLOBAL GRID SIZE
        flgreq = (nxgreq == nxglbl) .AND.(nygreq == nyglbl) .AND.(nzgreq == nzglbl)
        IF(.NOT.flgreq) THEN
            WRITE(ncrept,*)'Warning: INDATA: global domain size mismatch'
            WRITE(ncrept,*)'Required values in control file:'
            WRITE(ncrept,9200)nxgreq,nygreq,nzgreq
            WRITE(ncrept,*)'Global values take precedence'
            WRITE(ncrept,*)
        END IF

!       NUMBERS OF PROCESSORS
        WRITE(ncrept,*)'NXPROC,NYPROC,NZPROC,NPROC'
        WRITE(ncrept,9200)nxproc,nyproc,nzproc,nproc
        WRITE(ncrept,*)
        IF(nproc /= nxproc*nyproc*nzproc) THEN
            WRITE(ncrept,*) 'Warning: INDATA: mismatched total no. of processors'
            WRITE(ncrept,*)
        END IF

!       CHECK NUMBER OF PROCESSORS
        flgreq = (nxpreq == nxproc) .AND.(nypreq == nyproc) .AND.(nzpreq == nzproc)
        IF(.NOT.flgreq) THEN
            WRITE(ncrept,*)'Warning: INDATA: no. of processors mismatch'
            WRITE(ncrept,*)'Required values in control file:'
            WRITE(ncrept,9200)nxpreq,nypreq,nzpreq
            WRITE(ncrept,*)'Global values take precedence'
            WRITE(ncrept,*)
        END IF

!       LOCAL ARRAY SIZES
        WRITE(ncrept,*)'NXSIZE,NYSIZE,NZSIZE'
        WRITE(ncrept,9200)nxsize,nysize,nzsize
        WRITE(ncrept,*)

!       LOCAL GRID SIZES
        WRITE(ncrept,*)'NPMAPX'
        WRITE(ncrept,9200)(npmapx(ic),ic=0,nxprm1)
        WRITE(ncrept,*)'NPMAPY'
        WRITE(ncrept,9200)(npmapy(ic),ic=0,nyprm1)
        WRITE(ncrept,*)'NPMAPZ'
        WRITE(ncrept,9200)(npmapz(ic),ic=0,nzprm1)

        WRITE(ncrept,*)'--------------------------------'
        WRITE(ncrept,*)
        WRITE(ncrept,9000)cldash

!       =======================================================================

!       INITIAL DATA
        WRITE(ncrept,*)
        WRITE(ncrept,*)'Initial data:'
        WRITE(ncrept,*)'------------'
        WRITE(ncrept,*)

!       GLOBAL DOMAIN SIZE (x,y,z)
        WRITE(ncrept,*)'XGDLEN,YGDLEN,ZGDLEN'
        WRITE(ncrept,9400)xgdlen,ygdlen,zgdlen
        WRITE(ncrept,*)

!       NUMERICAL CONTROL
        WRITE(ncrept,*)'TSTEP,NTIME1,NTIME,NSTPSW,NTDUMP,NTREPT,NTSTAT'
        WRITE(ncrept,9300)tstep,ntime1,ntime,nstpsw,ntdump,ntrept,ntstat
        WRITE(ncrept,*)

!       COLD START SWITCH
        WRITE(ncrept,*)'NCDMPI'
        WRITE(ncrept,9200)ncdmpi

!       INITIAL TURBULENCE GENERATOR
        WRITE(ncrept,*)'INTURB, INSEED, SPARAM'
        WRITE(ncrept,9270)inturb,inseed,(sparam(ic),ic=1,nsparm)

!       FLAME START OPTION
        WRITE(ncrept,*)'INFLAM'
        WRITE(ncrept,9200)inflam

!       DEFAULT INITIAL DATA
        WRITE(ncrept,*)
        WRITE(ncrept,*)'PRIN,TRIN,URIN,VRIN,WRIN'
        WRITE(ncrept,9400)prin,trin,urin,vrin,wrin

!       CHECK NUMBER OF SPECIES
        IF(nspreq /= nspcmx) THEN
            WRITE(ncrept,*)'Warning: INDATA: no. of species mismatch'
            WRITE(ncrept,*)'Global value:'
            WRITE(ncrept,9200)nspcmx
            WRITE(ncrept,*)'Required value in control file:'
            WRITE(ncrept,9200)nspreq
            WRITE(ncrept,*)'Global value takes precedence'
            WRITE(ncrept,*)
        END IF

!       INITIAL SPECIES MASS FRACTIONS
        WRITE(ncrept,*)'ISPEC,YRIN'
        DO ispec = 1,nspec
            WRITE(ncrept,9450)ispec,yrin(ispec)
        END DO

!       CHECK INITIAL SPECIES MASS FRACTIONS
        fornow = zero
        DO ispec = 1,nspec
            fornow = fornow + yrin(ispec)
        END DO
        IF(ABS(fornow-one) > ytoler) THEN
            WRITE(ncrept,*)  &
                        'Warning: INDATA: initial mass fractions do not sum to unity'
            WRITE(ncrept,*)
        END IF

!       GLOBAL BOUNDARY CONDITION TYPES
        WRITE(ncrept,*)
        WRITE(ncrept,*)'NGBCXL,NGBCXR,NGBCYL,NGBCYR,NGBCZL,NGBCZR'
        WRITE(ncrept,9250)ngbcxl,ngbcxr,ngbcyl,ngbcyr,ngbczl,ngbczr

        WRITE(ncrept,*)
        WRITE(ncrept,*)'End of initial data'
        WRITE(ncrept,*)'-------------------'
        WRITE(ncrept,*)
        WRITE(ncrept,9000)cldash

    END IF

!   =========================================================================

!   GET THE CHEMICAL DATA
!   =====================
    call chemin

!   =========================================================================

    IF(iproc == 0) THEN

!       REPORT THE CHEMICAL DATA
!       ========================

        WRITE(ncrept,*)
        WRITE(ncrept,*)'Chemical data:'
        WRITE(ncrept,*)'-------------'
        WRITE(ncrept,*)

!       SPECIES LIST
        WRITE(ncrept,*)'Number of species:'
        WRITE(ncrept,'(I5)')nspec

!       CHECK NUMBER OF SPECIES
        IF(nspec /= nspcmx) THEN
            WRITE(ncrept,*)'Warning: INDATA: no. of species mismatch'
            WRITE(ncrept,*)'Global value:'
            WRITE(ncrept,9200)nspcmx
            WRITE(ncrept,*)'Required value in chemical data file:'
            WRITE(ncrept,9200)nspec
            WRITE(ncrept,*)'Chemical data file value takes precedence'
            WRITE(ncrept,*)
        END IF
        WRITE(ncrept,*)'List of species:'
        WRITE(ncrept,*)'  No.  Symbol         Mol. Mass'
        DO ispec = 1,nspec
            WRITE(ncrept,'(I5,3X,A10,3X,1PE12.4)') ispec,spcsym(ispec),wmolar(ispec)
        END DO
        WRITE(ncrept,*)

!       THERMODYNAMIC DATA
        WRITE(ncrept,*)'Species thermodynamic data:'
        WRITE(ncrept,'("  Reference pressure: ",1PE12.4)')prefgb
        WRITE(ncrept,*)' Spec.  No of T intervals   ',  &
                    'Interval   T low       T high       No of coeffs'
        DO ispec = 1,nspec
            ic = 1
            WRITE(ncrept,'(I5,6X,I5,9X,I5,8X,2(1PE12.4),I8)') ispec,ntint(ispec),  &
                                 ic,tintlo(ic,ispec),tinthi(ic,ispec),ncofcp(ic,ispec)
            DO ic = 2,ntint(ispec)
                WRITE(ncrept,'(25X,I5,8X,2(1PE12.4),I8)')  &
                                ic,tintlo(ic,ispec),tinthi(ic,ispec),ncofcp(ic,ispec)

            END DO
        END DO
        WRITE(ncrept,*)'Cp coeffs by mass'
        WRITE(ncrept,*)' Spec.  T int.  Coeff no.  Coeff.'
        DO ispec = 1,nspec
            DO ic = 1,ntint(ispec)
                jc = 1
                WRITE(ncrept,'(I5,2X,I5,5X,I5,4X,1PE15.7)')  &
                                ispec,ic,jc,amascp(jc,ic,ispec)
                DO jc = 2,ncofcp(ic,ispec)
                    WRITE(ncrept,'(17X,I5,4X,1PE15.7)')jc,amascp(jc,ic,ispec)
                END DO
            END DO
        END DO

        WRITE(ncrept,*)
        WRITE(ncrept,*)'Mass-specific Cp, Enthalpy, Entropy;', '  Molar Gibbs fn.:'
        WRITE(ncrept,*)' Spec.  T int.  Temp.       Cp',  &
                    '          Enthalpy    Entropy     Molar Gibbs'
        DO ispec = 1,nspec
            DO ic = 1,ntint(ispec)
                DO icp = 1,2

!                   TEMPERATURE
                    IF(icp == 1) ttemp(1) = tintlo(ic,ispec)
                    IF(icp == 2) ttemp(1) = tinthi(ic,ispec)

!                   MASS-SPECIFIC CP
                    fornow = amascp(ncpoly(ic,ispec),ic,ispec)
                    DO jc = ncpom1(ic,ispec),1,-1
                        fornow = fornow*ttemp(1) + amascp(jc,ic,ispec)
                    END DO
                    ttemp(2) = fornow

!                   MASS-SPECIFIC ENTHALPY
                    fornow = amasch(ncpoly(ic,ispec),ic,ispec)
                    DO jc = ncpom1(ic,ispec),1,-1
                        fornow = fornow*ttemp(1) + amasch(jc,ic,ispec)
                    END DO
                    fornow = amasch(ncenth(ic,ispec),ic,ispec) + fornow*ttemp(1)
                    ttemp(3) = fornow

!                   MASS-SPECIFIC ENTROPY
                    fornow = amascs(ncpoly(ic,ispec),ic,ispec)
                    DO jc = ncpom1(ic,ispec),2,-1
                        fornow = fornow*ttemp(1) + amascs(jc,ic,ispec)
                    END DO
                    fornow = amascs(ncenpy(ic,ispec),ic,ispec) + fornow*ttemp(1)  &
                           + amascs(1,ic,ispec)*LOG(ttemp(1))
                    ttemp(4) = fornow

!                   MOLAR GIBBS FUNCTION
!                   ACTUALLY GIBBS/(R^0 T) WITH PRESSURE TERM
                    fornow = amolgb(ncpoly(ic,ispec),ic,ispec)
                    DO jc = ncpom1(ic,ispec),1,-1
                        fornow = amolgb(jc,ic,ispec) + fornow*ttemp(1)
                    END DO
                    fornow = amolgb(ncenth(ic,ispec),ic,ispec) /ttemp(1)  &
                           - amolgb(ncenpy(ic,ispec),ic,ispec) *LOG(ttemp(1))  &
                           - fornow
                    ttemp(5) = fornow

                    IF(icp == 1) THEN
                        IF(ic == 1) THEN
                            WRITE(ncrept,'(I5,2X,I5,X,"l",X,5(1PE12.4))')  &
                                            ispec,ic,(ttemp(jc),jc=1,5)
                        ELSE
                            WRITE(ncrept,'(7X,I5,X,"l",X,5(1PE12.4))') ic,(ttemp(jc),jc=1,5)
                            DO jc = 1,5
                                IF(ABS(ttold(jc)-ttemp(jc)) > ABS(tttol*ttemp(jc))) THEN
                                    WRITE(ncrept,*)'Warning: INDATA: Mismatched thermo data'
                                    WRITE(ncrept,'(I7,1PE12.4)')jc,ttemp(jc)
                                END IF
                            END DO
                        END IF
                    ELSE
                        WRITE(ncrept,'(7X,I5,X,"h",X,5(1PE12.4))') ic,(ttemp(jc),jc=1,5)
                    END IF

                    IF(icp == 2) THEN
                        DO jc = 1,5
                            ttold(jc) = ttemp(jc)
                        END DO
                    END IF

                END DO
            END DO
        END DO
        WRITE(ncrept,*)

!       REACTION MECHANISM
!       STEP LIST
        WRITE(ncrept,*)'Reaction mechanism:'
        WRITE(ncrept,*)'  Number of steps:'
        WRITE(ncrept,'(I5)')nstep
        IF(nstep /= nstpmx) THEN
            WRITE(ncrept,*)'Warning: mismatch in number of steps:'
            WRITE(ncrept,9600)nstep,nstpmx
        END IF
        DO istep = 1,nstep

            WRITE(strstp,'(I5)')istep
            icoef1 = 1
            icoef2 = LEN(strstp) + 3
            strout(icoef1:icoef2) = strstp//'   '
            DO ispec = 1,nrslen(istep)
                streac = spcsym(nrspec(ispec,istep))
                icoef1 = icoef2 + 1
                icoef2 = icoef1 + LEN(streac)
                strout(icoef1:icoef2) = streac
                icoef1 = icoef2 + 1
                icoef2 = icoef1
                strout(icoef1:icoef2) = '+'
            END DO
            IF(mblist(istep) > 0) THEN
                streac = bdysym(mblist(istep))
                icoef1 = icoef2 + 1
                icoef2 = icoef1 + LEN(streac)
                strout(icoef1:icoef2) = streac
            ELSE
                icoef2 = icoef2 - 1
            END IF

            straro = ' => '
            IF(mglist(istep) > 0) straro = ' == '
            icoef1 = icoef2 + 1
            icoef2 = icoef1 + LEN(straro)
            strout(icoef1:icoef2) = straro
            DO ispec = 1,nsslen(istep)
                icp = 0
                DO jspec = 1,nrslen(istep)
                    IF(nrspec(jspec,istep) == nsspec(ispec,istep)) THEN
                        icp = icp + 1
                    END IF
                END DO
                icp = nint(diffmu(ispec,istep)) + icp
                IF(icp > 0) THEN
                    IF(icp > 1) THEN
                        WRITE(strcof,'(I1)')icp
                        icoef1 = icoef2 + 1
                        icoef2 = icoef1 + LEN(strcof)
                        strout(icoef1:icoef2) = strcof
                    END IF
                    streac = spcsym(nsspec(ispec,istep))
                    icoef1 = icoef2 + 1
                    icoef2 = icoef1 + LEN(streac)
                    strout(icoef1:icoef2) = streac
                    icoef1 = icoef2 + 1
                    icoef2 = icoef1
                    strout(icoef1:icoef2) = '+'
                END IF
            END DO

            IF(mblist(istep) > 0) THEN
                streac = bdysym(mblist(istep))
                icoef1 = icoef2 + 1
                icoef2 = icoef1 + LEN(streac)
                strout(icoef1:icoef2) = streac
            ELSE
                icoef2 = icoef2 - 1
            END IF

            WRITE(ncrept,'(A)')strout(1:icoef2)

        END DO

        WRITE(ncrept,*)
        WRITE(ncrept,*)'Reaction parameters A,n,E:'
        DO istep = 1,nstep
            WRITE(ncrept,'(I5,3(2X,1PE12.4))') istep,(rparam(icp,istep),icp=1,3)
        END DO

        IF(nlind > 0) THEN
            WRITE(ncrept,*)
            WRITE(ncrept,*)'Lindemann parameters:'
            DO istep = 1,nstep
                IF(mllist(istep) /= 0) THEN
                    WRITE(ncrept,'(I5,3(2X,1PE12.4))')  &
                                    istep,(rclind(icp,mllist(istep)),icp=1,3)
                END IF
            END DO
        END IF

        WRITE(ncrept,*)
        WRITE(ncrept,*)'Third body efficiencies:'
        DO jspec = 1,nbody
            WRITE(ncrept,'(I5,3X,A10)')jspec,bdysym(jspec)
            DO ispec = 1,nspec
                WRITE(ncrept,'(I5,3X,A10,1PE12.4)')  &
                                ispec,spcsym(ispec),effy3b(ispec,jspec)
            END DO
        END DO

        WRITE(ncrept,*)
        WRITE(ncrept,*)'End of chemical data'
        WRITE(ncrept,*)'--------------------'
        WRITE(ncrept,*)
        WRITE(ncrept,9000)cldash

!       =======================================================================

    END IF

!   =========================================================================

!   GET THE TRANSPORT DATA
!   ======================
!   RSC 06-JAN-2013 MIXTURE AVERAGED TRANSPORT
    call diffin

    IF(flmavt) THEN

!       =======================================================================

        IF(iproc == 0) THEN

!           REPORT THE TRANSPORT DATA
!           =========================

            WRITE(ncrept,*)
            WRITE(ncrept,*)'Transport data: mixture averaged transport'
            WRITE(ncrept,*)'--------------'
            WRITE(ncrept,*)

            WRITE(ncrept,'(2(1PE12.4))')pdifgb,tdifgb
            WRITE(ncrept,*)
            WRITE(ncrept,*)'Viscosity'
            DO ispec = 1,nspec
                WRITE(ncrept,'(2I5)')ispec,ncovis
                WRITE(ncrept,'(5(1PE15.7))') (viscco(icoeff,ispec),icoeff=1,ncovis)
            END DO
            WRITE(ncrept,*)
            WRITE(ncrept,*)'Thermal conductivity'
            DO ispec = 1,nspec
                WRITE(ncrept,'(2I5)')ispec,ncocon
                WRITE(ncrept,'(5(1PE15.7))') (condco(icoeff,ispec),icoeff=1,ncocon)
            END DO
            WRITE(ncrept,*)
            WRITE(ncrept,*)'Binary diffusion coefficient'
            DO ispec = 1,nspec
                WRITE(ncrept,'(2I5)')ispec,ncodif
                DO jspec = 1,ispec
                    WRITE(ncrept,'(I5,5(1PE15.7))')jspec,  &
                                (diffco(icoeff,jspec,ispec),icoeff=1,ncodif)
                END DO
            END DO
            WRITE(ncrept,*)
            WRITE(ncrept,*)'Thermal diffusion ratio'
            DO ispec = 1,nspec
                WRITE(ncrept,'(2I5)')ispec,ncotdr
                DO jspec = 1,ispec-1
                    WRITE(ncrept,'(I5,5(1PE15.7))')jspec,  &
                                (tdrcco(icoeff,jspec,ispec),icoeff=1,ncotdr)
                END DO
            END DO

!           EVALUATE TRANSPORT COEFFICIENTS FOR THE DEFAULT INITIAL STATE
            tmplog = LOG(trin/tdifgb)

!           INITIAL DENSITY
            rtin = zero
            DO ispec = 1,nspec
                rtin = rtin + rgspec(ispec)*yrin(ispec)
            END DO
            rtin = rtin*trin
            drin = prin/rtin

!           VISCOSITY FOR EACH SPECIES
            DO ispec = 1, nspec
                fornow = viscco(ncovis,ispec)
                DO icp = ncovm1,1,-1
                    fornow = fornow*tmplog + viscco(icp,ispec)
                END DO
                ctrans(ispec,1) = EXP(fornow)
            END DO

!           COMBINATION RULE FOR VISCOSITY
            combo1 = zero
            DO ispec = 1, nspec
                combo2 = zero
                DO jspec = 1, nspec
                    fornow = SQRT(ctrans(ispec,1)/ctrans(jspec,1))
                    fornow = one + fornow*wilko2(jspec,ispec)
                    fornow = wilko1(jspec,ispec)*fornow*fornow
                    combo2 = combo2 + yrin(jspec)*ovwmol(jspec)*fornow
                END DO
                fornow = ctrans(ispec,1)/combo2
                combo1 = combo1 + yrin(ispec)*ovwmol(ispec)*fornow

            END DO
            ttemp(1) = combo1

!           THERMAL CONDUCTIVITY FOR EACH SPECIES
            DO ispec = 1, nspec
                fornow = condco(ncocon,ispec)
                DO icp = ncocm1,1,-1
                    fornow = fornow*tmplog + condco(icp,ispec)
                END DO
                ctrans(ispec,2) = EXP(fornow)
            END DO

!           COMBINATION RULE FOR CONDUCTIVITY
            combo1 = zero
            combo2 = zero
            combo3 = zero
            DO ispec = 1, nspec
                fornow = yrin(ispec)*ovwmol(ispec)
                combo1 = combo1 + fornow*ctrans(ispec,2)
                combo2 = combo2 + fornow/ctrans(ispec,2)
                combo3 = combo3 + fornow
            END DO
            combo3 = one/combo3
            combo1 = combo1*combo3
            combo2 = combo2*combo3
            ttemp(2) = half*(combo1 + one/combo2)

!           MASS-SPECIFIC CP
            combo2 = zero
            ic = 1
            DO ispec = 1, nspec
                fornow = amascp(ncpoly(ic,ispec),ic,ispec)
                DO jc = ncpom1(ic,ispec),1,-1
                    fornow = fornow*trin + amascp(jc,ic,ispec)
                END DO
                combo1 = fornow
                combo2 = combo2 + yrin(ispec)*combo1
            END DO
            ttemp(3) = combo2

            WRITE(ncrept,*)
            WRITE(ncrept,*)'         Viscosity         Conductivity'
            DO ispec = 1, nspec
                WRITE(ncrept,'(I5,2(1PE18.7))')ispec,ctrans(ispec,1), ctrans(ispec,2)
            END DO
            WRITE(ncrept,*)
            WRITE(ncrept,*) '         Mix viscosity     Conductivity      Prandtl no'
            WRITE(ncrept,'(5X,3(1PE18.7))')ttemp(1),ttemp(2), ttemp(1)*combo2/ttemp(2)

!           MASS DIFFUSIVITY FOR EACH SPECIES
            DO ispec = 1, nspec

                DO jspec = 1, nspec
                    fornow = diffco(ncocon,jspec,ispec)
                    DO icp = ncodm1,1,-1
                        fornow = fornow*tmplog + diffco(icp,jspec,ispec)
                    END DO
                    dtrans(jspec) = EXP(fornow)*pdifgb/prin
                END DO

!               COMBINATION RULE FOR MASS DIFFUSIVITY
                combo1 = zero
                combo2 = zero
                combo3 = zero
                DO jspec = 1, nspec
                    fornow = yrin(jspec) + dfctol
                    combo1 = combo1 + fornow
                    combo2 = combo2 + fornow*ovwmol(jspec)/dtrans(jspec)
                    combo3 = combo3 + yrin(jspec)*ovwmol(jspec)
                END DO
                fornow = yrin(ispec) + dfctol
                combo1 = combo1 - fornow
                combo2 = combo2 - fornow*ovwmol(ispec)/dtrans(ispec)
                combo2 = combo2/combo3
                ctrans(ispec,3) = drin*combo1/combo2

            END DO

!           THERMAL DIFFUSION RATIO FOR EACH SPECIES
            tmplog = trin/tdifgb
            DO ispec = 1, nspec

                DO jspec = 1, nspec
                    fornow = tdrcco(ncocon,jspec,ispec)
                    DO icp = ncotm1,1,-1
                        fornow = fornow*tmplog + tdrcco(icp,jspec,ispec)
                    END DO
                    dtrans(jspec) = fornow
                END DO

!               COMBINATION RULE FOR THERMAL DIFFUSION RATIO
                combo1 = zero
                combo2 = zero
                DO jspec = 1, nspec
                    fornow = yrin(jspec)*ovwmol(jspec)
                    combo1 = combo1 + fornow*dtrans(jspec)
                    combo2 = combo2 + fornow
                END DO
                ctrans(ispec,4) = combo1/combo2

            END DO

            WRITE(ncrept,*)
            WRITE(ncrept,*)'         Diffusivity       Schmidt No',  &
                           '        Lewis No          Thermal diff'
            DO ispec = 1, nspec
                WRITE(ncrept,'(I5,4(1PE18.7))')ispec,ctrans(ispec,3),  &
                                ttemp(1)/ctrans(ispec,3), ttemp(2)/(ctrans(ispec,3)*ttemp(3)),  &
                                ctrans(ispec,4)
            END DO

            WRITE(ncrept,*)
            WRITE(ncrept,*)'End of transport data'
            WRITE(ncrept,*)'---------------------'
            WRITE(ncrept,*)
            WRITE(ncrept,9000)cldash

        END IF

!       =======================================================================

    ELSE

!       =======================================================================

        IF(iproc == 0) THEN

            WRITE(ncrept,*)
            WRITE(ncrept,*)'Transport data: constant Lewis numbers'
            WRITE(ncrept,*)'--------------'
            WRITE(ncrept,*)

            WRITE(ncrept,*)' constants A,r; ref T; Prandtl no.:'
            WRITE(ncrept,'(4(1PE12.4))')alamdc,rlamda,tlamda,prantl
            WRITE(ncrept,*)

            WRITE(ncrept,*)'Spec.No.  Symbol         Lewis no.'
            DO ispec = 1,nspec
                WRITE(ncrept,'(I8,3X,A10,3X,1PE12.4)') ispec,spcsym(ispec),clewis(ispec)
            END DO

            WRITE(ncrept,*)
            WRITE(ncrept,*)'End of transport data'
            WRITE(ncrept,*)'---------------------'
            WRITE(ncrept,*)
            WRITE(ncrept,9000)cldash

        END IF

!       =======================================================================

    END IF
!   MIXTURE AVERAGED TRANSPORT

!   =========================================================================

!   RADIATION TREATMENT
    call radcin

    IF(flradn) THEN

!       =======================================================================

        IF(iproc == 0) THEN

!           REPORT THE RADIATION DATA
!           =========================
            WRITE(ncrept,*)
            WRITE(ncrept,*)'Radiation data'
            WRITE(ncrept,*)'--------------'
            WRITE(ncrept,*)
            WRITE(ncrept,*)'Planck mean absorption coefficient data'
            DO ispec = 1, nsprad
                WRITE(ncrept,'(2I6)')nsprid(ispec),nkprad(ispec)
                WRITE(ncrept,'(5(1PE15.7))')(akprad(icp,ispec),icp=1,5)
                WRITE(ncrept,'(5(1PE15.7))')(akprad(icp,ispec), icp=6,nkprad(ispec))
            END DO

            WRITE(ncrept,*)
            WRITE(ncrept,*)'Reference ambient temperature:'
            WRITE(ncrept,'(1PE15.7)')trefrn

            WRITE(ncrept,*)
            WRITE(ncrept,*)'Planck mean absorption coefficients:'
            DO ispec = 1, nsprad
                fornow = akprad(nkprad(ispec),jspec)
                DO icp = nkprm1(ispec),1,-1
                    fornow = fornow*trefrn + akprad(icp,ispec)
                END DO
                jspec = nsprid(ispec)
                WRITE(ncrept,'(I6,5X,A10,1PE15.7)')jspec,spcsym(jspec),fornow
            END DO

            WRITE(ncrept,*)
            WRITE(ncrept,*)'End of radiation data'
            WRITE(ncrept,*)'---------------------'
            WRITE(ncrept,*)
            WRITE(ncrept,9000)cldash

        END IF

!       =======================================================================

    END IF

!   =========================================================================

!   END OF REPORTING OF INITIAL DATA FOR NOW
!   NOTE THAT THE REPORT FILE IS LEFT OPEN
!   FOR FURTHER REPORTING FROM OTHER INITIALISATION ROUTINES
!   AND IS CLOSED BEFORE EXIT FROM THIS SUBROUTINE

!   =========================================================================

!   INITIALISE ONLINE STATISTICS
!   ============================

!   STATISTICS ON ONE PROCESSOR ONLY
    IF(iproc == 0) THEN

!       INITIALISE THE STATISTICS FILE
        OPEN(UNIT=ncstat,FILE=fnstat,STATUS='REPLACE',FORM='FORMATTED')
        ENDFILE(ncstat)
        CLOSE(ncstat)

    END IF

!   ZERO THE STATISTICS STORAGE COUNTER
    itstat = 0

!   =========================================================================

!   CHECK AND INITIALISE DUMP FILES
!   -------------------------------
    INQUIRE(FILE=fndmpo(1),EXIST=fxdump)
    IF(.NOT.fxdump) THEN
        IF(ndofmt == 0) THEN
            OPEN(UNIT=ncdmpo,FILE=fndmpo(1),STATUS='REPLACE',FORM='UNFORMATTED')
        ELSE
            OPEN(UNIT=ncdmpo,FILE=fndmpo(1),STATUS='REPLACE',FORM='FORMATTED')
        END IF
        CLOSE(ncdmpo)
    END IF

    INQUIRE(FILE=fndmpo(2),EXIST=fxdump)
    IF(.NOT.fxdump) THEN
        IF(ndofmt == 0) THEN
            OPEN(UNIT=ncdmpo,FILE=fndmpo(2),STATUS='REPLACE',FORM='UNFORMATTED')
        ELSE
            OPEN(UNIT=ncdmpo,FILE=fndmpo(2),STATUS='REPLACE',FORM='FORMATTED')
        END IF
        CLOSE(ncdmpo)
    END IF

!   ==========================================================================

!   INITIALISE ADDITIONAL PARAMETERS
!   ================================

!   SET VALUE OF PI
!   ---------------
    pi = FOUR*ATAN(ONE)

!   SET VALUE OF LN(10)
!   -------------------
    clnten = LOG(10.0_8)

!   TIME STEPPING
!   -------------
    etime = zero
    ntime2 = ntime1+ntime-1
    itime = 0

!   =========================================================================

!   INITIALISE BOUNDARIES
!   =====================
    call bcinit

!   =========================================================================

!   INITIALISE THE VELOCITY FIELD
!   =============================
!   NOTE: ALL VARIABLES ARE INITIALISED ON STANDARD SIZE DOMAIN ONLY

!   PRE-INITIALISE THE VELOCITY FIELD TO ZERO
!   ---------------------------------
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!   -------------------------------------------------------------------------

!   INITIAL TURBULENCE GENERATOR
!   ----------------------------

!   GENERATE FRESH INITIAL TURBULENCE
    IF(inturb == 1) call turbin

!   COPY TURBULENT INFLOW VELOCITY FIELD INTO THE DOMAIN
    IF(inturb == 2) THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ))

        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ))

        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))

    END IF

!   -------------------------------------------------------------------------

!   ADD ON THE DEFAULT MEAN VELOCITY
!   --------------------------------
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqF, "A = A + var", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_gbl(urin, 1, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqF, "A = A + var", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_gbl(vrin, 1, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqF, "A = A + var", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_gbl(wrin, 1, "real(kind=8)", OPS_READ))

!   =========================================================================

!   INITIALISE THE SCALAR FIELD
!   ===========================
!   NOTE: ALL VARIABLES ARE INITIALISED ON STANDARD SIZE DOMAIN ONLY

!   PRE-INITIALISE SPECIES MASS FRACTIONS TO DEFAULT VALUES
!   -------------------------------------
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    DO ispec = 1,nspec
        call ops_par_loop(maths_kernel_eqE, "A = var(indx)", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_gbl(yrin, nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    END DO

!   PRE-INITIALISE PRESSURE AND TEMPERATURE TO DEFAULT VALUES
!   ---------------------------------------
!   BIGGER SIZE ARRAYS
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    call ops_par_loop(maths_kernel_eqD, "A = var", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(prin, 1, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqD, "A = var", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(trin, 1, "real(kind=8)", OPS_READ))

!   INITIAL DENSITY
!   ---------------
    rtin = zero
    DO ispec = 1,nspec
        rtin = rtin + rgspec(ispec)*yrin(ispec)
    END DO
    rtin = rtin*trin
    drin = prin/rtin

!   INITIAL INTERNAL ENERGY
!   -----------------------
    erin = zero
    DO ispec = 1,nspec

!       TEMPERATURE INTERVAL INDEX
        itint = 1
        DO WHILE (trin > tinthi(itint,ispec) .and. itint < ntint(ispec))
            itint = itint + 1
        END DO

        fornow = amasch(ncpoly(itint,ispec),itint,ispec)
        DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*trin + amasch(icp,itint,ispec)
        END DO
        fornow = amasch(ncenth(itint,ispec),itint,ispec) + fornow*trin
        erin = erin + fornow*yrin(ispec)

    END DO
    fornow = urin*urin + vrin*vrin + wrin*wrin
    erin = erin - rtin + half*fornow

!   -------------------------------------------------------------------------

!   INITIAL FLAME GENERATOR
!   -----------------------
    IF(inflam == 1) THEN

!       GENERATE AN INITIAL THERMOCHEMICAL FIELD
!       TEMPERATURE, PRESSURE AND MASS FRACTIONS
!       MODIFY VELOCITIES AS REQUIRED
        call flamin

    END IF

!   -------------------------------------------------------------------------

!   DENSITY FIELD
!   -------------
!   USE STORE1 AS AN ACCUMULATOR FOR MIXTURE GAS CONSTANT
!   INITIALISE TO ZERO
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_WRITE))


!   MIXTURE GAS CONSTANT
    DO ispec = 1,nspec
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(maths_kernel_eqK, "A = A + var(indx)*B", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    END DO

!   DENSITY
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqV, "A = A*B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqU, "A = B/C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                           STORE1 = T*MIX RG
!   -------------------------------------------------------------------------

!   INTERNAL ENERGY FIELD
!   ---------------------
!   INITIALISE TEMPERATURE INTERVAL INDEX
!   FOR ALL SPECIES LOCATE TEMPERATURE IN AN INTERVAL
!   BIGGER SIZE ARRAY
    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    DO iindex = 1,nintmx
        call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
    END DO

    DO ispec = 1, nspec
!       SET THE TEMPERATURE INDEX
        iindex = 1 + (ispec-1)/nspimx
        ipower = ispec - (iindex-1)*nspimx - 1

        call ops_par_loop(maths_kernel_eqBO, "INTERNAL ENERGY FIELD", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC), &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                    ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                    ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
    END DO

!   PRE-INITIALISE INTERNAL ENERGY TO ZERO
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!   INITIALISE INTERNAL ENERGY
    DO ispec = 1,nspec

!       TEMPERATURE INTERVAL INDEXING
        iindex = 1 + (ispec-1)/nspimx
        ipower = ispec - (iindex-1)*nspimx - 1
        icoef2 = ntbase**ipower
        icoef1 = icoef2*ntbase

!       USE THERMOCHEMICAL DATA
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(maths_kernel_eqBP, "INITIALISE INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ), &
                        ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    END DO

!   FINALISE INTERNAL ENERGY
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(maths_kernel_eqAU, "A = A - B + half*(C*C+D*D+E*E)", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_READ))

!                                                            ALL STORES CLEAR
!   =========================================================================

!   CONVERT VARIABLES TO CONSERVATIVE FORM
!   ======================================
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]

    call ops_par_loop(maths_kernel_eqV_fused, "A = A*B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_READ))

    DO ispec = 1,nspec
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(maths_kernel_eqV, "A = B*A", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_RW), &
                        ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO

!   =========================================================================

!   EXECUTE REQUIRED STARTUP
!   ========================

!   COLD START SWITCH
!   -----------------
    IF(ncdmpi == 1) THEN

!       =======================================================================
!       WARM START
!       =======================================================================

!       -----------------------------------------------------------------------
!       THIS BLOCK MAY BE MODIFIED AS REQUIRED
!       TO BLEND INITIAL VELOCITY AND SCALAR FIELDS
!       WITH PREVIOUSLY DUMPED DATA
!       -----------------------------------------------------------------------

!       Vishnu: Prevent overwriting TSTEP
        tsold=tstep

!       RESTART FROM FULL DUMP FILES
!       ----------------------------
!       READ THE DATA FROM DUMP INPUT FILE 1
!       NOTE THAT URUN,VRUN,WRUN,ERUN AND YRUN ARE ALL IN CONSERVATIVE FORM
!       RSC 11-JUL-2009 ADD A DUMP FORMAT SWITCH
        idflag = MOD(INT((ntime1-1)/ntdump), 2) + 1
        IF(ndifmt == 0) THEN

!           UNFORMATTED DUMP INPUT
            OPEN(UNIT=ncdmpi,FILE=fndmpo(idflag),STATUS='OLD', FORM='UNFORMATTED')

            IF (ops_is_root() == 1) THEN
                WRITE(*,*) "Reading previous run information from file(unformatted): ", trim(fndmpo(idflag)), " idflag: ", idflag
            END IF

            READ(ncdmpi)nxdmax,nydmax,nzdmax,ndspec,&
                        etime,tstep,errold,errldr

!           SIZE ERROR CHECK
            IF(nxdmax /= nxglbl) WRITE(6,*)'Dump input size error: x'
            IF(nydmax /= nyglbl) WRITE(6,*)'Dump input size error: y'
            IF(nzdmax /= nzglbl) WRITE(6,*)'Dump input size error: z'
            IF(ndspec /= nspec)  WRITE(6,*)'Dump input size error: species'

            CLOSE(ncdmpi)
        ELSE

!           FORMATTED DUMP INPUT
            OPEN(UNIT=ncdmpi,FILE=fndmpo(idflag),STATUS='OLD',FORM='FORMATTED')

            IF (ops_is_root() == 1) THEN
                WRITE(*,*) "Reading previous run information from file(formatted): ", trim(fndmpo(idflag)), " idflag: ", idflag
            END IF

            READ(ncdmpi,*)nxdmax,nydmax,nzdmax,ndspec

!           SIZE ERROR CHECK
            IF(nxdmax /= nxglbl) WRITE(6,*)'Dump input size error: x'
            IF(nydmax /= nyglbl) WRITE(6,*)'Dump input size error: y'
            IF(nzdmax /= nzglbl) WRITE(6,*)'Dump input size error: z'
            IF(ndspec /= nspec)  WRITE(6,*)'Dump input size error: species'

            READ(ncdmpi,*)etime,tstep,errold,errldr

            CLOSE(ncdmpi)

        END IF

        IF (ops_is_root() == 1) THEN
            WRITE(*,*) "etime: ", etime, "  tstep: ", tstep, "  errold:", errold, "  errldr: ",errldr
            WRITE(*,*) "Warm Start: Copying dumped dats....."
        END IF

        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_drun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_urun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_vrun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_wrun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_erun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

        DO ispec = 1,nspcmx
            call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_yrun_dump(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))
        END DO

        call ops_free_dat(d_drun_dump)
        call ops_free_dat(d_urun_dump)
        call ops_free_dat(d_vrun_dump)
        call ops_free_dat(d_wrun_dump)
        call ops_free_dat(d_erun_dump)
        DO ispec = 1,nspcmx
            call ops_free_dat(d_yrun_dump(ispec))
        END DO

!        dtime=INT((ntime1-1)/ntdump)
!        WRITE(citime,'(I8.8)') dtime
!        fname = 'output/timestep_warmstart'//citime//pnxhdf
!        call ops_fetch_block_hdf5_file(senga_grid, trim(fname))

!        call ops_fetch_dat_hdf5_file(d_drun, trim(fname))
!        call ops_fetch_dat_hdf5_file(d_urun, trim(fname))
!        call ops_fetch_dat_hdf5_file(d_vrun, trim(fname))
!        call ops_fetch_dat_hdf5_file(d_wrun, trim(fname))
!        call ops_fetch_dat_hdf5_file(d_erun, trim(fname))

!        DO ispec = 1,nspcmx
!            call ops_fetch_dat_hdf5_file(d_yrun(ispec), trim(fname))
!        END DO

!       Vishnu: TSTEP not overwritten
        IF(nstpsw == 0) tstep=tsold

!       =======================================================================
!       WARM START COMPLETE
!       =======================================================================

    END IF

!   ==========================================================================

!   INITIALISE TIME-STEPPING RHS TERMS
!   ==================================
!   BIGGER SIZE ARRAYS

!   PRE-INITIALISE TO DEFAULT VALUES
    durin = drin*urin
    dvrin = drin*vrin
    dwrin = drin*wrin
    derin = drin*erin

    DO ispec = 1,nspec
        dyrin(ispec) = drin*yrin(ispec)
    END DO

    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    call ops_par_loop(maths_kernel_eqD, "A = var", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(drin, 1, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqD, "A = var", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(durin, 1, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqD, "A = var", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(dvrin, 1, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqD, "A = var", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(dwrin, 1, "real(kind=8)", OPS_READ))

    call ops_par_loop(maths_kernel_eqD, "A = var", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(derin, 1, "real(kind=8)", OPS_READ))

    rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    DO ispec = 1,nspec
        call ops_par_loop(maths_kernel_eqE, "A = var(indx)", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_gbl(dyrin, nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    END DO

!   NOTE THAT THESE ARE BIGGER SIZE ARRAYS
!   BUT THAT NO HALO DATA IS YET AVAILABLE
!   THEREFORE ONLY THE STANDARD DOMAIN SIZE IS INITIALISED
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_READ))

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    DO ispec = 1,nspec
        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO

!   =========================================================================

!   INITIAL PARALLEL DATA TRANSFER
!   ==============================
!   INITIALISES TIME-STEPPING RHS TERMS TO BIGGER DOMAIN SIZE
    call parfer

!   =========================================================================

!   INITIAL TEMPERATURE AND PRESSURE
!   ================================
    call tempin

!#ifdef OPS_LAZY
!    call ops_execute()
!#endif

!   =========================================================================

!   INITIALISE TIME STEPPING
!   ========================
    call dtinit

!   ==========================================================================
!
!   INITIALISE SPATIAL DIFFERENTIATORS
!   ==================================
    call dfinit

!   ==========================================================================

!   INITIALISE SPATIAL FILTERING
!   ============================
!    call flinit

!   ==========================================================================

!   CLOSE THE REPORT FILE
!   =====================
    IF(iproc == 0) THEN

        WRITE(ncrept,9000)cldash
        WRITE(ncrept,*)

        CLOSE(ncrept)

    END IF

!   ==========================================================================

!     FORMATS FOR OUTPUT OF INITIAL DATA
9000  FORMAT(a79)
9010  FORMAT('**',75X,'**')
9020  FORMAT('**',18X,'DIRECT NUMERICAL SIMULATION CODE SENGA2', 18X,'**')
9030  FORMAT('**',31X,'Stewart Cant',32X,'**')
9040  FORMAT('**',16X,'Cambridge University Engineering Department', 16X,'**')
9200  FORMAT(12I6)
9250  FORMAT(6I7)
9270  FORMAT(2I5,5(1PE12.4))
9300  FORMAT(1PE12.4,6I8)
9400  FORMAT(6(1PE12.4))
9450  FORMAT(i7,1PE15.7)
9600  FORMAT(' Actual (NSPEC) = ',i5,2X,'Max (NSPCMX) = ',i5)

END SUBROUTINE indata
