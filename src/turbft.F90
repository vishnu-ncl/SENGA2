SUBROUTINE turbft

!   *************************************************************************

!   TURBFT
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   28-MAR-1997:  CREATED
!   15-MAY-1999:  KWJ PARALLEL IMPLEMENTATION (INDEPENDENT)
!   20-OCT-1999:  KWJ PARALLEL FFT IMPLEMENTATION
!   24-MAY-2003:  RSC UPDATED FOR SENGA2

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   INITIAL TURBULENCE GENERATOR
!   INVERSE FFT WITH CONJUGATE ANTISYMMETRY

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
    use com_senga
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    real(kind=8) :: scalex,scaley,scalez
    integer(kind=4) :: ipen(npenmx),jpen(npenmx),kpen(npenmx)
    integer(kind=4) :: ic,jc,kc,ix,jx,kx,ipencl,icproc
    integer(kind=4) :: iim,iic,jjm,jjc,kkm,kkc
    integer(kind=4) :: igofst,jgofst,kgofst,igofm1,jgofm1,kgofm1
    integer(kind=4) :: nodblx,nodbly,nodblz
    integer(kind=4) :: ndosym,npencl,ncount

!   BEGIN
!   =====

!   =========================================================================

!   INDEXING
!   --------

!   SET ODDBALL WAVENUMBER INDICES
    nodblx = nxglbl/2
    nodbly = nyglbl/2
    nodblz = nzglbl/2

!   PHYSICAL-SPACE GLOBAL INDEX OFFSETS
    igofst = 0
    DO icproc = 0, ixproc-1
        igofst = igofst + npmapx(icproc)
    END DO
    igofm1 = igofst-1

    jgofst = 0
    DO icproc = 0, iyproc-1
        jgofst = jgofst + npmapy(icproc)
    END DO
    jgofm1 = jgofst-1

    kgofst = 0
    DO icproc = 0, izproc-1
        kgofst = kgofst + npmapz(icproc)
    END DO
    kgofm1 = kgofst-1

!   =========================================================================

!   SET SCALE FACTORS
    scalex = one
    scaley = one
    scalez = one

!   =========================================================================

!   CARRY OUT AN INVERSE FOURIER TRANSFORM
!   ======================================

!   CHECK FFT PARALLEL DATA SIZE AGAINST PARALLEL ARRAY SIZE
    IF(iproc == 0) THEN
        ncount = 2*nszmax*3*npenmx
        IF(ncount > nparay) THEN
            WRITE(ncrept,*)'Warning: TURBFT: parallel array size too small'
        END IF
    END IF

!   =========================================================================

    write(*, '(a)') "Using the arrays not allocated by OPS, &
                        Please implement the function in OPS first, turbft.F90: ID=105"
    STOP

!   X-DIRECTION
!   -----------

!   PENCIL COUNTER
    ipencl = 0

!   LOOP OVER X-LINES
    DO kc = 1,nzsize
!       FOURIER-SPACE GLOBAL INDEXING
        kx = kgofm1 + kc
        IF(kx > nodblz) kx = kx-nzglbl
        DO jc = 1,nysize
!           FOURIER-SPACE GLOBAL INDEXING
            jx = jgofm1 + jc
            IF(jx > nodbly) jx = jx-nyglbl

!           TRANSFORM ONLY FOR UPPER-CENTRAL WAVENUMBERS
            IF(kx > 0) THEN
!               STANDARD UPPER-CENTRAL X-LINES
!               PENCIL INDEXING
                ipencl = ipencl + 1
                jpen(ipencl) = jc
                kpen(ipencl) = kc

!               ASSEMBLE FFT DATA
                DO ic = 1,nxsize
                    iic = 2*ic
                    iim = iic-1
                    ftpart(iim,1,ipencl) = urun(ic,jc,kc)
                    ftpart(iic,1,ipencl) = utmp(ic,jc,kc)
                    ftpart(iim,2,ipencl) = vrun(ic,jc,kc)
                    ftpart(iic,2,ipencl) = vtmp(ic,jc,kc)
                    ftpart(iim,3,ipencl) = wrun(ic,jc,kc)
                    ftpart(iic,3,ipencl) = wtmp(ic,jc,kc)

                END DO

!               CHECK FOR MAX NO OF PENCILS
                IF(ipencl == npenmx) THEN
!                   DO THE FFT (NO SYMMETRY IMPOSED)
                    ndosym = 0
                    npencl = ipencl

                    CALL fftsym(ndosym,npencl,ixproc,nxproc,nxsize,nxglbl,nprocx)

!                   RESTORE TRANSFORMED DATA
                    DO ipencl = 1,npencl
                        DO ic = 1,nxsize
                            iic = 2*ic
                            iim = iic-1
                            urun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,1,ipencl)*scalex
                            utmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,1,ipencl)*scalex
                            vrun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,2,ipencl)*scalex
                            vtmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,2,ipencl)*scalex
                            wrun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,3,ipencl)*scalex
                            wtmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,3,ipencl)*scalex

                        END DO
                    END DO

                    ipencl = 0

                END IF

            ELSE IF(kx == 0) THEN
                IF(jx > 0) THEN
!                   KX = 0 PLANE UPPER-CENTRAL X-LINES
!                   PENCIL INDEXING
                    ipencl = ipencl + 1
                    jpen(ipencl) = jc
                    kpen(ipencl) = kc

!                   ASSEMBLE FFT DATA
                    DO ic = 1,nxsize
                        iic = 2*ic
                        iim = iic-1
                        ftpart(iim,1,ipencl) = urun(ic,jc,kc)
                        ftpart(iic,1,ipencl) = utmp(ic,jc,kc)
                        ftpart(iim,2,ipencl) = vrun(ic,jc,kc)
                        ftpart(iic,2,ipencl) = vtmp(ic,jc,kc)
                        ftpart(iim,3,ipencl) = wrun(ic,jc,kc)
                        ftpart(iic,3,ipencl) = wtmp(ic,jc,kc)
                    END DO

!                   CHECK FOR MAX NO OF PENCILS
                    IF(ipencl == npenmx) THEN
!                       DO THE FFT (NO SYMMETRY IMPOSED)
                        ndosym = 0
                        npencl = ipencl

                        CALL fftsym(ndosym,npencl,ixproc,nxproc,nxsize,nxglbl,nprocx)

!                       RESTORE TRANSFORMED DATA
                        DO ipencl = 1,npencl
                            DO ic = 1,nxsize
                                iic = 2*ic
                                iim = iic-1
                                urun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,1,ipencl)*scalex
                                utmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,1,ipencl)*scalex
                                vrun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,2,ipencl)*scalex
                                vtmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,2,ipencl)*scalex
                                wrun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,3,ipencl)*scalex
                                wtmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,3,ipencl)*scalex
                            END DO
                        END DO

                        ipencl = 0

                    END IF

                ELSE IF(jx == 0) THEN
!                   SPECIAL CASE OF KY=KZ=0
!                   CHECK FOR PRE-STORED PENCILS
                    IF(ipencl /= 0) THEN
!                       DO THE FFT (NO SYMMETRY IMPOSED)
                        ndosym = 0
                        npencl = ipencl

                        CALL fftsym(ndosym,npencl,ixproc,nxproc,nxsize,nxglbl,nprocx)

!                       RESTORE TRANSFORMED DATA
                        DO ipencl = 1,npencl
                            DO ic = 1,nxsize
                                iic = 2*ic
                                iim = iic-1
                                urun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,1,ipencl)*scalex
                                utmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,1,ipencl)*scalex
                                vrun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,2,ipencl)*scalex
                                vtmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,2,ipencl)*scalex
                                wrun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,3,ipencl)*scalex
                                wtmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,3,ipencl)*scalex
                            END DO
                        END DO

                    END IF

!                   ASSEMBLE FFT DATA FOR KY=KZ=0
                    ipencl = 1
                    DO ic = 1,nxsize
                        iic = 2*ic
                        iim = iic-1
                        ftpart(iim,1,ipencl) = urun(ic,jc,kc)
                        ftpart(iic,1,ipencl) = utmp(ic,jc,kc)
                        ftpart(iim,2,ipencl) = vrun(ic,jc,kc)
                        ftpart(iic,2,ipencl) = vtmp(ic,jc,kc)
                        ftpart(iim,3,ipencl) = wrun(ic,jc,kc)
                        ftpart(iic,3,ipencl) = wtmp(ic,jc,kc)
                    END DO

!                   DO THE FFT (SYMMETRY IMPOSED)
                    ndosym = 1
                    npencl = 1

                    CALL fftsym(ndosym,npencl,ixproc,nxproc,nxsize,nxglbl,nprocx)

!                   RESTORE TRANSFORMED DATA
                    DO ic = 1,nxsize
                        iic = 2*ic
                        iim = iic-1
                        urun(ic,jc,kc) = ftpart(iim,1,ipencl)*scalex
                        utmp(ic,jc,kc) = ftpart(iic,1,ipencl)*scalex
                        vrun(ic,jc,kc) = ftpart(iim,2,ipencl)*scalex
                        vtmp(ic,jc,kc) = ftpart(iic,2,ipencl)*scalex
                        wrun(ic,jc,kc) = ftpart(iim,3,ipencl)*scalex
                        wtmp(ic,jc,kc) = ftpart(iic,3,ipencl)*scalex
                    END DO

                    ipencl = 0

                ELSE
!                   KX = 0 JX < 0 LOWER-CENTRAL WAVENUMBERS - NO ACTION
                    CONTINUE

                END IF

            ELSE
!               KX < 0 LOWER-CENTRAL WAVENUMBERS - NO ACTION
                CONTINUE

            END IF

        END DO

    END DO
!   END OF LOOP OVER X-LINES

!   CHECK FOR ANY REMAINING STORED PENCILS
    IF(ipencl /= 0) THEN
!       DO THE FFT (NO SYMMETRY IMPOSED)
        ndosym = 0
        npencl = ipencl

        CALL fftsym(ndosym,npencl,ixproc,nxproc,nxsize,nxglbl,nprocx)

!       RESTORE TRANSFORMED DATA
        DO ipencl = 1,npencl
            DO ic = 1,nxsize
                iic = 2*ic
                iim = iic-1
                urun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,1,ipencl)*scalex
                utmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,1,ipencl)*scalex
                vrun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,2,ipencl)*scalex
                vtmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,2,ipencl)*scalex
                wrun(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iim,3,ipencl)*scalex
                wtmp(ic,jpen(ipencl),kpen(ipencl)) = ftpart(iic,3,ipencl)*scalex
            END DO
        END DO

    END IF

!   =========================================================================

!   Y-DIRECTION
!   -----------

!   PENCIL COUNTER
    ipencl = 0

    DO kc = 1,nzsize
!       FOURIER-SPACE GLOBAL INDEXING
        kx = kgofm1 + kc
        IF(kx > nodblz) kx = kx-nzglbl

!       TRANSFORM ONLY FOR UPPER-CENTRAL WAVENUMBERS
        IF(kx > 0) THEN

            DO ic = 1,nxsize
!               FOURIER-SPACE GLOBAL INDEXING
                ix = igofm1 + ic
                IF(ix > nodblx) ix = ix-nxglbl

!               PENCIL INDEXING
                ipencl = ipencl + 1
                ipen(ipencl) = ic
                kpen(ipencl) = kc

!               ASSEMBLE FFT DATA
                DO jc = 1,nysize
                    jjc = 2*jc
                    jjm = jjc-1
                    ftpart(jjm,1,ipencl) = urun(ic,jc,kc)
                    ftpart(jjc,1,ipencl) = utmp(ic,jc,kc)
                    ftpart(jjm,2,ipencl) = vrun(ic,jc,kc)
                    ftpart(jjc,2,ipencl) = vtmp(ic,jc,kc)
                    ftpart(jjm,3,ipencl) = wrun(ic,jc,kc)
                    ftpart(jjc,3,ipencl) = wtmp(ic,jc,kc)
                END DO

!               CHECK FOR MAX NO OF PENCILS
                IF(ipencl == npenmx) THEN
!                   DO THE FFT (NO SYMMETRY IMPOSED)
                    ndosym = 0
                    npencl = ipencl

                    CALL fftsym(ndosym,npencl,iyproc,nyproc,nysize,nyglbl,nprocy)

!                   RESTORE TRANSFORMED DATA
                    DO ipencl = 1,npencl
                        DO jc = 1,nysize
                            jjc = 2*jc
                            jjm = jjc-1
                            urun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,1,ipencl)*scaley
                            utmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,1,ipencl)*scaley
                            vrun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,2,ipencl)*scaley
                            vtmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,2,ipencl)*scaley
                            wrun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,3,ipencl)*scaley
                            wtmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,3,ipencl)*scaley
                        END DO
                    END DO

                    ipencl = 0

                END IF

            END DO

        ELSE IF(kx == 0) THEN
!           SPECIAL CASE OF KZ=0
!           CHECK FOR PRE-STORED PENCILS
            IF(ipencl /= 0) THEN
!               DO THE FFT (NO SYMMETRY IMPOSED)
                ndosym = 0
                npencl = ipencl

                CALL fftsym(ndosym,npencl,iyproc,nyproc,nysize,nyglbl,nprocy)

!               RESTORE TRANSFORMED DATA
                DO ipencl = 1,npencl
                    DO jc = 1,nysize
                        jjc = 2*jc
                        jjm = jjc-1
                        urun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,1,ipencl)*scaley
                        utmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,1,ipencl)*scaley
                        vrun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,2,ipencl)*scaley
                        vtmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,2,ipencl)*scaley
                        wrun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,3,ipencl)*scaley
                        wtmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,3,ipencl)*scaley
                    END DO
                END DO

            END IF

!           ASSEMBLE FFT DATA FOR KZ=0
            ipencl = 0

            DO ic = 1,nxsize
!               FOURIER-SPACE GLOBAL INDEXING
                ix = igofm1 + ic
                IF(ix > nodblx) ix = ix-nxglbl

!               PENCIL INDEXING
                ipencl = ipencl + 1
                ipen(ipencl) = ic
                kpen(ipencl) = kc

!               ASSEMBLE FFT DATA
                DO jc = 1,nysize
                    jjc = 2*jc
                    jjm = jjc-1
                    ftpart(jjm,1,ipencl) = urun(ic,jc,kc)
                    ftpart(jjc,1,ipencl) = utmp(ic,jc,kc)
                    ftpart(jjm,2,ipencl) = vrun(ic,jc,kc)
                    ftpart(jjc,2,ipencl) = vtmp(ic,jc,kc)
                    ftpart(jjm,3,ipencl) = wrun(ic,jc,kc)
                    ftpart(jjc,3,ipencl) = wtmp(ic,jc,kc)
                END DO

!               CHECK FOR MAX NO OF PENCILS
                IF((ipencl == npenmx) .or. (ic == nxsize)) THEN
!                   DO THE FFT (SYMMETRY IMPOSED)
                    ndosym = 1
                    npencl = ipencl

                    CALL fftsym(ndosym,npencl,iyproc,nyproc,nysize,nyglbl,nprocy)

!                   RESTORE TRANSFORMED DATA
                    DO ipencl = 1,npencl
                        DO jc = 1,nysize
                            jjc = 2*jc
                            jjm = jjc-1
                            urun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,1,ipencl)*scaley
                            utmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,1,ipencl)*scaley
                            vrun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,2,ipencl)*scaley
                            vtmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,2,ipencl)*scaley
                            wrun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,3,ipencl)*scaley
                            wtmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,3,ipencl)*scaley
                        END DO
                    END DO

                    ipencl = 0

                END IF

            END DO

        ELSE
!           LOWER-CENTRAL WAVENUMBERS - NO ACTION
            CONTINUE

        END IF

    END DO

!   CHECK FOR ANY REMAINING STORED PENCILS
    IF(ipencl /= 0)THEN
!       DO THE FFT (NO SYMMETRY IMPOSED)
        ndosym = 0
        npencl = ipencl

        CALL fftsym(ndosym,npencl,iyproc,nyproc,nysize,nyglbl,nprocy)

!       RESTORE TRANSFORMED DATA
        DO ipencl = 1,npencl
            DO jc = 1,nysize
                jjc = 2*jc
                jjm = jjc-1
                urun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,1,ipencl)*scaley
                utmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,1,ipencl)*scaley
                vrun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,2,ipencl)*scaley
                vtmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,2,ipencl)*scaley
                wrun(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjm,3,ipencl)*scaley
                wtmp(ipen(ipencl),jc,kpen(ipencl)) = ftpart(jjc,3,ipencl)*scaley
            END DO
        END DO

    END IF

!     =========================================================================

!     Z-DIRECTION
!     -----------

!   PENCIL COUNTER
    ipencl = 0

    DO jc = 1,nysize
!       FOURIER-SPACE GLOBAL INDEXING
        jx = jgofm1 + jc
        IF(jx > nodbly) jx = jx-nyglbl

        DO ic = 1,nxsize

!           FOURIER-SPACE GLOBAL INDEXING
            ix = igofm1 + ic
            IF(ix > nodblx) ix = ix-nxglbl

!           PENCIL INDEXING
            ipencl = ipencl + 1
            ipen(ipencl) = ic
            jpen(ipencl) = jc

!           ASSEMBLE FFT DATA
            DO kc = 1,nzsize
                kkc = 2*kc
                kkm = kkc-1
                ftpart(kkm,1,ipencl) = urun(ic,jc,kc)
                ftpart(kkc,1,ipencl) = utmp(ic,jc,kc)
                ftpart(kkm,2,ipencl) = vrun(ic,jc,kc)
                ftpart(kkc,2,ipencl) = vtmp(ic,jc,kc)
                ftpart(kkm,3,ipencl) = wrun(ic,jc,kc)
                ftpart(kkc,3,ipencl) = wtmp(ic,jc,kc)
            END DO

!           CHECK FOR MAX NO OF PENCILS
            IF(ipencl == npenmx) THEN
!               DO THE FFT (SYMMETRY IMPOSED)
                ndosym = 1
                npencl = ipencl

                CALL fftsym(ndosym,npencl,izproc,nzproc,nzsize,nzglbl,nprocz)

!               RESTORE TRANSFORMED DATA
                DO ipencl = 1,npencl
                    DO kc = 1,nzsize
                        kkc = 2*kc
                        kkm = kkc-1
                        urun(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkm,1,ipencl)*scalez
                        utmp(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkc,1,ipencl)*scalez
                        vrun(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkm,2,ipencl)*scalez
                        vtmp(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkc,2,ipencl)*scalez
                        wrun(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkm,3,ipencl)*scalez
                        wtmp(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkc,3,ipencl)*scalez
                    END DO
                END DO

                ipencl = 0

            END IF

        END DO

    END DO

!   CHECK FOR ANY REMAINING STORED PENCILS
    IF(ipencl /= 0)THEN
!       DO THE FFT (SYMMETRY IMPOSED)
        ndosym = 1
        npencl = ipencl

        CALL fftsym(ndosym,npencl,izproc,nzproc,nzsize,nzglbl,nprocz)

!       RESTORE TRANSFORMED DATA
        DO ipencl = 1,npencl
            DO kc = 1,nzsize
                kkc = 2*kc
                kkm = kkc-1
                urun(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkm,1,ipencl)*scalez
                utmp(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkc,1,ipencl)*scalez
                vrun(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkm,2,ipencl)*scalez
                vtmp(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkc,2,ipencl)*scalez
                wrun(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkm,3,ipencl)*scalez
                wtmp(ipen(ipencl),jpen(ipencl),kc) = ftpart(kkc,3,ipencl)*scalez
            END DO
        END DO

    END IF

!   =========================================================================

END SUBROUTINE turbft
