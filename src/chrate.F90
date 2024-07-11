SUBROUTINE chrate

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   CHRATE
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   04-JAN-2003:  CREATED
!   22-MAR-2008:  RSC BUG FIX LINDEMANN RATE EXPRESSION
!   04-SEP-2009:  RSC RECODE REACTANT AND PRODUCT COEFFICIENT LISTS
!   13-SEP-2009:  RSC REFORMULATE FOR LINDEMANN RATE EXPRESSIONS
!   04-OCT-2009:  RSC REFORMULATE FOR TROE/SRI RATE EXPRESSIONS
!   08-AUG-2012:  RSC/ZN BUG FIX PRESSURE-DEPENDENT RATES

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   COMPUTES CHEMICAL REACTION RATES
!   USING MULTI-SPECIES MULTI-STEP CHEMISTRY

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   PARAMETERS
!   ==========
    real(kind=8), parameter :: ysmall = 0.000000000000000000000000000001_8,ydenom = 0.000000000000001_8
    real(kind=8), target :: ysmall_ops=ysmall, ydenom_ops=ydenom

!   LOCAL DATA
!   ==========
    real(kind=8) :: racnst,rncnst,reovrr
    real(kind=8) :: ovwmas,gibbsp,scoef,preduc
    real(kind=8) :: ovtst1,tstar2,ovtst3,omalph
    real(kind=8) :: fbroad,ftcent,trats1,trats2,trats3,cfactr,enfact
    real(kind=8) :: acfsri,bcfsrm,ovcsrm,dcfsri,ecfsri
    real(kind=8) :: fornow
    real(kind=8) :: talpha,flcnst,cfcst1,cfcst2,encst1,encst2,dtcnst
    integer(kind=4) :: ispec,isspec
    integer(kind=4) :: istep,ibody
    integer(kind=4) :: iindex,ipower,icoef1,icoef2,itint,icp
    LOGICAL :: flthrd
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   GLOBAL DATA
!   -----------
!   NSPEC total no of species
!   NSTEP total no of steps
!   NBODY total no of third bodies
!   NGIBB total no of Gibbs steps
!   NLIND total no of Lindemann steps
!   NTROE total no of Troe steps
!   NSRIF total no of SRI steps
!   NSPCMX max no of species
!   NSTPMX max no of steps
!   NSSMAX max length of step species-list
!   NRSMAX max length of step reactant- and product-lists
!   NBDYMX max no of third bodies
!   NLLMAX max no of Lindemann steps
!   RPARAM(1,NSTEP) pre-exponential factor ln(A)
!   RPARAM(2,NSTEP) temperature exponent n
!   RPARAM(3,NSTEP) activation energy E/R0
!   NSSPEC(NSSMAX,NSTEP) is the step species-list
!   NRSPEC(NRSMAX,NSTEP) is the step reactant-list
!   NPSPEC(NRSMAX,NSTEP) is the step product-list
!   NRCPEC(NRSMAX,NSTEP) is the step coefficient reactant-list
!   NPCPEC(NRSMAX,NSTEP) is the step coefficient product-list
!   NSSLEN(NSTPMX) contains the length of the species-list for each step
!   NRSLEN(NSTPMX) contains the length of the reactant-list for each step
!   NPSLEN(NSTPMX) contains the length of the product-list for each step
!   NRCLEN(NSTPMX) contains the length of the reactant coefficient-list for
!                  each step
!   NPCLEN(NSTPMX) contains the length of the product coefficient-list for
!                  each step
!   WMOLAR(NSPCMX) molar mass of species
!   OVWMOL(NSPCMX) reciprocal of molar mass of species
!   DIFFMU(NSSMAX,NSTPMX) is the species delta-mu for each step
!   DIFFMW(NSSMAX,NSTPMX) is the species (delta-mu times molar mass)
!                         for each step
!   CRSPEC(NSSMAX,NSTPMX) is the reactant stoichiometric coefficient-list
!                         for each step
!   CPSPEC(NSSMAX,NSTPMX) is the product stoichiometric coefficient-list
!                         for each step
!   MBLIST(NSTPMX) contains 0 if no third body in this step
!                  else contains the index number of the third body
!   MGLIST(NSTPMX) contains 0 if no Gibbs function evaluation in this step
!                  else contains the index number of the Gibbs step
!   MLLIST(NSTPMX) contains 0 if no Lindemann rate in this step
!                  else contains the index number of the Lindemann step
!   MTLIST(NSTPMX) contains 0 if no Troe rate in this step
!                  else contains the index number of the Troe step
!   MSLIST(NSTPMX) contains 0 if no SRI rate in this step
!                  else contains the index number of the SRI step
!   EFFY3B(NSPEC,NBDYMX) contains the third-body efficiencies
!   RCLIND(1,MLLIST(NSTPMX)) Lindemann rate pre-exponential factor ln(A)
!   RCLIND(2,MLLIST(NSTPMX)) Lindemann rate temperature exponent n
!   RCLIND(3,MLLIST(NSTPMX)) Lindemann rate activation energy E/R0
!   RCLIND(4,MLLIST(NSTPMX)) Lindemann rate broadening factor
!   RCTROE(1,MTLIST(NSTPMX)) Troe form rate pre-exponential factor ln(A)
!   RCTROE(2,MTLIST(NSTPMX)) Troe form rate temperature exponent n
!   RCTROE(3,MTLIST(NSTPMX)) Troe form rate activation energy E/R0
!   RCTROE(4,MTLIST(NSTPMX)) Troe form rate temperature parameter alpha
!   RCTROE(5,MTLIST(NSTPMX)) Troe form rate temperature parameter no 1*
!   RCTROE(6,MTLIST(NSTPMX)) Troe form rate temperature parameter no 2*
!   RCTROE(7,MTLIST(NSTPMX)) Troe form rate temperature parameter no 3*
!   RCTROE(8,MTLIST(NSTPMX)) Troe form rate c-const no 1
!   RCTROE(9,MTLIST(NSTPMX)) Troe form rate c-const no 2
!   RCTROE(10,MTLIST(NSTPMX)) Troe form rate n-const no 1
!   RCTROE(11,MTLIST(NSTPMX)) Troe form rate n-const no 2
!   RCTROE(12,MTLIST(NSTPMX)) Troe form rate d-const
!   RCSRIF(1,MSLIST(NSTPMX)) SRI form rate pre-exponential factor ln(A)
!   RCSRIF(2,MSLIST(NSTPMX)) SRI form rate temperature exponent n
!   RCSRIF(3,MSLIST(NSTPMX)) SRI form rate activation energy E/R0
!   RCSRIF(4,MSLIST(NSTPMX)) SRI form rate parameter a
!   RCSRIF(5,MSLIST(NSTPMX)) SRI form rate parameter b
!   RCSRIF(6,MSLIST(NSTPMX)) SRI form rate parameter c
!   RCSRIF(7,MSLIST(NSTPMX)) SRI form rate parameter d
!   RCSRIF(8,MSLIST(NSTPMX)) SRI form rate parameter e
!   =========================================================================

!   ZERO THE REACTION RATE ACCUMULATOR ARRAYS
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    DO ispec = 1,nspec
        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))

    END DO

!   =========================================================================

!   RUN THROUGH ALL STEPS IN THE REACTION MECHANISM
!   -----------------------------------------------

    DO istep = 1, nstep

!       =======================================================================

!       INCLUDE THIRD BODY CONCENTRATIONS, IF ANY
!       -----------------------------------------

!       SET THIRD-BODY EVALUATION FLAG
!       RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
        flthrd = .true.

!       CHECK THIRD-BODY LIST
        IF(mblist(istep) /= 0) THEN

!           GET THE INDEX OF THE THIRD BODY
            ibody = mblist(istep)

!           ZERO THE THIRD-BODY CONCENTRATION ACCUMULATOR
            rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           USE THE THIRD-BODY SPECIES-LIST
!           AND THE THIRD-BODY EFFICIENCY LIST
!           TO GET THE THIRD-BODY CONCENTRATION
            DO ispec = 1,nspec

!               GET THE THIRD-BODY EFFICIENCY FROM THE LIST
!               AND COMBINE WITH THE RECIPROCAL MOLAR MASS
                ovwmas = ovwmol(ispec)*effy3b(ispec,ibody)

!               EVALUATE THE THIRD BODY CONCENTRATION
!               AND STORE IN TEMPORARY ARRAY
                rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(maths_kernel_eqJ, "EVALUATE THE THIRD BODY CONCENTRATION", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ovwmas, 1, "real(kind=8)", OPS_READ))

            END DO
!           THIRD-BODY SPECIES-LIST

!          DIAGNOSTICS
!          WRITE(6,*)'3body',ISTEP,IBODY,FLTHRD,STORE3(250,1,1)

        END IF
!                                             STORE3 = THIRD BODY CONCENTRATION
!       =======================================================================

!       EVALUATE THE SPECIFIC REACTION RATE CONSTANT FOR THIS STEP
!       ----------------------------------------------------------
        racnst = rparam(1,istep)
        rncnst = rparam(2,istep)
        reovrr = rparam(3,istep)

        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(maths_kernel_eqAW, "EVALUATE THE SPECIFIC REACTION RATE CONSTANT", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(racnst, 1, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(rncnst, 1, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(reovrr, 1, "real(kind=8)", OPS_READ))

!       DIAGNOSTICS
!       WRITE(6,'("kf",I5,5(1PE12.4))') &
!       ISTEP,RACNST,RNCNST,REOVRR,STORE1(250,1,1),STORE2(250,1,1)

!                                          STORE1 = FORWARD RATE CONSTANT
!                                          STORE2 = LN(FORWARD RATE CONSTANT)
!                                          STORE3 = THIRD BODY CONCENTRATION
!       =====================================================================

!       INCLUDE LINDEMANN RATE EVALUATION, IF ANY
!       -----------------------------------------
!       RSC 13-SEP-2009 REFORMULATE FOR LINDEMANN RATE EXPRESSIONS
        IF(mllist(istep) /= 0) THEN

            racnst = rclind(1,mllist(istep))
            rncnst = rclind(2,mllist(istep))
            reovrr = rclind(3,mllist(istep))
            flcnst = rclind(4,mllist(istep))

            rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(maths_kernel_eqAX, "REFORMULATE FOR LINDEMANN RATE EXPRESSIONS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(racnst, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(rncnst, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(reovrr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flcnst, 1, "real(kind=8)", OPS_READ))

!           RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
!           RESET THIRD BODY EVALUATION FLAG
            flthrd = .false.

        END IF
!       END OF LINDEMANN RATE EVALUATION
!                                            STORE1 = FORWARD RATE CONSTANT
!                                            STORE2 = LN(FORWARD RATE CONSTANT)
!                                            STORE3 = THIRD BODY CONCENTRATION
!       =======================================================================

!       INCLUDE TROE FORM RATE EVALUATION, IF ANY
!       -----------------------------------------
!       RSC 04-OCT-2009
        IF(mtlist(istep) /= 0) THEN

            racnst = rctroe(1,mtlist(istep))
            rncnst = rctroe(2,mtlist(istep))
            reovrr = rctroe(3,mtlist(istep))
            talpha = rctroe(4,mtlist(istep))
            ovtst1 = rctroe(5,mtlist(istep))
            tstar2 = rctroe(6,mtlist(istep))
            ovtst3 = rctroe(7,mtlist(istep))
            cfcst1 = rctroe(8,mtlist(istep))
            cfcst2 = rctroe(9,mtlist(istep))
            encst1 = rctroe(10,mtlist(istep))
            encst2 = rctroe(11,mtlist(istep))
            dtcnst = rctroe(12,mtlist(istep))

            omalph = one - talpha

            rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(maths_kernel_eqAY, "INCLUDE TROE FORM RATE EVALUATION", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(racnst, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(rncnst, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(reovrr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(talpha, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ovtst1, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(tstar2, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ovtst3, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cfcst1, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cfcst2, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(encst1, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(encst2, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(dtcnst, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(omalph, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(clnten, 1, "real(kind=8)", OPS_READ))

!           RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
!           RESET THIRD BODY EVALUATION FLAG
            flthrd = .false.

        END IF
!       END OF TROE FORM RATE EVALUATION
!                                            STORE1 = FORWARD RATE CONSTANT
!                                            STORE2 = LN(FORWARD RATE CONSTANT)
!                                            STORE3 = THIRD BODY CONCENTRATION
!       =======================================================================

!       INCLUDE SRI FORM RATE EVALUATION, IF ANY
!       ----------------------------------------
!       RSC 04-OCT-2009
        IF(mslist(istep) /= 0)THEN

            racnst = rcsrif(1,mslist(istep))
            rncnst = rcsrif(2,mslist(istep))
            reovrr = rcsrif(3,mslist(istep))
            acfsri = rcsrif(4,mslist(istep))
            bcfsrm = rcsrif(5,mslist(istep))
            ovcsrm = rcsrif(6,mslist(istep))
            dcfsri = rcsrif(7,mslist(istep))
            ecfsri = rcsrif(8,mslist(istep))

            rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(maths_kernel_eqAZ, "INCLUDE SRI FORM RATE EVALUATION", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(racnst, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(rncnst, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(reovrr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(acfsri, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(bcfsrm, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ovcsrm, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(dcfsri, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ecfsri, 1, "real(kind=8)", OPS_READ))

!           RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
!           RESET THIRD BODY EVALUATION FLAG
            flthrd = .false.

        END IF
!       END OF SRI FORM RATE EVALUATION
!                                            STORE1 = FORWARD RATE CONSTANT
!                                            STORE2 = LN(FORWARD RATE CONSTANT)
!                                            STORE3 = THIRD BODY CONCENTRATION
!       =======================================================================

!       USE STEP REACTANT-LIST
!       ----------------------
!       TO FIND REACTANT SPECIES FOR THIS STEP
!       HAVING POSITIVE integer STOICHIOMETRIC COEFFICIENTS

!       NOTE: THIRD BODIES ARE NOT INCLUDED IN THE STEP REACTANT-LIST
!       AND ARE TREATED SEPARATELY (SEE BELOW)

        DO isspec = 1, nrslen(istep)

!           GET THE SPECIES NUMBER FROM THE REACTANT-LIST
            ispec = nrspec(isspec,istep)

!           EVALUATE REACTANT CONCENTRATIONS AND MULTIPLY UP
!           CONCENTRATION IS RHO*Y/WMOLAR
!           YRHS CONTAINS RHO*Y
            rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(maths_kernel_eqP, "A = A*max(B*var    zero)", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

        END DO
!       END OF STEP REACTANT-LIST
!                                            STORE1 = FORWARD REACTION RATE
!                                            STORE2 = LN(FORWARD RATE CONSTANT)
!                                            STORE3 = THIRD BODY CONCENTRATION
!       =======================================================================

!       USE STEP REACTANT COEFFICIENT-LIST
!       ----------------------------------
!       TO FIND REACTANT SPECIES FOR THIS STEP
!       HAVING NON-POSITIVE-integer STOICHIOMETRIC COEFFICIENTS
!       AND OBTAIN THE COEFFICIENTS
!       RSC 04-SEP-2009 RECODE REACTANT AND PRODUCT COEFFICIENT LISTS

!       NOTE: THIRD BODIES ARE NOT INCLUDED IN THE STEP REACTANT-LIST
!       AND ARE TREATED SEPARATELY (SEE BELOW)

        IF(nrclen(istep) > 0) THEN

            DO isspec = 1, nrclen(istep)

!               GET THE SPECIES NUMBER FROM THE STEP COEFFICIENT REACTANT-LIST
                ispec = nrcpec(isspec,istep)

!               GET THE REACTANT STOICHIOMETRIC COEFFICIENT
                scoef = crspec(isspec,istep)

!               EVALUATE REACTANT CONCENTRATIONS AND MULTIPLY UP
!               CONCENTRATION IS RHO*Y/WMOLAR
!               YRHS CONTAINS RHO*Y
                rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(maths_kernel_eqBM, "EVALUATE REACTANT CONCENTRATIONS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(scoef, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ysmall_ops, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ydenom_ops, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

            END DO

        END IF
!       END OF STEP REACTANT COEFFICIENT-LIST
!                                            STORE1 = FORWARD REACTION RATE
!                                            STORE2 = LN(FORWARD RATE CONSTANT)
!                                            STORE3 = THIRD BODY CONCENTRATION
!       =======================================================================

!       INCLUDE BACKWARD RATE EVALUATION, IF ANY
!       ----------------------------------------
        IF(mglist(istep) /= 0) THEN

!           =====================================================================

!           GIBBS FUNCTION EVALUATION
!           -------------------------

!           RUN THROUGH ALL SPECIES FOR THIS STEP
            DO isspec = 1, nsslen(istep)

!               GET THE SPECIES NUMBER FROM THE STEP SPECIES-LIST
                ispec = nsspec(isspec,istep)

!               LOCATE THE TEMPERATURE IN AN INTERVAL
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

!               EVALUATE GIBBS FUNCTION FOR EACH SPECIES
                rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(maths_kernel_eqBN, "EVALUATE GIBBS FUNCTION FOR EACH SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(amolgb, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenpy, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(diffmu, nssmax*nstpmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(isspec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(istep, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

!           STEP SPECIES-LIST
!                                           STORE1 = FORWARD REACTION RATE
!                                           STORE2 = LN(BACKWARD RATE CONSTANT)
!                                           STORE3 = THIRD BODY CONCENTRATION
!           =====================================================================

!           FINALISE BACKWARD RATE COEFFICIENT
!           ----------------------------------
            rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(maths_kernel_eqB, "A = exp(A)", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_RW))

!                                               STORE1 = FORWARD REACTION RATE
!                                               STORE2 = BACKWARD RATE CONSTANT
!                                               STORE3 = THIRD BODY CONCENTRATION
!           =====================================================================

!           USE STEP PRODUCT-LIST
!           ---------------------
!           TO FIND PRODUCT SPECIES FOR THIS STEP
!           HAVING POSITIVE integer STOICHIOMETRIC COEFFICIENTS

!           NOTE: THIRD BODIES ARE NOT INCLUDED IN THE STEP PRODUCT-LIST
!           AND ARE TREATED SEPARATELY (SEE BELOW)

            DO isspec = 1, npslen(istep)

!               GET THE SPECIES NUMBER FROM THE PRODUCT-LIST
                ispec = npspec(isspec,istep)

!               EVALUATE PRODUCT CONCENTRATIONS AND MULTIPLY UP
!               CONCENTRATION IS RHO*Y/WMOLAR
!               YRHS CONTAINS RHO*Y
                rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(maths_kernel_eqP, "A = A*max(B*var    zero)", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

            END DO
!           END OF STEP PRODUCT-LIST
!                                               STORE1 = FORWARD REACTION RATE
!                                               STORE2 = BACKWARD REACTION RATE
!                                               STORE3 = THIRD BODY CONCENTRATION
!           =====================================================================

!           USE STEP PRODUCT COEFFICIENT-LIST
!           ---------------------------------
!           TO FIND PRODUCT SPECIES FOR THIS STEP
!           HAVING NON-POSITIVE-integer STOICHIOMETRIC COEFFICIENTS
!           AND OBTAIN THE COEFFICIENTS
!           RSC 04-SEP-2009 RECODE REACTANT AND PRODUCT COEFFICIENT LISTS

!           NOTE: THIRD BODIES ARE NOT INCLUDED IN THE STEP PRODUCT-LIST
!           AND ARE TREATED SEPARATELY (SEE BELOW)

            IF(npclen(istep) > 0) THEN

                DO isspec = 1, npclen(istep)

!                   GET THE SPECIES NUMBER FROM THE STEP COEFFICIENT PRODUCT-LIST
                    ispec = npcpec(isspec,istep)

!                   GET THE PRODUCT STOICHIOMETRIC COEFFICIENT
                    scoef = cpspec(isspec,istep)

!                   EVALUATE PRODUCT CONCENTRATIONS AND MULTIPLY UP
!                   CONCENTRATION IS RHO*Y/WMOLAR
!                   YRHS CONTAINS RHO*Y
                    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
                    call ops_par_loop(maths_kernel_eqBM, "EVALUATE PRODUCT CONCENTRATIONS", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_RW), &
                                    ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                                    ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
                                    ops_arg_gbl(scoef, 1, "real(kind=8)", OPS_READ), &
                                    ops_arg_gbl(ysmall_ops, 1, "real(kind=8)", OPS_READ), &
                                    ops_arg_gbl(ydenom_ops, 1, "real(kind=8)", OPS_READ), &
                                    ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

                END DO

            END IF
!           END OF STEP PRODUCT COEFFICIENT-LIST
!                                               STORE1 = FORWARD REACTION RATE
!                                               STORE2 = BACKWARD REACTION RATE
!                                               STORE3 = THIRD BODY CONCENTRATION
!           =====================================================================

!           EVALUATE NET FORWARD REACTION RATE
!           ----------------------------------
            rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(maths_kernel_eqS, "A = A-B", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

!           =====================================================================

        END IF
!       END OF BACKWARD RATE EVALUATION
!                                            STORE1 = NET FORWARD REACTION RATE
!                                            STORE3 = THIRD BODY CONCENTRATION
!       =======================================================================

!       INCLUDE THIRD BODY CONCENTRATIONS, IF ANY
!       -----------------------------------------

!       CHECK THIRD-BODY LIST
        IF(mblist(istep) /= 0) THEN

!           RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
!           CHECK THIRD-BODY EVALUATION FLAG
            IF(flthrd) THEN

!               INCLUDE THIRD BODY CONCENTRATION IN REACTION RATE FOR THE STEP
                rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(maths_kernel_eqV, "A = A*B", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

            END IF
!           END OF THIRD BODY EVALUATION FLAG

        END IF
!       END OF THIRD BODY TREATMENT
!                                            STORE1 = NET FORWARD REACTION RATE
!       =======================================================================

!       USE STEP SPECIES-LIST AND STEP SPECIES DELTA-LIST
!       -------------------------------------------------
!       TO GET DIFFERENCE OF STOICHIOMETRIC COEFFICIENTS FOR EACH SPECIES
!       (ACTUALLY DELTA-MU TIMES WMOLAR)
!       ADD STEP CONTRIBUTION TO TOTAL RATE FOR EACH SPECIES INVOLVED

        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        DO isspec = 1, nsslen(istep)

!           GET THE SPECIES NUMBER FROM THE STEP SPECIES-LIST
            ispec = nsspec(isspec,istep)
            fornow = diffmw(isspec,istep)
            call ops_par_loop(maths_kernel_eqJ, "A = A + var*B", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(fornow, 1, "real(kind=8)", OPS_READ))

        END DO
!       STEP SPECIES-LIST

!     =========================================================================

    END DO

!   END OF RUN THROUGH ALL STEPS IN THE REACTION MECHANISM
!   ------------------------------------------------------

!   =========================================================================

END SUBROUTINE chrate
