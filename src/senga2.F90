PROGRAM senga2

    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use omp_lib
    use com_ops_senga

!   *************************************************************************

!   SENGA2
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   12-NOV-2002:  CREATED
!   30-JUN-2006:  RSC SENGA2 VERSION 1.0 (beta)
!   29-DEC-2006:  RSC UPDATES TO INDEXING
!   24-SEP-2009:  RSC SENGA2 VERSION 1.1: UPDATES AND BUG FIXES
!   08-AUG-2012:  RSC SENGA2 VERSION 1.2: UPDATES AND BUG FIXES
!   08-MAY-2013:  RSC SENGA2 VERSION 2.0: UPDATES AND BUG FIXES
!   01-JUL-2015:  RSC SENGA2 VERSION 2.1: UPDATES AND BUG FIXES

!   BASED CLOSELY ON SENGA, CREATED 01-AUG-1996
!   01-AUG-1997:  RSC PARALLEL VERSION 1.0
!   1997-2002:    DR KARL JENKINS: PARALLEL DNS OF FLAME KERNELS
!   2001-2005:    DR NILAN CHAKRABORTY: INFLOW BOUNDARY CONDITIONS
!   2008-2009:    DR TOM DUNSTAN: BUG FIXES
!   2009-2012:    PRESSURE-DEPENDENT REACTION RATES
!   2012:         RSC/ZACK NICOLAOU BUG FIXES
!   2013:         RSC MIXTURE AVERAGED TRANSPORT
!   2015:         RSC/ALEX PHILPOTT UPDATED WALL BOUNDARY CONDITIONS

!   DESCRIPTION
!   -----------
!   THIS PROGRAM CARRIES OUT DIRECT NUMERICAL SIMULATION
!   OF COMPRESSSIBLE TURBULENT FLOW AND COMBUSTION
!   USING HIGH-ORDER METHODS IN TIME AND SPACE

!   PROVISION IS MADE FOR MULTIPLE SPECIES MASS FRACTIONS
!                         VARIABLE TRANSPORT PROPERTIES
!                         MULTIPLE STEP CHEMISTRY

!   REFERENCES
!   ==========
!   1) Kennedy, C.A., Carpenter, M.H., Lewis, R.M.: "Low-storage
!      Explicit Runge-Kutta Schemes for the Compressible Navier-Stokes
!      Equations", Appl. Numer. Math. 35 (3) 177-219, 2000.

!   2) Lele, S.K.:  "Compact Finite Difference Schemes with Spectral-like
!                    Resolution", J. Comput. Phys. 103, 16-29, 1992.

!   3) Poinsot, T.J., Lele, S.K.: "Boundary Conditions for Direct
!      Simulations of Compressible Viscous Flow", J. Comput. Phys. 101,
!      104-129, 1992.

!   4) Sutherland, J.C., Kennedy, C.A.: "Improved Boundary Conditions for
!      Viscous Reacting Compressible Flows", J. Comput. Phys. 2004.

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
!   RSC 29-DEC-2006 UPDATED INDEXING
    integer(kind=4) :: jtime,jrkstp

!   profiling
    real(kind=c_double) :: startTime = 0
    real(kind=c_double) :: endTime = 0

!   Timings
    DOUBLE PRECISION :: indata_stime, indata_etime, indata_ttime, &
                    parfer_stime, parfer_etime, parfer_ttime, &
                    rhscal_stime, rhscal_etime, rhscal_ttime, &
                    rhsvel_stime, rhsvel_etime, rhsvel_ttime, &
                    boundt_stime, boundt_etime, boundt_ttime, &
                    bounds_stime, bounds_etime, bounds_ttime, &
                    bountt_stime, bountt_etime, bountt_ttime, &
                    lincom_stime, lincom_etime, lincom_ttime, &
                    fincom_stime, fincom_etime, fincom_ttime, &
                    adaptt_stime, adaptt_etime, adaptt_ttime

!   BEGIN
!   =====
    indata_ttime = 0d0
    parfer_ttime = 0d0
    rhscal_ttime = 0d0
    rhsvel_ttime = 0d0
    boundt_ttime = 0d0
    bounds_ttime = 0d0
    bountt_ttime = 0d0
    lincom_ttime = 0d0
    fincom_ttime = 0d0
    adaptt_ttime = 0d0
    total_ttime  = 0d0

!   =========================================================================

!   INITIALISATION
!   ==============

    call ops_init(6)
    call ops_set_soa(1)

    call ops_timers ( startTime )

    call ops_data_init

!   PARALLEL DOMAIN DECOMPOSITION
    call pardom

!   INITIALISE THE DATA
    indata_stime = omp_get_wtime()
    call indata
    indata_etime = omp_get_wtime()
    indata_ttime = indata_ttime + (indata_etime-indata_stime)

!   RECORD INITIAL CONDITIONS
    call output

!   INITIALISE CONSTANT USED IN KERNEL FUNCTIONS FOR CUDA
    call cuda_const_init

#ifdef OPS_LAZY
    call ops_execute()
#endif

!   =========================================================================

!   TIME STEP LOOP
!   ==============

!   RSC 29-DEC-2006 UPDATED INDEXING
    DO jtime = ntime1,ntime2

        itime = jtime

!        IF(nxlprm(1) == 4) THEN
!            CALL inflow
!            CALL p_sync
!        END IF

!       =======================================================================

!       RUNGE-KUTTA SUBSTEPS
!       ====================

!       STANDARD SUBSTEPS
!       -----------------
!       RSC 29-DEC-2006 UPDATED INDEXING
        DO jrkstp = 1,nrksm1

            irkstp = jrkstp

!           APPLY BCS ON PRIMITIVE VARIABLES
            boundt_stime = omp_get_wtime()
            call boundt
            boundt_etime = omp_get_wtime()
            boundt_ttime = boundt_ttime + (boundt_etime-boundt_stime)

!           PARALLEL DATA TRANSFER
            parfer_stime = omp_get_wtime()
            call parfer
            parfer_etime = omp_get_wtime()
            parfer_ttime = parfer_ttime + (parfer_etime-parfer_stime)
#ifdef OPS_LAZY
    call ops_execute()
#endif

!           EVALUATE RHS FOR SCALARS
            rhscal_stime = omp_get_wtime()
            call rhscal
            rhscal_etime = omp_get_wtime()
            rhscal_ttime = rhscal_ttime + (rhscal_etime-rhscal_stime)

!           EVALUATE RHS FOR VELOCITIES
            rhsvel_stime = omp_get_wtime()
            call rhsvel
            rhsvel_etime = omp_get_wtime()
            rhsvel_ttime = rhsvel_ttime + (rhsvel_etime-rhsvel_stime)

!           APPLY BCS ON SOURCE TERMS
            bounds_stime = omp_get_wtime()
            call bounds
            bounds_etime = omp_get_wtime()
            bounds_ttime = bounds_ttime + (bounds_etime-bounds_stime)

!           RUNGE-KUTTA ADVANCEMENT
            lincom_stime = omp_get_wtime()
            call lincom
            lincom_etime = omp_get_wtime()
            lincom_ttime = lincom_ttime + (lincom_etime-lincom_stime)
#ifdef OPS_LAZY
    call ops_execute()
#endif
        END DO

!       =======================================================================

!       LAST SUBSTEP IS DIFFERENT
!       -------------------------
        irkstp = nrkstp

!       APPLY BCS ON PRIMITIVE VARIABLES
        boundt_stime = omp_get_wtime()
        call boundt
        boundt_etime = omp_get_wtime()
        boundt_ttime = boundt_ttime + (boundt_etime-boundt_stime)

!       PARALLEL DATA TRANSFER
        parfer_stime = omp_get_wtime()
        call parfer
        parfer_etime = omp_get_wtime()
        parfer_ttime = parfer_ttime + (parfer_etime-parfer_stime)
#ifdef OPS_LAZY
    call ops_execute()
#endif

!       EVALUATE RHS FOR SCALARS
        rhscal_stime = omp_get_wtime()
        call rhscal
        rhscal_etime = omp_get_wtime()
        rhscal_ttime = rhscal_ttime + (rhscal_etime-rhscal_stime)

!       EVALUATE RHS FOR VELOCITIES
        rhsvel_stime = omp_get_wtime()
        call rhsvel
        rhsvel_etime = omp_get_wtime()
        rhsvel_ttime = rhsvel_ttime + (rhsvel_etime-rhsvel_stime)

!       APPLY BCS ON SOURCE TERMS
        bounds_stime = omp_get_wtime()
        call bounds
        bounds_etime = omp_get_wtime()
        bounds_ttime = bounds_ttime + (bounds_etime-bounds_stime)

!       RUNGE-KUTTA ADVANCEMENT
        fincom_stime = omp_get_wtime()
        call fincom
        fincom_etime = omp_get_wtime()
        fincom_ttime = fincom_ttime + (fincom_etime-fincom_stime)

!       =======================================================================

!       UPDATE THE ELAPSED TIME
!       =======================
        etime = etime + tstep

!       =======================================================================

!       SYNCHRONISE THE TIME-DEPENDENT BCS
!       ==================================
        bountt_stime = omp_get_wtime()
        call bountt
        bountt_etime = omp_get_wtime()
        bountt_ttime = bountt_ttime + (bountt_etime-bountt_stime)

!       =======================================================================

!       ADJUST THE TIME STEP
!       ====================
        adaptt_stime = omp_get_wtime()
        call adaptt
        adaptt_etime = omp_get_wtime()
        adaptt_ttime = adaptt_ttime + (adaptt_etime-adaptt_stime)

!       =======================================================================

!       PROCESS THE RESULTS
!       ===================
        call output

!        IF (itime == 1 .or. itime == 500 .or. itime == 1000) THEN
!            call print_dats()
!            IF (itime == 1000) STOP
!        END IF
#ifdef OPS_LAZY
    call ops_execute()
#endif

!       =======================================================================
    END DO
!   END OF TIME STEP LOOP

!   =========================================================================

!   TERMINATION
!   ===========
    call ops_timers ( endTime )
    call ops_timing_output( )
    IF (ops_is_root() == 1) THEN
        write (*,'(a,i5,a,f16.7,a)') 'IPROC: ', iproc,' indata_ttime  =', indata_ttime,' seconds'
        write (*,'(a,i5,a,f16.7,a)') 'IPROC: ', iproc,' parfer_ttime  =', parfer_ttime,' seconds'
        write (*,'(a,i5,a,f16.7,a)') 'IPROC: ', iproc,' rhscal_ttime  =', rhscal_ttime,' seconds'
        write (*,'(a,i5,a,f16.7,a)') 'IPROC: ', iproc,' rhsvel_ttime  =', rhsvel_ttime,' seconds'
        write (*,'(a,i5,a,f16.7,a)') 'IPROC: ', iproc,' boundt_ttime  =', boundt_ttime,' seconds'
        write (*,'(a,i5,a,f16.7,a)') 'IPROC: ', iproc,' bounds_ttime  =', bounds_ttime,' seconds'
        write (*,'(a,i5,a,f16.7,a)') 'IPROC: ', iproc,' bountt_ttime  =', bountt_ttime,' seconds'
        write (*,'(a,i5,a,f16.7,a)') 'IPROC: ', iproc,' lincom_ttime  =', lincom_ttime,' seconds'
        write (*,'(a,i5,a,f16.7,a)') 'IPROC: ', iproc,' fincom_ttime  =', fincom_ttime,' seconds'
        write (*,'(a,i5,a,f16.7,a)') 'IPROC: ', iproc,' adaptt_ttime  =', adaptt_ttime,' seconds'
        write (*,'(a,i5,a,f16.7,a)') 'IPROC: ', iproc,' total_ttime   =', endTime - startTime,' seconds'
    END IF

!   TERMINATE THE PROGRAM
    call finish

!   =========================================================================

    call ops_exit( )

END PROGRAM senga2
