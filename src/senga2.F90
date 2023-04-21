PROGRAM senga2
 
    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    use com_senga
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
    INTEGER :: jtime,jrkstp

!   profiling
    real(kind=c_double) :: startTime = 0
    real(kind=c_double) :: endTime = 0

!   BEGIN
!   =====

!   =========================================================================

!   INITIALISATION
!   ==============

    call ops_init(6)
    call ops_set_soa(0)

    call ops_timers ( startTime )
    
    call ops_data_init

!   PARALLEL DOMAIN DECOMPOSITION
    call pardom

!   INITIALISE THE DATA
    call indata

!   RECORD INITIAL CONDITIONS
    call output

!   =========================================================================

!   TIME STEP LOOP
!   ==============

!   RSC 29-DEC-2006 UPDATED INDEXING
    DO jtime = ntime1,ntime2

        itime = jtime

!       =======================================================================

!       RUNGE-KUTTA SUBSTEPS
!       ====================

!       STANDARD SUBSTEPS
!       -----------------
!       RSC 29-DEC-2006 UPDATED INDEXING
        DO jrkstp = 1,nrksm1

            irkstp = jrkstp

!           APPLY BCS ON PRIMITIVE VARIABLES
            call boundt

!           PARALLEL DATA TRANSFER
            call parfer

!           EVALUATE RHS FOR SCALARS
            call rhscal
#ifdef OPS_LAZY
    call ops_execute()
#endif
!           EVALUATE RHS FOR VELOCITIES
            call rhsvel

!           APPLY BCS ON SOURCE TERMS
            call bounds

!           RUNGE-KUTTA ADVANCEMENT
            call lincom

        END DO

!       =======================================================================

!       LAST SUBSTEP IS DIFFERENT
!       -------------------------
        irkstp = nrkstp

!       APPLY BCS ON PRIMITIVE VARIABLES
        call boundt

!       PARALLEL DATA TRANSFER
        call parfer

!       EVALUATE RHS FOR SCALARS
        call rhscal
#ifdef OPS_LAZY
    call ops_execute()
#endif
!       EVALUATE RHS FOR VELOCITIES
        call rhsvel

!       APPLY BCS ON SOURCE TERMS
        call bounds

!       RUNGE-KUTTA ADVANCEMENT
        call fincom

!       =======================================================================

!       UPDATE THE ELAPSED TIME
!       =======================
        etime = etime + tstep

!       =======================================================================

!       SYNCHRONISE THE TIME-DEPENDENT BCS
!       ==================================
        call bountt

!       =======================================================================

!       ADJUST THE TIME STEP
!       ====================
        call adaptt

!       =======================================================================

!       PROCESS THE RESULTS
!       ===================
        call output

!        IF (itime == 2100) THEN
!            call print_dats()
!        END IF

!       =======================================================================
    END DO
!   END OF TIME STEP LOOP

!   =========================================================================

!   TERMINATION
!   ===========
    call ops_timers ( endTime )
    call ops_timing_output( )
    IF (ops_is_root() .eq. 1) THEN
        write (*,'(a,f16.7,a)') 'Max total runtime =', endTime - startTime,' seconds'
    END IF

!   TERMINATE THE PROGRAM
    call finish

!   =========================================================================

    call ops_exit( )

END PROGRAM senga2
