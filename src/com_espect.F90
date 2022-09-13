MODULE com_espect
!     ESPCOM-------------------------------------------------------------------
 

!     COMMON DATA FOR ENERGY SPECTRUM FUNCTION
!     ----------------------------------------

!     GENERAL SPECTRUM PARAMETERS
!     ---------------------------
!     NUMBER OF SPECTRUM PARAMETERS
use data_types

INTEGER :: nsparm
PARAMETER(nsparm = 4)

!     SPECTRUM PARAMETERS
REAL(KIND=dp) :: sparam(nsparm)


!     DATA SPECIFIC TO THE SPECTRUM FUNCTION
!     --------------------------------------
!     BATCHELOR-TOWNSEND SPECTRUM
REAL(KIND=dp) :: const0,ck0,ovk0,covk0


COMMON/espcom/sparam, const0,ck0,ovk0,covk0

!     ESPCOM-------------------------------------------------------------------
END MODULE com_espect
