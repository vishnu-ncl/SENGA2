! This file is specific to the chemical mechanism used
!
! Chemical mechanism: NC7H16 36 species, 205 reactions
!
module chmech

use decl_var
implicit none

type (tag),dimension(:),allocatable :: wrtspc
type (tag),dimension(:),allocatable :: wrtrrt

end module chmech

subroutine set_mech_tags
use decl_var
use chmech

implicit none

!---------------------------------------------
! HDF5 Tags for the chemical mechanism
!---------------------------------------------

wrtspc(1)%wrtdst = "H2"
wrtspc(2)%wrtdst = "O2"
wrtspc(3)%wrtdst = "H2O" 
wrtspc(4)%wrtdst = "O"
wrtspc(5)%wrtdst = "OH"
wrtspc(6)%wrtdst = "H"
wrtspc(7)%wrtdst = "HO2"
wrtspc(8)%wrtdst = "H2O2"
wrtspc(9)%wrtdst = "N2"
                                       
wrtrrt(1)%wrtdst = "RRTE_H2"
wrtrrt(2)%wrtdst = "RRTE_O2"
wrtrrt(3)%wrtdst = "RRTE_H2O"
wrtrrt(4)%wrtdst = "RRTE_O"
wrtrrt(5)%wrtdst = "RRTE_OH"
wrtrrt(6)%wrtdst = "RRTE_H"
wrtrrt(7)%wrtdst = "RRTE_HO2"
wrtrrt(8)%wrtdst = "RRTE_H2O2"
wrtrrt(9)%wrtdst = "RRTE_N2"

end subroutine set_mech_tags

