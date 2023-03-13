! This file is specific to the chemical mechanism used
!
! Chemical mechanism: NC7H16 36 species, 205 reactions
!
module chmech

use dtyps
implicit none


integer, parameter :: nspec=16

type (tag) :: wrtspc(nspec)
type (tag) :: wrtrrt(nspec)

end module chmech

subroutine set_mech_tags
use dtyps
use chmech

implicit none

!---------------------------------------------
! HDF5 Tags for the chemical mechanism
!---------------------------------------------

wrtspc(1)%wrtdst = "CH4"
wrtspc(2)%wrtdst = "O2"
wrtspc(3)%wrtdst = "CO2"
wrtspc(4)%wrtdst = "H2O"
wrtspc(5)%wrtdst = "H2"
wrtspc(6)%wrtdst = "O"
wrtspc(7)%wrtdst = "OH"
wrtspc(8)%wrtdst = "H"
wrtspc(9)%wrtdst = "HO2"
wrtspc(10)%wrtdst = "H2O2"
wrtspc(11)%wrtdst = "CO"
wrtspc(12)%wrtdst = "CH2O"
wrtspc(13)%wrtdst = "HCO"
wrtspc(14)%wrtdst = "CH3"
wrtspc(15)%wrtdst = "CH3O"
wrtspc(16)%wrtdst = "N2"
                                       
wrtrrt(1)%wrtdst = "RRTE_CH4"
wrtrrt(2)%wrtdst = "RRTE_O2"
wrtrrt(3)%wrtdst = "RRTE_CO2"
wrtrrt(4)%wrtdst = "RRTE_H2O"
wrtrrt(5)%wrtdst = "RRTE_H2"
wrtrrt(6)%wrtdst = "RRTE_O"
wrtrrt(7)%wrtdst = "RRTE_OH"
wrtrrt(8)%wrtdst = "RRTE_H"
wrtrrt(9)%wrtdst = "RRTE_HO2"
wrtrrt(10)%wrtdst = "RRTE_H2O2"
wrtrrt(11)%wrtdst = "RRTE_CO"
wrtrrt(12)%wrtdst = "RRTE_CH2O"
wrtrrt(13)%wrtdst = "RRTE_HCO"
wrtrrt(14)%wrtdst = "RRTE_CH3"
wrtrrt(15)%wrtdst = "RRTE_CH3O"
wrtrrt(16)%wrtdst = "RRTE_N2"

end subroutine set_mech_tags

