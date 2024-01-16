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

wrtspc(1)%wrtdst = "R"
wrtspc(2)%wrtdst = "P"
                                       
wrtrrt(1)%wrtdst = "RRTE_R"
wrtrrt(2)%wrtdst = "RRTE_P"

end subroutine set_mech_tags

