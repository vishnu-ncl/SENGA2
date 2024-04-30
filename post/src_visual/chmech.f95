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
wrtspc(3)%wrtdst = "Y1"
wrtspc(4)%wrtdst = "Y2"
wrtspc(5)%wrtdst = "Y3"
wrtspc(6)%wrtdst = "Y4"
wrtspc(7)%wrtdst = "Y5"
wrtspc(8)%wrtdst = "Y6"
wrtspc(9)%wrtdst = "Y7"

wrtrrt(1)%wrtdst = "RRTE_R"
wrtrrt(2)%wrtdst = "RRTE_P"
wrtrrt(3)%wrtdst = "RRTE_Y1"
wrtrrt(4)%wrtdst = "RRTE_Y2"
wrtrrt(5)%wrtdst = "RRTE_Y3"
wrtrrt(6)%wrtdst = "RRTE_Y4"
wrtrrt(7)%wrtdst = "RRTE_Y5"
wrtrrt(8)%wrtdst = "RRTE_Y6"
wrtrrt(9)%wrtdst = "RRTE_Y7"

end subroutine set_mech_tags

