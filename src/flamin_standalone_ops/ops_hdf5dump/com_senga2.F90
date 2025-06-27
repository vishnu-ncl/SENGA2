MODULE com_senga

    implicit none

!   GLOBAL GRID SIZE
    integer(kind=4), parameter :: nxglbl=3000, nyglbl=1, nzglbl=1

!   SIZE OF HALO
    integer(kind=4), parameter :: nhalox=5,nhaloy=0,nhaloz=0

    integer(kind=4), parameter :: nspcmx=9

END MODULE
