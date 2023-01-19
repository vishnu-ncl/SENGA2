SUBROUTINE d2fdy2_kernel_null(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE d2fdy2_kernel_null

SUBROUTINE d2fdy2_kernel_interior(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep
    real(kind=8) :: fdifam,fdifbm,fdifcm,fdifdm,fdifem

    fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0)) 
    fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
    fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))
    fdifdp = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,0,0))
    fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-4,0))
    fdifep = functn(OPS_ACC1(0,5,0)) - functn(OPS_ACC1(0,0,0))
    fdifem = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-5,0))

    fderiv(OPS_ACC2(0,0,0)) = acofsy*(fdifap-fdifam) + bcofsy*(fdifbp-fdifbm)  &
          + ccofsy*(fdifcp-fdifcm) + dcofsy*(fdifdp-fdifdm)  &
          + ecofsy*(fdifep-fdifem)

END SUBROUTINE d2fdy2_kernel_interior

SUBROUTINE d2fdy2_kernel_lhpoint_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0)) 
    fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(0,5,0)) - functn(OPS_ACC1(0,0,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs1y*fdifap + bcfs1y*fdifbp  &
          + ccfs1y*fdifcp + dcfs1y*fdifdp  &
          + ecfs1y*fdifep

END SUBROUTINE d2fdy2_kernel_lhpoint_4th_onesided    

SUBROUTINE d2fdy2_kernel_lhpoint_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(0,-1,0)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
    fdifcp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs2y*fdifap + bcfs2y*fdifbp  &
          + ccfs2y*fdifcp + dcfs2y*fdifdp  &
          + ecfs2y*fdifep

END SUBROUTINE d2fdy2_kernel_lhpoint_4th_mixed

SUBROUTINE d2fdy2_kernel_lhpoint_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp
    real(kind=8) :: fdifam,fdifbm
    
    fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
    fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0)) 
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs3y*(fdifap-fdifam) + bcfs3y*(fdifbp-fdifbm)

END SUBROUTINE d2fdy2_kernel_lhpoint_4th_centered 

SUBROUTINE d2fdy2_kernel_lhpoint_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp
    real(kind=8) :: fdifam,fdifbm,fdifcm
    
    fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
    fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
    fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs4y*(fdifap-fdifam)  &
          + bcfs4y*(fdifbp-fdifbm) + ccfs4y*(fdifcp-fdifcm)

END SUBROUTINE d2fdy2_kernel_lhpoint_6th_centered

SUBROUTINE d2fdy2_kernel_lhpoint_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp
    real(kind=8) :: fdifam,fdifbm,fdifcm,fdifdm

    fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
    fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
    fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))
    fdifdp = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,0,0))
    fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-4,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs5y*(fdifap-fdifam)  &
          + bcfs5y*(fdifbp-fdifbm) + ccfs5y*(fdifcp-fdifcm)  &
          + dcfs5y*(fdifdp-fdifdm)

END SUBROUTINE d2fdy2_kernel_lhpoint_8th_centered

SUBROUTINE d2fdy2_kernel_rhpoint_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp
    real(kind=8) :: fdifam,fdifbm,fdifcm,fdifdm

    fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
    fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
    fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))
    fdifdp = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,0,0))
    fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-4,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs5y*(fdifap-fdifam)  &
          + bcfs5y*(fdifbp-fdifbm) + ccfs5y*(fdifcp-fdifcm)  &
          + dcfs5y*(fdifdp-fdifdm)

END SUBROUTINE d2fdy2_kernel_rhpoint_8th_centered

SUBROUTINE d2fdy2_kernel_rhpoint_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp
    real(kind=8) :: fdifam,fdifbm,fdifcm

    fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
    fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
    fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs4y*(fdifap-fdifam)  &
          + bcfs4y*(fdifbp-fdifbm) + ccfs4y*(fdifcp-fdifcm)

END SUBROUTINE d2fdy2_kernel_rhpoint_6th_centered

SUBROUTINE d2fdy2_kernel_rhpoint_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp
    real(kind=8) :: fdifam,fdifbm

    fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
    fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs3y*(fdifap-fdifam) + bcfs3y*(fdifbp-fdifbm)

END SUBROUTINE d2fdy2_kernel_rhpoint_4th_centered

SUBROUTINE d2fdy2_kernel_rhpoint_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(0,-1,0)) - functn(OPS_ACC1(0,0,0))
    fdifcp = functn(OPS_ACC1(0,-2,0)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(0,-3,0)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(0,-4,0)) - functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs2y*fdifap + bcfs2y*fdifbp  &
          + ccfs2y*fdifcp + dcfs2y*fdifdp  &
          + ecfs2y*fdifep

END SUBROUTINE d2fdy2_kernel_rhpoint_4th_mixed

SUBROUTINE d2fdy2_kernel_rhpoint_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(0,-1,0)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(0,-2,0)) - functn(OPS_ACC1(0,0,0))
    fdifcp = functn(OPS_ACC1(0,-3,0)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(0,-4,0)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(0,-5,0)) - functn(OPS_ACC1(0,0,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs1y*fdifap + bcfs1y*fdifbp  &
          + ccfs1y*fdifcp + dcfs1y*fdifdp  &
          + ecfs1y*fdifep

END SUBROUTINE d2fdy2_kernel_rhpoint_4th_onesided

SUBROUTINE d2fdy2_kernel_scaling(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv
        
    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdly2

END SUBROUTINE d2fdy2_kernel_scaling
