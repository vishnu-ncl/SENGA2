SUBROUTINE d2fdx2_kernel_interior(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv
    real(kind=dp) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep
    real(kind=dp) :: fdifam,fdifbm,fdifcm,fdifdm,fdifem

    fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0)) 
    fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))
    fdifdp = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-4,0,0))
    fdifep = functn(OPS_ACC1(5,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifem = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-5,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acofsx*(fdifap-fdifam) + bcofsx*(fdifbp-fdifbm)  &
          + ccofsx*(fdifcp-fdifcm) + dcofsx*(fdifdp-fdifdm)  &
          + ecofsx*(fdifep-fdifem)

END SUBROUTINE d2fdx2_kernel_interior

SUBROUTINE d2fdx2_kernel_lhpoint_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv
    real(kind=dp) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0)) 
    fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(5,0,0)) - functn(OPS_ACC1(0,0,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs1x*fdifap + bcfs1x*fdifbp  &
          + ccfs1x*fdifcp + dcfs1x*fdifdp  &
          + ecfs1x*fdifep

END SUBROUTINE d2fdx2_kernel_lhpoint_4th_onesided    

SUBROUTINE d2fdx2_kernel_lhpoint_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv
    real(kind=dp) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(-1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifcp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs2x*fdifap + bcfs2x*fdifbp  &
          + ccfs2x*fdifcp + dcfs2x*fdifdp  &
          + ecfs2x*fdifep

END SUBROUTINE d2fdx2_kernel_lhpoint_4th_mixed

SUBROUTINE d2fdx2_kernel_lhpoint_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv
    real(kind=dp) :: fdifap,fdifbp
    real(kind=dp) :: fdifam,fdifbm
    
    fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0)) 
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs3x*(fdifap-fdifam) + bcfs3x*(fdifbp-fdifbm)

END SUBROUTINE d2fdx2_kernel_lhpoint_4th_centered 

SUBROUTINE d2fdx2_kernel_lhpoint_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv
    real(kind=dp) :: fdifap,fdifbp,fdifcp
    real(kind=dp) :: fdifam,fdifbm,fdifcm
    
    fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs4x*(fdifap-fdifam)  &
          + bcfs4x*(fdifbp-fdifbm) + ccfs4x*(fdifcp-fdifcm)

END SUBROUTINE d2fdx2_kernel_lhpoint_6th_centered

SUBROUTINE d2fdx2_kernel_lhpoint_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv
    real(kind=dp) :: fdifap,fdifbp,fdifcp,fdifdp
    real(kind=dp) :: fdifam,fdifbm,fdifcm,fdifdm

    fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))
    fdifdp = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-4,0,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs5x*(fdifap-fdifam)  &
          + bcfs5x*(fdifbp-fdifbm) + ccfs5x*(fdifcp-fdifcm)  &
          + dcfs5x*(fdifdp-fdifdm)

END SUBROUTINE d2fdx2_kernel_lhpoint_8th_centered

SUBROUTINE d2fdx2_kernel_rhpoint_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv
    real(kind=dp) :: fdifap,fdifbp,fdifcp,fdifdp
    real(kind=dp) :: fdifam,fdifbm,fdifcm,fdifdm

    fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))
    fdifdp = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-4,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs5x*(fdifap-fdifam)  &
          + bcfs5x*(fdifbp-fdifbm) + ccfs5x*(fdifcp-fdifcm)  &
          + dcfs5x*(fdifdp-fdifdm)

END SUBROUTINE d2fdx2_kernel_rhpoint_8th_centered

SUBROUTINE d2fdx2_kernel_rhpoint_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv
    real(kind=dp) :: fdifap,fdifbp,fdifcp
    real(kind=dp) :: fdifam,fdifbm,fdifcm

    fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs4x*(fdifap-fdifam)  &
          + bcfs4x*(fdifbp-fdifbm) + ccfs4x*(fdifcp-fdifcm)

END SUBROUTINE d2fdx2_kernel_rhpoint_6th_centered

SUBROUTINE d2fdx2_kernel_rhpoint_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv
    real(kind=dp) :: fdifap,fdifbp
    real(kind=dp) :: fdifam,fdifbm

    fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs3x*(fdifap-fdifam) + bcfs3x*(fdifbp-fdifbm)

END SUBROUTINE d2fdx2_kernel_rhpoint_4th_centered

SUBROUTINE d2fdx2_kernel_rhpoint_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv
    real(kind=dp) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(-1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifcp = functn(OPS_ACC1(-2,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(-3,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(-4,0,0)) - functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs2x*fdifap + bcfs2x*fdifbp  &
          + ccfs2x*fdifcp + dcfs2x*fdifdp  &
          + ecfs2x*fdifep

END SUBROUTINE d2fdx2_kernel_rhpoint_4th_mixed

SUBROUTINE d2fdx2_kernel_rhpoint_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv
    real(kind=dp) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(-1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(-2,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifcp = functn(OPS_ACC1(-3,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(-4,0,0)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(-5,0,0)) - functn(OPS_ACC1(0,0,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs1x*fdifap + bcfs1x*fdifbp  &
          + ccfs1x*fdifcp + dcfs1x*fdifdp  &
          + ecfs1x*fdifep

END SUBROUTINE d2fdx2_kernel_rhpoint_4th_onesided

SUBROUTINE d2fdx2_kernel_scaling(fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: fderiv
        
    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdlx2

END SUBROUTINE d2fdx2_kernel_scaling
