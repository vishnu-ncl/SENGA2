SUBROUTINE d2fdz2_kernel_null(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE d2fdz2_kernel_null

SUBROUTINE d2fdz2_kernel_interior(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep
    real(kind=8) :: fdifam,fdifbm,fdifcm,fdifdm,fdifem

    fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1)) 
    fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
    fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
    fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
    fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-4))
    fdifep = functn(OPS_ACC1(0,0,5)) - functn(OPS_ACC1(0,0,0))
    fdifem = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-5))

    fderiv(OPS_ACC2(0,0,0)) = acofsz*(fdifap-fdifam) + bcofsz*(fdifbp-fdifbm)  &
          + ccofsz*(fdifcp-fdifcm) + dcofsz*(fdifdp-fdifdm)  &
          + ecofsz*(fdifep-fdifem)

END SUBROUTINE d2fdz2_kernel_interior

SUBROUTINE d2fdz2_kernel_lhpoint_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0)) 
    fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(0,0,5)) - functn(OPS_ACC1(0,0,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs1z*fdifap + bcfs1z*fdifbp  &
          + ccfs1z*fdifcp + dcfs1z*fdifdp  &
          + ecfs1z*fdifep

END SUBROUTINE d2fdz2_kernel_lhpoint_4th_onesided    

SUBROUTINE d2fdz2_kernel_lhpoint_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(0,0,-1)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdifcp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs2z*fdifap + bcfs2z*fdifbp  &
          + ccfs2z*fdifcp + dcfs2z*fdifdp  &
          + ecfs2z*fdifep

END SUBROUTINE d2fdz2_kernel_lhpoint_4th_mixed

SUBROUTINE d2fdz2_kernel_lhpoint_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp
    real(kind=8) :: fdifam,fdifbm
    
    fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
    fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0)) 
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs3z*(fdifap-fdifam) + bcfs3z*(fdifbp-fdifbm)

END SUBROUTINE d2fdz2_kernel_lhpoint_4th_centered 

SUBROUTINE d2fdz2_kernel_lhpoint_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp
    real(kind=8) :: fdifam,fdifbm,fdifcm
    
    fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
    fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
    fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs4z*(fdifap-fdifam)  &
          + bcfs4z*(fdifbp-fdifbm) + ccfs4z*(fdifcp-fdifcm)

END SUBROUTINE d2fdz2_kernel_lhpoint_6th_centered

SUBROUTINE d2fdz2_kernel_lhpoint_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp
    real(kind=8) :: fdifam,fdifbm,fdifcm,fdifdm

    fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
    fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
    fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
    fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
    fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-4))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs5z*(fdifap-fdifam)  &
          + bcfs5z*(fdifbp-fdifbm) + ccfs5z*(fdifcp-fdifcm)  &
          + dcfs5z*(fdifdp-fdifdm)

END SUBROUTINE d2fdz2_kernel_lhpoint_8th_centered

SUBROUTINE d2fdz2_kernel_rhpoint_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp
    real(kind=8) :: fdifam,fdifbm,fdifcm,fdifdm

    fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
    fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
    fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
    fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
    fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-4))

    fderiv(OPS_ACC2(0,0,0)) = acfs5z*(fdifap-fdifam)  &
          + bcfs5z*(fdifbp-fdifbm) + ccfs5z*(fdifcp-fdifcm)  &
          + dcfs5z*(fdifdp-fdifdm)

END SUBROUTINE d2fdz2_kernel_rhpoint_8th_centered

SUBROUTINE d2fdz2_kernel_rhpoint_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp
    real(kind=8) :: fdifam,fdifbm,fdifcm

    fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
    fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
    fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
    fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))

    fderiv(OPS_ACC2(0,0,0)) = acfs4z*(fdifap-fdifam)  &
          + bcfs4z*(fdifbp-fdifbm) + ccfs4z*(fdifcp-fdifcm)

END SUBROUTINE d2fdz2_kernel_rhpoint_6th_centered

SUBROUTINE d2fdz2_kernel_rhpoint_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp
    real(kind=8) :: fdifam,fdifbm

    fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
    fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
    fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))

    fderiv(OPS_ACC2(0,0,0)) = acfs3z*(fdifap-fdifam) + bcfs3z*(fdifbp-fdifbm)

END SUBROUTINE d2fdz2_kernel_rhpoint_4th_centered

SUBROUTINE d2fdz2_kernel_rhpoint_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(0,0,-1)) - functn(OPS_ACC1(0,0,0))
    fdifcp = functn(OPS_ACC1(0,0,-2)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(0,0,-3)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(0,0,-4)) - functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acfs2z*fdifap + bcfs2z*fdifbp  &
          + ccfs2z*fdifcp + dcfs2z*fdifdp  &
          + ecfs2z*fdifep

END SUBROUTINE d2fdz2_kernel_rhpoint_4th_mixed

SUBROUTINE d2fdz2_kernel_rhpoint_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep

    fdifap = functn(OPS_ACC1(0,0,-1)) - functn(OPS_ACC1(0,0,0))
    fdifbp = functn(OPS_ACC1(0,0,-2)) - functn(OPS_ACC1(0,0,0))
    fdifcp = functn(OPS_ACC1(0,0,-3)) - functn(OPS_ACC1(0,0,0))
    fdifdp = functn(OPS_ACC1(0,0,-4)) - functn(OPS_ACC1(0,0,0))
    fdifep = functn(OPS_ACC1(0,0,-5)) - functn(OPS_ACC1(0,0,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acfs1z*fdifap + bcfs1z*fdifbp  &
          + ccfs1z*fdifcp + dcfs1z*fdifdp  &
          + ecfs1z*fdifep

END SUBROUTINE d2fdz2_kernel_rhpoint_4th_onesided

SUBROUTINE d2fdz2_kernel_scaling(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv
        
    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdlz2

END SUBROUTINE d2fdz2_kernel_scaling
