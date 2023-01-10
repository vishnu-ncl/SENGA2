      subroutine print_param
        use com_senga
        
        WRITE(*,'(a,1PE14.7)') 'test_param_xgdlen:  ',xgdlen
        WRITE(*,'(a,1PE14.7)') 'test_param_ygdlen:  ',ygdlen
        WRITE(*,'(a,1PE14.7)') 'test_param_zgdlen:  ',zgdlen
        WRITE(*,'(a,1PE14.7)') 'test_param_deltax:  ',deltax
        WRITE(*,'(a,1PE14.7)') 'test_param_ovdelx:  ',ovdelx
        WRITE(*,'(a,1PE14.7)') 'test_param_ovdlx2:  ',ovdlx2
        WRITE(*,'(a,1PE14.7)') 'test_param_deltay:  ',deltay
        WRITE(*,'(a,1PE14.7)') 'test_param_ovdely:  ',ovdely
        WRITE(*,'(a,1PE14.7)') 'test_param_ovdly2:  ',ovdly2
        WRITE(*,'(a,1PE14.7)') 'test_param_deltaz:  ',deltaz
        WRITE(*,'(a,1PE14.7)') 'test_param_ovdelz:  ',ovdelz
        WRITE(*,'(a,1PE14.7)') 'test_param_ovdlz2:  ',ovdlz2
        WRITE(*,'(a,1PE14.7)') 'test_param_one:  ',one
        WRITE(*,'(a,1PE14.7)') 'test_param_acoeff: ',acoeff
        WRITE(*,'(a,1PE14.7)') 'test_param_bcoeff: ',bcoeff
        WRITE(*,'(a,1PE14.7)') 'test_param_ccoeff: ',ccoeff
        WRITE(*,'(a,1PE14.7)') 'test_param_dcoeff: ',dcoeff
        WRITE(*,'(a,1PE14.7)') 'test_param_ecoeff: ',ecoeff

        WRITE(*,'(a,1PE14.7)') 'test_param_acoffx: ',acoffx
        WRITE(*,'(a,1PE14.7)') 'test_param_bcoffx: ',bcoffx
        WRITE(*,'(a,1PE14.7)') 'test_param_ccoffx: ',ccoffx
        WRITE(*,'(a,1PE14.7)') 'test_param_dcoffx: ',dcoffx
        WRITE(*,'(a,1PE14.7)') 'test_param_ecoffx: ',ecoffx

        WRITE(*,'(a,1PE14.7)') 'test_param_acoffy: ',acoffy
        WRITE(*,'(a,1PE14.7)') 'test_param_bcoffy: ',bcoffy
        WRITE(*,'(a,1PE14.7)') 'test_param_ccoffy: ',ccoffy
        WRITE(*,'(a,1PE14.7)') 'test_param_dcoffy: ',dcoffy
        WRITE(*,'(a,1PE14.7)') 'test_param_ecoffy: ',ecoffy

        WRITE(*,'(a,1PE14.7)') 'test_param_acoffz: ',acoffz
        WRITE(*,'(a,1PE14.7)') 'test_param_bcoffz: ',bcoffz
        WRITE(*,'(a,1PE14.7)') 'test_param_ccoffz: ',ccoffz
        WRITE(*,'(a,1PE14.7)') 'test_param_dcoffz: ',dcoffz
        WRITE(*,'(a,1PE14.7)') 'test_param_ecoffz: ',ecoffz

        WRITE(*,'(a,1PE14.7)') 'test_param_acoef2: ',acoef2
        WRITE(*,'(a,1PE14.7)') 'test_param_bcoef2: ',bcoef2
        WRITE(*,'(a,1PE14.7)') 'test_param_ccoef2: ',ccoef2
        WRITE(*,'(a,1PE14.7)') 'test_param_dcoef2: ',dcoef2

        WRITE(*,'(a,1PE14.7)') 'test_param_acof1x: ',acof1x
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof1x: ',bcof1x
        WRITE(*,'(a,1PE14.7)') 'test_param_ccof1x: ',ccof1x
        WRITE(*,'(a,1PE14.7)') 'test_param_dcof1x: ',dcof1x

        WRITE(*,'(a,1PE14.7)') 'test_param_acof1y: ',acof1y
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof1y: ',bcof1y
        WRITE(*,'(a,1PE14.7)') 'test_param_ccof1y: ',ccof1y
        WRITE(*,'(a,1PE14.7)') 'test_param_dcof1y: ',dcof1y

        WRITE(*,'(a,1PE14.7)') 'test_param_acof1z: ',acof1z
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof1z: ',bcof1z
        WRITE(*,'(a,1PE14.7)') 'test_param_ccof1z: ',ccof1z
        WRITE(*,'(a,1PE14.7)') 'test_param_dcof1z: ',dcof1z

        WRITE(*,'(a,1PE14.7)') 'test_param_acoef3: ',acoef3
        WRITE(*,'(a,1PE14.7)') 'test_param_bcoef3: ',bcoef3

        WRITE(*,'(a,1PE14.7)') 'test_param_acof3x: ',acof3x
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof3x: ',bcof3x

        WRITE(*,'(a,1PE14.7)') 'test_param_acof3y: ',acof3y
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof3y: ',bcof3y

        WRITE(*,'(a,1PE14.7)') 'test_param_acof3z: ',acof3z
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof3z: ',bcof3z

        WRITE(*,'(a,1PE14.7)') 'test_param_acoef4: ',acoef4
        WRITE(*,'(a,1PE14.7)') 'test_param_bcoef4: ',bcoef4
        WRITE(*,'(a,1PE14.7)') 'test_param_ccoef4: ',ccoef4

        WRITE(*,'(a,1PE14.7)') 'test_param_acof4x: ',acof4x
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof4x: ',bcof4x
        WRITE(*,'(a,1PE14.7)') 'test_param_ccof4x: ',ccof4x

        WRITE(*,'(a,1PE14.7)') 'test_param_acof4y: ',acof4y
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof4y: ',bcof4y
        WRITE(*,'(a,1PE14.7)') 'test_param_ccof4y: ',ccof4y

        WRITE(*,'(a,1PE14.7)') 'test_param_acof4z: ',acof4z
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof4z: ',bcof4z
        WRITE(*,'(a,1PE14.7)') 'test_param_ccof4z: ',ccof4z

        WRITE(*,'(a,1PE14.7)') 'test_param_acoef5: ',acoef5
        WRITE(*,'(a,1PE14.7)') 'test_param_bcoef5: ',bcoef5
        WRITE(*,'(a,1PE14.7)') 'test_param_ccoef5: ',ccoef5
        WRITE(*,'(a,1PE14.7)') 'test_param_dcoef5: ',dcoef5

        WRITE(*,'(a,1PE14.7)') 'test_param_acof5x: ',acof5x
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof5x: ',bcof5x
        WRITE(*,'(a,1PE14.7)') 'test_param_ccof5x: ',ccof5x
        WRITE(*,'(a,1PE14.7)') 'test_param_dcof5x: ',dcof5x

        WRITE(*,'(a,1PE14.7)') 'test_param_acof5y: ',acof5y
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof5y: ',bcof5y
        WRITE(*,'(a,1PE14.7)') 'test_param_ccof5y: ',ccof5y
        WRITE(*,'(a,1PE14.7)') 'test_param_dcof5y: ',dcof5y

        WRITE(*,'(a,1PE14.7)') 'test_param_acof5z: ',acof5z
        WRITE(*,'(a,1PE14.7)') 'test_param_bcof5z: ',bcof5z
        WRITE(*,'(a,1PE14.7)') 'test_param_ccof5z: ',ccof5z
        WRITE(*,'(a,1PE14.7)') 'test_param_dcof5z: ',dcof5z
      end subroutine print_param
