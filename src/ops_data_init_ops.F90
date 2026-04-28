
! Auto-generated at 2026-04-28 18:43:09.017470 by ops-translator
SUBROUTINE ops_data_init
  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE set_zero_kernel_module
  USE set_zero_kernel_xdir_module
  USE set_zero_kernel_ydir_module
  USE set_zero_kernel_zdir_module
  USE set_zero_kernel_int_module
  USE OPS_Fortran_hdf5_Declarations
  USE OPS_CONSTANTS

  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  INTEGER(KIND = 4) :: d_size(3)
  INTEGER(KIND = 4) :: d_base(3) = [1, 1, 1]
  !array indexing - start from 1
  INTEGER(KIND = 4) :: d_p(3)
  !max boundary depths for the dat in the possitive direction
  INTEGER(KIND = 4) :: d_m(3)
  !max boundary depths for the dat in the negative direction

  REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: temp_real_null
  INTEGER(KIND = 4), DIMENSION(:), ALLOCATABLE :: temp_int_null

  INTEGER(KIND = 4) :: ispec, jspec, iindex, line, status, ic
  INTEGER(KIND = 4) :: rangexyz(6)

  CHARACTER(LEN = 20) :: buf

  INTEGER(KIND = 4) :: halo_idx
  INTEGER(KIND = 4) :: iter_size(3), base_from(3), base_to(3), dir_from(3), dir_to(3)

  INTEGER(KIND = 4) :: a3d_000(3) = [0, 0, 0]

  INTEGER(KIND = 4) :: stride3d_x(3) = [1, 0, 0]
  INTEGER(KIND = 4) :: stride3d_y(3) = [0, 1, 0]
  INTEGER(KIND = 4) :: stride3d_z(3) = [0, 0, 1]

  INTEGER(KIND = 4) :: stride3d_xy(3) = [1, 1, 0]
  INTEGER(KIND = 4) :: stride3d_xz(3) = [1, 0, 1]
  INTEGER(KIND = 4) :: stride3d_yz(3) = [0, 1, 1]

  CHARACTER(LEN = 10) :: pncont
  CHARACTER(LEN = 4) :: pnxdat
  PARAMETER(pncont = 'input/cont', pnxdat = '.dat')

  CHARACTER(LEN = 60) :: fname
  CHARACTER(LEN = 3) :: pnxhdf
  PARAMETER(pnxhdf = '.h5')
  INTEGER(KIND = 4) :: dtime
  CHARACTER(LEN = 8) :: citime

  !------------------------------------------------------------------------------------------------------------

  INTEGER(KIND = 4) :: a3d_000_to_p400_x(15) = [0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0]
  INTEGER(KIND = 4) :: a3d_000_to_m400_x(15) = [0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0, 0]

  INTEGER(KIND = 4) :: a3d_p100_to_p400_x(12) = [1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0]
  INTEGER(KIND = 4) :: a3d_m100_to_m400_x(12) = [- 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0, 0]

  INTEGER(KIND = 4) :: a3d_p500_to_m500_x(33) = [5, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0, 0, - 5, 0, 0]

  !------------------------------------------------------------------------------------------------------------
  INTEGER(KIND = 4) :: a3d_000_to_p040_y(15) = [0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0]
  INTEGER(KIND = 4) :: a3d_000_to_m040_y(15) = [0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0]

  INTEGER(KIND = 4) :: a3d_p010_to_p040_y(12) = [0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0]
  INTEGER(KIND = 4) :: a3d_m010_to_m040_y(12) = [0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0]

  INTEGER(KIND = 4) :: a3d_p050_to_m050_y(33) = [0, 5, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0, 0, - 5, 0]

  !------------------------------------------------------------------------------------------------------------
  INTEGER(KIND = 4) :: a3d_000_to_p004_z(15) = [0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4]
  INTEGER(KIND = 4) :: a3d_000_to_m004_z(15) = [0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4]

  INTEGER(KIND = 4) :: a3d_p001_to_p004_z(12) = [0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4]
  INTEGER(KIND = 4) :: a3d_m001_to_m004_z(12) = [0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4]

  INTEGER(KIND = 4) :: a3d_p005_to_m005_z(33) = [0, 0, 5, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0, 0, - 5]

  !------------------------------------------------------------------------------------------------------------

  INTEGER(KIND = 4) :: a3d_p000_m440_mixed_xy(39) = [0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0, - 1, 0, 0, - 1, - 1, 0, - 2, 0, 0, - 2, - 2, 0, - 3, 0, 0, - 3, - 3, 0, - 4, 0, 0, - 4, - 4, 0]
  INTEGER(KIND = 4) :: a3d_p010_m430_mixed_xy(75) = [0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, - 1, 1, 0, - 1, 0, 0, - 1, - 1, 0, - 1, - 2, 0, - 1, - 3, 0, - 2, 1, 0, - 2, 0, 0, - 2, - 1, 0, - 2, - 2, 0, - 2, - 3, 0, - 3, 1, 0, - 3, 0, 0, - 3, - 1, 0, - 3, - 2, 0, - 3, - 3, 0, - 4, 1, 0, - 4, 0, 0, - 4, - 1, 0, - 4, - 2, 0, - 4, - 3, 0]
  INTEGER(KIND = 4) :: a3d_p020_m420_mixed_xy(60) = [0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2, 0, - 1, 2, 0, - 1, 1, 0, - 1, - 1, 0, - 1, - 2, 0, - 2, 2, 0, - 2, 1, 0, - 2, - 1, 0, - 2, - 2, 0, - 3, 2, 0, - 3, 1, 0, - 3, - 1, 0, - 3, - 2, 0, - 4, 2, 0, - 4, 1, 0, - 4, - 1, 0, - 4, - 2, 0]
  INTEGER(KIND = 4) :: a3d_p030_m410_mixed_xy(75) = [0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, - 1, 3, 0, - 1, 2, 0, - 1, 1, 0, - 1, 0, 0, - 1, - 1, 0, - 2, 3, 0, - 2, 2, 0, - 2, 1, 0, - 2, 0, 0, - 2, - 1, 0, - 3, 3, 0, - 3, 2, 0, - 3, 1, 0, - 3, 0, 0, - 3, - 1, 0, - 4, 3, 0, - 4, 2, 0, - 4, 1, 0, - 4, 0, 0, - 4, - 1, 0]
  INTEGER(KIND = 4) :: a3d_p040_m400_mixed_xy(39) = [0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, - 1, 1, 0, - 1, 0, 0, - 2, 2, 0, - 2, 0, 0, - 3, 3, 0, - 3, 0, 0, - 4, 4, 0, - 4, 0, 0]

  INTEGER(KIND = 4) :: a3d_p100_m300_mixed_xy(39) = [1, 0, 0, 1, - 1, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, - 1, 1, 0, - 1, 0, 0, - 2, 2, 0, - 2, 0, 0, - 3, 3, 0, - 3, 0, 0]
  INTEGER(KIND = 4) :: a3d_p100_m340_mixed_xy(75) = [1, 0, 0, 1, - 1, 0, 1, - 2, 0, 1, - 3, 0, 1, - 4, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0, - 1, 0, 0, - 1, - 1, 0, - 1, - 2, 0, - 1, - 3, 0, - 1, - 4, 0, - 2, 0, 0, - 2, - 1, 0, - 2, - 2, 0, - 2, - 3, 0, - 2, - 4, 0, - 3, 0, 0, - 3, - 1, 0, - 3, - 2, 0, - 3, - 3, 0, - 3, - 4, 0]
  INTEGER(KIND = 4) :: a3d_p110_m330_mixed_xy(39) = [1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, - 1, 0, 0, - 1, - 1, 0, - 2, 0, 0, - 2, - 2, 0, - 3, 0, 0, - 3, - 3, 0]
  INTEGER(KIND = 4) :: a3d_p120_m320_mixed_xy(60) = [1, 2, 0, 1, 1, 0, 1, - 1, 0, 1, - 2, 0, 0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2, 0, - 1, 2, 0, - 1, 1, 0, - 1, - 1, 0, - 1, - 2, 0, - 2, 2, 0, - 2, 1, 0, - 2, - 1, 0, - 2, - 2, 0, - 3, 2, 0, - 3, 1, 0, - 3, - 1, 0, - 3, - 2, 0]
  INTEGER(KIND = 4) :: a3d_p140_m300_mixed_xy(75) = [1, 4, 0, 1, 3, 0, 1, 2, 0, 1, 1, 0, 1, 0, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, - 1, 4, 0, - 1, 3, 0, - 1, 2, 0, - 1, 1, 0, - 1, 0, 0, - 2, 4, 0, - 2, 3, 0, - 2, 2, 0, - 2, 1, 0, - 2, 0, 0, - 3, 4, 0, - 3, 3, 0, - 3, 2, 0, - 3, 1, 0, - 3, 0, 0]

  INTEGER(KIND = 4) :: a3d_p200_m240_mixed_xy(60) = [2, 0, 0, 2, - 1, 0, 2, - 2, 0, 2, - 3, 0, 2, - 4, 0, 1, 0, 0, 1, - 1, 0, 1, - 2, 0, 1, - 3, 0, 1, - 4, 0, - 1, 0, 0, - 1, - 1, 0, - 1, - 2, 0, - 1, - 3, 0, - 1, - 4, 0, - 2, 0, 0, - 2, - 1, 0, - 2, - 2, 0, - 2, - 3, 0, - 2, - 4, 0]
  INTEGER(KIND = 4) :: a3d_p210_m230_mixed_xy(60) = [2, 1, 0, 2, 0, 0, 2, - 1, 0, 2, - 2, 0, 2, - 3, 0, 1, 1, 0, 1, 0, 0, 1, - 1, 0, 1, - 2, 0, 1, - 3, 0, - 1, 1, 0, - 1, 0, 0, - 1, - 1, 0, - 1, - 2, 0, - 1, - 3, 0, - 2, 1, 0, - 2, 0, 0, - 2, - 1, 0, - 2, - 2, 0, - 2, - 3, 0]
  INTEGER(KIND = 4) :: a3d_p220_m220_mixed_xy(24) = [2, 2, 0, 2, - 2, 0, 1, 1, 0, 1, - 1, 0, - 1, 1, 0, - 1, - 1, 0, - 2, 2, 0, - 2, - 2, 0]
  INTEGER(KIND = 4) :: a3d_p230_m210_mixed_xy(60) = [2, 3, 0, 2, 2, 0, 2, 1, 0, 2, 0, 0, 2, - 1, 0, 1, 3, 0, 1, 2, 0, 1, 1, 0, 1, 0, 0, 1, - 1, 0, - 1, 3, 0, - 1, 2, 0, - 1, 1, 0, - 1, 0, 0, - 1, - 1, 0, - 2, 3, 0, - 2, 2, 0, - 2, 1, 0, - 2, 0, 0, - 2, - 1, 0]
  INTEGER(KIND = 4) :: a3d_p240_m200_mixed_xy(60) = [2, 4, 0, 2, 3, 0, 2, 2, 0, 2, 1, 0, 2, 0, 0, 1, 4, 0, 1, 3, 0, 1, 2, 0, 1, 1, 0, 1, 0, 0, - 1, 4, 0, - 1, 3, 0, - 1, 2, 0, - 1, 1, 0, - 1, 0, 0, - 2, 4, 0, - 2, 3, 0, - 2, 2, 0, - 2, 1, 0, - 2, 0, 0]

  INTEGER(KIND = 4) :: a3d_p300_m100_mixed_xy(39) = [3, 0, 0, 3, - 3, 0, 2, 0, 0, 2, - 2, 0, 1, 0, 0, 1, - 1, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, - 1, 1, 0, - 1, 0, 0]
  INTEGER(KIND = 4) :: a3d_p300_m140_mixed_xy(75) = [3, 0, 0, 3, - 1, 0, 3, - 2, 0, 3, - 3, 0, 3, - 4, 0, 2, 0, 0, 2, - 1, 0, 2, - 2, 0, 2, - 3, 0, 2, - 4, 0, 1, 0, 0, 1, - 1, 0, 1, - 2, 0, 1, - 3, 0, 1, - 4, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0, - 1, 0, 0, - 1, - 1, 0, - 1, - 2, 0, - 1, - 3, 0, - 1, - 4, 0]
  INTEGER(KIND = 4) :: a3d_p320_m120_mixed_xy(60) = [3, 2, 0, 3, 1, 0, 3, - 1, 0, 3, - 2, 0, 2, 2, 0, 2, 1, 0, 2, - 1, 0, 2, - 2, 0, 1, 2, 0, 1, 1, 0, 1, - 1, 0, 1, - 2, 0, 0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2, 0, - 1, 2, 0, - 1, 1, 0, - 1, - 1, 0, - 1, - 2, 0]
  INTEGER(KIND = 4) :: a3d_p330_m110_mixed_xy(39) = [3, 3, 0, 3, 0, 0, 2, 2, 0, 2, 0, 0, 1, 1, 0, 1, 0, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, - 1, 0, 0, - 1, - 1, 0]
  INTEGER(KIND = 4) :: a3d_p330_m330_mixed_xy(36) = [3, 3, 0, 3, - 3, 0, 2, 2, 0, 2, - 2, 0, 1, 1, 0, 1, - 1, 0, - 1, 1, 0, - 1, - 1, 0, - 2, 2, 0, - 2, - 2, 0, - 3, 3, 0, - 3, - 3, 0]
  INTEGER(KIND = 4) :: a3d_p340_m100_mixed_xy(75) = [3, 4, 0, 3, 3, 0, 3, 2, 0, 3, 1, 0, 3, 0, 0, 2, 4, 0, 2, 3, 0, 2, 2, 0, 2, 1, 0, 2, 0, 0, 1, 4, 0, 1, 3, 0, 1, 2, 0, 1, 1, 0, 1, 0, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, - 1, 4, 0, - 1, 3, 0, - 1, 2, 0, - 1, 1, 0, - 1, 0, 0]

  INTEGER(KIND = 4) :: a3d_p400_p040_mixed_xy(39) = [4, 0, 0, 4, - 4, 0, 3, 0, 0, 3, - 3, 0, 2, 0, 0, 2, - 2, 0, 1, 0, 0, 1, - 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0]
  INTEGER(KIND = 4) :: a3d_p410_p030_mixed_xy(75) = [4, 1, 0, 4, 0, 0, 4, - 1, 0, 4, - 2, 0, 4, - 3, 0, 3, 1, 0, 3, 0, 0, 3, - 1, 0, 3, - 2, 0, 3, - 3, 0, 2, 1, 0, 2, 0, 0, 2, - 1, 0, 2, - 2, 0, 2, - 3, 0, 1, 1, 0, 1, 0, 0, 1, - 1, 0, 1, - 2, 0, 1, - 3, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0]
  INTEGER(KIND = 4) :: a3d_p420_m020_mixed_xy(60) = [4, 2, 0, 4, 1, 0, 4, - 1, 0, 4, - 2, 0, 3, 2, 0, 3, 1, 0, 3, - 1, 0, 3, - 2, 0, 2, 2, 0, 2, 1, 0, 2, - 1, 0, 2, - 2, 0, 1, 2, 0, 1, 1, 0, 1, - 1, 0, 1, - 2, 0, 0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2, 0]
  INTEGER(KIND = 4) :: a3d_p430_m010_mixed_xy(75) = [4, 3, 0, 4, 2, 0, 4, 1, 0, 4, 0, 0, 4, - 1, 0, 3, 3, 0, 3, 2, 0, 3, 1, 0, 3, 0, 0, 3, - 1, 0, 2, 3, 0, 2, 2, 0, 2, 1, 0, 2, 0, 0, 2, - 1, 0, 1, 3, 0, 1, 2, 0, 1, 1, 0, 1, 0, 0, 1, - 1, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0]
  INTEGER(KIND = 4) :: a3d_p440_p000_mixed_xy(39) = [4, 4, 0, 4, 0, 0, 3, 3, 0, 3, 0, 0, 2, 2, 0, 2, 0, 0, 1, 1, 0, 1, 0, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0]
  INTEGER(KIND = 4) :: a3d_p440_m440_mixed_xy(48) = [4, 4, 0, 4, - 4, 0, 3, 3, 0, 3, - 3, 0, 2, 2, 0, 2, - 2, 0, 1, 1, 0, 1, - 1, 0, - 1, 1, 0, - 1, - 1, 0, - 2, 2, 0, - 2, - 2, 0, - 3, 3, 0, - 3, - 3, 0, - 4, 4, 0, - 4, - 4, 0]

  INTEGER(KIND = 4) :: a3d_p550_m550_mixed_xy(60) = [5, 5, 0, 5, - 5, 0, 4, 4, 0, 4, - 4, 0, 3, 3, 0, 3, - 3, 0, 2, 2, 0, 2, - 2, 0, 1, 1, 0, 1, - 1, 0, - 1, 1, 0, - 1, - 1, 0, - 2, 2, 0, - 2, - 2, 0, - 3, 3, 0, - 3, - 3, 0, - 4, 4, 0, - 4, - 4, 0, - 5, 5, 0, - 5, - 5, 0]

  !------------------------------------------------------------------------------------------------------------

  INTEGER(KIND = 4) :: a3d_p000_m404_mixed_xz(39) = [0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, - 1, 0, 0, - 1, 0, - 1, - 2, 0, 0, - 2, 0, - 2, - 3, 0, 0, - 3, 0, - 3, - 4, 0, 0, - 4, 0, - 4]
  INTEGER(KIND = 4) :: a3d_p001_m403_mixed_xz(75) = [0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, - 1, 0, 1, - 1, 0, 0, - 1, 0, - 1, - 1, 0, - 2, - 1, 0, - 3, - 2, 0, 1, - 2, 0, 0, - 2, 0, - 1, - 2, 0, - 2, - 2, 0, - 3, - 3, 0, 1, - 3, 0, 0, - 3, 0, - 1, - 3, 0, - 2, - 3, 0, - 3, - 4, 0, 1, - 4, 0, 0, - 4, 0, - 1, - 4, 0, - 2, - 4, 0, - 3]
  INTEGER(KIND = 4) :: a3d_p002_m402_mixed_xz(60) = [0, 0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2, - 1, 0, 2, - 1, 0, 1, - 1, 0, - 1, - 1, 0, - 2, - 2, 0, 2, - 2, 0, 1, - 2, 0, - 1, - 2, 0, - 2, - 3, 0, 2, - 3, 0, 1, - 3, 0, - 1, - 3, 0, - 2, - 4, 0, 2, - 4, 0, 1, - 4, 0, - 1, - 4, 0, - 2]
  INTEGER(KIND = 4) :: a3d_p003_m401_mixed_xz(75) = [0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, - 1, 0, 3, - 1, 0, 2, - 1, 0, 1, - 1, 0, 0, - 1, 0, - 1, - 2, 0, 3, - 2, 0, 2, - 2, 0, 1, - 2, 0, 0, - 2, 0, - 1, - 3, 0, 3, - 3, 0, 2, - 3, 0, 1, - 3, 0, 0, - 3, 0, - 1, - 4, 0, 3, - 4, 0, 2, - 4, 0, 1, - 4, 0, 0, - 4, 0, - 1]
  INTEGER(KIND = 4) :: a3d_p004_m400_mixed_xz(39) = [0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, - 1, 0, 1, - 1, 0, 0, - 2, 0, 2, - 2, 0, 0, - 3, 0, 3, - 3, 0, 0, - 4, 0, 4, - 4, 0, 0]

  INTEGER(KIND = 4) :: a3d_p100_m300_mixed_xz(39) = [1, 0, 0, 1, 0, - 1, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, - 1, 0, 1, - 1, 0, 0, - 2, 0, 2, - 2, 0, 0, - 3, 0, 3, - 3, 0, 0]
  INTEGER(KIND = 4) :: a3d_p100_m304_mixed_xz(75) = [1, 0, 0, 1, 0, - 1, 1, 0, - 2, 1, 0, - 3, 1, 0, - 4, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, - 1, 0, 0, - 1, 0, - 1, - 1, 0, - 2, - 1, 0, - 3, - 1, 0, - 4, - 2, 0, 0, - 2, 0, - 1, - 2, 0, - 2, - 2, 0, - 3, - 2, 0, - 4, - 3, 0, 0, - 3, 0, - 1, - 3, 0, - 2, - 3, 0, - 3, - 3, 0, - 4]
  INTEGER(KIND = 4) :: a3d_p101_m303_mixed_xz(39) = [1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, - 1, 0, 0, - 1, 0, - 1, - 2, 0, 0, - 2, 0, - 2, - 3, 0, 0, - 3, 0, - 3]
  INTEGER(KIND = 4) :: a3d_p102_m302_mixed_xz(60) = [1, 0, 2, 1, 0, 1, 1, 0, - 1, 1, 0, - 2, 0, 0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2, - 1, 0, 2, - 1, 0, 1, - 1, 0, - 1, - 1, 0, - 2, - 2, 0, 2, - 2, 0, 1, - 2, 0, - 1, - 2, 0, - 2, - 3, 0, 2, - 3, 0, 1, - 3, 0, - 1, - 3, 0, - 2]
  INTEGER(KIND = 4) :: a3d_p104_m300_mixed_xz(75) = [1, 0, 4, 1, 0, 3, 1, 0, 2, 1, 0, 1, 1, 0, 0, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, - 1, 0, 4, - 1, 0, 3, - 1, 0, 2, - 1, 0, 1, - 1, 0, 0, - 2, 0, 4, - 2, 0, 3, - 2, 0, 2, - 2, 0, 1, - 2, 0, 0, - 3, 0, 4, - 3, 0, 3, - 3, 0, 2, - 3, 0, 1, - 3, 0, 0]

  INTEGER(KIND = 4) :: a3d_p200_m204_mixed_xz(60) = [2, 0, 0, 2, 0, - 1, 2, 0, - 2, 2, 0, - 3, 2, 0, - 4, 1, 0, 0, 1, 0, - 1, 1, 0, - 2, 1, 0, - 3, 1, 0, - 4, - 1, 0, 0, - 1, 0, - 1, - 1, 0, - 2, - 1, 0, - 3, - 1, 0, - 4, - 2, 0, 0, - 2, 0, - 1, - 2, 0, - 2, - 2, 0, - 3, - 2, 0, - 4]
  INTEGER(KIND = 4) :: a3d_p201_m203_mixed_xz(60) = [2, 0, 1, 2, 0, 0, 2, 0, - 1, 2, 0, - 2, 2, 0, - 3, 1, 0, 1, 1, 0, 0, 1, 0, - 1, 1, 0, - 2, 1, 0, - 3, - 1, 0, 1, - 1, 0, 0, - 1, 0, - 1, - 1, 0, - 2, - 1, 0, - 3, - 2, 0, 1, - 2, 0, 0, - 2, 0, - 1, - 2, 0, - 2, - 2, 0, - 3]
  INTEGER(KIND = 4) :: a3d_p202_m202_mixed_xz(24) = [2, 0, 2, 2, 0, - 2, 1, 0, 1, 1, 0, - 1, - 1, 0, 1, - 1, 0, - 1, - 2, 0, 2, - 2, 0, - 2]
  INTEGER(KIND = 4) :: a3d_p203_m201_mixed_xz(60) = [2, 0, 3, 2, 0, 2, 2, 0, 1, 2, 0, 0, 2, 0, - 1, 1, 0, 3, 1, 0, 2, 1, 0, 1, 1, 0, 0, 1, 0, - 1, - 1, 0, 3, - 1, 0, 2, - 1, 0, 1, - 1, 0, 0, - 1, 0, - 1, - 2, 0, 3, - 2, 0, 2, - 2, 0, 1, - 2, 0, 0, - 2, 0, - 1]
  INTEGER(KIND = 4) :: a3d_p204_m200_mixed_xz(60) = [2, 0, 4, 2, 0, 3, 2, 0, 2, 2, 0, 1, 2, 0, 0, 1, 0, 4, 1, 0, 3, 1, 0, 2, 1, 0, 1, 1, 0, 0, - 1, 0, 4, - 1, 0, 3, - 1, 0, 2, - 1, 0, 1, - 1, 0, 0, - 2, 0, 4, - 2, 0, 3, - 2, 0, 2, - 2, 0, 1, - 2, 0, 0]

  INTEGER(KIND = 4) :: a3d_p300_m100_mixed_xz(39) = [3, 0, 0, 3, 0, - 3, 2, 0, 0, 2, 0, - 2, 1, 0, 0, 1, 0, - 1, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, - 1, 0, 1, - 1, 0, 0]
  INTEGER(KIND = 4) :: a3d_p300_m104_mixed_xz(75) = [3, 0, 0, 3, 0, - 1, 3, 0, - 2, 3, 0, - 3, 3, 0, - 4, 2, 0, 0, 2, 0, - 1, 2, 0, - 2, 2, 0, - 3, 2, 0, - 4, 1, 0, 0, 1, 0, - 1, 1, 0, - 2, 1, 0, - 3, 1, 0, - 4, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, - 1, 0, 0, - 1, 0, - 1, - 1, 0, - 2, - 1, 0, - 3, - 1, 0, - 4]
  INTEGER(KIND = 4) :: a3d_p302_m102_mixed_xz(60) = [3, 0, 2, 3, 0, 1, 3, 0, - 1, 3, 0, - 2, 2, 0, 2, 2, 0, 1, 2, 0, - 1, 2, 0, - 2, 1, 0, 2, 1, 0, 1, 1, 0, - 1, 1, 0, - 2, 0, 0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2, - 1, 0, 2, - 1, 0, 1, - 1, 0, - 1, - 1, 0, - 2]
  INTEGER(KIND = 4) :: a3d_p303_m101_mixed_xz(39) = [3, 0, 3, 3, 0, 0, 2, 0, 2, 2, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, - 1, 0, 0, - 1, 0, - 1]
  INTEGER(KIND = 4) :: a3d_p303_m303_mixed_xz(36) = [3, 0, 3, 3, 0, - 3, 2, 0, 2, 2, 0, - 2, 1, 0, 1, 1, 0, - 1, - 1, 0, 1, - 1, 0, - 1, - 2, 0, 2, - 2, 0, - 2, - 3, 0, 3, - 3, 0, - 3]
  INTEGER(KIND = 4) :: a3d_p304_m100_mixed_xz(75) = [3, 0, 4, 3, 0, 3, 3, 0, 2, 3, 0, 1, 3, 0, 0, 2, 0, 4, 2, 0, 3, 2, 0, 2, 2, 0, 1, 2, 0, 0, 1, 0, 4, 1, 0, 3, 1, 0, 2, 1, 0, 1, 1, 0, 0, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, - 1, 0, 4, - 1, 0, 3, - 1, 0, 2, - 1, 0, 1, - 1, 0, 0]

  INTEGER(KIND = 4) :: a3d_p400_p004_mixed_xz(39) = [4, 0, 0, 4, 0, - 4, 3, 0, 0, 3, 0, - 3, 2, 0, 0, 2, 0, - 2, 1, 0, 0, 1, 0, - 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4]
  INTEGER(KIND = 4) :: a3d_p401_p003_mixed_xz(75) = [4, 0, 1, 4, 0, 0, 4, 0, - 1, 4, 0, - 2, 4, 0, - 3, 3, 0, 1, 3, 0, 0, 3, 0, - 1, 3, 0, - 2, 3, 0, - 3, 2, 0, 1, 2, 0, 0, 2, 0, - 1, 2, 0, - 2, 2, 0, - 3, 1, 0, 1, 1, 0, 0, 1, 0, - 1, 1, 0, - 2, 1, 0, - 3, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3]
  INTEGER(KIND = 4) :: a3d_p402_m002_mixed_xz(60) = [4, 0, 2, 4, 0, 1, 4, 0, - 1, 4, 0, - 2, 3, 0, 2, 3, 0, 1, 3, 0, - 1, 3, 0, - 2, 2, 0, 2, 2, 0, 1, 2, 0, - 1, 2, 0, - 2, 1, 0, 2, 1, 0, 1, 1, 0, - 1, 1, 0, - 2, 0, 0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2]
  INTEGER(KIND = 4) :: a3d_p403_m001_mixed_xz(75) = [4, 0, 3, 4, 0, 2, 4, 0, 1, 4, 0, 0, 4, 0, - 1, 3, 0, 3, 3, 0, 2, 3, 0, 1, 3, 0, 0, 3, 0, - 1, 2, 0, 3, 2, 0, 2, 2, 0, 1, 2, 0, 0, 2, 0, - 1, 1, 0, 3, 1, 0, 2, 1, 0, 1, 1, 0, 0, 1, 0, - 1, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1]
  INTEGER(KIND = 4) :: a3d_p404_p000_mixed_xz(39) = [4, 0, 4, 4, 0, 0, 3, 0, 3, 3, 0, 0, 2, 0, 2, 2, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0]
  INTEGER(KIND = 4) :: a3d_p404_m404_mixed_xz(48) = [4, 0, 4, 4, 0, - 4, 3, 0, 3, 3, 0, - 3, 2, 0, 2, 2, 0, - 2, 1, 0, 1, 1, 0, - 1, - 1, 0, 1, - 1, 0, - 1, - 2, 0, 2, - 2, 0, - 2, - 3, 0, 3, - 3, 0, - 3, - 4, 0, 4, - 4, 0, - 4]

  INTEGER(KIND = 4) :: a3d_p505_m505_mixed_xz(60) = [5, 0, 5, 5, 0, - 5, 4, 0, 4, 4, 0, - 4, 3, 0, 3, 3, 0, - 3, 2, 0, 2, 2, 0, - 2, 1, 0, 1, 1, 0, - 1, - 1, 0, 1, - 1, 0, - 1, - 2, 0, 2, - 2, 0, - 2, - 3, 0, 3, - 3, 0, - 3, - 4, 0, 4, - 4, 0, - 4, - 5, 0, 5, - 5, 0, - 5]

  !------------------------------------------------------------------------------------------------------------

  INTEGER(KIND = 4) :: a3d_p000_m044_mixed_yz(39) = [0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0, - 1, 0, 0, - 1, - 1, 0, - 2, 0, 0, - 2, - 2, 0, - 3, 0, 0, - 3, - 3, 0, - 4, 0, 0, - 4, - 4]
  INTEGER(KIND = 4) :: a3d_p001_m043_mixed_yz(75) = [0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, - 1, 1, 0, - 1, 0, 0, - 1, - 1, 0, - 1, - 2, 0, - 1, - 3, 0, - 2, 1, 0, - 2, 0, 0, - 2, - 1, 0, - 2, - 2, 0, - 2, - 3, 0, - 3, 1, 0, - 3, 0, 0, - 3, - 1, 0, - 3, - 2, 0, - 3, - 3, 0, - 4, 1, 0, - 4, 0, 0, - 4, - 1, 0, - 4, - 2, 0, - 4, - 3]
  INTEGER(KIND = 4) :: a3d_p002_m042_mixed_yz(60) = [0, 0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2, 0, - 1, 2, 0, - 1, 1, 0, - 1, - 1, 0, - 1, - 2, 0, - 2, 2, 0, - 2, 1, 0, - 2, - 1, 0, - 2, - 2, 0, - 3, 2, 0, - 3, 1, 0, - 3, - 1, 0, - 3, - 2, 0, - 4, 2, 0, - 4, 1, 0, - 4, - 1, 0, - 4, - 2]
  INTEGER(KIND = 4) :: a3d_p003_m041_mixed_yz(75) = [0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, - 1, 3, 0, - 1, 2, 0, - 1, 1, 0, - 1, 0, 0, - 1, - 1, 0, - 2, 3, 0, - 2, 2, 0, - 2, 1, 0, - 2, 0, 0, - 2, - 1, 0, - 3, 3, 0, - 3, 2, 0, - 3, 1, 0, - 3, 0, 0, - 3, - 1, 0, - 4, 3, 0, - 4, 2, 0, - 4, 1, 0, - 4, 0, 0, - 4, - 1]
  INTEGER(KIND = 4) :: a3d_p004_m040_mixed_yz(39) = [0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, - 1, 1, 0, - 1, 0, 0, - 2, 2, 0, - 2, 0, 0, - 3, 3, 0, - 3, 0, 0, - 4, 4, 0, - 4, 0]

  INTEGER(KIND = 4) :: a3d_p010_m030_mixed_yz(39) = [0, 1, 0, 0, 1, - 1, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, - 1, 1, 0, - 1, 0, 0, - 2, 2, 0, - 2, 0, 0, - 3, 3, 0, - 3, 0]
  INTEGER(KIND = 4) :: a3d_p010_m034_mixed_yz(75) = [0, 1, 0, 0, 1, - 1, 0, 1, - 2, 0, 1, - 3, 0, 1, - 4, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0, - 1, 0, 0, - 1, - 1, 0, - 1, - 2, 0, - 1, - 3, 0, - 1, - 4, 0, - 2, 0, 0, - 2, - 1, 0, - 2, - 2, 0, - 2, - 3, 0, - 2, - 4, 0, - 3, 0, 0, - 3, - 1, 0, - 3, - 2, 0, - 3, - 3, 0, - 3, - 4]
  INTEGER(KIND = 4) :: a3d_p011_m033_mixed_yz(39) = [0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, - 1, 0, 0, - 1, - 1, 0, - 2, 0, 0, - 2, - 2, 0, - 3, 0, 0, - 3, - 3]
  INTEGER(KIND = 4) :: a3d_p012_m032_mixed_yz(60) = [0, 1, 2, 0, 1, 1, 0, 1, - 1, 0, 1, - 2, 0, 0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2, 0, - 1, 2, 0, - 1, 1, 0, - 1, - 1, 0, - 1, - 2, 0, - 2, 2, 0, - 2, 1, 0, - 2, - 1, 0, - 2, - 2, 0, - 3, 2, 0, - 3, 1, 0, - 3, - 1, 0, - 3, - 2]
  INTEGER(KIND = 4) :: a3d_p014_m030_mixed_yz(75) = [0, 1, 4, 0, 1, 3, 0, 1, 2, 0, 1, 1, 0, 1, 0, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, - 1, 4, 0, - 1, 3, 0, - 1, 2, 0, - 1, 1, 0, - 1, 0, 0, - 2, 4, 0, - 2, 3, 0, - 2, 2, 0, - 2, 1, 0, - 2, 0, 0, - 3, 4, 0, - 3, 3, 0, - 3, 2, 0, - 3, 1, 0, - 3, 0]

  INTEGER(KIND = 4) :: a3d_p020_m024_mixed_yz(60) = [0, 2, 0, 0, 2, - 1, 0, 2, - 2, 0, 2, - 3, 0, 2, - 4, 0, 1, 0, 0, 1, - 1, 0, 1, - 2, 0, 1, - 3, 0, 1, - 4, 0, - 1, 0, 0, - 1, - 1, 0, - 1, - 2, 0, - 1, - 3, 0, - 1, - 4, 0, - 2, 0, 0, - 2, - 1, 0, - 2, - 2, 0, - 2, - 3, 0, - 2, - 4]
  INTEGER(KIND = 4) :: a3d_p021_m023_mixed_yz(60) = [0, 2, 1, 0, 2, 0, 0, 2, - 1, 0, 2, - 2, 0, 2, - 3, 0, 1, 1, 0, 1, 0, 0, 1, - 1, 0, 1, - 2, 0, 1, - 3, 0, - 1, 1, 0, - 1, 0, 0, - 1, - 1, 0, - 1, - 2, 0, - 1, - 3, 0, - 2, 1, 0, - 2, 0, 0, - 2, - 1, 0, - 2, - 2, 0, - 2, - 3]
  INTEGER(KIND = 4) :: a3d_p022_m022_mixed_yz(24) = [0, 2, 2, 0, 2, - 2, 0, 1, 1, 0, 1, - 1, 0, - 1, 1, 0, - 1, - 1, 0, - 2, 2, 0, - 2, - 2]
  INTEGER(KIND = 4) :: a3d_p023_m021_mixed_yz(60) = [0, 2, 3, 0, 2, 2, 0, 2, 1, 0, 2, 0, 0, 2, - 1, 0, 1, 3, 0, 1, 2, 0, 1, 1, 0, 1, 0, 0, 1, - 1, 0, - 1, 3, 0, - 1, 2, 0, - 1, 1, 0, - 1, 0, 0, - 1, - 1, 0, - 2, 3, 0, - 2, 2, 0, - 2, 1, 0, - 2, 0, 0, - 2, - 1]
  INTEGER(KIND = 4) :: a3d_p024_m020_mixed_yz(60) = [0, 2, 4, 0, 2, 3, 0, 2, 2, 0, 2, 1, 0, 2, 0, 0, 1, 4, 0, 1, 3, 0, 1, 2, 0, 1, 1, 0, 1, 0, 0, - 1, 4, 0, - 1, 3, 0, - 1, 2, 0, - 1, 1, 0, - 1, 0, 0, - 2, 4, 0, - 2, 3, 0, - 2, 2, 0, - 2, 1, 0, - 2, 0]

  INTEGER(KIND = 4) :: a3d_p030_m010_mixed_yz(39) = [0, 3, 0, 0, 3, - 3, 0, 2, 0, 0, 2, - 2, 0, 1, 0, 0, 1, - 1, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, - 1, 1, 0, - 1, 0]
  INTEGER(KIND = 4) :: a3d_p030_m014_mixed_yz(75) = [0, 3, 0, 0, 3, - 1, 0, 3, - 2, 0, 3, - 3, 0, 3, - 4, 0, 2, 0, 0, 2, - 1, 0, 2, - 2, 0, 2, - 3, 0, 2, - 4, 0, 1, 0, 0, 1, - 1, 0, 1, - 2, 0, 1, - 3, 0, 1, - 4, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4, 0, - 1, 0, 0, - 1, - 1, 0, - 1, - 2, 0, - 1, - 3, 0, - 1, - 4]
  INTEGER(KIND = 4) :: a3d_p032_m012_mixed_yz(60) = [0, 3, 2, 0, 3, 1, 0, 3, - 1, 0, 3, - 2, 0, 2, 2, 0, 2, 1, 0, 2, - 1, 0, 2, - 2, 0, 1, 2, 0, 1, 1, 0, 1, - 1, 0, 1, - 2, 0, 0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2, 0, - 1, 2, 0, - 1, 1, 0, - 1, - 1, 0, - 1, - 2]
  INTEGER(KIND = 4) :: a3d_p033_m011_mixed_yz(39) = [0, 3, 3, 0, 3, 0, 0, 2, 2, 0, 2, 0, 0, 1, 1, 0, 1, 0, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, - 1, 0, 0, - 1, - 1]
  INTEGER(KIND = 4) :: a3d_p033_m033_mixed_yz(36) = [0, 3, 3, 0, 3, - 3, 0, 2, 2, 0, 2, - 2, 0, 1, 1, 0, 1, - 1, 0, - 1, 1, 0, - 1, - 1, 0, - 2, 2, 0, - 2, - 2, 0, - 3, 3, 0, - 3, - 3]
  INTEGER(KIND = 4) :: a3d_p034_m010_mixed_yz(75) = [0, 3, 4, 0, 3, 3, 0, 3, 2, 0, 3, 1, 0, 3, 0, 0, 2, 4, 0, 2, 3, 0, 2, 2, 0, 2, 1, 0, 2, 0, 0, 1, 4, 0, 1, 3, 0, 1, 2, 0, 1, 1, 0, 1, 0, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, - 1, 4, 0, - 1, 3, 0, - 1, 2, 0, - 1, 1, 0, - 1, 0]

  INTEGER(KIND = 4) :: a3d_p040_p004_mixed_yz(39) = [0, 4, 0, 0, 4, - 4, 0, 3, 0, 0, 3, - 3, 0, 2, 0, 0, 2, - 2, 0, 1, 0, 0, 1, - 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3, 0, 0, - 4]
  INTEGER(KIND = 4) :: a3d_p041_p003_mixed_yz(75) = [0, 4, 1, 0, 4, 0, 0, 4, - 1, 0, 4, - 2, 0, 4, - 3, 0, 3, 1, 0, 3, 0, 0, 3, - 1, 0, 3, - 2, 0, 3, - 3, 0, 2, 1, 0, 2, 0, 0, 2, - 1, 0, 2, - 2, 0, 2, - 3, 0, 1, 1, 0, 1, 0, 0, 1, - 1, 0, 1, - 2, 0, 1, - 3, 0, 0, 1, 0, 0, 0, 0, 0, - 1, 0, 0, - 2, 0, 0, - 3]
  INTEGER(KIND = 4) :: a3d_p042_m002_mixed_yz(60) = [0, 4, 2, 0, 4, 1, 0, 4, - 1, 0, 4, - 2, 0, 3, 2, 0, 3, 1, 0, 3, - 1, 0, 3, - 2, 0, 2, 2, 0, 2, 1, 0, 2, - 1, 0, 2, - 2, 0, 1, 2, 0, 1, 1, 0, 1, - 1, 0, 1, - 2, 0, 0, 2, 0, 0, 1, 0, 0, - 1, 0, 0, - 2]
  INTEGER(KIND = 4) :: a3d_p043_m001_mixed_yz(75) = [0, 4, 3, 0, 4, 2, 0, 4, 1, 0, 4, 0, 0, 4, - 1, 0, 3, 3, 0, 3, 2, 0, 3, 1, 0, 3, 0, 0, 3, - 1, 0, 2, 3, 0, 2, 2, 0, 2, 1, 0, 2, 0, 0, 2, - 1, 0, 1, 3, 0, 1, 2, 0, 1, 1, 0, 1, 0, 0, 1, - 1, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, - 1]
  INTEGER(KIND = 4) :: a3d_p044_p000_mixed_yz(39) = [0, 4, 4, 0, 4, 0, 0, 3, 3, 0, 3, 0, 0, 2, 2, 0, 2, 0, 0, 1, 1, 0, 1, 0, 0, 0, 4, 0, 0, 3, 0, 0, 2, 0, 0, 1, 0, 0, 0]
  INTEGER(KIND = 4) :: a3d_p044_m044_mixed_yz(48) = [0, 4, 4, 0, 4, - 4, 0, 3, 3, 0, 3, - 3, 0, 2, 2, 0, 2, - 2, 0, 1, 1, 0, 1, - 1, 0, - 1, 1, 0, - 1, - 1, 0, - 2, 2, 0, - 2, - 2, 0, - 3, 3, 0, - 3, - 3, 0, - 4, 4, 0, - 4, - 4]

  INTEGER(KIND = 4) :: a3d_p055_m055_mixed_yz(60) = [0, 5, 5, 0, 5, - 5, 0, 4, 4, 0, 4, - 4, 0, 3, 3, 0, 3, - 3, 0, 2, 2, 0, 2, - 2, 0, 1, 1, 0, 1, - 1, 0, - 1, 1, 0, - 1, - 1, 0, - 2, 2, 0, - 2, - 2, 0, - 3, 3, 0, - 3, - 3, 0, - 4, 4, 0, - 4, - 4, 0, - 5, 5, 0, - 5, - 5]

  !------------------------------------------------------------------------------------------------------------

  !   *-----------------------------------------OPS Declarations-----------------------------------------------*

  !   Declare OPS Block
  CALL ops_decl_block(3, senga_grid, "SENGA_GRID")

  !   Declare OPS Dats
  d_size = [nxglbl, nyglbl, nzglbl]
  d_m = [0, 0, 0]
  d_p = [0, 0, 0]

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store1, "real(kind=8)", "STORE1")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store2, "real(kind=8)", "STORE2")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store3, "real(kind=8)", "STORE3")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store4, "real(kind=8)", "STORE4")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store5, "real(kind=8)", "STORE5")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store6, "real(kind=8)", "STORE6")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_divm, "real(kind=8)", "DIVM")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ucor, "real(kind=8)", "UCOR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vcor, "real(kind=8)", "VCOR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wcor, "real(kind=8)", "WCOR")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd1x, "real(kind=8)", "WD1X")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd1x, "real(kind=8)", "PD1X")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td1x, "real(kind=8)", "TD1X")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd1y, "real(kind=8)", "WD1Y")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd1y, "real(kind=8)", "PD1Y")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td1y, "real(kind=8)", "TD1Y")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd1z, "real(kind=8)", "WD1Z")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd1z, "real(kind=8)", "PD1Z")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td1z, "real(kind=8)", "TD1Z")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd2x, "real(kind=8)", "WD2X")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd2x, "real(kind=8)", "PD2X")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td2x, "real(kind=8)", "TD2X")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd2y, "real(kind=8)", "WD2Y")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd2y, "real(kind=8)", "PD2Y")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td2y, "real(kind=8)", "TD2Y")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wd2z, "real(kind=8)", "WD2Z")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_pd2z, "real(kind=8)", "PD2Z")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_td2z, "real(kind=8)", "TD2Z")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ufxl, "real(kind=8)", "UFXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vfxl, "real(kind=8)", "VFXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wfxl, "real(kind=8)", "WFXL")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_drun, "real(kind=8)", "DRUN")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_urun, "real(kind=8)", "URUN")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vrun, "real(kind=8)", "VRUN")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wrun, "real(kind=8)", "WRUN")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_erun, "real(kind=8)", "ERUN")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_derr, "real(kind=8)", "DERR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_uerr, "real(kind=8)", "UERR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_verr, "real(kind=8)", "VERR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_werr, "real(kind=8)", "WERR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_eerr, "real(kind=8)", "EERR")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d2prun, "real(kind=8)", "PRN2")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d2trun, "real(kind=8)", "TRN2")

  !---------------------------------------MULTI-DIM DAT--------------------------------------------------------

  d_size = [nxglbl, nyglbl, nzglbl]
  d_m = [0, 0, 0]
  d_p = [0, 0, 0]
  DO ispec = 1, nspcmx
    WRITE(buf, "(A4,I2.2)") "YRUN", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yrun(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A4,I2.2)") "YERR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yerr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A4,I2.2)") "RATE", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_rate(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A4,I2.2)") "RRTE", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_rrte(ispec), "real(kind=8)", TRIM(buf))
  END DO

  d_size = [nxglbl, nyglbl, nzglbl]
  d_m = [- nhalox, - nhaloy, - nhaloz]
  d_p = [nhalox, nhaloy, nhaloz]
  DO iindex = 1, nintmx
    WRITE(buf, "(A6,I2.2)") "ITNDEX", iindex
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_int_null, d_itndex(iindex), "integer(kind=4)", TRIM(buf))
  END DO
  CALL ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, temp_real_null, d_yrhs_mdim, "real(kind=8)", "YRHS-MDIM")

    DO ispec = 1, nspcmx
    WRITE(buf, "(A4,I2.2)") "YRHS", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yrhs(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "CTRANS", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ctrans(ispec), "real(kind=8)", TRIM(buf))
  END DO

  CALL ops_decl_dat(senga_grid, nctmax + 1, d_size, d_base, d_m, d_p, temp_real_null, d_tcoeff, "real(kind=8)", "TCOEFF")
  CALL ops_decl_dat(senga_grid, nctmax, d_size, d_base, d_m, d_p, temp_real_null, d_tderiv, "real(kind=8)", "TDERIV")

  d_size = [1, nyglbl, nzglbl]
  d_m = [0, 0, 0]
  d_p = [0, 0, 0]
  DO ispec = 1, nspcmx
    WRITE(buf, "(A6,I2.2)") "BCLYXL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyxl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "BCLYXR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyxr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRYXL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryxl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRYXR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryxr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "DYDTXL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtxl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "DYDTXR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtxr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "RATEXL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ratexl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "RATEXR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ratexr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRHXL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhxl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRHXR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhxr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "T6BXL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6bxl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "T6BXR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6bxr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "TT6XL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6xl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "TT6XR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6xr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "YINF1", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yinf1(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STYLXL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stlyxl(ispec), "real(kind=8)", TRIM(buf))
  END DO

  d_size = [nxglbl, 1, nzglbl]
  d_m = [0, 0, 0]
  d_p = [0, 0, 0]
  DO ispec = 1, nspcmx
    WRITE(buf, "(A6,I2.2)") "BCLYYL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyyl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "BCLYYR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyyr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRYYL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryyl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRYYR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryyr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "DYDTYL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtyl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "DYDTYR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtyr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "RATEYL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_rateyl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "RATEYR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_rateyr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRHYL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhyl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRHYR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhyr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "T6BYL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6byl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "T6BYR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6byr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "TT6YL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6yl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "TT6YR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6yr(ispec), "real(kind=8)", TRIM(buf))
  END DO

  d_size = [nxglbl, nyglbl, 1]
  d_m = [0, 0, 0]
  d_p = [0, 0, 0]
  DO ispec = 1, nspcmx
    WRITE(buf, "(A6,I2.2)") "BCLYZL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyzl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "BCLYZR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bclyzr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRYZL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryzl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRYZR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stryzr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "DYDTZL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtzl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "DYDTZR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dydtzr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "RATEZL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ratezl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "RATEZR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ratezr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRHZL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhzl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A6,I2.2)") "STRHZR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strhzr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "T6BZL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6bzl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "T6BZR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t6bzr(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "TT6ZL", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6zl(ispec), "real(kind=8)", TRIM(buf))
    WRITE(buf, "(A5,I2.2)") "TT6ZR", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt6zr(ispec), "real(kind=8)", TRIM(buf))
  END DO

  !---------------------------------------WITH HALOS-----------------------------------------------------------

  d_size = [nxglbl, nyglbl, nzglbl]
  d_m = [- nhalox, - nhaloy, - nhaloz]
  d_p = [nhalox, nhaloy, nhaloz]
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_drhs, "real(kind=8)", "DRHS")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_urhs, "real(kind=8)", "URHS")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vrhs, "real(kind=8)", "VRHS")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wrhs, "real(kind=8)", "WRHS")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_erhs, "real(kind=8)", "ERHS")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_utmp, "real(kind=8)", "UTMP")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vtmp, "real(kind=8)", "VTMP")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wtmp, "real(kind=8)", "WTMP")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_prun, "real(kind=8)", "PRUN")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_trun, "real(kind=8)", "TRUN")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_transp, "real(kind=8)", "TRANSP")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_store7, "real(kind=8)", "STORE7")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wmomix, "real(kind=8)", "WMOMIX")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_difmix, "real(kind=8)", "DIFMIX")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tdrmix, "real(kind=8)", "TDRMIX")

  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_combo1, "real(kind=8)", "COMBO1")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_combo2, "real(kind=8)", "COMBO2")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_combo3, "real(kind=8)", "COMBO3")

  !-----------------------------------------Boundary YZ--------------------------------------------------------

  d_size = [1, nyglbl, nzglbl]
  d_m = [0, 0, 0]
  d_p = [0, 0, 0]
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_totyxl, "real(kind=8)", "TOTYXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1xl, "real(kind=8)", "BCL1XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1xr, "real(kind=8)", "BCL1XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2xl, "real(kind=8)", "BCL2XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2xr, "real(kind=8)", "BCL2XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3xl, "real(kind=8)", "BCL3XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3xr, "real(kind=8)", "BCL3XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4xl, "real(kind=8)", "BCL4XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4xr, "real(kind=8)", "BCL4XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5xl, "real(kind=8)", "BCL5XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5xr, "real(kind=8)", "BCL5XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltxl, "real(kind=8)", "BCLTXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltxr, "real(kind=8)", "BCLTXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struxl, "real(kind=8)", "STRUXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struxr, "real(kind=8)", "STRUXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvxl, "real(kind=8)", "STRVXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvxr, "real(kind=8)", "STRVXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwxl, "real(kind=8)", "STRWXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwxr, "real(kind=8)", "STRWXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpxl, "real(kind=8)", "STRPXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpxr, "real(kind=8)", "STRPXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdxl, "real(kind=8)", "STRDXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdxr, "real(kind=8)", "STRDXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtxl, "real(kind=8)", "STRTXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtxr, "real(kind=8)", "STRTXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strexl, "real(kind=8)", "STREXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strexr, "real(kind=8)", "STREXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgxl, "real(kind=8)", "STRGXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgxr, "real(kind=8)", "STRGXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrxl, "real(kind=8)", "STRRXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrxr, "real(kind=8)", "STRRXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtxl, "real(kind=8)", "DUDTXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtxr, "real(kind=8)", "DUDTXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtxl, "real(kind=8)", "DVDTXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtxr, "real(kind=8)", "DVDTXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtxl, "real(kind=8)", "DWDTXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtxr, "real(kind=8)", "DWDTXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtxl, "real(kind=8)", "DTDTXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtxr, "real(kind=8)", "DTDTXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtxl, "real(kind=8)", "DDDTXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtxr, "real(kind=8)", "DDDTXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouxl, "real(kind=8)", "ACOUXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouxr, "real(kind=8)", "ACOUXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2xl, "real(kind=8)", "OVA2XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2xr, "real(kind=8)", "OVA2XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1xl, "real(kind=8)", "GAM1XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1xr, "real(kind=8)", "GAM1XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmxl, "real(kind=8)", "OVGMXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmxr, "real(kind=8)", "OVGMXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtxl, "real(kind=8)", "SYDTXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtxr, "real(kind=8)", "SYDTXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpxl, "real(kind=8)", "SORPXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpxr, "real(kind=8)", "SORPXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1bxl, "real(kind=8)", "T1BXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1bxr, "real(kind=8)", "T1BXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2bxl, "real(kind=8)", "T2BXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2bxr, "real(kind=8)", "T2BXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3bxl, "real(kind=8)", "T3BXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3bxr, "real(kind=8)", "T3BXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4bxl, "real(kind=8)", "T4BXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4bxr, "real(kind=8)", "T4BXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51bxl, "real(kind=8)", "T51BXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51bxr, "real(kind=8)", "T51BXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52bxl, "real(kind=8)", "T52BXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52bxr, "real(kind=8)", "T52BXR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_uinf1, "real(kind=8)", "UINF1")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vinf1, "real(kind=8)", "VINF1")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_winf1, "real(kind=8)", "WINF1")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ustead, "real(kind=8)", "USTEAD")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1xl, "real(kind=8)", "TT1XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1xr, "real(kind=8)", "TT1XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2xl, "real(kind=8)", "TT2XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2xr, "real(kind=8)", "TT2XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3xl, "real(kind=8)", "TT3XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3xr, "real(kind=8)", "TT3XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4xl, "real(kind=8)", "TT4XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4xr, "real(kind=8)", "TT4XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5xl, "real(kind=8)", "TT5XL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5xr, "real(kind=8)", "TT5XR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stluxl, "real(kind=8)", "STLUXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stlvxl, "real(kind=8)", "STLVXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stlwxl, "real(kind=8)", "STLWXL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_stltxl, "real(kind=8)", "STLTXL")

  !-----------------------------------------Boundary XZ--------------------------------------------------------

  d_size = [nxglbl, 1, nzglbl]
  d_m = [0, 0, 0]
  d_p = [0, 0, 0]
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1yl, "real(kind=8)", "BCL1YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1yr, "real(kind=8)", "BCL1YR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2yl, "real(kind=8)", "BCL2YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2yr, "real(kind=8)", "BCL2YR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3yl, "real(kind=8)", "BCL3YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3yr, "real(kind=8)", "BCL3YR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4yl, "real(kind=8)", "BCL4YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4yr, "real(kind=8)", "BCL4YR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5yl, "real(kind=8)", "BCL5YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5yr, "real(kind=8)", "BCL5YR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltyl, "real(kind=8)", "BCLTYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltyr, "real(kind=8)", "BCLTYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struyl, "real(kind=8)", "STRUYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struyr, "real(kind=8)", "STRUYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvyl, "real(kind=8)", "STRVYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvyr, "real(kind=8)", "STRVYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwyl, "real(kind=8)", "STRWYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwyr, "real(kind=8)", "STRWYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpyl, "real(kind=8)", "STRPYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpyr, "real(kind=8)", "STRPYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdyl, "real(kind=8)", "STRDYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdyr, "real(kind=8)", "STRDYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtyl, "real(kind=8)", "STRTYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtyr, "real(kind=8)", "STRTYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_streyl, "real(kind=8)", "STREYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_streyr, "real(kind=8)", "STREYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgyl, "real(kind=8)", "STRGYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgyr, "real(kind=8)", "STRGYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strryl, "real(kind=8)", "STRRYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strryr, "real(kind=8)", "STRRYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtyl, "real(kind=8)", "DUDTYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtyr, "real(kind=8)", "DUDTYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtyl, "real(kind=8)", "DVDTYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtyr, "real(kind=8)", "DVDTYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtyl, "real(kind=8)", "DWDTYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtyr, "real(kind=8)", "DWDTYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtyl, "real(kind=8)", "DTDTYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtyr, "real(kind=8)", "DTDTYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtyl, "real(kind=8)", "DDDTYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtyr, "real(kind=8)", "DDDTYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouyl, "real(kind=8)", "ACOUYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouyr, "real(kind=8)", "ACOUYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2yl, "real(kind=8)", "OVA2YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2yr, "real(kind=8)", "OVA2YR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1yl, "real(kind=8)", "GAM1YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1yr, "real(kind=8)", "GAM1YR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmyl, "real(kind=8)", "OVGMYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmyr, "real(kind=8)", "OVGMYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtyl, "real(kind=8)", "SYDTYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtyr, "real(kind=8)", "SYDTYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpyl, "real(kind=8)", "SORPYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpyr, "real(kind=8)", "SORPYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1byl, "real(kind=8)", "T1BYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1byr, "real(kind=8)", "T1BYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2byl, "real(kind=8)", "T2BYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2byr, "real(kind=8)", "T2BYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3byl, "real(kind=8)", "T3BYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3byr, "real(kind=8)", "T3BYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4byl, "real(kind=8)", "T4BYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4byr, "real(kind=8)", "T4BYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51byl, "real(kind=8)", "T51BYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51byr, "real(kind=8)", "T51BYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52byl, "real(kind=8)", "T52BYL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52byr, "real(kind=8)", "T52BYR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1yl, "real(kind=8)", "TT1YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1yr, "real(kind=8)", "TT1YR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2yl, "real(kind=8)", "TT2YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2yr, "real(kind=8)", "TT2YR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3yl, "real(kind=8)", "TT3YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3yr, "real(kind=8)", "TT3YR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4yl, "real(kind=8)", "TT4YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4yr, "real(kind=8)", "TT4YR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5yl, "real(kind=8)", "TT5YL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5yr, "real(kind=8)", "TT5YR")

  !-----------------------------------------Boundary XY--------------------------------------------------------

  d_size = [nxglbl, nyglbl, 1]
  d_m = [0, 0, 0]
  d_p = [0, 0, 0]
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1zl, "real(kind=8)", "BCL1ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl1zr, "real(kind=8)", "BCL1ZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2zl, "real(kind=8)", "BCL2ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl2zr, "real(kind=8)", "BCL2ZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3zl, "real(kind=8)", "BCL3ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl3zr, "real(kind=8)", "BCL3ZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4zl, "real(kind=8)", "BCL4ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl4zr, "real(kind=8)", "BCL4ZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5zl, "real(kind=8)", "BCL5ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcl5zr, "real(kind=8)", "BCL5ZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltzl, "real(kind=8)", "BCLTZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_bcltzr, "real(kind=8)", "BCLTZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struzl, "real(kind=8)", "STRUZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_struzr, "real(kind=8)", "STRUZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvzl, "real(kind=8)", "STRVZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strvzr, "real(kind=8)", "STRVZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwzl, "real(kind=8)", "STRWZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strwzr, "real(kind=8)", "STRWZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpzl, "real(kind=8)", "STRPZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strpzr, "real(kind=8)", "STRPZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdzl, "real(kind=8)", "STRDZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strdzr, "real(kind=8)", "STRDZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtzl, "real(kind=8)", "STRTZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strtzr, "real(kind=8)", "STRTZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strezl, "real(kind=8)", "STREZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strezr, "real(kind=8)", "STREZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgzl, "real(kind=8)", "STRGZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strgzr, "real(kind=8)", "STRGZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrzl, "real(kind=8)", "STRRZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_strrzr, "real(kind=8)", "STRRZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtzl, "real(kind=8)", "DUDTZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dudtzr, "real(kind=8)", "DUDTZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtzl, "real(kind=8)", "DVDTZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dvdtzr, "real(kind=8)", "DVDTZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtzl, "real(kind=8)", "DWDTZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dwdtzr, "real(kind=8)", "DWDTZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtzl, "real(kind=8)", "DTDTZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dtdtzr, "real(kind=8)", "DTDTZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtzl, "real(kind=8)", "DDDTZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_dddtzr, "real(kind=8)", "DDDTZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouzl, "real(kind=8)", "ACOUZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_acouzr, "real(kind=8)", "ACOUZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2zl, "real(kind=8)", "OVA2ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ova2zr, "real(kind=8)", "OVA2ZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1zl, "real(kind=8)", "GAM1ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_gam1zr, "real(kind=8)", "GAM1ZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmzl, "real(kind=8)", "OVGMZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ovgmzr, "real(kind=8)", "OVGMZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtzl, "real(kind=8)", "SYDTZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sydtzr, "real(kind=8)", "SYDTZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpzl, "real(kind=8)", "SORPZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_sorpzr, "real(kind=8)", "SORPZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1bzl, "real(kind=8)", "T1BZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t1bzr, "real(kind=8)", "T1BZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2bzl, "real(kind=8)", "T2BZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t2bzr, "real(kind=8)", "T2BZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3bzl, "real(kind=8)", "T3BZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t3bzr, "real(kind=8)", "T3BZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4bzl, "real(kind=8)", "T4BZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t4bzr, "real(kind=8)", "T4BZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51bzl, "real(kind=8)", "T51BZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t51bzr, "real(kind=8)", "T51BZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52bzl, "real(kind=8)", "T52BZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_t52bzr, "real(kind=8)", "T52BZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1zl, "real(kind=8)", "TT1ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt1zr, "real(kind=8)", "TT1ZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2zl, "real(kind=8)", "TT2ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt2zr, "real(kind=8)", "TT2ZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3zl, "real(kind=8)", "TT3ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt3zr, "real(kind=8)", "TT3ZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4zl, "real(kind=8)", "TT4ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt4zr, "real(kind=8)", "TT4ZR")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5zl, "real(kind=8)", "TT5ZL")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_tt5zr, "real(kind=8)", "TT5ZR")

  !------------------------------------Only X-direction--------------------------------------------------------

  d_size = [nxglbl, 1, 1]
  d_m = [0, 0, 0]
  d_p = [0, 0, 0]
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_crin, "real(kind=8)", "CRIN")

  !------------------------------------Inflow Filter YZ--------------------------------------------------------

  d_size = [1, nyglbl, nzglbl]
  d_m = [0, 0, 0]
  d_p = [0, 0, 0]
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ufilt, "real(kind=8)", "UFILT")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vfilt, "real(kind=8)", "VFILT")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wfilt, "real(kind=8)", "WFILT")
  DO ispec = 1, nspcmx
    WRITE(buf, "(A5,I2.2)") "YFILT", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yfilt(ispec), "real(kind=8)", TRIM(buf))
  END DO


  d_size = [1, nyglbl, nzglbl]
  d_m = [0, 0, - nfz]
  d_p = [0, 0, nfz]
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_ufold, "real(kind=8)", "UFOLD")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vfold, "real(kind=8)", "VFOLD")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wfold, "real(kind=8)", "WFOLD")
  DO ispec = 1, nspcmx
    WRITE(buf, "(A5,I2.2)") "YFOLD", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yfold(ispec), "real(kind=8)", TRIM(buf))
  END DO


  d_size = [1, nyglbl, nzglbl]
  d_m = [0, - nfy, - nfz]
  d_p = [0, nfy, nfz]
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_urand, "real(kind=8)", "URAND")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vrand, "real(kind=8)", "VRAND")
  CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_wrand, "real(kind=8)", "WRAND")
  DO ispec = 1, nspcmx
    WRITE(buf, "(A5,I2.2)") "YRAND", ispec
    CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yrand(ispec), "real(kind=8)", TRIM(buf))
  END DO

  !------------------------------------Check for COLD/WARM start-----------------------------------------------

  fncont = pncont // pnxdat
  !   Open the input file and find if its cold start/restart
  OPEN(UNIT = nccont, FILE = fncont, STATUS = 'OLD', FORM = 'FORMATTED')
  DO line = 1, 16
    READ(nccont, *)
  END DO
  READ(nccont, *) tstep, ntime1, ntime, nstpsw
  READ(nccont, *)
  READ(nccont, *)
  READ(nccont, *) ntdump, ntrept, ntstat, ndofmt
  READ(nccont, *)
  READ(nccont, *)
  !   COLD START SWITCH (0=COLD START, 1=RESTART), DUMP INPUT FORMAT
  READ(nccont, *) ncdmpi, ndifmt
  DO line = 1, 6
    READ(nccont, *)
  END DO
  READ(nccont, *) inflam
  DO line = 1, 5
    READ(nccont, *)
  END DO
  READ(nccont, *) nspreq
  DO ispec = 1, nspreq
    READ(nccont, *) jspec, yrin(ispec)
  END DO
  DO line = 1, 5
    READ(nccont, *)
  END DO
  READ(nccont, *) ngbcxl, (nxlprm(ic), ic = 1, nbcpri), (rxlprm(ic), ic = 1, nbcprr)

  CLOSE(UNIT = nccont)

    !   ====================
    !   WARM START
    !   ====================
    IF (ncdmpi == 1) THEN

    dtime = INT((ntime1 - 1) / ntdump)
    WRITE(citime, '(I8.8)') dtime

    fname = 'output/timestep' // citime // pnxhdf

      !       -----------------------------------------------------------------------
      !       THIS BLOCK MAY BE MODIFIED AS REQUIRED
      !       TO BLEND INITIAL VELOCITY AND SCALAR FIELDS
      !       WITH PREVIOUSLY DUMPED DATA
      !       -----------------------------------------------------------------------

      !       RESTART FROM FULL DUMP FILES
      !       ----------------------------
      !       READ THE DATA FROM DUMP INPUT FILE
      !       NOTE THAT URUN,VRUN,WRUN,ERUN AND YRUN ARE ALL IN CONSERVATIVE FORM
      IF (ops_is_root() == 1) THEN
      WRITE(*, *) "Warm Start: start step -> ", ntime1
      WRITE(*, *) "Reading from dumped file: ", TRIM(fname)
    END IF

    CALL ops_decl_dat_hdf5(d_drun_dump, senga_grid, 1, "real(kind=8)", "DRUN", TRIM(fname), status)
    CALL ops_decl_dat_hdf5(d_urun_dump, senga_grid, 1, "real(kind=8)", "URUN", TRIM(fname), status)
    CALL ops_decl_dat_hdf5(d_vrun_dump, senga_grid, 1, "real(kind=8)", "VRUN", TRIM(fname), status)
    CALL ops_decl_dat_hdf5(d_wrun_dump, senga_grid, 1, "real(kind=8)", "WRUN", TRIM(fname), status)
    CALL ops_decl_dat_hdf5(d_erun_dump, senga_grid, 1, "real(kind=8)", "ERUN", TRIM(fname), status)

      DO ispec = 1, nspcmx
      WRITE(buf, "(A4,I2.2)") "YRUN", ispec
      CALL ops_decl_dat_hdf5(d_yrun_dump(ispec), senga_grid, 1, "real(kind=8)", TRIM(buf), TRIM(fname), status)
    END DO

  END IF

    IF (nxlprm(1) == 4) THEN
    IF (ncdmpi == 1) THEN

      fname = 'output/inflow' // pnxhdf

      CALL ops_decl_dat_hdf5(d_uinf2, senga_grid, 1, "real(kind=8)", "UINF2", TRIM(fname), status)
      CALL ops_decl_dat_hdf5(d_vinf2, senga_grid, 1, "real(kind=8)", "VINF2", TRIM(fname), status)
      CALL ops_decl_dat_hdf5(d_winf2, senga_grid, 1, "real(kind=8)", "WINF2", TRIM(fname), status)

        DO ispec = 1, nspcmx
        WRITE(buf, "(A5,I2.2)") "YINF2", ispec
        CALL ops_decl_dat_hdf5(d_yinf2(ispec), senga_grid, 1, "real(kind=8)", TRIM(buf), TRIM(fname), status)
      END DO

    ELSE
      d_size = [1, nyglbl, nzglbl]
      d_m = [0, 0, 0]
      d_p = [0, 0, 0]

      CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_uinf2, "real(kind=8)", "UINF2")
      CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_vinf2, "real(kind=8)", "VINF2")
      CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_winf2, "real(kind=8)", "WINF2")
      DO ispec = 1, nspcmx
        WRITE(buf, "(A5,I2.2)") "YINF2", ispec
        CALL ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_yinf2(ispec), "real(kind=8)", TRIM(buf))
      END DO

    END IF
  END IF

    IF (inflam == 1) THEN

    fname = 'input/flame_init' // pnxhdf

    CALL ops_decl_dat_hdf5(d_drun_dump, senga_grid, 1, "real(kind=8)", "DRUN", TRIM(fname), status)
    CALL ops_decl_dat_hdf5(d_urun_dump, senga_grid, 1, "real(kind=8)", "URUN", TRIM(fname), status)
    CALL ops_decl_dat_hdf5(d_vrun_dump, senga_grid, 1, "real(kind=8)", "VRUN", TRIM(fname), status)
    CALL ops_decl_dat_hdf5(d_wrun_dump, senga_grid, 1, "real(kind=8)", "WRUN", TRIM(fname), status)
    CALL ops_decl_dat_hdf5(d_trun_dump, senga_grid, 1, "real(kind=8)", "TRUN", TRIM(fname), status)

      DO ispec = 1, nspcmx
      WRITE(buf, "(A4,I2.2)") "YRUN", ispec
      CALL ops_decl_dat_hdf5(d_yrun_dump(ispec), senga_grid, 1, "real(kind=8)", TRIM(buf), TRIM(fname), status)
    END DO

  END IF

  !------------------------------------OPS Reduction Handles---------------------------------------------------

  CALL ops_decl_reduction_handle(8, h_erdtot, "real(kind=8)", "erdtot")
  CALL ops_decl_reduction_handle(8, h_erutot, "real(kind=8)", "erutot")
  CALL ops_decl_reduction_handle(8, h_ervtot, "real(kind=8)", "ervtot")
  CALL ops_decl_reduction_handle(8, h_erwtot, "real(kind=8)", "erwtot")
  CALL ops_decl_reduction_handle(8, h_eretot, "real(kind=8)", "eretot")
  CALL ops_decl_reduction_handle(8, h_erytot, "real(kind=8)", "erytot")
  CALL ops_decl_reduction_handle(8, h_tket, "real(kind=8)", "tket")
  CALL ops_decl_reduction_handle(8, h_tkes, "real(kind=8)", "tkes")
  CALL ops_decl_reduction_handle(8, h_enstro, "real(kind=8)", "enstro")
  CALL ops_decl_reduction_handle(8, h_ubart, "real(kind=8)", "ubart")
  CALL ops_decl_reduction_handle(8, h_vbart, "real(kind=8)", "vbart")
  CALL ops_decl_reduction_handle(8, h_wbart, "real(kind=8)", "wbart")
  CALL ops_decl_reduction_handle(8, h_uvart, "real(kind=8)", "uvart")
  CALL ops_decl_reduction_handle(8, h_vvart, "real(kind=8)", "vvart")
  CALL ops_decl_reduction_handle(8, h_wvart, "real(kind=8)", "wvart")
  CALL ops_decl_reduction_handle(8, h_umean, "real(kind=8)", "umean")
  CALL ops_decl_reduction_handle(8, h_denom, "real(kind=8)", "denom")

  !------------------------------------OPS Stencil-------------------------------------------------------------

  CALL ops_decl_stencil(3, 1, a3d_000, s3d_000, "0,0,0")

  CALL ops_decl_strided_stencil(3, 1, a3d_000, stride3d_x, s3d_000_strid3d_x, "stride 3D X dir")
  CALL ops_decl_strided_stencil(3, 1, a3d_000, stride3d_y, s3d_000_strid3d_y, "stride 3D Y dir")
  CALL ops_decl_strided_stencil(3, 1, a3d_000, stride3d_z, s3d_000_strid3d_z, "stride 3D Z dir")

  CALL ops_decl_strided_stencil(3, 1, a3d_000, stride3d_xy, s3d_000_strid3d_xy, "stride 3D XY dir")
  CALL ops_decl_strided_stencil(3, 1, a3d_000, stride3d_xz, s3d_000_strid3d_xz, "stride 3D XZ dir")
  CALL ops_decl_strided_stencil(3, 1, a3d_000, stride3d_yz, s3d_000_strid3d_yz, "stride 3D YZ dir")

  !------------------------------------------------------------------------------------------------------------

  CALL ops_decl_stencil(3, 5, a3d_000_to_p400_x, s3d_000_to_p400_x, "0,0,0 to 4,0,0")
  CALL ops_decl_stencil(3, 5, a3d_000_to_m400_x, s3d_000_to_m400_x, "0,0,0 to -4,0,0")

  CALL ops_decl_stencil(3, 4, a3d_p100_to_p400_x, s3d_p100_to_p400_x, "1,0,0 to 4,0,0")
  CALL ops_decl_stencil(3, 4, a3d_m100_to_m400_x, s3d_m100_to_m400_x, "-1,0,0 to -4,0,0")

  CALL ops_decl_stencil(3, 11, a3d_p500_to_m500_x, s3d_p500_to_m500_x, "5,0,0 to -5,0,0")

  !------------------------------------------------------------------------------------------------------------

  CALL ops_decl_stencil(3, 5, a3d_000_to_p040_y, s3d_000_to_p040_y, "0,0,0 to 0,4,0")
  CALL ops_decl_stencil(3, 5, a3d_000_to_m040_y, s3d_000_to_m040_y, "0,0,0 to  0,-4,0")

  CALL ops_decl_stencil(3, 4, a3d_p010_to_p040_y, s3d_p010_to_p040_y, "0,1,0 to 0,4,0")
  CALL ops_decl_stencil(3, 4, a3d_m010_to_m040_y, s3d_m010_to_m040_y, "0,-1,0 to 0,-4,0")

  CALL ops_decl_stencil(3, 11, a3d_p050_to_m050_y, s3d_p050_to_m050_y, "0,5,0 to  0,-5,0")

  !------------------------------------------------------------------------------------------------------------

  CALL ops_decl_stencil(3, 5, a3d_000_to_p004_z, s3d_000_to_p004_z, "0,0,0 to 0,0,4")
  CALL ops_decl_stencil(3, 5, a3d_000_to_m004_z, s3d_000_to_m004_z, "0,0,0 to  0,0,-4")

  CALL ops_decl_stencil(3, 4, a3d_p001_to_p004_z, s3d_p001_to_p004_z, "0,0,1 to 0,0,4")
  CALL ops_decl_stencil(3, 4, a3d_m001_to_m004_z, s3d_m001_to_m004_z, "0,0,-1 to 0,0,-4")

  CALL ops_decl_stencil(3, 11, a3d_p005_to_m005_z, s3d_p005_to_m005_z, "0,0,5 to  0,0,-5")

  !------------------------------------------------------------------------------------------------------------

  CALL ops_decl_stencil(3, 13, a3d_p000_m440_mixed_xy, s3d_p000_m440_mixed_xy, "0,0,0 to -4,-4,0")
  CALL ops_decl_stencil(3, 25, a3d_p010_m430_mixed_xy, s3d_p010_m430_mixed_xy, "0,1,0 to -4,-3,0")
  CALL ops_decl_stencil(3, 20, a3d_p020_m420_mixed_xy, s3d_p020_m420_mixed_xy, "0,2,0 to -4,-2,0")
  CALL ops_decl_stencil(3, 25, a3d_p030_m410_mixed_xy, s3d_p030_m410_mixed_xy, "0,3,0 to -4,-1,0")
  CALL ops_decl_stencil(3, 13, a3d_p040_m400_mixed_xy, s3d_p040_m400_mixed_xy, "0,4,0 to -4,0,0")

  CALL ops_decl_stencil(3, 13, a3d_p100_m300_mixed_xy, s3d_p100_m300_mixed_xy, "1,0,0 to -3,0,0")
  CALL ops_decl_stencil(3, 25, a3d_p100_m340_mixed_xy, s3d_p100_m340_mixed_xy, "1,0,0 to -3,-4,0")
  CALL ops_decl_stencil(3, 13, a3d_p110_m330_mixed_xy, s3d_p110_m330_mixed_xy, "1,1,0 to -3,-3,0")
  CALL ops_decl_stencil(3, 20, a3d_p120_m320_mixed_xy, s3d_p120_m320_mixed_xy, "1,2,0 to -3,-2,0")
  CALL ops_decl_stencil(3, 25, a3d_p140_m300_mixed_xy, s3d_p140_m300_mixed_xy, "1,4,0 to -3,0,0")

  CALL ops_decl_stencil(3, 20, a3d_p200_m240_mixed_xy, s3d_p200_m240_mixed_xy, "2,0,0 to -2,-4,0")
  CALL ops_decl_stencil(3, 20, a3d_p210_m230_mixed_xy, s3d_p210_m230_mixed_xy, "2,1,0 to -2,-3,0")
  CALL ops_decl_stencil(3, 8, a3d_p220_m220_mixed_xy, s3d_p220_m220_mixed_xy, "2,2,0 to -2,-2,0")
  CALL ops_decl_stencil(3, 20, a3d_p230_m210_mixed_xy, s3d_p230_m210_mixed_xy, "2,3,0 to -2,-1,0")
  CALL ops_decl_stencil(3, 20, a3d_p240_m200_mixed_xy, s3d_p240_m200_mixed_xy, "2,4,0 to -2,0,0")

  CALL ops_decl_stencil(3, 13, a3d_p300_m100_mixed_xy, s3d_p300_m100_mixed_xy, "3,0,0 to -1,0,0")
  CALL ops_decl_stencil(3, 25, a3d_p300_m140_mixed_xy, s3d_p300_m140_mixed_xy, "3,0,0 to -1,-4,0")
  CALL ops_decl_stencil(3, 20, a3d_p320_m120_mixed_xy, s3d_p320_m120_mixed_xy, "3,2,0 to -1,-2,0")
  CALL ops_decl_stencil(3, 13, a3d_p330_m110_mixed_xy, s3d_p330_m110_mixed_xy, "3,3,0 to -1,-1,0")
  CALL ops_decl_stencil(3, 12, a3d_p330_m330_mixed_xy, s3d_p330_m330_mixed_xy, "3,3,0 to -3,-3,0")
  CALL ops_decl_stencil(3, 25, a3d_p340_m100_mixed_xy, s3d_p340_m100_mixed_xy, "3,4,0 to -1,0,0")

  CALL ops_decl_stencil(3, 13, a3d_p400_p040_mixed_xy, s3d_p400_p040_mixed_xy, "4,0,0 to 0,-4,0")
  CALL ops_decl_stencil(3, 25, a3d_p410_p030_mixed_xy, s3d_p410_p030_mixed_xy, "4,1,0 to 0,-3,0")
  CALL ops_decl_stencil(3, 20, a3d_p420_m020_mixed_xy, s3d_p420_m020_mixed_xy, "4,2,0 to 0,-2,0")
  CALL ops_decl_stencil(3, 25, a3d_p430_m010_mixed_xy, s3d_p430_m010_mixed_xy, "4,3,0 to 0,-1,0")
  CALL ops_decl_stencil(3, 13, a3d_p440_p000_mixed_xy, s3d_p440_p000_mixed_xy, "4,4,0 to 0,0,0")
  CALL ops_decl_stencil(3, 16, a3d_p440_m440_mixed_xy, s3d_p440_m440_mixed_xy, "4,4,0 to -4,-4,0")

  CALL ops_decl_stencil(3, 20, a3d_p550_m550_mixed_xy, s3d_p550_m550_mixed_xy, "5,5,0 to -5,-5,0")

  !------------------------------------------------------------------------------------------------------------

  CALL ops_decl_stencil(3, 13, a3d_p000_m404_mixed_xz, s3d_p000_m404_mixed_xz, "0,0,0 to -4,0,-4")
  CALL ops_decl_stencil(3, 25, a3d_p001_m403_mixed_xz, s3d_p001_m403_mixed_xz, "0,0,1 to -4,0,-3")
  CALL ops_decl_stencil(3, 20, a3d_p002_m402_mixed_xz, s3d_p002_m402_mixed_xz, "0,0,2 to -4,0,-2")
  CALL ops_decl_stencil(3, 25, a3d_p003_m401_mixed_xz, s3d_p003_m401_mixed_xz, "0,0,3 to -4,0,-1")
  CALL ops_decl_stencil(3, 13, a3d_p004_m400_mixed_xz, s3d_p004_m400_mixed_xz, "0,0,4 to -4,0,0")

  CALL ops_decl_stencil(3, 13, a3d_p100_m300_mixed_xz, s3d_p100_m300_mixed_xz, "1,0,0 to -3,0,0")
  CALL ops_decl_stencil(3, 25, a3d_p100_m304_mixed_xz, s3d_p100_m304_mixed_xz, "1,0,0 to -3,0,-4")
  CALL ops_decl_stencil(3, 13, a3d_p101_m303_mixed_xz, s3d_p101_m303_mixed_xz, "1,0,1 to -3,0,-3")
  CALL ops_decl_stencil(3, 20, a3d_p102_m302_mixed_xz, s3d_p102_m302_mixed_xz, "1,0,2 to -3,0,-2")
  CALL ops_decl_stencil(3, 25, a3d_p104_m300_mixed_xz, s3d_p104_m300_mixed_xz, "1,0,4 to -3,0,0")

  CALL ops_decl_stencil(3, 20, a3d_p200_m204_mixed_xz, s3d_p200_m204_mixed_xz, "2,0,0 to -2,0,-4")
  CALL ops_decl_stencil(3, 20, a3d_p201_m203_mixed_xz, s3d_p201_m203_mixed_xz, "2,0,1 to -2,0,-3")
  CALL ops_decl_stencil(3, 8, a3d_p202_m202_mixed_xz, s3d_p202_m202_mixed_xz, "2,0,2 to -2,0,-2")
  CALL ops_decl_stencil(3, 20, a3d_p203_m201_mixed_xz, s3d_p203_m201_mixed_xz, "2,0,3 to -2,0,-1")
  CALL ops_decl_stencil(3, 20, a3d_p204_m200_mixed_xz, s3d_p204_m200_mixed_xz, "2,0,4 to -2,0,0")

  CALL ops_decl_stencil(3, 13, a3d_p300_m100_mixed_xz, s3d_p300_m100_mixed_xz, "3,0,0 to -1,0,0")
  CALL ops_decl_stencil(3, 25, a3d_p300_m104_mixed_xz, s3d_p300_m104_mixed_xz, "3,0,0 to -1,0,-4")
  CALL ops_decl_stencil(3, 20, a3d_p302_m102_mixed_xz, s3d_p302_m102_mixed_xz, "3,0,2 to -1,0,-2")
  CALL ops_decl_stencil(3, 13, a3d_p303_m101_mixed_xz, s3d_p303_m101_mixed_xz, "3,0,3 to -1,0,-1")
  CALL ops_decl_stencil(3, 12, a3d_p303_m303_mixed_xz, s3d_p303_m303_mixed_xz, "3,0,3 to -3,0,-3")
  CALL ops_decl_stencil(3, 25, a3d_p304_m100_mixed_xz, s3d_p304_m100_mixed_xz, "3,0,4 to -1,0,0")

  CALL ops_decl_stencil(3, 13, a3d_p400_p004_mixed_xz, s3d_p400_p004_mixed_xz, "4,0,0 to 0,0,-4")
  CALL ops_decl_stencil(3, 25, a3d_p401_p003_mixed_xz, s3d_p401_p003_mixed_xz, "4,0,1 to 0,0,-3")
  CALL ops_decl_stencil(3, 20, a3d_p402_m002_mixed_xz, s3d_p402_m002_mixed_xz, "4,0,2 to 0,0,-2")
  CALL ops_decl_stencil(3, 25, a3d_p403_m001_mixed_xz, s3d_p403_m001_mixed_xz, "4,0,3 to 0,0,-1")
  CALL ops_decl_stencil(3, 13, a3d_p404_p000_mixed_xz, s3d_p404_p000_mixed_xz, "4,0,4 to 0,0,0")
  CALL ops_decl_stencil(3, 16, a3d_p404_m404_mixed_xz, s3d_p404_m404_mixed_xz, "4,0,4 to -4,0,-4")

  CALL ops_decl_stencil(3, 20, a3d_p505_m505_mixed_xz, s3d_p505_m505_mixed_xz, "5,0,5 to -5,0,-5")

  !------------------------------------------------------------------------------------------------------------

  CALL ops_decl_stencil(3, 13, a3d_p000_m044_mixed_yz, s3d_p000_m044_mixed_yz, "0,0,0 to 0,-4,-4")
  CALL ops_decl_stencil(3, 25, a3d_p001_m043_mixed_yz, s3d_p001_m043_mixed_yz, "0,0,1 to 0,-4,-3")
  CALL ops_decl_stencil(3, 20, a3d_p002_m042_mixed_yz, s3d_p002_m042_mixed_yz, "0,0,2 to 0,-4,-2")
  CALL ops_decl_stencil(3, 25, a3d_p003_m041_mixed_yz, s3d_p003_m041_mixed_yz, "0,0,3 to 0,-4,-1")
  CALL ops_decl_stencil(3, 13, a3d_p004_m040_mixed_yz, s3d_p004_m040_mixed_yz, "0,0,4 to 0,-4,0")

  CALL ops_decl_stencil(3, 13, a3d_p010_m030_mixed_yz, s3d_p010_m030_mixed_yz, "0,1,0 to 0,-3,0")
  CALL ops_decl_stencil(3, 25, a3d_p010_m034_mixed_yz, s3d_p010_m034_mixed_yz, "0,1,0 to 0,-3,-4")
  CALL ops_decl_stencil(3, 13, a3d_p011_m033_mixed_yz, s3d_p011_m033_mixed_yz, "0,1,1 to 0,-3,-3")
  CALL ops_decl_stencil(3, 20, a3d_p012_m032_mixed_yz, s3d_p012_m032_mixed_yz, "0,1,2 to 0,-3,-2")
  CALL ops_decl_stencil(3, 25, a3d_p014_m030_mixed_yz, s3d_p014_m030_mixed_yz, "0,1,4 to 0,-3,0")

  CALL ops_decl_stencil(3, 20, a3d_p020_m024_mixed_yz, s3d_p020_m024_mixed_yz, "0,2,0 to 0,-2,-4")
  CALL ops_decl_stencil(3, 20, a3d_p021_m023_mixed_yz, s3d_p021_m023_mixed_yz, "0,2,1 to 0,-2,-3")
  CALL ops_decl_stencil(3, 8, a3d_p022_m022_mixed_yz, s3d_p022_m022_mixed_yz, "0,2,2 to 0,-2,-2")
  CALL ops_decl_stencil(3, 20, a3d_p023_m021_mixed_yz, s3d_p023_m021_mixed_yz, "0,2,3 to 0,-2,-1")
  CALL ops_decl_stencil(3, 20, a3d_p024_m020_mixed_yz, s3d_p024_m020_mixed_yz, "0,2,4 to 0,-2,0")

  CALL ops_decl_stencil(3, 13, a3d_p030_m010_mixed_yz, s3d_p030_m010_mixed_yz, "0,3,0 to 0,-1,0")
  CALL ops_decl_stencil(3, 25, a3d_p030_m014_mixed_yz, s3d_p030_m014_mixed_yz, "0,3,0 to 0,-1,-4")
  CALL ops_decl_stencil(3, 20, a3d_p032_m012_mixed_yz, s3d_p032_m012_mixed_yz, "0,3,2 to 0,-1,-2")
  CALL ops_decl_stencil(3, 13, a3d_p033_m011_mixed_yz, s3d_p033_m011_mixed_yz, "0,3,3 to 0,-1,-1")
  CALL ops_decl_stencil(3, 12, a3d_p033_m033_mixed_yz, s3d_p033_m033_mixed_yz, "0,3,3 to 0,-3,-3")
  CALL ops_decl_stencil(3, 25, a3d_p034_m010_mixed_yz, s3d_p034_m010_mixed_yz, "0,3,4 to 0,-1,0")

  CALL ops_decl_stencil(3, 13, a3d_p040_p004_mixed_yz, s3d_p040_p004_mixed_yz, "0,4,0 to 0,0,-4")
  CALL ops_decl_stencil(3, 25, a3d_p041_p003_mixed_yz, s3d_p041_p003_mixed_yz, "0,4,1 to 0,0,-3")
  CALL ops_decl_stencil(3, 20, a3d_p042_m002_mixed_yz, s3d_p042_m002_mixed_yz, "0,4,2 to 0,0,-2")
  CALL ops_decl_stencil(3, 25, a3d_p043_m001_mixed_yz, s3d_p043_m001_mixed_yz, "0,4,3 to 0,0,-1")
  CALL ops_decl_stencil(3, 13, a3d_p044_p000_mixed_yz, s3d_p044_p000_mixed_yz, "0,4,4 to 0,0,0")
  CALL ops_decl_stencil(3, 16, a3d_p044_m044_mixed_yz, s3d_p044_m044_mixed_yz, "0,4,4 to 0,-4,-4")

  CALL ops_decl_stencil(3, 20, a3d_p055_m055_mixed_yz, s3d_p055_m055_mixed_yz, "0,5,5 to 0,-5,-5")

  !------------------------------------------------------------------------------------------------------------

  !---------------------OPS DECL HALO for periodic transfer on single processor--------------------------------

  dir_from = [1, 2, 3]
  dir_to = [1, 2, 3]

  !   X-DIRECTION : RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
  iter_size = [nhalox, nyglbl, nzglbl]
  base_from = [1, 1, 1]
  base_to = [nxglbl + 1, 1, 1]

  halo_idx = 0

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    DO ispec = 1, nspcmx
    halo_idx = halo_idx + 1
    CALL ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))
  END DO

  !   X-DIRECTION : LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
  iter_size = [nhalox, nyglbl, nzglbl]
  base_from = [nxglbl - nhalox + 1, 1, 1]
  base_to = [1 - nhalox, 1, 1]

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))

    DO ispec = 1, nspcmx
    halo_idx = halo_idx + 1
    CALL ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_x(halo_idx))
  END DO

  CALL ops_decl_halo_group(halo_idx, halos_x, halos_grp_x)

  !------------------------------------------------------------------------------------------------------------

  !   Y-DIRECTION : RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
  iter_size = [nxglbl + 2 * nhalox, nhaloy, nzglbl]
  base_from = [1 - nhalox, 1, 1]
  base_to = [1 - nhalox, nyglbl + 1, 1]

  halo_idx = 0

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    DO ispec = 1, nspcmx
    halo_idx = halo_idx + 1
    CALL ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))
  END DO

  !   Y-DIRECTION : LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
  iter_size = [nxglbl + 2 * nhalox, nhaloy, nzglbl]
  base_from = [1 - nhalox, nyglbl - nhaloy + 1, 1]
  base_to = [1 - nhalox, 1 - nhaloy, 1]

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))

    DO ispec = 1, nspcmx
    halo_idx = halo_idx + 1
    CALL ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_y(halo_idx))
  END DO

  CALL ops_decl_halo_group(halo_idx, halos_y, halos_grp_y)

  !------------------------------------------------------------------------------------------------------------

  !   Z-DIRECTION : RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
  iter_size = [nxglbl + 2 * nhalox, nyglbl + 2 * nhaloy, nhaloz]
  base_from = [1 - nhalox, 1 - nhaloy, 1]
  base_to = [1 - nhalox, 1 - nhaloy, nzglbl + 1]

  halo_idx = 0

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    DO ispec = 1, nspcmx
    halo_idx = halo_idx + 1
    CALL ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))
  END DO

  !   Z-DIRECTION : LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
  iter_size = [nxglbl + 2 * nhalox, nyglbl + 2 * nhaloy, nhaloz]
  base_from = [1 - nhalox, 1 - nhaloy, nzglbl - nhaloz + 1]
  base_to = [1 - nhalox, 1 - nhaloy, 1 - nhaloz]

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_drhs, d_drhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_urhs, d_urhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_vrhs, d_vrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_wrhs, d_wrhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

  halo_idx = halo_idx + 1
  CALL ops_decl_halo(d_erhs, d_erhs, iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))

    DO ispec = 1, nspcmx
    halo_idx = halo_idx + 1
    CALL ops_decl_halo(d_yrhs(ispec), d_yrhs(ispec), iter_size, base_from, base_to, dir_from, dir_to, halos_z(halo_idx))
  END DO

  CALL ops_decl_halo_group(halo_idx, halos_z, halos_grp_z)

  !------------------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------------------
  CALL ops_partition(" ")
  !------------------------------------------------------------------------------------------------------------

  !---------------------------------First touch - OPS DATS without halos---------------------------------------
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_divm, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wd1x, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_pd1x, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_td1x, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wd1y, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_pd1y, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_td1y, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wd1z, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_pd1z, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_td1z, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wd2x, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_pd2x, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_td2x, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wd2y, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_pd2y, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_td2y, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wd2z, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_pd2z, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_td2z, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ufxl, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_vfxl, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wfxl, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  !-------------------------------First touch - OPS DATS without halos-----------------------------------------
  rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wmomix, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_difmix, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  !------------------------------------First touch - OPS DATS Boundary YZ--------------------------------------
  rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]

  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcltxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strexl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strgxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dtdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dddtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ovgmxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sydtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t1bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t2bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t3bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t4bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t51bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t52bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

  rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]

  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl5xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcltxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strexr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strgxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strrxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dudtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dvdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dwdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dtdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dddtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ovgmxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sydtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t1bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t2bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t3bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t4bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t51bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t52bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

  !--------------------------First touch - OPS DATS Boundary XZ------------------------------------------------
  rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]

  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcltyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_streyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strgyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dudtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dvdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dwdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dtdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dddtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ovgmyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sydtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))


  rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]

  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl5yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcltyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_streyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strgyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strryr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dudtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dvdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dwdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dtdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dddtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ovgmyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sydtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

  !------------------------------First touch - OPS DATS Boundary XY--------------------------------------------
  rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]

  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcltzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strezl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strgzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dudtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dvdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dwdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dtdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dddtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ovgmzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sydtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))


  rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]

  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcl5zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bcltzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strezr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strgzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strrzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dudtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dvdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dwdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dtdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dddtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ovgmzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sydtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

  !--------------------------------First touch - MULTI-DIM DAT-------------------------------------------------
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  DO ispec = 1, nspcmx
    CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rrte(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))

  END DO

  rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
  DO iindex = 1, nintmx
    CALL set_zero_kernel_int_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
  END DO

    DO ispec = 1, nspcmx
    CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ctrans(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
  END DO

  rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
  DO ispec = 1, nspcmx
    CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bclyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dydtxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ratexl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  END DO

  rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
  DO ispec = 1, nspcmx
    CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bclyxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_stryxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dydtxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ratexr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
  END DO

  rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
  DO ispec = 1, nspcmx
    CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bclyyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_stryyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dydtyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rateyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  END DO

  rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
  DO ispec = 1, nspcmx
    CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bclyyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_stryyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dydtyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rateyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_ydir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
  END DO

  rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
  DO ispec = 1, nspcmx
    CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bclyzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_stryzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dydtzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ratezl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  END DO

  rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
  DO ispec = 1, nspcmx
    CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_bclyzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_stryzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_dydtzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ratezr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_zdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
  END DO

  !------------------------------------------------------------------------------------------------------------
#ifdef OPS_WITH_CUDAFOR

  ncofmx_opsconstant = ncofmx
  ntinmx_opsconstant = ntinmx
  nspcmx_opsconstant = nspcmx
  nssmax_opsconstant = nssmax
  nstpmx_opsconstant = nstpmx
  ndcfmx_opsconstant = ndcfmx
  nvcfmx_opsconstant = nvcfmx
  nccfmx_opsconstant = nccfmx
  nrkmax_opsconstant = nrkmax
  ncbcsz_opsconstant = ncbcsz
  nbcprr_opsconstant = nbcprr
  nspimx_opsconstant = nspimx
  ntbase_opsconstant = ntbase
  nintmx_opsconstant = nintmx
  nctmax_opsconstant = nctmax
  nctmm1_opsconstant = nctmm1
  nrsmax_opsconstant = nrsmax
  nbcpri_opsconstant = nbcpri
  ncfrmx_opsconstant = ncfrmx
#endif
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ncofmx", 1, "integer(kind=4)", ncofmx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ntinmx", 1, "integer(kind=4)", ntinmx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nspcmx", 1, "integer(kind=4)", nspcmx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nssmax", 1, "integer(kind=4)", nssmax)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nstpmx", 1, "integer(kind=4)", nstpmx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ndcfmx", 1, "integer(kind=4)", ndcfmx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nvcfmx", 1, "integer(kind=4)", nvcfmx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nccfmx", 1, "integer(kind=4)", nccfmx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nrkmax", 1, "integer(kind=4)", nrkmax)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ncbcsz", 1, "integer(kind=4)", ncbcsz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nbcprr", 1, "integer(kind=4)", nbcprr)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nspimx", 1, "integer(kind=4)", nspimx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ntbase", 1, "integer(kind=4)", ntbase)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nintmx", 1, "integer(kind=4)", nintmx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nctmax", 1, "integer(kind=4)", nctmax)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nctmm1", 1, "integer(kind=4)", nctmm1)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nrsmax", 1, "integer(kind=4)", nrsmax)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("nbcpri", 1, "integer(kind=4)", nbcpri)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ncfrmx", 1, "integer(kind=4)", ncfrmx)
#endif
#ifdef OPS_WITH_OMPOFFLOADFOR
#endif

END SUBROUTINE ops_data_init