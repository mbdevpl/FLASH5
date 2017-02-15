!!***h* source/physics/Hydro/HydroMain/unsplit/UHD.h
!!
!! This is the internal header file for the Unsplit MHD/Hydro units.
!!
!!***

#include "Flash.h"
#include "constants.h"

#define DIR_X 1
#define DIR_Y 2
#define DIR_Z 3




!!--------------------------------------------------------------!!
!! FOR UNSPLIT HYDRO IMPLEMENTATION ----------------------------!!
!!--------------------------------------------------------------!!
 
#define ROE   1
#define LLF   2
#define HLL   3

!! PRIMITIVE VARIABLES FOR PURE HYDRO
#define HY_DENS 1
#define HY_VELX 2
#define HY_VELY 3
#define HY_VELZ 4
#define HY_PRES 5
#define HY_GAMC 6
#define HY_GAME 7
#define HY_EINT 8

#define HY_GRAV 9
#define HY_TEMP 10
#define HY_ABAR 11
#define HY_ZBAR 12

!! CONSERVATIVE VARIABLES FOR PURE HYDRO
#define HY_XMOM 2
#define HY_YMOM 3
#define HY_ZMOM 4
#define HY_ENER 5

!! FLUX VARIABLES FOR PURE HYDRO
#define HY_DENS_FLUX 1
#define HY_XMOM_FLUX 2
#define HY_YMOM_FLUX 3
#define HY_ZMOM_FLUX 4
#define HY_ENER_FLUX 5
#define HY_EINT_FLUX 6
#define HY_PRES_FLUX 7

!! WAVE STRUCTURES FOR PURE HYDRO
#define HY_FASTLEFT 1
#define HY_SLOWLEFT 2
#define HY_ENTROPY  3
#define HY_SLOWRGHT 4
#define HY_FASTRGHT 5

!! EXTRA PARAMETERS FOR PURE HYDRO
#define HY_VARINUM HY_PRES

#define HY_SCRATCH_NUM 8



!!--------------------------------------------------------------!!
!! END OF HYDRO DEFINITIONS ------------------------------------!!
!!--------------------------------------------------------------!!


!! DEFINE TOTAL NUMBERS OF VARIABLES AND WAVES
#define HY_WAVENUM  HY_FASTRGHT
#define HY_VARINUM2 HY_VARINUM+2
#define HY_VARINUM3 HY_VARINUM+3
#define HY_VARINUM4 HY_VARINUM+4
#define HY_VARINUM7 HY_VARINUM+7


#define HY_VARINUMMAX HY_VARINUM4

#define HY_END_VARS HY_EINT
#define HY_END_FLUX HY_PRES_FLUX
