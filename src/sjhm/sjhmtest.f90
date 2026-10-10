! Module to test HSHMv4 with a set of standard problems in DISORT

module Disotest_vars

  use Mod_Gdrtes, only: Kreal, Kreal2, Kint

  implicit none 

  LOGICAL :: DOPROB(17) = .true.
  !DATA DOPROB / 17*.TRUE. /
  LOGICAL  USRANG, USRTAU, ONLYFL, PRNT(5), &
           PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE
  data PRNT / .TRUE., .FALSE., .FALSE., .FALSE., .TRUE. /
  INTEGER(Kint)::  IBCND, NMOM, NLYR, NUMU, NSTR, NPHI, NTAU
  INTEGER(Kint)::  NUMU_O
  LOGICAL::  DEBUG
  real(Kreal) :: ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, &
                  PHI0, TEMIS, TTEMP, WVNMLO, WVNMHI, UMU0 
  real(Kreal),parameter :: EARTH_RADIUS = 6371.0     

  INTEGER(Kint)::       BRDF_TYPE
  REAL(Kreal)::         BRDF_ARG(4)
  LOGICAL::             DO_SHADOW
  REAL(Kreal)::         WIND_SPD, REFRAC_INDEX
  REAL(Kreal)::         B0, HH, W
  REAL(Kreal)::         K_VOL, K_ISO, K_GEO
  REAL(Kreal)::         RHO_0, KAPPA, G, H0
  REAL(Kreal)::         FLUX_UP, DFDTAU
  INTEGER(Kint)::       NMUG
  REAL(Kreal), external::         BDREF

  real(Kreal),dimension(:),allocatable:: DTAUC, PHI, SSALB, TEMPER, UMU, UTAU
  real(Kreal),dimension(:,:),allocatable   :: PMOM
  real(Kreal),dimension(:,:,:),allocatable :: RHOQ, RHOU 
  real(Kreal),dimension(:),allocatable     :: EMUST, BEMST
  real(Kreal),dimension(:,:),allocatable   :: RHO_ACCURATE
  real(Kreal),dimension(:),allocatable:: RFLDIR, RFLDN, FLUP, DFDT, UAVG, &
                                       & ALBMED, TRNMED
  real(Kreal),dimension(:,:,:),allocatable :: UU
  real(Kreal),dimension(:),allocatable     :: H_LYR

contains

!-------------------------------------------------------------------------------

subroutine allocate_disort_allocatable_arrays(NLYR, NMOM, NSTR, NUMU, NPHI, NTAU)
implicit none
integer,intent(in) :: NLYR, NMOM, NSTR, NUMU, NPHI, NTAU 

allocate( DTAUC(NLYR), SSALB(NLYR), PMOM(0:NMOM, NLYR), TEMPER(0:NLYR), &
&         UTAU(NTAU), UMU(NUMU), PHI(NPHI), H_LYR(0:NLYR) )  
allocate( RHOQ(NSTR/2, 0:NSTR/2, 0:(NSTR-1)), &
&         RHOU(NUMU, 0:NSTR/2, 0:(NSTR-1)), &
&         EMUST(NUMU), BEMST(NSTR/2), RHO_ACCURATE(NUMU, NPHI) )
allocate( RFLDIR(NTAU), RFLDN(NTAU), FLUP(NTAU), DFDT(NTAU), UAVG(NTAU), &
&         ALBMED(NUMU), TRNMED(NUMU), UU(NUMU, NTAU, NPHI) )

DTAUC = 0.0; SSALB = 0.0; PMOM = 0.0; TEMPER = 0.0; 
UTAU = 0.0; UMU = 0.0; PHI = 0.0;

H_LYR = 0.0; RHOQ = 0.0; RHOU = 0.0; EMUST = 0.0; BEMST = 0.0; 
RHO_ACCURATE = 0.0;

RFLDIR = 0.0; RFLDN = 0.0; FLUP = 0.0; DFDT = 0.0; UAVG = 0.0; UU = 0.0;
ALBMED = 0.0; TRNMED = 0.0; 

end subroutine allocate_disort_allocatable_arrays

!-------------------------------------------------------------------------------

subroutine deallocate_disort_allocatable_arrays()

deallocate( DTAUC, SSALB, PMOM, TEMPER, UTAU, UMU, PHI, H_LYR )  
deallocate( RHOQ, RHOU, EMUST, BEMST, RHO_ACCURATE )                
deallocate( RFLDIR, RFLDN, FLUP, DFDT, UAVG, ALBMED, TRNMED, UU )  

end subroutine deallocate_disort_allocatable_arrays

!-------------------------------------------------------------------------------

end module Disotest_vars

!===============================================================================
module Disotestaux
  use Mod_Gdrtes, only: Kreal, Kint

  integer, parameter:: MAXCLY = 13, MAXMOM = 999, MAXPHI = 5, MAXULV = 5, &
     &            MAXCMU = 48, MAXUMU = 90
  integer, parameter:: MXPROB = 16, MXCASE = 8, MXTAU = MAXULV, &
     &            MXMU = MAXUMU, MXPHI = MAXPHI
!     ..
!     .. Local Scalars ..

!  INTEGER::   I, J, K
!     ..
!     .. Common blocks ..

  REAL(Kreal)::        TSTDFD( MXTAU, MXCASE, MXPROB ), &
     &          TSTFDN( MXTAU, MXCASE, MXPROB ), &
     &          TSTFIR( MXTAU, MXCASE, MXPROB ), &
     &          TSTFUP( MXTAU, MXCASE, MXPROB ), &
     &          TSTUU( MXTAU, MXMU, MXPHI, MXCASE, MXPROB )
! ********************* Test Case 1A *********************************
      DATA (TSTFIR(I,1,1), I = 1, 2) / 3.14159E+00, 2.29844E+00 /
      DATA (TSTFDN(I,1,1), I = 1, 2) / 0.0, 7.94108E-02 /
      DATA (TSTFUP(I,1,1), I = 1, 2) / 7.99451E-02, 0.0 /
      DATA (TSTDFD(I,1,1), I = 1, 2) / 2.54067E+01, 1.86531E+01 /
      DATA ((TSTUU(I,J,1,1,1), J = 1, 6), I = 1, 2) &
     &  / 3*0.0, 1.17771E-01, 2.64170E-02, 1.34041E-02, &
     &    1.33826E-02, 2.63324E-02, 1.15898E-01, 3*0.0 /

! ********************* Test Case 1B *********************************

      DATA (TSTFIR(I,2,1), I = 1, 2) / 3.14159E+00, 2.29844E+00 /
      DATA (TSTFDN(I,2,1), I = 1, 2) / 0.0, 4.20233E-01 /
      DATA (TSTFUP(I,2,1), I = 1, 2) / 4.22922E-01, 0.0 /
      DATA (TSTDFD(I,2,1), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,1), J = 1, 6), I = 1, 2) &
     &  / 3*0.0, 6.22884E-01, 1.39763E-01, 7.09192E-02, &
     &    7.08109E-02, 1.39337E-01, 6.13458E-01, 3*0.0 /

! ********************* Test Case 1C *********************************

      DATA (TSTFIR(I,3,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,3,1), I = 1, 2) / 3.14159E+00, 3.04897E+00 /
      DATA (TSTFUP(I,3,1), I = 1, 2) / 9.06556E-02, 0.0 /
      DATA (TSTDFD(I,3,1), I = 1, 2) / 6.66870E-02, 5.88936E-02 /
      DATA ((TSTUU(I,J,1,3,1), J = 1, 6), I = 1, 2) &
     &  / 3*1.0, 1.33177E-01, 2.99879E-02, 1.52233E-02, &
     &    9.84447E-01, 9.69363E-01, 8.63946E-01, 3*0.0 /

! ********************* Test Case 1D *********************************

      DATA (TSTFIR(I,4,1), I = 1, 2) / 3.14159E+00, 0.00000E+00 /
      DATA (TSTFDN(I,4,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFUP(I,4,1), I = 1, 2) / 2.59686E-01, 0.0 /
      DATA (TSTDFD(I,4,1), I = 1, 2) / 2.57766E+01, 0.0 /
      DATA ((TSTUU(I,J,1,4,1), J = 1, 6), I = 1, 2) &
     &  / 3*0.0, 2.62972E-01, 9.06967E-02, 5.02853E-02, &
     &    1.22980E-15, 1.30698E-17, 6.88840E-18, 3*0.0 /

! ********************* Test Case 1E *********************************

      DATA (TSTFIR(I,5,1), I = 1, 2) / 3.14159E+00, 0.00000E+00 /
      DATA (TSTFDN(I,5,1), I = 1, 2) / 0.0, 6.76954E-02 /
      DATA (TSTFUP(I,5,1), I = 1, 2) / 3.07390E+00, 0.0 /
      DATA (TSTDFD(I,5,1), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,5,1), J = 1, 6), I = 1, 2) &
     &  / 3*0.0, 1.93321E+00, 1.02732E+00, 7.97199E-01, &
     &    2.71316E-02, 1.87805E-02, 1.16385E-02, 3*0.0 /

! ********************* Test Case 1F *********************************

      DATA (TSTFIR(I,6,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,6,1), I = 1, 2) / 3.14159E+00, 4.60048E-03 /
      DATA (TSTFUP(I,6,1), I = 1, 2) / 2.49618E+00, 0.0 /
      DATA (TSTDFD(I,6,1), I = 1, 2) / 1.14239E-01, 7.93633E-05 /
      DATA ((TSTUU(I,J,1,6,1), J = 1, 6), I = 1, 2) &
     &  / 3*1.0, 8.77510E-01, 8.15136E-01, 7.52715E-01, &
     &    1.86840E-03, 1.26492E-03, 7.79280E-04, 3*0.0 /


! ********************* Test Case 2A *********************************

      DATA (TSTFIR(I,1,2), I = 1, 2) / 2.52716E-01, 2.10311E-02 /
      DATA (TSTFDN(I,1,2), I = 1, 2) / 0.0, 4.41791E-02 /
      DATA (TSTFUP(I,1,2), I = 1, 2) / 5.35063E-02, 0.0 /
      DATA (TSTDFD(I,1,2), I = 1, 2) / 1.66570E+00, 1.89848E-01 /
      DATA ((TSTUU(I,J,1,1,2), J = 1, 6), I = 1, 2) &
     &  / 3*0.0, 1.61796E-01, 2.11501E-02, 7.86713E-03, &
     &    7.71897E-03, 2.00778E-02, 2.57685E-02, 3*0.0 /

! ********************* Test Case 2B *********************************

      DATA (TSTFIR(I,2,2), I = 1, 2) / 2.52716E-01, 2.10311E-02 /
      DATA (TSTFDN(I,2,2), I = 1, 2) / 0.0, 1.06123E-01 /
      DATA (TSTFUP(I,2,2), I = 1, 2) / 1.25561E-01, 0.0 /
      DATA (TSTDFD(I,2,2), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,2), J = 1, 6), I = 1, 2) &
     &  / 3*0.0, 3.47678E-01, 4.87120E-02, 1.89387E-02, &
     &    1.86027E-02, 4.64061E-02, 6.77603E-02, 3*0.0 /

! ********************* Test Case 2C *********************************

      DATA (TSTFIR(I,3,2), I = 1, 2) / 2.52716E-01, 2.56077E-28 /
      DATA (TSTFDN(I,3,2), I = 1, 2) / 0.0, 2.51683E-04 /
      DATA (TSTFUP(I,3,2), I = 1, 2) / 6.24730E-02, 0.0 /
      DATA (TSTDFD(I,3,2), I = 1, 2) / 1.67462E+00, 1.75464E-04 /
      DATA ((TSTUU(I,J,1,3,2), J = 1, 6), I = 1, 2) &
     &  / 3*0.0, 1.62566E-01, 2.45786E-02, 1.01498E-02, &
     &    1.70004E-04, 3.97168E-05, 1.32472E-05, 3*0.0 /

! ********************* Test Case 2D *********************************

      DATA (TSTFIR(I,4,2), I = 1, 2) / 2.52716E-01, 0.0 /
      DATA (TSTFDN(I,4,2), I = 1, 2) / 0.0, 2.68008E-02 /
      DATA (TSTFUP(I,4,2), I = 1, 2) / 2.25915E-01, 0.0 /
      DATA (TSTDFD(I,4,2), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,4,2), J = 1, 6), I = 1, 2) &
     &  / 3*0.0, 3.64010E-01, 8.26993E-02, 4.92370E-02, &
     &    1.05950E-02, 7.69337E-03, 3.79276E-03, 3*0.0 /

! ********************* Test Case 3A *********************************

      DATA (TSTFIR(I,1,3), I = 1, 2) / 3.14159E+00, 1.15573E+00 /
      DATA (TSTFDN(I,1,3), I = 1, 2) / 0.0, 1.73849E+00 /
      DATA (TSTFUP(I,1,3), I = 1, 2) / 2.47374E-01, 0.0 /
      DATA (TSTDFD(I,1,3), I = 1, 2) / 0.0, 0.0 /
      DATA ((TSTUU(I,J,1,1,3), J = 1, 6), I = 1, 2) / &
     &    3*0.0, 1.51159E-01, 1.01103E-01, 3.95460E-02, &
     &    3.05855E+00, 2.66648E-01, 2.13750E-01, 3*0.0 /

! ********************* Test Case 3B *********************************

      DATA (TSTFIR(I,2,3), I = 1, 2) / 3.14159E+00, 1.05389E-03 /
      DATA (TSTFDN(I,2,3), I = 1, 2) / 0.0, 1.54958E+00 /
      DATA (TSTFUP(I,2,3), I = 1, 2) / 1.59096E+00, 0.0 /
      DATA (TSTDFD(I,2,3), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,3), J = 1, 6), I = 1, 2) / &
     &    3*0.0, 3.79740E-01, 5.19598E-01, 4.93302E-01, &
     &    6.69581E-01, 4.22350E-01, 2.36362E-01, 3*0.0 /

! ********************* Test Case 4A *********************************

      DATA (TSTFIR(I,1,4), I = 1, 3) &
     &  / 3.14159E+00, 1.90547E+00, 1.15573E+00 /
      DATA (TSTFDN(I,1,4), I = 1, 3) &
     &  / 0.0, 1.17401E+00, 1.81264E+00 /
      DATA (TSTFUP(I,1,4), I = 1, 3) &
     &  / 1.73223E-01, 1.11113E-01, 0.0 /
      DATA (TSTDFD(I,1,4), I = 1, 3) / 3*0.0 /
      DATA ((TSTUU(I,J,1,1,4), J = 1, 6), I = 1, 3) &
     &  / 3*0.0, 9.26837E-02, &
     &    6.59569E-02, 3.64755E-02, 2.51608E+00, 1.19287E-01, &
     &    1.34962E-01, 1.23887E-01, 4.02058E-02, 1.77746E-02, &
     &    3.37302E+00, 2.19835E-01, 1.56893E-01, 3*0.0 /

! ********************* Test Case 4B *********************************

      DATA (TSTFIR(I,2,4), I = 1, 3) &
     &  / 3.14159E+00, 1.90547E+00, 1.15573E+00 /
      DATA (TSTFDN(I,2,4), I = 1, 3) &
     &  / 0.0, 1.01517E+00, 1.51554E+00 /
      DATA (TSTFUP(I,2,4), I = 1, 3) &
     &  / 1.23665E-01, 7.88690E-02, 0.0 /
      DATA (TSTDFD(I,2,4), I = 1, 3) &
     &  / 3.43724E-01, 3.52390E-01, 3.19450E-01 /
      DATA ((TSTUU(I,J,1,2,4), J = 1, 6), I = 1, 3) &
     &  / 3*0.0, 6.53056E-02, &
     &    4.55144E-02, 2.82693E-02, 2.24258E+00, 9.66049E-02, &
     &    9.61335E-02, 8.43278E-02, 2.79473E-02, 1.38835E-02, &
     &    2.97057E+00, 1.67698E-01, 1.08115E-01, 3*0.0 /

! ********************* Test Case 4C *********************************

      DATA (TSTFIR(I,3,4), I = 1, 3) &
     &  / 1.57080E+00, 5.77864E-01, 2.12584E-01 /
      DATA (TSTFDN(I,3,4), I = 1, 3) &
     &  / 0.0, 7.02764E-01, 8.03294E-01 /
      DATA (TSTFUP(I,3,4), I = 1, 3) &
     &  / 2.25487E-01, 1.23848E-01, 0.0 /
      DATA (TSTDFD(I,3,4), I = 1, 3) &
     &  / 3.85003E-01, 3.37317E-01, 2.16403E-01 /
      DATA (((TSTUU(I,J,K,3,4), J = 1, 6), I = 1, 3), K = 1, 3) &
     &  / 3*0.0, 8.70812E-01, &
     &    2.24960E-01, 2.27572E-02, 4.77016E-02, 3.02631E+00, &
     &    1.41195E+00, 6.97692E-01, 1.09130E-01, 9.32861E-03, &
     &    8.38488E-02, 2.70538E+00, 8.76523E-01, 6*0.0, &
     &    8.88117E-02, 5.77411E-02, 2.27572E-02, &
     &    4.77016E-02, 5.80971E-02, 1.04502E-01, 9.16071E-02, &
     &    2.95842E-02, 9.32861E-03, 8.38488E-02, 9.42187E-02, &
     &    8.95457E-02, 6*0.0, 6.98247E-02, &
     &    5.02877E-02, 2.27572E-02, 4.77016E-02, 2.58544E-02, &
     &    6.25954E-02, 5.91273E-02, 2.47702E-02, 9.32861E-03, &
     &    8.38488E-02, 3.99383E-02, 4.67155E-02, 3*0.0 /

! ********************* Test Case 5A *********************************

      DATA (TSTFIR(I,1,5), I = 1, 3) &
     &  / 3.14159E+00, 3.97856E-14, 5.03852E-28 /
      DATA (TSTFDN(I,1,5), I = 1, 3) &
     &  / 0.0, 2.24768E+00, 4.79851E-01 /
      DATA (TSTFUP(I,1,5), I = 1, 3) &
     &  / 2.66174E+00, 1.76783E+00, 0.0 /
      DATA (TSTDFD(I,1,5), I = 1, 3) / 3*0.0 /
      DATA ((TSTUU(I,J,1,1,5), J = 1, 6), I = 1, 3) &
     &  / 3*0.0, 4.58927E-01, &
     &    7.72983E-01, 1.07196E+00, 7.53662E-01, 6.96362E-01, &
     &    6.50541E-01, 6.27631E-01, 5.81809E-01, 5.24532E-01, &
     &    1.95230E-01, 1.31990E-01, 7.20655E-02, 3*0.0 /

! ********************* Test Case 5B *********************************

      DATA (TSTFIR(I,2,5), I = 1, 3) &
     &  / 1.28058E-01, 8.67322E-06, 4.47729E-21 /
      DATA (TSTFDN(I,2,5), I = 1, 3) &
     &  / 1.74767E+00, 2.33975E-01, 6.38345E-05 /
      DATA (TSTFUP(I,2,5), I = 1, 3) &
     &  / 2.70485E-01, 3.74252E-02, 1.02904E-05 /
      DATA (TSTDFD(I,2,5), I = 1, 3) &
     &  / 3.10129E-01, 4.52671E-02, 1.25021E-05 /
      DATA ((TSTUU(I,J,1,2,5), J = 1, 6), I = 1, 3) &
     &  / 6.79623E+01, 2.21027E-01, 1.36619E-01, 1.14084E-01, &
     &    8.73870E-02, 8.81626E-02, 2.05706E-01, 4.92736E-02, &
     &    2.65449E-02, 2.02154E-02, 1.29661E-02, 9.51334E-03, &
     &    3.41286E-05, 1.39916E-05, 7.47039E-06, 5.65602E-06, &
     &    3.58245E-06, 2.57858E-06 /

! ********************* Test Case 6A *********************************

      DATA (TSTFIR(I,1,6), I = 1, 2) / 2*100.0 /
      DATA (TSTFDN(I,1,6), I = 1, 2) / 2*0.0 /
      DATA (TSTFUP(I,1,6), I = 1, 2) / 2*0.0 /
      DATA (TSTDFD(I,1,6), I = 1, 2) / 2*200.0 /
      DATA ((TSTUU(I,J,1,1,6), J = 1, 4), I = 1, 2) / 8*0.0 /

! ********************* Test Case 6B *********************************

      DATA (TSTFIR(I,2,6), I = 1, 3) &
     &  / 1.000000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,2,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,2,6), I = 1, 3) / 3*0.0 /
      DATA (TSTDFD(I,2,6), I = 1, 3) &
     &  / 2.00000E+02, 7.35759E+01, 2.70671E+01 /
      DATA ((TSTUU(I,J,1,2,6), J = 1, 4), I = 1, 3) / 12*0.0 /

! ********************* Test Case 6C *********************************

      DATA (TSTFIR(I,3,6), I = 1, 3) &
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,3,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,3,6), I = 1, 3) &
     &  / 1.48450E+00, 2.99914E+00, 6.76676E+00 /
      DATA (TSTDFD(I,3,6), I = 1, 3) &
     &  / 2.02010E+02, 7.79962E+01, 4.06006E+01 /
      DATA ((TSTUU(I,J,1,3,6), J = 1, 4), I = 1, 3) &
     &  / 2*0.0, 9.77882E-05, 7.92386E-01, &
     &    2*0.0, 1.45131E-02, 1.30642E+00, &
     &    2*0.0, 2.15393E+00, 2.15393E+00 /

! ********************* Test Case 6D *********************************

      DATA (TSTFIR(I,4,6), I = 1, 3) &
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,4,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,4,6), I = 1, 3) &
     &  / 6.70783E-01, 1.39084E+00, 3.31655E+00 /
      DATA (TSTDFD(I,4,6), I = 1, 3) &
     &  / 2.00936E+02, 7.57187E+01, 3.45317E+01 /
      DATA ((TSTUU(I,J,1,4,6), J = 1, 4), I = 1, 3) &
     &  / 2*0.0, 6.80068E-05, 3.15441E-01, &
     &    2*0.0, 1.00931E-02, 5.20074E-01, &
     &    2*0.0, 1.49795E+00, 8.57458E-01 /

! ********************* Test Case 6E *********************************

      DATA (TSTFIR(I,5,6), I = 1, 3) &
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,5,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,5,6), I = 1, 3) &
!     &  / 7.95458E+01, 1.59902E+02, 3.56410E+02 /
     &  / 9.44651E+01, 1.90616E+02, 4.28806E+02 /
      DATA (TSTDFD(I,5,6), I = 1, 3) &
!     &  / 3.07079E+02, 3.07108E+02, 7.17467E+02 /
     &  / 3.27710E+02, 3.53897E+02, 8.78110E+02 /
      DATA ((TSTUU(I,J,1,5,6), J = 1, 4), I = 1, 3) &
!     &  / 2*0.0, 4.53789E-03, 4.33773E+01,
!     &    2*0.0, 6.73483E-01, 7.15170E+01,
!     &    2*0.0, 9.99537E+01, 1.17912E+02 /
     &  / 2*0.0, 6.01546E-03, 5.06858E+01, &
     &    2*0.0, 8.92774E-01, 8.35668E+01, &
     &    2*0.0, 1.32499E+02, 1.37778E+02 /

! ********************* Test Case 6F *********************************

      DATA (TSTFIR(I,6,6), I = 1, 3) &
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,6,6), I = 1, 3) &
     &  / 3.21497E+02, 1.42493E+02, 7.05305E+01 /
      DATA (TSTFUP(I,6,6), I = 1, 3) &
!    &  / 8.27917E+01, 1.66532E+02, 3.71743E+02 /
     &  / 9.77110E+01, 1.97247E+02, 4.41139E+02 /
      DATA (TSTDFD(I,6,6), I = 1, 3) &
!    &  / 9.54523E+02, 5.27085E+02, 8.45341E+02 /
     &  / 9.75154E+02, 5.73874E+02, 1.00598E+03 /
      DATA ((TSTUU(I,J,1,6,6), J = 1, 4), I = 1, 3) &
!    &  / 1.02336E+02, 1.02336E+02, 4.80531E-03, 4.50168E+01,
!    &    6.20697E+01, 6.89532E-01, 7.13172E-01, 7.42191E+01,
!    &    3.76472E+01, 4.64603E-03, 1.05844E+02, 1.22368E+02 /
     &  / 1.02336E+02, 1.02336E+02, 6.28289E-03, 5.23253E+01, &
     &    6.20697E+01, 6.89532E-01, 9.32463E-01, 8.62699E+01, &
     &    3.76472E+01, 4.64603E-03, 1.38390E+02, 1.42235E+02 /

! ********************* Test Case 6G *********************************

      DATA (TSTFIR(I,7,6), I = 1, 3) &
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,7,6), I = 1, 3) &
     &  / 3.21497E+02, 3.04775E+02, 3.63632E+02 /
      DATA (TSTFUP(I,7,6), I = 1, 3) &
!    &  / 3.35292E+02, 4.12540E+02, 4.41125E+02 /
     &  / 3.50211E+02, 4.43254E+02, 5.13521E+02 /
      DATA (TSTDFD(I,7,6), I = 1, 3) &
!    &  / 5.80394E+02, 1.27117E+02, -1.68003E+02 /
     &  / 6.01025E+02, 1.73906E+02, -7.36023E+00 /
      DATA ((TSTUU(I,J,1,7,6), J = 1, 4), I = 1, 3) &
!    &  / 1.02336E+02, 1.02336E+02, 7.80733E+01, 1.16430E+02,
!    &    9.78748E+01, 1.01048E+02, 1.15819E+02, 1.34966E+02,
!    &    1.10061E+02, 1.38631E+02, 1.38695E+02, 1.40974E+02 /
     &  / 1.02336E+02, 1.02336E+02, 7.80748E+01, 1.23739E+02, &
     &    9.78748E+01, 1.01048E+02, 1.16039E+02, 1.47015E+02, &
     &    1.10061E+02, 1.38631E+02, 1.71240E+02, 1.60840E+02 /

! ********************* Test Case 6H *********************************

      DATA (TSTFIR(I,8,6), I = 1, 3) &
     &  / 1.00000E+02, 1.35335E+01, 2.06115E-07 /
      DATA ( TSTFDN(I,8,6), I = 1, 3) &
     &  / 3.21497E+02, 2.55455E+02, 4.43444E+02 /
      DATA (TSTFUP(I,8,6), I = 1, 3) &
!     &  / 2.37350E+02, 2.61130E+02, 4.55861E+02 /
     &  / 2.37351E+02, 2.61131E+02, 5.28258E+02 /
      DATA (TSTDFD(I,8,6), I = 1, 3) &
!     &  / 4.23780E+02, 6.19828E+01, -3.11658E+01 /
     &  / 4.23780E+02, 6.19842E+01, 1.29477E+02  /
      DATA ((TSTUU(I,J,1,8,6), J = 1, 4), I = 1, 3) &
!     &  / 1.02336E+02, 1.02336E+02, 7.12616E+01, 7.80736E+01,
!     &    8.49992E+01, 7.73186E+01, 7.88310E+01, 8.56423E+01,
!     &    1.38631E+02, 1.45441E+02, 1.44792E+02, 1.45163E+02 /
     &  / 1.02336E+02, 1.02336E+02, 7.12616E+01, 7.80745E+01, &
     &    8.49992E+01, 7.73186E+01, 7.88310E+01, 8.56448E+01, &
     &    1.38631E+02, 1.45441E+02, 1.77337E+02, 1.65030E+02 /

! ********************* Test Case 7A *********************************

      DATA (TSTFIR(I,1,7), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,1,7), I = 1, 2) / 0.0, 1.21204E+02 /
      DATA (TSTFUP(I,1,7), I = 1, 2) / 8.62936E+01,  0.0 /
      DATA (TSTDFD(I,1,7), I = 1, 2) /-5.13731E+01,-5.41036E+02 /

! ********************* Test Case 7B *********************************

      DATA (TSTFIR(I,2,7), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,2,7), I = 1, 2) / 0.0, 2.07786E-05 /
      DATA (TSTFUP(I,2,7), I = 1, 2) / 1.10949E-06,  0.0 /
      DATA (TSTDFD(I,2,7), I = 1, 2) / 8.23219E-08, -5.06461E-06 /
      DATA ((TSTUU(I,J,1,2,7), J = 1, 2), I = 1, 2) / &
     &    0.00000E+00, 4.65744E-07, 7.52311E-06, 0.00000E+00 /

! ********************* Test Case 7C *********************************

      DATA (TSTFIR(I,3,7), I = 1, 3) &
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,3,7), I = 1, 3) &
     &  / 3.19830E+02, 3.54099E+02, 3.01334E+02 /
      DATA (TSTFUP(I,3,7), I = 1, 3) &
     &  / 4.29572E+02, 4.47018E+02, 5.94576E+02 /
      DATA (TSTDFD(I,3,7), I = 1, 3) &
     &  /-8.04270E+01, 2.51589E+02, 7.15964E+02 /
      DATA (((TSTUU(I,J,K,3,7), J = 1, 4), I = 1, 3), K = 1, 2) &
     &  / 1.01805E+02, 1.01805E+02, 1.46775E+02, 1.49033E+02, &
     &    1.06583E+02, 1.28565E+02, 1.04464E+02, 1.59054E+02, &
     &    9.66519E+01, 8.65854E+01, 1.89259E+02, 1.89259E+02, &
     &    1.01805E+02, 1.01805E+02, 1.29641E+02, 1.49033E+02, &
     &    1.06583E+02, 1.06408E+02, 9.48418E+01, 1.59054E+02, &
     &    9.66519E+01, 7.49310E+01, 1.89259E+02, 1.89259E+02 /

! ********************* Test Case 7D *********************************

      DATA (TSTFIR(I,4,7), I = 1, 3) &
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,4,7), I = 1, 3) &
     &  / 3.19830E+02, 3.50555E+02, 2.92063E+02 /
      DATA (TSTFUP(I,4,7), I = 1, 3) &
     &  / 3.12563E+02, 2.68126E+02, 3.05596E+02 /
      DATA (TSTDFD(I,4,7), I = 1, 3) &
     &  /-1.68356E+02, 1.01251E+02, 4.09326E+02 /
      DATA (((TSTUU(I,J,K,4,7), J = 1, 4), I = 1, 3), K = 1, 2) &
     &  / 1.01805E+02, 1.01805E+02, 1.40977E+02, 9.62764E+01, &
     &    1.06203E+02, 1.23126E+02, 9.19545E+01, 8.89528E+01, &
     &    9.56010E+01, 7.25576E+01, 9.72743E+01, 9.72743E+01, &
     &    1.01805E+02, 1.01805E+02, 1.23843E+02, 9.62764E+01, &
     &    1.06203E+02, 1.00969E+02, 8.23318E+01, 8.89528E+01, &
     &    9.56010E+01, 6.09031E+01, 9.72743E+01, 9.72743E+01 /

! ********************* Test Case 7E *********************************

      DATA (TSTFIR(I,5,7), I = 1, 3) &
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,5,7), I = 1, 3) &
!    &  / 3.19830E+02, 3.53275E+02, 2.99002E+02 /
     &  / 3.19830E+02, 3.54468E+02, 3.02366E+02 /
      DATA (TSTFUP(I,5,7), I = 1, 3) &
!    &  / 4.04300E+02, 4.07843E+02, 5.29248E+02 /
     &  / 4.40940E+02, 4.64624E+02, 6.23842E+02 /
      DATA (TSTDFD(I,5,7), I = 1, 3) &
!    &  /-9.98568E+01, 2.17387E+02, 6.38461E+02 /
     &  /-7.16977E+01, 2.66918E+02, 7.50170E+02 /
      DATA (((TSTUU(I,J,K,5,7), J = 1, 4), I = 1, 3), K = 1, 2) &
!    &  / 1.01805E+02, 1.01805E+02, 1.45448E+02, 1.38554E+02,
!    &    1.06496E+02, 1.27296E+02, 1.01395E+02, 1.45229E+02,
!    &    9.63993E+01, 8.29009E+01, 1.60734E+02, 1.71307E+02,
!    &    1.01805E+02, 1.01805E+02, 1.28281E+02, 1.38554E+02,
!    &    1.06496E+02, 1.05111E+02, 9.16726E+01, 1.45229E+02,
!    &    9.63993E+01, 7.11248E+01, 1.59286E+02, 1.71307E+02 /
     &  / 1.01805E+02, 1.01805E+02, 1.47357E+02, 1.53713E+02, &
     &    1.06621E+02, 1.29120E+02, 1.05792E+02, 1.65216E+02, &
     &    9.67644E+01, 8.81547E+01, 2.00886E+02, 1.97240E+02, &
     &    1.01805E+02, 1.01805E+02, 1.30246E+02, 1.53713E+02, &
     &    1.06621E+02, 1.06984E+02, 9.62371E+01, 1.65216E+02, &
     &    9.67644E+01, 7.65889E+01, 2.01922E+02, 1.97240E+02 /

! ********************* Test Case 8A *********************************

      DATA (TSTFIR(I,1,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,1,8), I = 1, 3) / &
     &  1.00000E+00, 7.22235E-01, 5.13132E-01 /
      DATA (TSTFUP(I,1,8), I = 1, 3) / 9.29633E-02, 2.78952E-02, 0.0 /
      DATA (TSTDFD(I,1,8), I = 1, 3) / &
     &  1.12474E+00, 6.51821E-01, 5.63361E-01 /
      DATA ((TSTUU(I,J,1,1,8), J = 1, 4), I = 1, 3) / &
     &  2*3.18310E-01, 5.62566E-02, 1.94423E-02, &
     &  2.62711E-01, 1.36952E-01, 1.84909E-02, 5.52188E-03, &
     &  2.10014E-01, 5.60376E-02, 2*0.0 /

! ********************* Test Case 8B *********************************

      DATA (TSTFIR(I,2,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,2,8), I = 1, 3) / &
     &  1.00000E+00, 7.95332E-01, 6.50417E-01 /
      DATA (TSTFUP(I,2,8), I = 1, 3) / 2.25136E-01, 1.26349E-01, 0.0 /
      DATA (TSTDFD(I,2,8), I = 1, 3) / &
     &  5.12692E-01, 3.56655E-01, 5.68095E-02 /
      DATA ((TSTUU(I,J,1,2,8), J = 1, 4), I = 1, 3) / &
     &  2*3.18310E-01, 1.23687E-01, 4.95581E-02, &
     &  2.77499E-01, 1.83950E-01, 8.35695E-02, 2.50575E-02, &
     &  2.40731E-01, 1.29291E-01, 2*0.0 /

! ********************* Test Case 8C *********************************

      DATA (TSTFIR(I,3,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,3,8), I = 1, 3) / &
     &  1.00000E+00, 4.86157E-01, 1.59984E-01 /
      DATA (TSTFUP(I,3,8), I = 1, 3) / 3.78578E-01, 2.43397E-01, 0.0 /
      DATA (TSTDFD(I,3,8), I = 1, 3) / &
     &  5.65095E-01, 2.76697E-01, 1.35679E-02 /
      DATA ((TSTUU(I,J,1,3,8), J = 1, 4), I = 1, 3) / &
     &  2*3.18310E-01, 1.49335E-01, 1.04766E-01, &
     &  1.89020E-01, 9.88158E-02, 9.65192E-02, 6.54445E-02, &
     &  6.84762E-02, 2.96698E-02, 2*0.0 /

! ********************* Test Case 9A *********************************

      DATA (TSTFIR(I,1,9), I = 1, 5) / 5*0.0 /
      DATA (TSTFDN(I,1,9), I = 1, 5) / &
     &  1.00000E+00, 3.55151E-01, 1.44265E-01, 6.71445E-03, 6.16968E-07/
      DATA (TSTFUP(I,1,9), I = 1, 5) / &
     &  2.27973E-01, 8.75098E-02, 3.61819E-02, 2.19291E-03, 0.0 /
      DATA (TSTDFD(I,1,9), I = 1, 5) / &
     &  8.82116E-01, 2.32366E-01, 9.33443E-02, 3.92782E-03, 1.02500E-07/
      DATA ((TSTUU(I,J,1,1,9), J = 1, 4), I = 1, 5) / &
     &  2*3.18310E-01, 9.98915E-02, 5.91345E-02, &
     &  1.53507E-01, 5.09531E-02, 3.67006E-02, 2.31903E-02, &
     &  7.06614E-02, 2.09119E-02, 1.48545E-02, 9.72307E-03, &
     &  3.72784E-03, 1.08815E-03, 8.83316E-04, 5.94743E-04, &
     &  2.87656E-07, 1.05921E-07, 2*0.0 /

! ********************* Test Case 9B *********************************

      DATA (TSTFIR(I,2,9), I = 1, 5) / 5*0.0 /
      DATA (TSTFDN(I,2,9), I = 1, 5) / &
     &  1.00000E+00, 4.52357E-01, 2.36473E-01, 2.76475E-02, 7.41853E-05/
      DATA (TSTFUP(I,2,9), I = 1, 5) / &
     &  1.00079E-01, 4.52014E-02, 2.41941E-02, 4.16016E-03, 0.0 /
      DATA (TSTDFD(I,2,9), I = 1, 5) / &
     &  8.04577E-01, 2.55330E-01, 1.30976E-01, 1.36227E-02, 1.22022E-05/
      DATA ((TSTUU(I,J,1,2,9), J = 1, 4), I = 1, 5) / &
     &  2*3.18310E-01, 7.39198E-02, 1.32768E-02, &
     &  1.96609E-01, 5.92369E-02, 3.00230E-02, 7.05566E-03, &
     &  1.15478E-01, 3.01809E-02, 1.52672E-02, 4.06932E-03, &
     &  1.46177E-02, 3.85590E-03, 2.38301E-03, 7.77890E-04, &
     &  3.37742E-05, 1.20858E-05, 2*0.0 /

! ********************* Test Case 9C *********************************

      DATA (TSTFIR(I,3,9), I = 1, 5) / &
     &  1.57080E+00, 1.92354E-01, 2.35550E-02, 9.65131E-06, 9.03133E-19/
      DATA (TSTFDN(I,3,9), I = 1, 5 ) / &
     &  6.09217E+00, 4.97279E+00, 4.46616E+00, 4.22731E+00, 4.73767E+00/
      DATA (TSTFUP(I,3,9), I = 1, 5) / &
     &  4.68414E+00, 4.24381E+00, 4.16941E+00, 4.30667E+00, 5.11524E+00/
      DATA (TSTDFD(I,3,9), I = 1, 5) / &
     &  3.49563E+00, 8.81206E-01, 3.50053E-01, 1.93471E-02, 7.15349E-02/
      DATA (((TSTUU(I,J,K,3,9), J = 1, 4), I = 1, 5), K = 1, 3) &
     & / 1.93920E+00, 1.93920E+00, 1.61855E+00, 1.43872E+00, &
     &   1.66764E+00, 1.44453E+00, 1.38339E+00, 1.33890E+00, &
     &   1.48511E+00, 1.35009E+00, 1.33079E+00, 1.32794E+00, &
     &   1.34514E+00, 1.35131E+00, 1.35980E+00, 1.37918E+00, &
     &   1.48927E+00, 1.54270E+00, 1.62823E+00, 1.62823E+00, &
     &   1.93920E+00, 1.93920E+00, 1.57895E+00, 1.43872E+00, &
     &   1.66764E+00, 1.42925E+00, 1.37317E+00, 1.33890E+00, &
     &   1.48511E+00, 1.34587E+00, 1.32921E+00, 1.32794E+00, &
     &   1.34514E+00, 1.35129E+00, 1.35979E+00, 1.37918E+00, &
     &   1.48927E+00, 1.54270E+00, 1.62823E+00, 1.62823E+00, &
     &   1.93920E+00, 1.93920E+00, 1.56559E+00, 1.43872E+00, &
     &   1.66764E+00, 1.42444E+00, 1.37034E+00, 1.33890E+00, &
     &   1.48511E+00, 1.34469E+00, 1.32873E+00, 1.32794E+00, &
     &   1.34514E+00, 1.35128E+00, 1.35979E+00, 1.37918E+00, &
     &   1.48927E+00, 1.54270E+00, 2*1.62823E+00 /

contains

!-------------------------------------------------------------------------------
subroutine Getmom(Iphas, GG, Nmom, Pmom)
  implicit none
  integer(Kint), intent(in) :: Iphas, Nmom
  real(Kreal), intent(in) :: GG
  real(Kreal), intent(out) :: Pmom(0:Nmom)

  ! Haze L from Garcia and Siewert (1985)
  real(Kreal):: Hazelm(0:82) = (/1.0, &
      2.41260,3.23047,3.37296,3.23150,2.89350,2.49594,2.11361, & ! 7
      1.74812,1.44692,1.17714,0.96643,0.78237,0.64114,0.51966, & !14
      0.42563,0.34688,0.28351,0.23317,0.18963,0.15788,0.12739, & !21
      0.10762,0.08597,0.07381,0.05828,0.05089,0.03971,0.03524, & !28
      0.02720,0.02451,0.01874,0.01711,0.01298,0.01198,0.00904, & !35
      0.00841,0.00634,0.00592,0.00446,0.00418,0.00316,0.00296, & !42
      0.00225,0.00210,0.00160,0.00150,0.00115,0.00107,0.00082, & !49
      0.00077,0.00059,0.00055,0.00043,0.00040,0.00031,0.00029, & !56
      0.00023,0.00021,0.00017,0.00015,0.00012,0.00011,0.00009, & !63
      0.00008,0.00006,0.00006,0.00005,0.00004,0.00004,0.00003, & !70
      0.00003,0.00002,0.00002,0.00002,0.00001,0.00001,0.00001, & !77
      0.00001,0.00001,0.00001,0.00001,0.00001/)

  ! Cloud C1 from Garcia and Siewert (1985)
  real(Kreal):: Cldmom(0:299) = (/1.0, &
      2.544, 3.883, 4.568, 5.235, 5.887, 6.457, 7.177, 7.859, 8.494, 9.286, &
      9.856,10.615,11.229,11.851,12.503,13.058,13.626,14.209,14.660,15.231, &
     15.641,16.126,16.539,16.934,17.325,17.673,17.999,18.329,18.588,18.885, &
     19.103,19.345,19.537,19.721,19.884,20.024,20.145,20.251,20.330,20.401, &
     20.444,20.477,20.489,20.483,20.467,20.427,20.382,20.310,20.236,20.136, &
     20.036,19.909,19.785,19.632,19.486,19.311,19.145,18.949,18.764,18.551, &
     18.348,18.119,17.901,17.659,17.428,17.174,16.931,16.668,16.415,16.144, &
     15.883,15.606,15.338,15.058,14.784,14.501,14.225,13.941,13.662,13.378, &
     13.098,12.816,12.536,12.257,11.978,11.703,11.427,11.156,10.884,10.618, &
     10.350,10.090, 9.827, 9.574, 9.318, 9.072, 8.822, 8.584, 8.340, 8.110, &
      7.874, 7.652, 7.424, 7.211, 6.990, 6.785, 6.573, 6.377, 6.173, 5.986, &
      5.790, 5.612, 5.424, 5.255, 5.075, 4.915, 4.744, 4.592, 4.429, 4.285, &
      4.130, 3.994, 3.847, 3.719, 3.580, 3.459, 3.327, 3.214, 3.090, 2.983, &
      2.866, 2.766, 2.656, 2.562, 2.459, 2.372, 2.274, 2.193, 2.102, 2.025, &
      1.940, 1.869, 1.790, 1.723, 1.649, 1.588, 1.518, 1.461, 1.397, 1.344, &
      1.284, 1.235, 1.179, 1.134, 1.082, 1.040, 0.992, 0.954, 0.909, 0.873, &
      0.832, 0.799, 0.762, 0.731, 0.696, 0.668, 0.636, 0.610, 0.581, 0.557, &
      0.530, 0.508, 0.483, 0.463, 0.440, 0.422, 0.401, 0.384, 0.364, 0.349, &
      0.331, 0.317, 0.301, 0.288, 0.273, 0.262, 0.248, 0.238, 0.225, 0.215, &
      0.204, 0.195, 0.185, 0.177, 0.167, 0.160, 0.151, 0.145, 0.137, 0.131, &
      0.124, 0.118, 0.112, 0.107, 0.101, 0.097, 0.091, 0.087, 0.082, 0.079, &
      0.074, 0.071, 0.067, 0.064, 0.060, 0.057, 0.054, 0.052, 0.049, 0.047, &
      0.044, 0.042, 0.039, 0.038, 0.035, 0.034, 0.032, 0.030, 0.029, 0.027, &
      0.026, 0.024, 0.023, 0.022, 0.021, 0.020, 0.018, 0.018, 0.017, 0.016, &
      0.015, 0.014, 0.013, 0.013, 0.012, 0.011, 0.011, 0.010, 0.009, 0.009, &
      0.008, 0.008, 0.008, 0.007, 0.007, 0.006, 0.006, 0.006, 0.005, 0.005, &
      0.005, 0.005, 0.004, 0.004, 0.004, 0.004, 0.003, 0.003, 0.003, 0.003, &
      0.003, 0.003, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, &
      0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, &
      0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001 /)


  integer(Kint):: k


  Pmom(0) = 1.0_Kreal
  Pmom(1:) = 0.0_Kreal

  if (Iphas == 1) then
    ! Isotropic phase function
  elseif (Iphas == 2) then
    Pmom(2) = 0.1     ! Rayleigh phase function
  elseif (Iphas == 3) then
    do k = 1, Nmom
      Pmom(k) = GG**k
    enddo
  elseif (Iphas == 4) then
    do k = 1, min(82,Nmom)
      Pmom(k) = Hazelm(k) / (2*k+1)
    enddo
  elseif (Iphas == 5) then
    do k = 1, min(298,Nmom)
      Pmom(k) = Cldmom(k) / (2*k+1)
    enddo

  else
    print*, 'Bad Iphas in Getmom: ', Iphas
    stop
  endif

end subroutine Getmom
!-------------------------------------------------------------------------------
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! $Rev: 55 $ $Date: 2014-12-31 12:16:59 -0500 (Wed, 31 Dec 2014) $
! FORTRAN 77
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE  ErrMsg( MESSAG, FATAL )

!        Print out a warning or error message;  abort if error

      LOGICAL       FATAL, MsgLim
      CHARACTER*(*) MESSAG
      INTEGER       MaxMsg, NumMsg
      SAVE          MaxMsg, NumMsg, MsgLim
      DATA NumMsg / 0 /,  MaxMsg / 100 /,  MsgLim / .FALSE. /


      IF ( FATAL )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg.LE.MaxMsg )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG
      ELSE
         WRITE ( *,99 )
         MsgLim = .True.
      ENDIF

      RETURN

   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ', &
     &   'They will no longer be printed  <<<<<<<', // )
      ENDsubroutine ErrMsg

! ---------------------------------------------------------------------
      REAL FUNCTION RATIO( A, B )

!        Calculate ratio  A/B  with over- and under-flow protection
!        (thanks to Prof. Jeff Dozier for some suggestions here).
!        Since this routine takes two logs, it is no speed demon,
!        but it is invaluable for comparing results from two runs
!        of a program under development.
!
!        NOTE:  In Fortran90, built-in functions TINY and HUGE
!               can replace the R1MACH calls.
!
!   Called by- DISORT
!   Calls- R1MACH
! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..

      REAL(Kreal)      A, B
!     ..
!     .. Local Scalars ..

      LOGICAL   PASS1
      REAL(Kreal)      ABSA, ABSB, HUGE0, POWA, POWB, POWMAX, POWMIN, TINY0
!     ..
!     .. External Functions ..

      !REAL      R1MACH
      !EXTERNAL  R1MACH
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC ABS, LOG10, SIGN, TINY, HUGE
!     ..
!     .. Save statement ..

      SAVE      PASS1, TINY0, HUGE0, POWMAX, POWMIN
!     ..
!     .. Data statements ..

      DATA      PASS1 / .TRUE. /
!     ..


      IF( PASS1 ) THEN

         TINY0   = TINY(1.0) != R1MACH( 1 )
         HUGE0   = HUGE(1.0) != R1MACH( 2 )
         POWMAX = LOG10( HUGE0 )
         POWMIN = LOG10( TINY0 )
         PASS1  = .FALSE.

      END IF


      IF( A.EQ.0.0 ) THEN

         IF( B.EQ.0.0 ) THEN

            RATIO  = 1.0

         ELSE

            RATIO  = 0.0

         END IF


      ELSE IF( B.EQ.0.0 ) THEN

         RATIO  = SIGN( HUGE0, A )

      ELSE

         ABSA   = ABS( A )
         ABSB   = ABS( B )
         POWA   = LOG10( ABSA )
         POWB   = LOG10( ABSB )

         IF( ABSA.LT.TINY0 .AND. ABSB.LT.TINY0 ) THEN

            RATIO  = 1.0

         ELSE IF( POWA - POWB.GE.POWMAX ) THEN

            RATIO  = HUGE0

         ELSE IF( POWA - POWB.LE.POWMIN ) THEN

            RATIO  = TINY0

         ELSE

            RATIO  = ABSA / ABSB

         END IF
!                      ** DONT use old trick of determining sign
!                      ** from A*B because A*B may (over/under)flow

         IF( ( A.GT.0.0 .AND. B.LT.0.0 ) .OR. &
     &       ( A.LT.0.0 .AND. B.GT.0.0 ) ) RATIO = -RATIO

      END IF


      RETURN
    ENDfunction RATIO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!-------------------------------------------------------------------------------

      SUBROUTINE  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV, &
     &                    MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, &
     &                    UU, TSTFIR, TSTFDN, TSTFUP, TSTDFD, TSTUU, &
     &                    MXTAU, MXMU, MXPHI, NTEST, NPASS )

!        Print DISORT results and, directly beneath them, their
!        ratios to the correct answers;  print number of non-unit
!        ratios that occur but try to count just the cases where
!        there is a real disagreement and not those where flux or
!        intensity are down at their noise level (defined as 10^(-6)
!        times their maximum value).  d(flux)/d(tau) is treated the
!        same as fluxes in this noise estimation even though it
!        is a different type of quantity (although with flux units).

!     INPUT :   TSTFIR  correct direct flux
!               TSTFDN  correct diffuse down flux
!               TSTFUP  correct diffuse up flux
!               TSTDFD  correct d(flux)/d(optical depth)
!               TSTUU   correct intensity
!               (remaining input = DISORT I/O variables)

! --------------------------------------------------------------------

!     .. Parameters ..

      INTEGER, parameter:: MAXRAT=90
!     ..
!     .. Scalar Arguments ..

      LOGICAL   ONLYFL
      INTEGER   MXPHI, MXMU, MXTAU, MAXULV, MAXUMU, NPHI, NTAU, NUMU
!     ..
!     .. Array Arguments ..

      REAL(Kreal)      DFDT( MXTAU ), FLUP( MXTAU ), PHI( MXPHI ), &
     &          RFLDIR( MXTAU ), RFLDN( MXTAU ), &
     &          TSTDFD( MXTAU ), TSTFDN( MXTAU ), &
     &          TSTFIR( MXTAU ), TSTFUP( MXTAU ), &
     &          TSTUU( MXTAU, MXMU, MXPHI ), UMU( MXMU ), UTAU( MXTAU ), &
     &          UU( MAXUMU, MAXULV, * )
!     ..
!     .. Local Scalars ..

      INTEGER  IU, J, LU, NUMBAD
      REAL(Kreal)     FLXMAX, FNOISE, RAT, RAT1, RAT2, RAT3, RAT4, UMAX, UNOISE
!     ..
!     .. Local Arrays ..

      REAL(Kreal)      RATV( NPHI )
!     ..
!     .. External Functions ..

      !REAL      RATIO
      !EXTERNAL  RATIO
!     ..
!     .. External Subroutines ..

      !EXTERNAL  ERRMSG
!     ..
!     .. Statement Functions ..

      LOGICAL   BADRAT
!     ..
!     .. Statement Function definitions ..

      BADRAT( RAT ) = (RAT.LT.0.99) .OR. (RAT.GT.1.01)
!     ..

      IF ( NTAU.GT.MXTAU .OR. NUMU.GT.MXMU .OR. NPHI.GT.MXPHI )  CALL &
     &   ERRMSG( 'PRTFIN--out of bounds in comparator arrays', .TRUE.)

      FLXMAX = 0.0
      DO 5  LU = 1, NTAU
         FLXMAX = MAX( FLXMAX, TSTFIR(LU), TSTFDN(LU), TSTFUP(LU) )
 5    CONTINUE
      FNOISE = 1.E-6 * FLXMAX
      IF( FLXMAX.LE.0.0 ) &
     &    CALL ERRMSG( 'PRTFIN--all fluxes zero or negative', .FALSE.)
      IF( FNOISE.LE.0.0 ) &
     &    CALL ERRMSG( 'PRTFIN--all fluxes near underflowing', .FALSE.) 

      NUMBAD = 0

      WRITE(*,'(//,A,/,A,/,A)') &
     &  '                  <-------------- FLUXES -------------->', &
     &  '    Optical       Downward       Downward         Upward'// &
     &  '    d(Net Flux)', &
     &  '      Depth         Direct        Diffuse        Diffuse'// &
     &  '    / d(Op Dep)'

      DO 10  LU = 1, NTAU

!         WRITE( *,'(0P,F11.4,1P,4E15.4)')  UTAU(LU), RFLDIR(LU),
         WRITE( *,'(0P,F11.5,1P,4E15.5)')  UTAU(LU), RFLDIR(LU), &
     &          RFLDN(LU), FLUP(LU), DFDT(LU)
         RAT1 = RATIO( RFLDIR(LU), TSTFIR(LU) )
         RAT2 = RATIO( RFLDN(LU),  TSTFDN(LU) )
         RAT3 = RATIO(  FLUP(LU),  TSTFUP(LU) )
         RAT4 = RATIO(  DFDT(LU),  TSTDFD(LU) )
!         WRITE( *,'(11X,4( ''    ('',F9.4,'')''))')
!         WRITE( *,'(11X,4( ''   ('',G11.5,'')''))')
         WRITE( *,'(11X,4(''  ('',F11.5,'')''))') &
     &          RAT1, RAT2, RAT3, RAT4

         IF( BADRAT(RAT1) .AND. ABS(RFLDIR(LU)).GT.FNOISE &
     &      .AND. TSTFIR(LU) .NE. 0. ) &
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT2) .AND. ABS(RFLDN(LU)).GT.FNOISE &
     &      .AND. TSTFDN(LU) .NE. 0. ) &
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT3) .AND. ABS(FLUP(LU)).GT.FNOISE &
     &      .AND. TSTFUP(LU) .NE. 0. ) &
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT4) .AND. ABS(DFDT(LU)).GT.FNOISE &
     &      .AND. TSTDFD(LU) .NE. 0. ) &
     &        NUMBAD = NUMBAD+1


10    CONTINUE


      IF ( ONLYFL )  GO TO 100


!                                       ** Print intensities

      IF ( NPHI.GT.8 ) CALL ERRMSG &
     &      ( 'PRTFIN--intensity FORMATs inadequate',.FALSE.)

      UMAX = 0.0
      DO 36  LU = 1, NTAU
         DO 35  IU = 1, NUMU
            DO 34  J = 1, NPHI
               UMAX = MAX( UMAX, TSTUU(LU,IU,J) )
 34         CONTINUE
 35      CONTINUE
 36   CONTINUE
      UNOISE = 1.E-6 * UMAX
      IF( UMAX.LE.0.0 )  CALL ERRMSG &
     &     ( 'PRTFIN--all intensities zero or negative',.FALSE.)
      IF( UNOISE.LE.0.0 ) CALL ERRMSG &
     &     ( 'PRTFIN--all intensities near underflowing',.FALSE.)

      WRITE( *,'(//,A,//,A,/,A,/,A,8(F10.1,4X))' ) &
     &            ' ********  I N T E N S I T I E S  *********', &
     &            '             Polar   Azimuthal Angles (Degrees)', &
     &            '   Optical   Angle', &
     &            '     Depth  Cosine', ( PHI(J), J = 1, NPHI )

      DO 60  LU = 1, NTAU

         DO 50  IU = 1, NUMU

!            IF( IU.EQ.1 ) WRITE( *,'(/,0P,F10.3,F8.3,1P,8E14.4)')
            IF( IU.EQ.1 ) WRITE( *,'(/,0P,F10.3,F8.3,1P,8E16.5)') &
     &          UTAU(LU), UMU(IU), ( UU( IU,LU,J ), J = 1, NPHI )

!            IF( IU.GT.1 ) WRITE( *,'(10X,0P,F8.3, 1P,8E14.4)')
            IF( IU.GT.1 ) WRITE( *,'(10X,0P,F8.3, 1P,8E16.5)') &
     &                    UMU(IU), ( UU( IU,LU,J ), J = 1, NPHI )

            DO 40  J = 1, NPHI
               RATV(J) = RATIO( UU(IU,LU,J), TSTUU(LU,IU,J) )
               IF( BADRAT(RATV(J)) .AND. ABS(UU(IU,LU,J)).GT.UNOISE ) &
     &                  NUMBAD = NUMBAD + 1
   40       CONTINUE

!            WRITE( *,'(18X, 8(:,''   ('',F9.4,'')''))')
            WRITE( *,'(18X, 8(:,''   ('',F11.5,'')''))') &
     &           ( RATV(J), J = 1, NPHI )

   50    CONTINUE
   60 CONTINUE

  100 CONTINUE
      IF( NUMBAD.GT.0 ) THEN
          WRITE( *,300)  ' ====  ', NUMBAD, &
     &    '  SERIOUSLY NON-UNIT RATIOS    ===='
      ELSE
          NPASS = NPASS+1
      ENDIF
      NTEST = NTEST+1
      RETURN

300   FORMAT( //,1X,45('='),/,A,I4,A,/,1X,45('=') )
      ENDsubroutine PRTFIN

!-------------------------------------------------------------------------------
      SUBROUTINE PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER, &
     &                   WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU, &
     &                   NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, &
     &                   LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, &
     &                   PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT, &
     &                   OPRIM, TAUC, TAUCPR, MAXMOM, PRTMOM, &
     &                   DO_PSEUDO_SPHERE, H_LYR, DELTAMPLUS )

!        Print values of input variables
!
!   Called by- DISORT
! --------------------------------------------------------------------

!     .. Scalar Arguments ..

      LOGICAL   CORINT, DELTAM, LAMBER, LYRCUT, ONLYFL, PLANK, PRTMOM
      LOGICAL   DO_PSEUDO_SPHERE, DELTAMPLUS
      INTEGER   IBCND, MAXMOM, NLYR, NMOM, NPHI, NSTR, NTAU, NUMU
      REAL(Kreal) ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP, &
     &          UMU0, WVNMHI, WVNMLO
!     ..
!     .. Array Arguments ..

      REAL(Kreal)      DTAUC( * ), DTAUCP( * ), FLYR( * ), OPRIM( * ), &
!     &          PHI( * ), PMOM( 0:MAXMOM, * ), SSALB( * ), TAUC( 0:* ), &
     &          PHI( * ), PMOM( 0:, : ), SSALB( * ), TAUC( 0:* ), &
     &          TAUCPR( 0:* ), TEMPER( 0:* ), UMU( * ), UTAU( * ), &
     &          H_LYR( * ) 
!     ..
!     .. Local Scalars ..

      INTEGER   IU, J, K, LC, LU
      REAL(Kreal)      YESSCT
!     ..


      WRITE( *, '(/,A,I4,A,I4)' ) ' No. streams =', NSTR, &
     &       '     No. computational layers =', NLYR

      IF( IBCND.NE.1 ) WRITE( *, '(I4,A,10F10.4,/,(26X,10F10.4))' ) &
     &    NTAU, ' User optical depths :', ( UTAU(LU), LU = 1, NTAU )

      IF( .NOT.ONLYFL ) WRITE( *, '(I4,A,10F9.5,/,(31X,10F9.5))' ) &
     &    NUMU, ' User polar angle cosines :', ( UMU(IU), IU = 1, NUMU )

      IF( .NOT.ONLYFL .AND. IBCND.NE.1 ) &
     &    WRITE( *, '(I4,A,10F9.2,/,(28X,10F9.2))' ) &
     &           NPHI,' User azimuthal angles :',( PHI(J), J = 1, NPHI )

      IF( .NOT.PLANK .OR. IBCND.EQ.1 ) &
     &    WRITE( *, '(A)' ) ' No thermal emission'


      WRITE( *, '(A,I2)' ) ' Boundary condition flag: IBCND =', IBCND

      IF( IBCND.EQ.0 ) THEN

         WRITE( *, '(A,1P,E11.3,A,0P,F8.5,A,F7.2,/,A,1P,E11.3)' ) &
     &          '    Incident beam with intensity =', FBEAM, &
     &          ' and polar angle cosine = ', UMU0, &
     &          '  and azimuth angle =', PHI0, &
     &          '    plus isotropic incident intensity =', FISOT

         IF( LAMBER ) WRITE( *, '(A,0P,F8.4)' ) &
     &                '    Bottom albedo (Lambertian) =', ALBEDO

         IF( .NOT.LAMBER ) WRITE( *, '(A)' ) &
     &       '    Bidirectional reflectivity at bottom'

         IF( PLANK ) WRITE( *, '(A,2F14.4,/,A,F10.2,A,F10.2,A,F8.4)' ) &
     &       '    Thermal emission in wavenumber interval :', WVNMLO, &
     &       WVNMHI, &
     &       '    Bottom temperature =', BTEMP, &
     &       '    Top temperature =', TTEMP, &
     &       '    Top emissivity =', TEMIS

      ELSE IF( IBCND.EQ.1 ) THEN

         WRITE( *, '(A)' ) &
     &          '    Isotropic illumination from top and bottom'
         WRITE( *, '(A,0P,F8.4)' ) &
     &          '    Bottom albedo (Lambertian) =', ALBEDO

      END IF


      IF( DELTAM ) WRITE( *, '(A)' ) ' Uses delta-M method'
      IF( DELTAMPLUS ) WRITE( *, '(A)' ) ' Uses New-Delta-M+ method'
      IF( .NOT.DELTAM .AND. .NOT. DELTAMPLUS ) &
     & WRITE( *, '(A)' ) ' Does not use delta-M / delta-M Plus method'

      IF( CORINT ) WRITE( *, '(A)' ) ' Uses TMS/IMS method'
      IF( .NOT.CORINT ) WRITE( *,'(A)' ) ' Does not use TMS/IMS method'

      IF(DO_PSEUDO_SPHERE) THEN
        WRITE(*,'(A)') ' Uses pseudo spherical method'
        WRITE(*,'(A15, 13(F5.2,2x))')' Layer height: ', &
     &                               (H_LYR(LC),LC=1,NLYR)
      ENDIF
      IF(.NOT.DO_PSEUDO_SPHERE) WRITE(*, '(A)' ) ' Uses plane'// &
     & ' parallel method'



      IF( IBCND.EQ.1 ) THEN

         WRITE( *, '(A)' ) ' Calculate albedo and transmissivity of'// &
     &                     ' medium vs. incident beam angle'

      ELSE IF( ONLYFL ) THEN

         WRITE( *, '(A)' ) &
     &          ' Calculate fluxes only'

      ELSE

         WRITE( *, '(A)' ) ' Calculate fluxes and intensities'

      END IF

      WRITE( *, '(A,1P,E11.2)' ) &
     &       ' Relative convergence criterion for azimuth series =', &
     &       ACCUR

      IF( LYRCUT ) WRITE( *, '(A)' ) &
     &    ' Sets radiation = 0 below absorption optical depth 10'


!                                    ** Print layer variables
!                                    ** (to read, skip every other line)

      IF( PLANK ) WRITE( *, '(/,37X,A,3(/,2A))' ) &
     &  '<------------- Delta-M --------------->', &
     &  '                   Total    Single                           ', &
     &  'Total    Single', &
     &  '       Optical   Optical   Scatter   Separated   ', &
     &  'Optical   Optical   Scatter    Asymm', &
     &  '         Depth     Depth    Albedo    Fraction     ', &
     &  'Depth     Depth    Albedo   Factor   Temperature'

      IF( .NOT.PLANK ) WRITE( *, '(/,37X,A,3(/,2A))' ) &
     &  '<------------- Delta-M --------------->', &
     &  '                   Total    Single                           ', &
     &  'Total    Single', &
     &  '       Optical   Optical   Scatter   Separated   ', &
     &  'Optical   Optical   Scatter    Asymm', &
     &  '         Depth     Depth    Albedo    Fraction     ', &
     &  'Depth     Depth    Albedo   Factor'


      YESSCT = 0.0

      DO 10 LC = 1, NLYR

         YESSCT = YESSCT + SSALB( LC )
!                                       ** f90 nonadvancing I/O would
!                                       ** simplify this a lot (also the
!                                       ** two WRITEs above)
         IF( PLANK ) &
     &       WRITE( *,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4,F14.3)') &
     &             LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ), &
     &             DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), &
     &             PMOM( 1,LC ), TEMPER( LC-1 )

         IF( .NOT.PLANK ) &
     &       WRITE( *,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4)' ) &
     &             LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ), &
     &             DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM( 1,LC )
   10 CONTINUE

      IF( PLANK ) WRITE( *, '(85X,F14.3)' ) TEMPER( NLYR )


      IF( PRTMOM .AND. YESSCT.GT.0.0 ) THEN

         WRITE( *, '(/,A,I5)' ) ' Number of Phase Function Moments = ', &
     &        NMOM + 1
         WRITE( *, '(A)' ) ' Layer   Phase Function Moments'

         DO 20 LC = 1, NLYR

            IF( SSALB( LC ).GT.0.0 ) &
     &          WRITE( *, '(I6,10F11.6,/,(6X,10F11.6))' ) &
     &                 LC, ( PMOM( K, LC ), K = 0, NMOM )
   20    CONTINUE

      END IF


      RETURN
      END subroutine PRTINP


end module Disotestaux  
!===============================================================================

program Sjhmtest

  use Disotest_vars
  use Mod_Gdrtes, only: Kreal, Kint
  use Sjhm_module, only: Nm, Ng, Nh, Lup, Ldn, Numu1=>Numu, Nphi1=>Nphi, Umu1=>Umu, Phi1=>Phi, &
    Sjhminit, Sjhmfin, Sjhm
  use Disotestaux, only: Getmom, PRTFIN, TSTFIR, TSTFDN, TSTFUP, TSTDFD, TSTUU, Mxmu, Mxphi, Mxtau, PRTINP, MAXMOM

  implicit none 

  !logical::  DOPROB(17)
  integer(Kint):: ICAS, NTEST, NPASS, NPROB, Lev
  CHARACTER:: HEADER*127
  CHARACTER::  ABC(18)*1, TITLE*100, BLANKS*3
  DATA ABC / 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', &
            'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r' /, &
       BLANKS / '   ' /

  real(Kreal),parameter :: PI = 2.*ASIN(1.0)
  real(Kreal):: Gs(0:0)

  real(Kreal), allocatable:: DTAUCP(:), FLYR(:), OPRIM(:), TAUC(:), TAUCPR(:)
  logical:: DELTAM, CORINT, LYRCUT=.false.
  integer:: LU, LENTIT, K, LC

  real(Kreal), allocatable:: Rad(:,:,:)

  !DOPROB(:) = .false.
  !DOPROB(2) = .true.
  ACCUR = 0.0

  if (DOPROB(1)) then

  !c **********************************************************************
  !c ****  Test Problem 1:  Isotropic Scattering                       ****
  !c ****  (Compare to Ref. VH1, Table 12)                             ****
  !c **********************************************************************

    USRTAU    = .TRUE.
    USRANG    = .TRUE.
    LAMBER    = .TRUE.
    PLANK     = .FALSE.
    ONLYFL    = .FALSE.
    DO_PSEUDO_SPHERE = .FALSE.
    DELTAMPLUS = .FALSE.

    NSTR = 16; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
    NLYR = 1
    NMOM = NSTR 
    NTAU      = 2; IF(.NOT.USRTAU) NTAU = NLYR + 1
    NUMU      = 6; IF(.NOT.USRANG) NUMU = NSTR
    NPHI      = 1; 
    IBCND     = 0

    ! Initialize HSHM parameters
    Nm = NSTR
    Ng = NMOM+1
    Nh = 1
    Lup = 0 !-Nstr/2
    Ldn = 0 !-Nstr/2

    call Sjhminit(.false.)

    ! for PRTINP
    allocate(DTAUCP(NLYR),FLYR(NLYR),OPRIM(NLYR),TAUC(0:NLYR),TAUCPR(0:NLYR))
  
    do ICAS = 1, 6

      if (ICAS.eq.1) then
  
        call allocate_disort_allocatable_arrays(NLYR, NMOM, NSTR, NUMU, NPHI, NTAU)

        UMU0      = 0.1;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0; Gs(0)     = ALBEdO/PI
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 0.03125
        CALL  GETMOM( 1, 0.0_Kreal, NMOM, PMOM(:,1) )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 0.2
        FBEAM      = PI / UMU0; FISOT      = 0.0

      ELSE IF ( ICAS.EQ.2 ) THEN

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 0.1;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 0.03125
        CALL  GETMOM( 1, 0.0_Kreal, NMOM, PMOM(:,1) )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 1.0
        FBEAM      = PI / UMU0; FISOT      = 0.0

      ELSE IF ( ICAS.EQ.3 ) THEN

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )
  
        UMU0      = 0.1;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 0.03125
        CALL  GETMOM( 1, 0.0_Kreal, NMOM, PMOM(:,1) )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 0.99
        FBEAM      = 0.0;       FISOT      = 1.0

      ELSE IF ( ICAS.EQ.4 ) THEN

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 0.1;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 32.0
        CALL  GETMOM( 1, 0.0_Kreal, NMOM, PMOM(:,1) )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 0.2
        FBEAM      = PI / UMU0; FISOT      = 0.0

      ELSE IF ( ICAS.EQ.5 ) THEN

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 0.1;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 32.0
        CALL  GETMOM( 1, 0.0_Kreal, NMOM, PMOM(:,1) )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 1.0
        FBEAM      = PI / UMU0; FISOT      = 0.0

      ELSE IF ( ICAS.EQ.6 ) THEN

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 0.1;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 32.0
        CALL  GETMOM( 1, 0.0_Kreal, NMOM, PMOM(:,1) )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 0.99
        FBEAM      = 0.0;       FISOT      = 1.0

      endif

      Numu1 = NUMU
      Nphi1 = NPHI
      Umu1 = UMU
      Phi1 = PHI

      GS(:) = 0.0
      GS(0) = ALBEDO / PI

      if (USRTAU) then
        allocate(Rad(Numu1,Nphi1,0:NTAU-1))
      else
        allocate(Rad(Numu1,Nphi1,0:Nlyr))
      endif

      WRITE( HEADER,'(3A,F9.5,A,F5.2)') 'Test Case No. 1',ABC(ICAS), &
        ':  Isotropic Scattering, Ref. VH1, Table 12:  b =', &
        UTAU(2), ', a =', SSALB(1)

      WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' ) ' DISORT: '// trim(HEADER)
      
      ! print input parameters
      DELTAM     = .FALSE.
      CORINT     = .FALSE.
      DTAUCP = DTAUC
      FLYR = PMOM(NMOM,1:NLYR)
      OPRIM = SSALB
      TAUC(0) = 0.0
      TAUCPR(0) = 0.0
      do LU = 1, NLYR
        TAUC(LU) = TAUC(LU-1) + DTAUC(LU)
        TAUCPR(LU) = TAUCPR(LU-1) + DTAUCP(LU)
      enddo

      CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER, &
     &               WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU, &
     &               NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, &
     &               LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, &
     &               PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT, &
     &               OPRIM, TAUC, TAUCPR, MAXMOM, PRNT( 5 ), &
     &               DO_PSEUDO_SPHERE, H_LYR, DELTAMPLUS)


      call Sjhm(Nlay=NLYR, Od0=DTAUC, Ssa0=SSALB, Pf0=PMOM, &
      & I0=FBEAM, Mu0=UMU0, T0=TEMPER, Tsurf=BTEMP, GS=Gs, &
      & Thermal=PLANK, Deltam=DELTAM, Wnumlo=WVNMLO, Wnumhi=WVNMHI, &
      & Beam0=RFLDIR, Fisot=FISOT, Ttemp=TTEMP, Temis=TEMIS, &
      & Usertau=USRTAU, Nutau=NTAU, Utau=UTAU, &
      & Rad=Rad, Fup=FLUP, Fdn=RFLDN, Act=UAVG, Dfdt=DFDT, Crt=ACCUR) 

      if (USRTAU) then
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Ntau-1
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      else
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Nlyr
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      endif


      NPROB = 1
      CALL PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, NTAU,            &
                 NUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,       &
                 TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),        &
                 TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),        &
                 TSTUU(1,1,1,ICAS,NPROB), MXTAU, MXMU, MXPHI,       &
                 NTEST, NPASS )

      deallocate(Rad)
      call deallocate_disort_allocatable_arrays()

      deallocate(DTAUCP,FLYR,OPRIM)
    
    enddo

    call Sjhmfin()

  endif

  IF( DOPROB(2) )  THEN
  !c **********************************************************************
  !c ****  Test Problem 2:  Rayleigh Scattering, Beam Source           ****
  !c ****  (Compare To Ref. SW, Table 1)                               ****
  !c **********************************************************************

    USRTAU    = .TRUE.
    USRANG    = .TRUE.
    LAMBER    = .TRUE.
    PLANK     = .FALSE.
    ONLYFL    = .FALSE.
    DO_PSEUDO_SPHERE = .FALSE.
    DELTAMPLUS = .FALSE.

    NSTR = 16; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
    NLYR = 1; 
    NMOM = NSTR 
    NTAU      = 2; IF(.NOT.USRTAU) NTAU = NLYR + 1
    NUMU      = 6; IF(.NOT.USRANG) NUMU = NSTR
    NPHI      = 1; 
    IBCND     = 0

    ! Initialize HSHM parameters
    Nm = NSTR
    Ng = NMOM+1
    Nh = 1
    Lup = 0 !-Nstr/2
    Ldn = 0 !-Nstr/2
    call Sjhminit(.false.)

    DO  ICAS = 1, 4

      IF ( ICAS.EQ.1 ) THEN
 
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 0.080442;  PHI0      = 0.0
        UMU( 1 )  = -0.981986; UMU( 2 )  = -0.538263; UMU( 3 )  = -0.018014
        UMU( 4 )  =  0.018014; UMU( 5 )  =  0.538263; UMU( 6 )  =  0.981986
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 0.2
        CALL  GETMOM( 2, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 0.5
        FBEAM      = PI;        FISOT      = 0.0

      ELSEIF ( ICAS .EQ. 2 ) THEN 
       
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 0.080442;  PHI0      = 0.0
        UMU( 1 )  = -0.981986; UMU( 2 )  = -0.538263; UMU( 3 )  = -0.018014
        UMU( 4 )  =  0.018014; UMU( 5 )  =  0.538263; UMU( 6 )  =  0.981986
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 0.2
        CALL  GETMOM( 2, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 1.0
        FBEAM      = PI;        FISOT      = 0.0
              
      ELSEIF ( ICAS .EQ. 3 ) THEN 
       
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 0.080442;  PHI0      = 0.0
        UMU( 1 )  = -0.981986; UMU( 2 )  = -0.538263; UMU( 3 )  = -0.018014
        UMU( 4 )  =  0.018014; UMU( 5 )  =  0.538263; UMU( 6 )  =  0.981986
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 5.0
        CALL  GETMOM( 2, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 0.5
        FBEAM      = PI;        FISOT      = 0.0

      ELSEIF ( ICAS .EQ. 4 ) THEN 
       
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 0.080442;  PHI0      = 0.0
        UMU( 1 )  = -0.981986; UMU( 2 )  = -0.538263; UMU( 3 )  = -0.018014
        UMU( 4 )  =  0.018014; UMU( 5 )  =  0.538263; UMU( 6 )  =  0.981986
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 5.0
        CALL  GETMOM( 2, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 1.0
        FBEAM      = PI;        FISOT      = 0.0                     
       
      ENDIF

      Numu1 = NUMU
      Nphi1 = NPHI
      Umu1 = UMU
      Phi1 = PHI

      GS(:) = 0.0
      GS(0) = ALBEDO / PI

      if (USRTAU) then
        allocate(Rad(Numu1,Nphi1,0:max(Nlyr,NTAU-1)))
      else
        allocate(Rad(Numu1,Nphi1,0:Nlyr))
      endif

      WRITE( HEADER, '(3A,F5.2,A,F9.6,A,F4.2)') &
           'Test Case No. 2', ABC(ICAS),      &
           ', Rayleigh Scattering, Ref. SW, Table 1:  tau =', &
           UTAU(2), ', mu0 =', UMU0, ', ss-albedo =', SSALB(1)          

      WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' ) ' DISORT: '// trim(HEADER)
      
      ! print input parameters
      DELTAM     = .FALSE.
      CORINT     = .FALSE.
      DTAUCP = DTAUC
      FLYR = PMOM(NMOM,1:NLYR)
      OPRIM = SSALB
      TAUC(0) = 0.0
      TAUCPR(0) = 0.0
      do LU = 1, NLYR
        TAUC(LU) = TAUC(LU-1) + DTAUC(LU)
        TAUCPR(LU) = TAUCPR(LU-1) + DTAUCP(LU)
      enddo

      CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER, &
     &               WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU, &
     &               NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, &
     &               LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, &
     &               PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT, &
     &               OPRIM, TAUC, TAUCPR, MAXMOM, PRNT( 5 ), &
     &               DO_PSEUDO_SPHERE, H_LYR, DELTAMPLUS)


      call Sjhm(Nlay=NLYR, Od0=DTAUC, Ssa0=SSALB, Pf0=PMOM, &
      & I0=FBEAM, Mu0=UMU0, T0=TEMPER, Tsurf=BTEMP, GS=Gs, &
      & Thermal=PLANK, Deltam=DELTAM, Wnumlo=WVNMLO, Wnumhi=WVNMHI, &
      & Beam0=RFLDIR, Fisot=FISOT, Ttemp=TTEMP, Temis=TEMIS, &
      & Usertau=USRTAU, Nutau=NTAU, Utau=UTAU, &
      & Rad=Rad, Fup=FLUP, Fdn=RFLDN, Act=UAVG, Dfdt=DFDT, Crt=ACCUR)

      if (USRTAU) then
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Ntau-1
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      else
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Nlyr
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      endif

      NPROB = 2
      CALL PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, NTAU,            &
                 NUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,       &
                 TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),        &
                 TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),        &
                 TSTUU(1,1,1,ICAS,NPROB), MXTAU, MXMU, MXPHI,       &
                 NTEST, NPASS )

      deallocate(Rad)
      call deallocate_disort_allocatable_arrays()

      deallocate(DTAUCP,FLYR,OPRIM)

    enddo

    call Sjhmfin()

  endif

  IF( DOPROB(3) )  THEN
  ! c **********************************************************************
  ! c ****  Test Problem 3:  Henyey-Greenstein Scattering               ****
  ! c ****  (Compare To Ref. VH2, Table 37)                             ****
  ! c **********************************************************************

    USRTAU    = .TRUE.
    USRANG    = .TRUE.
    LAMBER    = .TRUE.
    PLANK     = .FALSE.
    ONLYFL    = .FALSE.
    DO_PSEUDO_SPHERE = .FALSE.
    DELTAMPLUS = .FALSE.

    NSTR = 16; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
    NLYR = 1; 
    NMOM = 32
    NTAU      = 2; IF(.NOT.USRTAU) NTAU = NLYR + 1
    NUMU      = 6; IF(.NOT.USRANG) NUMU = NSTR
    NPHI      = 1; 
    IBCND     = 0

    ! Initialize HSHM parameters
    Nm = NSTR
    Ng = NMOM+1
    Nh = 1
    Lup = 0 !-Nstr/2
    Ldn = 0 !-Nstr/2
    call Sjhminit(.false.)

    DO ICAS = 1, 2

      IF ( ICAS.EQ.1 ) THEN

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 1.0;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 1.0
        CALL  GETMOM( 3, 0.75, NMOM, PMOM )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 1.0
        FBEAM      = PI / UMU0; FISOT      = 0.0

      ELSEIF ( ICAS .EQ. 2 ) THEN 

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 1.0;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 )  = 8.0
        CALL  GETMOM( 3, 0.75, NMOM, PMOM )
        DTAUC( 1 ) = UTAU( 2 ); SSALB( 1 ) = 1.0
        FBEAM      = PI / UMU0; FISOT      = 0.0

      END IF

      Numu1 = NUMU
      Nphi1 = NPHI
      Umu1 = UMU
      Phi1 = PHI

      GS(:) = 0.0
      GS(0) = ALBEDO / PI


      if (USRTAU) then
        allocate(Rad(Numu1,Nphi1,0:max(Nlyr,NTAU-1)))
      else
        allocate(Rad(Numu1,Nphi1,0:Nlyr))
      endif

      WRITE( HEADER, '(3A,F9.5,A,F5.2)') 'Test Case No. 3',ABC(ICAS),  &
            ', Henyey-Greenstein Scattering, Ref. VH2, Table 37,'      &
            //' g = 0.75, b =', UTAU(2), ', a =', SSALB(1)

      WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' ) ' DISORT: '// trim(HEADER)
      
      ! print input parameters
      DELTAM     = .FALSE.
      CORINT     = .FALSE.
      DTAUCP = DTAUC
      FLYR = PMOM(NMOM,1:NLYR)
      OPRIM = SSALB
      TAUC(0) = 0.0
      TAUCPR(0) = 0.0
      do LU = 1, NLYR
        TAUC(LU) = TAUC(LU-1) + DTAUC(LU)
        TAUCPR(LU) = TAUCPR(LU-1) + DTAUCP(LU)
      enddo

      CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER, &
     &               WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU, &
     &               NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, &
     &               LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, &
     &               PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT, &
     &               OPRIM, TAUC, TAUCPR, MAXMOM, PRNT( 5 ), &
     &               DO_PSEUDO_SPHERE, H_LYR, DELTAMPLUS)

      call Sjhm(Nlay=NLYR, Od0=DTAUC, Ssa0=SSALB, Pf0=PMOM, &
      & I0=FBEAM, Mu0=UMU0, T0=TEMPER, Tsurf=BTEMP, GS=Gs, &
      & Thermal=PLANK, Deltam=DELTAM, Wnumlo=WVNMLO, Wnumhi=WVNMHI, &
      & Beam0=RFLDIR, Fisot=FISOT, Ttemp=TTEMP, Temis=TEMIS, &
      & Usertau=USRTAU, Nutau=NTAU, Utau=UTAU, &
      & Rad=Rad, Fup=FLUP, Fdn=RFLDN, Act=UAVG, Dfdt=DFDT, Crt=ACCUR) 

      if (USRTAU) then
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Ntau-1
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      else
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Nlyr
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      endif

      NPROB = 3
      !IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL &
      !    ERRMSG( 'Out of bounds in exact-answer arrays', .FALSE. )

      CALL PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, NTAU,            &
                   NUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,       &
                   TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),        &
                   TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),        &
                   TSTUU(1,1,1,ICAS,NPROB), MXTAU, MXMU, MXPHI,       &
                   NTEST, NPASS )

      call deallocate_disort_allocatable_arrays()

      deallocate(Rad)
      deallocate(DTAUCP,FLYR,OPRIM)

    enddo         

    call Sjhmfin()
      
  ENDIF

  IF( DOPROB(4) )  THEN

  ! c **********************************************************************
  ! c ****  Test Problem 4:  Haze-L Scattering, Beam Source             ****
  ! c ****  (Compare to Ref. GS, Tables 12-16)                          ****
  ! c **********************************************************************

    USRTAU    = .TRUE.
    USRANG    = .TRUE.
    LAMBER    = .TRUE.
    PLANK     = .FALSE.
    ONLYFL    = .FALSE.
    DO_PSEUDO_SPHERE = .FALSE.
    DELTAMPLUS = .FALSE.

    NSTR = 32; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
    NLYR = 1; 
    NMOM = 32 !NSTR 
    NTAU      = 3; IF(.NOT.USRTAU) NTAU = NLYR + 1
    NUMU      = 6; IF(.NOT.USRANG) NUMU = NSTR
    ! NPHI      = 3; 
    IBCND     = 0

        ! Initialize HSHM parameters
    Nm = NSTR
    Ng = NMOM+1
    Nh = 1
    Lup = 0 !-Nstr/2
    Ldn = 0 !-Nstr/2
    !call Sjhminit(.false.)

    DO 40  ICAS = 1, 3

      WRITE( TITLE, '(3A)' ) 'Test Case No. 4', ABC(ICAS),   &
            ', Haze-L Scattering, Ref. GS, Table '           
            LENTIT = INDEX( TITLE,BLANKS )

      IF ( ICAS.EQ.1 ) THEN

        NPHI = 1
        call Sjhminit(.false.)
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 1.0;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0; 
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 ) = 0.5; UTAU( 3 )  = 1.0
        CALL  GETMOM( 4, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = 1.0; SSALB( 1 ) = 1.0
        FBEAM      = PI ; FISOT      = 0.0

        HEADER = TITLE(1:LENTIT) // ' 12'

      ELSE IF ( ICAS.EQ.2 ) THEN

        NPHI = 1
        call Sjhminit(.false.)
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 1.0;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0; 
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 ) = 0.5; UTAU( 3 )  = 1.0
        CALL  GETMOM( 4, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = 1.0; SSALB( 1 ) = 0.9
        FBEAM      = PI; FISOT      = 0.0

        HEADER = TITLE(1:LENTIT) // ' 13'

      ELSE IF ( ICAS.EQ.3 ) THEN

        NPHI = 3
        call Sjhminit(.false.)
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0;  PHI( 2 )  = 90.0; PHI( 3 )  = 180.0;
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 ) = 0.5; UTAU( 3 )  = 1.0
        CALL  GETMOM( 4, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = 1.0; SSALB( 1 ) = 0.9
        FBEAM      = PI ; FISOT      = 0.0

        HEADER = TITLE(1:LENTIT) // ' 14 - 16'

      ENDIF

      Numu1 = NUMU
      Nphi1 = NPHI
      Umu1 = UMU
      Phi1 = PHI/180.*PI

      GS(:) = 0.0
      GS(0) = ALBEDO / PI


      if (USRTAU) then
        allocate(Rad(Numu1,Nphi1,0:max(Nlyr,NTAU-1)))
      else
        allocate(Rad(Numu1,Nphi1,0:Nlyr))
      endif

      WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' ) ' DISORT: '// trim(HEADER)

      ! print input parameters
      DELTAM     = .FALSE.
      CORINT     = .FALSE.
      DTAUCP = DTAUC
      FLYR = PMOM(NMOM,1:NLYR)
      OPRIM = SSALB
      TAUC(0) = 0.0
      TAUCPR(0) = 0.0
      do LU = 1, NLYR
        TAUC(LU) = TAUC(LU-1) + DTAUC(LU)
        TAUCPR(LU) = TAUCPR(LU-1) + DTAUCP(LU)
      enddo

      CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER, &
     &               WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU, &
     &               NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, &
     &               LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, &
     &               PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT, &
     &               OPRIM, TAUC, TAUCPR, MAXMOM, PRNT( 5 ), &
     &               DO_PSEUDO_SPHERE, H_LYR, DELTAMPLUS)

      call Sjhm(Nlay=NLYR, Od0=DTAUC, Ssa0=SSALB, Pf0=PMOM, &
      & I0=FBEAM, Mu0=UMU0, T0=TEMPER, Tsurf=BTEMP, GS=Gs, &
      & Thermal=PLANK, Deltam=DELTAM, Wnumlo=WVNMLO, Wnumhi=WVNMHI, &
      & Beam0=RFLDIR, Fisot=FISOT, Ttemp=TTEMP, Temis=TEMIS, &
      & Usertau=USRTAU, Nutau=NTAU, Utau=UTAU, &
      & Rad=Rad, Fup=FLUP, Fdn=RFLDN, Act=UAVG, Dfdt=DFDT, Crt=ACCUR)

      if (USRTAU) then
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Ntau-1
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      else
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Nlyr
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      endif

      NPROB = 4
      !IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL &
      !    ERRMSG( 'Out of bounds in exact-answer arrays', .FALSE. )

      CALL PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, NTAU,            &
                   NUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,       &
                   TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),        &
                   TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),        &
                   TSTUU(1,1,1,ICAS,NPROB), MXTAU, MXMU, MXPHI,       &
                   NTEST, NPASS )

      call deallocate_disort_allocatable_arrays()   

      deallocate(Rad)
      deallocate(DTAUCP,FLYR,OPRIM)

      call Sjhmfin()

    40 CONTINUE 

  ENDIF

  IF( DOPROB(5) )  THEN

  !c **********************************************************************
  !c ****  Test Problem 5:  Cloud C.1 Scattering, Beam Source          ****
  !c ****  (Compare to Ref. GS, Tables 19-20)                          ****
  !c **********************************************************************

    USRTAU    = .TRUE.
    USRANG    = .TRUE.
    LAMBER    = .TRUE.
    PLANK     = .FALSE.
    ONLYFL    = .FALSE.
    DO_PSEUDO_SPHERE = .FALSE.
    DELTAMPLUS = .FALSE.

    NSTR = 48; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
    NLYR = 1; 
    NMOM = 299 
    NTAU      = 3; IF(.NOT.USRTAU) NTAU = NLYR + 1
    NUMU      = 6; IF(.NOT.USRANG) NUMU = NSTR
    NPHI      = 1; 
    IBCND     = 0

    ! Initialize HSHM parameters
    Nm = NSTR
    Ng = NMOM+1
    Nh = 1
    Lup = 0 !-Nstr/2
    Ldn = 0 !-Nstr/2


    DO 50  ICAS = 1, 2

      WRITE( TITLE, '(3A)' ) 'Test Case No. 5', ABC(ICAS), &
          ', Cloud C.1 Scattering, Ref. GS, Table '
      LENTIT = INDEX( TITLE,BLANKS )

      IF ( ICAS .EQ. 1 ) THEN

        call Sjhminit(.false.)
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 1.0;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0; 
        ALBEDO    = 0.0
        UTAU( 1 ) = 0.0; UTAU( 2 ) = 32.0; UTAU( 3 )  = 64.0
        CALL  GETMOM( 5, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = 64.0; SSALB( 1 ) = 1.0
        FBEAM      = PI; FISOT      = 0.0       
        HEADER = TITLE(1:LENTIT) // ' 19'

      ELSEIF ( ICAS .EQ. 2 ) THEN

        call Sjhminit(.false.)
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU+Nlyr+1 )

        UMU0      = 1.0;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.5; UMU( 3 )  = -0.1
        UMU( 4 )  =  0.1; UMU( 5 )  =  0.5; UMU( 6 )  =  1.0
        PHI( 1 )  = 0.0; 
        ALBEDO    = 0.0
        UTAU( 1 ) = 3.2; UTAU( 2 ) = 12.8; UTAU( 3 )  = 48.0
        CALL  GETMOM( 5, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = 64.0; SSALB( 1 ) = 0.9
        FBEAM      = PI; FISOT      = 0.0              
        HEADER = TITLE(1:LENTIT) // ' 20'

      ENDIF

      Numu1 = NUMU
      Nphi1 = NPHI
      Umu1 = UMU
      Phi1 = PHI/180.*PI

      GS = 0.0
      GS(0) = ALBEDO / PI

      if (USRTAU) then
        allocate(Rad(Numu1,Nphi1,0:NTAU+Nlyr))
      else
        allocate(Rad(Numu1,Nphi1,0:Nlyr))
      endif

      WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' ) ' DISORT: '// trim(HEADER)

      ! print input parameters
      DELTAM     = .FALSE.
      CORINT     = .FALSE.
      DTAUCP = DTAUC
      FLYR = PMOM(NMOM,1:NLYR)
      OPRIM = SSALB
      TAUC(0) = 0.0
      TAUCPR(0) = 0.0
      do LU = 1, NLYR
        TAUC(LU) = TAUC(LU-1) + DTAUC(LU)
        TAUCPR(LU) = TAUCPR(LU-1) + DTAUCP(LU)
      enddo

      CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER, &
     &               WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU, &
     &               NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, &
     &               LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, &
     &               PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT, &
     &               OPRIM, TAUC, TAUCPR, MAXMOM, PRNT( 5 ), &
     &               DO_PSEUDO_SPHERE, H_LYR, DELTAMPLUS)

      call Sjhm(Nlay=NLYR, Od0=DTAUC, Ssa0=SSALB, Pf0=PMOM, &
      & I0=FBEAM, Mu0=UMU0, T0=TEMPER, Tsurf=BTEMP, GS=Gs, &
      & Thermal=PLANK, Deltam=DELTAM, Wnumlo=WVNMLO, Wnumhi=WVNMHI, &
      & Beam0=RFLDIR, Fisot=FISOT, Ttemp=TTEMP, Temis=TEMIS, &
      & Usertau=USRTAU, Nutau=NTAU, Utau=UTAU, &
      & Rad=Rad, Fup=FLUP, Fdn=RFLDN, Act=UAVG, Dfdt=DFDT, Crt=ACCUR) 

      if (USRTAU) then
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Ntau-1
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      else
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Nlyr
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      endif

      NPROB = 5
      !IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL &
      !    ERRMSG( 'Out of bounds in exact-answer arrays', .FALSE. )

      CALL PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, NTAU,            &
                   NUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,       &
                   TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),        &
                   TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),        &
                   TSTUU(1,1,1,ICAS,NPROB), MXTAU, MXMU, MXPHI,       &
                   NTEST, NPASS )

      call deallocate_disort_allocatable_arrays()

      deallocate(Rad)
      deallocate(DTAUCP,FLYR,OPRIM)

      call Sjhmfin()

    50 CONTINUE 

  ENDIF

  IF( DOPROB(6) )  THEN

  !c **********************************************************************
  !c ****  Test Problem 6:  No Scattering, Increasingly Complex Sources****
  !c **********************************************************************
  ! So far only first 3 case used.

    USRTAU    = .TRUE.
    USRANG    = .TRUE.
    LAMBER    = .TRUE.
    PLANK     = .FALSE.
    ONLYFL    = .FALSE.
    DO_PSEUDO_SPHERE = .FALSE.
    DELTAMPLUS = .FALSE.

    NSTR = 16; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
    NLYR = 1; 
    NMOM = NSTR 
    !NTAU      = 1; IF(.NOT.USRTAU) NTAU = NLYR + 1
    NUMU      = 4; IF(.NOT.USRANG) NUMU = NSTR
    NPHI      = 1; 
    IBCND     = 0

    ! Initialize HSHM parameters
    Nm = NSTR
    Ng = NMOM+1
    Nh = 1
    Lup = 0 !-Nstr/2
    Ldn = 0 !-Nstr/2

    DO 60  ICAS = 1, 3 !8

      WRITE( TITLE, '(3A)' ) 'Test Case No. 6', ABC(ICAS),  &
             ': No Scattering; Source = Beam'
      LENTIT = INDEX( TITLE, BLANKS )
    
      IF ( ICAS.EQ.1 ) THEN
      !c                                    ** Transparent medium, beam source   
        NTAU = 2
        call Sjhminit(.false.)
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU+Nlyr+1 )

        UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.1
        UMU( 3 )  =  0.1; UMU( 4 )  =  1.0
        PHI( 1 )  = 90.0; 
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 0.0; 
        !CALL  GETMOM( 5, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = 0.0;   SSALB( 1 ) = 0.0
        FBEAM      = 200.0; FISOT      = 0.0     
        ALBEDO    = 0.0         
        HEADER = TITLE(1:LENTIT) // '; Bottom Albedo = 0'

      ELSE IF ( ICAS.EQ.2 ) THEN
      !c                                    ** Add some optical depth       
        NTAU = 3 

        call Sjhminit(.false.)
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU+Nlyr+1 )

        UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.1
        UMU( 3 )  =  0.1; UMU( 4 )  =  1.0
        PHI( 1 )  = 90.0; 
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 0.5;  UTAU( 3 ) = 1.0; 
        !CALL  GETMOM( 5, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = 1.0;   SSALB( 1 ) = 0.0
        FBEAM      = 200.0; FISOT      = 0.0       
        ALBEDO    = 0.0       
        HEADER = TITLE(1:LENTIT) // '; Bottom Albedo = 0'       

      ELSE IF ( ICAS.EQ.3 ) THEN
      !c                                   ** Add some isotropic reflection       
        NTAU = 3 

        call Sjhminit(.false.)
        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU+Nlyr+1 )

        UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.1
        UMU( 3 )  =  0.1; UMU( 4 )  =  1.0
        PHI( 1 )  = 90.0; 
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 0.5;  UTAU( 3 ) = 1.0; 
        !CALL  GETMOM( 5, 0.0, NMOM, PMOM )
        DTAUC( 1 ) = 1.0;   SSALB( 1 ) = 0.0
        FBEAM      = 200.0; FISOT      = 0.0     
        ALBEDO    = 0.5         
        HEADER = TITLE(1:LENTIT) // '; Bottom Albedo = 0.5 lambert'   

      endif

      Numu1 = NUMU
      Nphi1 = NPHI
      Umu1 = UMU
      Phi1 = PHI/180.*PI

      if (USRTAU) then
        allocate(Rad(Numu1,Nphi1,0:NTAU+Nlyr))
      else
        allocate(Rad(Numu1,Nphi1,0:Nlyr))
      endif

      WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' ) ' DISORT: '// trim(HEADER)

      ! print input parameters
      DELTAM     = .FALSE.
      CORINT     = .FALSE.
      DTAUCP = DTAUC
      FLYR = PMOM(NMOM,1:NLYR)
      OPRIM = SSALB
      TAUC(0) = 0.0
      TAUCPR(0) = 0.0
      do LU = 1, NLYR
        TAUC(LU) = TAUC(LU-1) + DTAUC(LU)
        TAUCPR(LU) = TAUCPR(LU-1) + DTAUCP(LU)
      enddo

      GS(:) = 0.0
      GS(0) = ALBEDO / PI

      CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER, &
     &               WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU, &
     &               NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, &
     &               LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, &
     &               PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT, &
     &               OPRIM, TAUC, TAUCPR, MAXMOM, PRNT( 5 ), &
     &               DO_PSEUDO_SPHERE, H_LYR, DELTAMPLUS)

      call Sjhm(Nlay=NLYR, Od0=DTAUC, Ssa0=SSALB, Pf0=PMOM, &
      & I0=FBEAM, Mu0=UMU0, T0=TEMPER, Tsurf=BTEMP, GS=Gs, &
      & Thermal=PLANK, Deltam=DELTAM, Wnumlo=WVNMLO, Wnumhi=WVNMHI, &
      & Beam0=RFLDIR, Fisot=FISOT, Ttemp=TTEMP, Temis=TEMIS, &
      & Usertau=USRTAU, Nutau=NTAU, Utau=UTAU, &
      & Rad=Rad, Fup=FLUP, Fdn=RFLDN, Act=UAVG, Dfdt=DFDT, Crt=ACCUR) 

      if (USRTAU) then
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Ntau-1
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
        !DFDT = UAVG * (1.0 - SSALB(1))
      else
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Nlyr
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
        !DFDT = UAVG * (1.0 - SSALB(1))
      endif

      NPROB = 6
      !IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL &
      !  ERRMSG( 'Out of bounds in exact-answer arrays', .FALSE. )

      CALL PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, NTAU,            &
                 NUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,       &
                 TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),        &
                 TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),        &
                 TSTUU(1,1,1,ICAS,NPROB), MXTAU, MXMU, MXPHI,       &
                 NTEST, NPASS )
                    
      call deallocate_disort_allocatable_arrays()

      deallocate(Rad)
      deallocate(DTAUCP,FLYR,OPRIM)

      call Sjhmfin()

    60 CONTINUE 

  endif

  IF( DOPROB(7) )  THEN

  !c **********************************************************************
  !c ****  Test Problem 7:  Absorption + Scattering + All Possible     ****
  !c ****  Sources, Various Surface Reflectivities ( One Layer )       ****
  !c **** (Compare 7a,f Fluxes and Intensities to Ref. KS, Tables I-II ****
  !c **********************************************************************
  ! So far only first 4 cases used

    DO 70 ICAS = 1, 4 !5

      WRITE( TITLE, '(2A)' ) 'Test Case No. 7', ABC(ICAS)
      LENTIT = INDEX( TITLE, BLANKS )

      IF ( ICAS.EQ.1 ) THEN

        USRTAU    = .TRUE.
        USRANG    = .TRUE.
        LAMBER    = .TRUE.
        PLANK     = .FALSE.
        ONLYFL    = .TRUE.
        DO_PSEUDO_SPHERE = .FALSE.
        DELTAMPLUS = .FALSE.

        NSTR = 16; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
        NLYR = 1; 
        NMOM = NSTR 
        NTAU      = 2; IF(.NOT.USRTAU) NTAU = NLYR + 1
        NUMU      = 2; IF(.NOT.USRANG) NUMU = NSTR
        NPHI      = 1; 
        IBCND     = 0     

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU+Nlyr+1 )
        Numu1 = NUMU
        Nphi1 = NPHI
        Umu1 = UMU
        Phi1 = PHI/180.*PI

        !UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = 1.0
        PHI( 1 )  = 0.0; 
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 1.0; 
        CALL  GETMOM( 3, 0.05, NMOM, PMOM )
        DTAUC( 1 ) = 1.0;   SSALB( 1 ) = 0.1
        FBEAM      = 0.0;   FISOT      = 0.0
        PLANK      = .TRUE.; BTEMP      = 0.0; TTEMP  =  0.0; TEMIS      = 1.0
        TEMPER( 0 ) = 200.0;  TEMPER( 1 ) = 300.0;  WVNMLO = 300.0; WVNMHI = 800.0      
        LAMBER     = .TRUE.; ALBEDO   = 0.0; 
        HEADER = TITLE(1:LENTIT) // ': Absorption + Scattering, '//  &
                 'Internal Thermal Sources; Ref. KS, Table I, '//    &
                 'tau = 1.0, a = 0.1, g = 0.05'

      ELSEIF ( ICAS.EQ.2 ) THEN

        USRTAU    = .TRUE.
        USRANG    = .TRUE.
        LAMBER    = .TRUE.
        PLANK     = .FALSE.
        ONLYFL    = .FALSE.
        DO_PSEUDO_SPHERE = .FALSE.
        DELTAMPLUS = .FALSE.

        NSTR = 16; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
        NLYR = 1; 
        NMOM = NSTR 
        NTAU      = 2; IF(.NOT.USRTAU) NTAU = NLYR + 1
        NUMU      = 2; IF(.NOT.USRANG) NUMU = NSTR
        NPHI      = 1;
        IBCND     = 0     

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU+Nlyr+1 )

        !UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = 1.0
        PHI( 1 )  = 0.0; 
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 100.0; 
        CALL  GETMOM( 3, 0.75, NMOM, PMOM )
        DTAUC( 1 ) = 100.0;   SSALB( 1 ) = 0.95
        FBEAM      = 0.0;   FISOT      = 0.0
        PLANK      = .TRUE.; BTEMP      = 0.0; TTEMP  =  0.0; TEMIS      = 1.0
        TEMPER( 0 ) = 200.0;  TEMPER( 1 ) = 300.0;  WVNMLO = 2702.99; WVNMHI = 2703.01      
        LAMBER     = .TRUE.; ALBEDO   = 0.0; 
        HEADER = TITLE(1:LENTIT) // ': Absorption + Scattering, '//  &
                 'Internal Thermal Sources; Ref. KS, Table I, '//    &
                 'tau = 100.0, a = 0.95, g = 0.75'

      ELSEIF ( ICAS.EQ.3 ) THEN

        USRTAU    = .TRUE.
        USRANG    = .TRUE.
        LAMBER    = .TRUE.
        PLANK     = .FALSE.
        ONLYFL    = .FALSE.
        DO_PSEUDO_SPHERE = .FALSE.
        DELTAMPLUS = .FALSE.

        NSTR = 12
        IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
        NLYR = 1; 
        NMOM = NSTR 
        NTAU      = 3; IF(.NOT.USRTAU) NTAU = NLYR + 1
        NUMU      = 4; IF(.NOT.USRANG) NUMU = NSTR
        NPHI      = 2; 
        IBCND     = 0     

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU)

        UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.1; UMU( 3 )  = 0.1; UMU( 4 )  = 1.0; 
        PHI( 1 )  = 0.0; PHI( 2 )  = 90.0; 
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 0.5;  UTAU( 3 ) = 1.0; 
        CALL  GETMOM( 3, 0.8, NMOM, PMOM )
        DTAUC( 1 ) = 1.0;   SSALB( 1 ) = 0.5
        FBEAM      = 200.0
        FISOT      = 100.0
        PLANK      = .TRUE.
        BTEMP      = 320.0
        TTEMP  =  100.0; TEMIS      = 1.0
        TEMPER( 0 ) = 300.0;  TEMPER( 1 ) = 200.0;  WVNMLO = 0.0; WVNMHI = 80000.0      
        LAMBER     = .TRUE.; ALBEDO   = 0.0;                 
        HEADER = TITLE(1:LENTIT) // ': Absorption + '//            &
                 'Henyey-Greenstein Scattering, All Sources, '//   &
                 'Bottom Albedo = 0'

      ELSEIF ( ICAS.EQ.4 ) THEN

        USRTAU    = .TRUE.
        USRANG    = .TRUE.
        LAMBER    = .TRUE.
        PLANK     = .FALSE.
        ONLYFL    = .FALSE.
        DO_PSEUDO_SPHERE = .FALSE.
        DELTAMPLUS = .FALSE.

        NSTR = 12; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
        NLYR = 1; 
        NMOM = NSTR 
        NTAU      = 3; IF(.NOT.USRTAU) NTAU = NLYR + 1
        NUMU      = 4; IF(.NOT.USRANG) NUMU = NSTR
        NPHI      = 2; 
        IBCND     = 0     

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.1; UMU( 3 )  = 0.1; UMU( 4 )  = 1.0; 
        PHI( 1 )  = 0.0; PHI( 2 )  = 90.0; 
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 0.5;  UTAU( 3 ) = 1.0; 
        CALL  GETMOM( 3, 0.8, NMOM, PMOM )
        DTAUC( 1 ) = 1.0;   SSALB( 1 ) = 0.5
        FBEAM      = 200.0;   FISOT      = 100.0
        PLANK      = .TRUE.; BTEMP      = 320.0; TTEMP  =  100.0; TEMIS      = 1.0
        TEMPER( 0 ) = 300.0;  TEMPER( 1 ) = 200.0;  WVNMLO = 0.0; WVNMHI = 80000.0      
        LAMBER     = .TRUE.; ALBEDO   = 1.0;                 
        HEADER = TITLE(1:LENTIT) // ': Absorption + '//            &
                 'Henyey-Greenstein Scattering, All Sources, '//   &
                 'Bottom Albedo = 1'

      endif

      ! Initialize HSHM parameters
      Nm = NSTR
      Ng = NMOM+1
      Nh = 1
      Lup = 0 !-Nstr/2
      Ldn = 0 !-Nstr/2

      Numu1 = NUMU
      Nphi1 = NPHI

      call Sjhminit(.false.)

      Umu1 = UMU
      Phi1 = PHI/180.*PI

      if (USRTAU) then
        allocate(Rad(Numu1,Nphi1,0:NTAU+Nlyr))
      else
        allocate(Rad(Numu1,Nphi1,0:Nlyr))
      endif

      WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' ) ' DISORT: '// trim(HEADER)

      ! print input parameters
      DELTAM     = .FALSE.
      CORINT     = .FALSE.
      DTAUCP = DTAUC
      FLYR = PMOM(NMOM,1:NLYR)
      OPRIM = SSALB
      TAUC(0) = 0.0
      TAUCPR(0) = 0.0
      do LU = 1, NLYR
        TAUC(LU) = TAUC(LU-1) + DTAUC(LU)
        TAUCPR(LU) = TAUCPR(LU-1) + DTAUCP(LU)
      enddo

      GS(:) = 0.0
      GS(0) = ALBEDO / PI

      !PLANK = .false.
      !FBEAM = 0.; FISOT=0.0; TTEMP = 0.0
      !TEMPER = 0.0
      !BTEMP = 0.0

      CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER, &
     &               WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU, &
     &               NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, &
     &               LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, &
     &               PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT, &
     &               OPRIM, TAUC, TAUCPR, MAXMOM, PRNT( 5 ), &
     &               DO_PSEUDO_SPHERE, H_LYR, DELTAMPLUS)

      call Sjhm(Nlay=NLYR, Od0=DTAUC, Ssa0=SSALB, Pf0=PMOM, &
      & I0=FBEAM, Mu0=UMU0, T0=TEMPER, Tsurf=BTEMP, GS=Gs, &
      & Thermal=PLANK, Deltam=DELTAM, Wnumlo=WVNMLO, Wnumhi=WVNMHI, &
      & Beam0=RFLDIR, Fisot=FISOT, Ttemp=TTEMP, Temis=TEMIS, &
      & Usertau=USRTAU, Nutau=NTAU, Utau=UTAU, &
      & Rad=Rad, Fup=FLUP, Fdn=RFLDN, Act=UAVG, Dfdt=DFDT, Crt=ACCUR)

      if (USRTAU) then
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Ntau-1
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      else
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Nlyr
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      endif

      NPROB = 7
      !IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL &
      !    ERRMSG( 'Out of bounds in exact-answer arrays', .FALSE. )

      CALL PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, NTAU,            &
                   NUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,       &
                   TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),        &
                   TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),        &
                   TSTUU(1,1,1,ICAS,NPROB), MXTAU, MXMU, MXPHI,       &
                   NTEST, NPASS )

      call deallocate_disort_allocatable_arrays()

      deallocate(Rad)
      deallocate(DTAUCP,FLYR,OPRIM)

      call Sjhmfin()

    70 enddo

    deallocate(TAUC,TAUCPR)

  ENDIF

  IF( DOPROB(8) )  THEN

  ! c **********************************************************************
  ! c ****  Test Problem 8:  Absorbing/Isotropic-Scattering Medium      ****
  ! c ****  With Two Computational Layers                               ****
  ! c **** (Compare Fluxes To Ref. OS, Table 1)                         ****
  ! c **********************************************************************

    USRTAU    = .TRUE.
    USRANG    = .TRUE.
    LAMBER    = .TRUE.
    PLANK     = .FALSE.
    ONLYFL    = .FALSE.
    DO_PSEUDO_SPHERE = .FALSE.
    DELTAMPLUS = .FALSE.

    NSTR = 8; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
    NLYR = 2; 
    NMOM = NSTR 
    NTAU      = 3; IF(.NOT.USRTAU) NTAU = NLYR + 1
    NUMU      = 4; IF(.NOT.USRANG) NUMU = NSTR
    NPHI      = 1; 
    IBCND     = 0  

    DO 80 ICAS = 1, 3

      WRITE( TITLE, '(2A)' ) 'Test Case No. 7', ABC(ICAS)
      LENTIT = INDEX( TITLE, BLANKS )

      IF ( ICAS.EQ.1 ) THEN

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        !UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.2; UMU( 3 )  = 0.2; UMU( 4 )  = 1.0; 
        PHI( 1 )  = 60.0;
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 0.25;  UTAU( 3 ) = 0.5; 
        CALL  GETMOM( 1, 0.0, NMOM, PMOM(:,1) )
        CALL  GETMOM( 1, 0.0, NMOM, PMOM(:,2) )
        DTAUC( 1 ) = 0.25;  DTAUC( 2 ) = 0.25;   
        SSALB( 1 ) = 0.5;   SSALB( 2 ) = 0.3;
        FBEAM      = 0.0;   FISOT      = 1.0/PI
        PLANK      = .FALSE.;   
        LAMBER     = .TRUE.; ALBEDO   = 0.0;                 

        HEADER = 'Test Case No. 8a:  Ref. OS, Table 1,'   &
                  // ' Line 4 (Two Inhomogeneous Layers)'

      ELSEIF ( ICAS.EQ.2 ) THEN

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )
        !UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.2; UMU( 3 )  = 0.2; UMU( 4 )  = 1.0; 
        PHI( 1 )  = 60.0;
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 0.25;  UTAU( 3 ) = 0.5; 
        CALL  GETMOM( 1, 0.0, NMOM, PMOM(:,1) )
        CALL  GETMOM( 1, 0.0, NMOM, PMOM(:,2) )
        DTAUC( 1 ) = 0.25;  DTAUC( 2 ) = 0.25;   
        SSALB( 1 ) = 0.8;   SSALB( 2 ) = 0.95;
        FBEAM      = 0.0;   FISOT      = 1.0/PI
        PLANK      = .FALSE.;   
        LAMBER     = .TRUE.; ALBEDO   = 0.0;                 

        HEADER = 'Test Case No. 8b:  Ref. OS, Table 1,'   &
                  // ' Line 1 (Two Inhomogeneous Layers)'

      ELSEIF ( ICAS.EQ.3 ) THEN

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

        !UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.2; UMU( 3 )  = 0.2; UMU( 4 )  = 1.0; 
        PHI( 1 )  = 60.0;
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 1.0;  UTAU( 3 ) = 3.0; 
        CALL  GETMOM( 1, 0.0, NMOM, PMOM(:,1) )
        CALL  GETMOM( 1, 0.0, NMOM, PMOM(:,2) )
        DTAUC( 1 ) = 1.0;   DTAUC( 2 ) = 2.0;   
        SSALB( 1 ) = 0.8;   SSALB( 2 ) = 0.95;
        FBEAM      = 0.0;   FISOT      = 1.0/PI
        PLANK      = .FALSE.;   
        LAMBER     = .TRUE.; ALBEDO   = 0.0;                 

        HEADER = 'Test Case No. 8c:  Ref. OS, Table 1,'   &
                  // ' Line 13 (Two Inhomogeneous Layers)'

      ENDIF

      ! Initialize HSHM parameters
      Nm = NSTR
      Ng = NMOM+1
      Nh = 1
      Lup = 0 !-Nstr/2
      Ldn = 0 !-Nstr/2

      Numu1 = NUMU
      Nphi1 = NPHI

      call Sjhminit(.false.)

      Umu1 = UMU
      Phi1 = PHI/180.*PI

      if (USRTAU) then
        allocate(Rad(Numu1,Nphi1,0:NTAU+Nlyr))
      else
        allocate(Rad(Numu1,Nphi1,0:Nlyr))
      endif

      WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' ) ' DISORT: '// trim(HEADER)

      ! print input parameters

      allocate(DTAUCP(NLYR),FLYR(NLYR),OPRIM(NLYR),TAUC(0:NLYR),TAUCPR(0:NLYR))

      DELTAM     = .FALSE.
      CORINT     = .FALSE.
      DTAUCP = DTAUC
      FLYR = PMOM(NMOM,1:NLYR)
      OPRIM = SSALB
      TAUC(0) = 0.0
      TAUCPR(0) = 0.0
      do LU = 1, NLYR
        TAUC(LU) = TAUC(LU-1) + DTAUC(LU)
        TAUCPR(LU) = TAUCPR(LU-1) + DTAUCP(LU)
      enddo

      GS(:) = 0.0
      GS(0) = ALBEDO / PI

      !PLANK = .false.
      !FBEAM = 0.; FISOT=0.0; TTEMP = 0.0
      !TEMPER = 0.0
      !BTEMP = 0.0

      CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER, &
     &               WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU, &
     &               NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, &
     &               LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, &
     &               PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT, &
     &               OPRIM, TAUC, TAUCPR, MAXMOM, PRNT( 5 ), &
     &               DO_PSEUDO_SPHERE, H_LYR, DELTAMPLUS)

      call Sjhm(Nlay=NLYR, Od0=DTAUC, Ssa0=SSALB, Pf0=PMOM, &
      & I0=FBEAM, Mu0=UMU0, T0=TEMPER, Tsurf=BTEMP, GS=Gs, &
      & Thermal=PLANK, Deltam=DELTAM, Wnumlo=WVNMLO, Wnumhi=WVNMHI, &
      & Beam0=RFLDIR, Fisot=FISOT, Ttemp=TTEMP, Temis=TEMIS, &
      & Usertau=USRTAU, Nutau=NTAU, Utau=UTAU, &
      & Rad=Rad, Fup=FLUP, Fdn=RFLDN, Act=UAVG, Dfdt=DFDT, Crt=ACCUR)

      if (USRTAU) then
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Ntau-1
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      else
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Nlyr
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      endif

      NPROB = 8
      !IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL &
      !    ERRMSG( 'Out of bounds in exact-answer arrays', .FALSE. )

      CALL PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, NTAU,            &
                   NUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,       &
                   TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),        &
                   TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),        &
                   TSTUU(1,1,1,ICAS,NPROB), MXTAU, MXMU, MXPHI,       &
                   NTEST, NPASS )

      call deallocate_disort_allocatable_arrays() 

      deallocate(Rad)
      deallocate(DTAUCP,FLYR,OPRIM,TAUC,TAUCPR)

      call Sjhmfin()

    80 CONTINUE

  endif

  IF( DOPROB(9) )  THEN

  ! c **********************************************************************
  ! c ****  Test Problem 9:  General Emitting/Absorbing/Scattering      ****
  ! c ****  Medium with Every Computational Layer Different.            ****
  ! c **** (Compare 9a,b Fluxes to Ref. DGIS, Tables VI-VII, beta = 0)  ****
  ! c **********************************************************************

    DO 90 ICAS = 1, 3     

      IF ( ICAS.EQ.1 ) THEN

        USRTAU    = .TRUE.
        USRANG    = .TRUE.
        LAMBER    = .TRUE.
        PLANK     = .FALSE.
        ONLYFL    = .FALSE.
        DO_PSEUDO_SPHERE = .FALSE.
        DELTAMPLUS = .FALSE.

        NSTR = 8; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
        NLYR = 6; 
        NMOM = NSTR 
        NTAU      = 5; IF(.NOT.USRTAU) NTAU = NLYR + 1
        NUMU      = 4; IF(.NOT.USRANG) NUMU = NSTR
        NPHI      = 1; 
        IBCND     = 0     

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU+NLYR )

        !UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.2; UMU( 3 ) = 0.2; UMU( 4 ) = 1.0;
        PHI( 1 )  = 60.0; 
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 1.05; UTAU( 3 ) = 2.1; UTAU( 4 ) = 6.0; UTAU( 5 ) = 21.0; 
        DO LC = 1, NLYR
          CALL  GETMOM( 1, 0.0, NMOM, PMOM(:,LC) )
        ENDDO       
        DO LC = 1, NLYR
          DTAUC( LC ) = LC
          SSALB( LC ) = 0.6 + LC*0.05
        ENDDO
        FBEAM      = 0.0;   FISOT      = 1.0/PI
        PLANK      = .FALSE.;     
        LAMBER     = .TRUE.; ALBEDO   = 0.0; 
        HEADER = 'Test Case No. 9a:  Ref. DGIS, Tables VI-VII,'  &
                 // ' beta=l=0 (multiple inhomogeneous layers)'

      ELSEIF ( ICAS .EQ. 2 ) THEN

        USRTAU    = .TRUE.
        USRANG    = .TRUE.
        LAMBER    = .TRUE.
        PLANK     = .FALSE.
        ONLYFL    = .FALSE.
        DO_PSEUDO_SPHERE = .FALSE.
        DELTAMPLUS = .FALSE.

        NSTR = 8; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
        NLYR = 6; 
        NMOM = NSTR 
        NTAU      = 5; IF(.NOT.USRTAU) NTAU = NLYR + 1
        NUMU      = 4; IF(.NOT.USRANG) NUMU = NSTR
        NPHI      = 1; 
        IBCND     = 0     

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU+NLYR )

        !UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.2; UMU( 3 ) = 0.2; UMU( 4 ) = 1.0;
        PHI( 1 )  = 60.0; 
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 1.05; UTAU( 3 ) = 2.1; UTAU( 4 ) = 6.0; UTAU( 5 ) = 21.0; 
        PMOM(0,1) = 1.0
        PMOM(1,1) = 2.00916/3.
        PMOM(2,1) = 1.56339/5.
        PMOM(3,1) = 0.67407/7.
        PMOM(4,1) = 0.22215/9.
        PMOM(5,1) = 0.04725/11.
        PMOM(6,1) = 0.00671/13.
        PMOM(7,1) = 0.00068/15.
        PMOM(8,1) = 0.00005/17.
        DO LC = 2, NLYR
          DO K = 0, 8
            PMOM(K,LC) = PMOM(K,1)
          ENDDO
        ENDDO   
        DO LC = 1, NLYR
          DTAUC( LC ) = LC
          SSALB( LC ) = 0.6 + LC*0.05
        ENDDO
        FBEAM      = 0.0;   FISOT      = 1.0/PI
        PLANK      = .FALSE.;     
        LAMBER     = .TRUE.; ALBEDO   = 0.0; 
        HEADER = 'Test Case No. 9b:  Ref. DGIS, Tables VI-VII,'    &
                 // ' beta=0,l=8 (multiple inhomogeneous layers)'

      ELSEIF ( ICAS .EQ. 3 ) THEN

        USRTAU    = .TRUE.
        USRANG    = .TRUE.
        LAMBER    = .TRUE.
        PLANK     = .FALSE.
        ONLYFL    = .FALSE.
        DO_PSEUDO_SPHERE = .FALSE.
        DELTAMPLUS = .FALSE.

        NSTR = 8; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
        NLYR = 6; 
        NMOM = NSTR 
        NTAU      = 5; IF(.NOT.USRTAU) NTAU = NLYR + 1
        NUMU      = 4; IF(.NOT.USRANG) NUMU = NSTR
        NPHI      = 3; 
        IBCND     = 0     

        call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU+NLYR )

        UMU0      = 0.5;  PHI0      = 0.0
        UMU( 1 )  = -1.0; UMU( 2 )  = -0.2; UMU( 3 ) = 0.2; UMU( 4 ) = 1.0;
        PHI( 1 )  = 60.0; PHI( 2 ) = 120.0; PHI( 3 ) = 180.0
        UTAU( 1 ) = 0.0;  UTAU( 2 ) = 1.05; UTAU( 3 ) = 2.1; UTAU( 4 ) = 6.0; UTAU( 5 ) = 21.0; 
        DO LC = 1, NLYR 
          CALL  GETMOM( 3, FLOAT(LC)/7.0, NMOM, PMOM(:,LC) )
        ENDDO 
        DO LC = 1, NLYR
          DTAUC( LC ) = LC
          SSALB( LC ) = 0.6 + LC*0.05
        ENDDO
        FBEAM      = PI;   FISOT      = 1.0
        PLANK      = .TRUE.; BTEMP      = 700.0; TTEMP  =  550.0; TEMIS      = 1.0
        WVNMLO = 999.0; WVNMHI = 1000.0  
        TEMPER( 0 ) = 600.0; 
        DO LC = 1, NLYR
          TEMPER( LC ) = 600.0 + LC*10.0
        ENDDO    
        LAMBER     = .TRUE.; ALBEDO   = 0.5; 
        HEADER = 'Test Case No. 9c:  Generalization of 9A '//    & 
                 'to include all possible complexity'

      ENDIF

      ! Initialize HSHM parameters
      Nm = NSTR
      Ng = NMOM+1
      Nh = 1
      Lup = 0 !-Nstr/2
      Ldn = 0 !-Nstr/2

      Numu1 = NUMU
      Nphi1 = NPHI

      call Sjhminit(.false.)

      Umu1 = UMU
      Phi1 = PHI/180.*PI

      if (USRTAU) then
        allocate(Rad(Numu1,Nphi1,0:NTAU+Nlyr))
      else
        allocate(Rad(Numu1,Nphi1,0:Nlyr))
      endif

      WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' ) ' DISORT: '// trim(HEADER)

      ! print input parameters
      allocate(DTAUCP(NLYR),FLYR(NLYR),OPRIM(NLYR),TAUC(0:NLYR),TAUCPR(0:NLYR))

      DELTAM     = .FALSE.
      CORINT     = .FALSE.
      DTAUCP = DTAUC
      FLYR = PMOM(NMOM,1:NLYR)
      OPRIM = SSALB
      TAUC(0) = 0.0
      TAUCPR(0) = 0.0
      do LU = 1, NLYR
        TAUC(LU) = TAUC(LU-1) + DTAUC(LU)
        TAUCPR(LU) = TAUCPR(LU-1) + DTAUCP(LU)
      enddo

      GS(:) = 0.0
      GS(0) = ALBEDO / PI

      !PLANK = .false.
      !FBEAM = 0.; FISOT=0.0; TTEMP = 0.0
      !TEMPER = 0.0
      !BTEMP = 0.0

      CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, NMOM, PMOM, TEMPER, &
     &               WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU, &
     &               NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, &
     &               LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, &
     &               PLANK, ONLYFL, CORINT, ACCUR, FLYR, LYRCUT, &
     &               OPRIM, TAUC, TAUCPR, MAXMOM, PRNT( 5 ), &
     &               DO_PSEUDO_SPHERE, H_LYR, DELTAMPLUS)

      call Sjhm(Nlay=NLYR, Od0=DTAUC, Ssa0=SSALB, Pf0=PMOM, &
      & I0=FBEAM, Mu0=UMU0, T0=TEMPER, Tsurf=BTEMP, GS=Gs, &
      & Thermal=PLANK, Deltam=DELTAM, Wnumlo=WVNMLO, Wnumhi=WVNMHI, &
      & Beam0=RFLDIR, Fisot=FISOT, Ttemp=TTEMP, Temis=TEMIS, &
      & Usertau=USRTAU, Nutau=NTAU, Utau=UTAU, &
      & Rad=Rad, Fup=FLUP, Fdn=RFLDN, Act=UAVG, Dfdt=DFDT, Crt=ACCUR)

      deallocate(UU)
      allocate(UU(NUMU, NTAU, NPHI))

      if (USRTAU) then
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Ntau-1
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      else
        RFLDIR = RFLDIR * UMU0
        do Lev = 0, Nlyr
          UU(1:Numu1,Lev+1,1:Nphi1) = Rad(1:Numu1,1:Nphi1,Lev)
        enddo
      endif

      NPROB = 9
      !IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL &
      !    ERRMSG( 'Out of bounds in exact-answer arrays', .FALSE. )

      CALL PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, NTAU,            &
                   NUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,       &
                   TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),        &
                   TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),        &
                   TSTUU(1,1,1,ICAS,NPROB), MXTAU, MXMU, MXPHI,       &
                   NTEST, NPASS )

      call deallocate_disort_allocatable_arrays()

      deallocate(Rad)
      deallocate(DTAUCP,FLYR,OPRIM,TAUC,TAUCPR)

      call Sjhmfin()



    90 CONTINUE

  ENDIF


end program Sjhmtest