module Gdrtesv0

  implicit none

  private
  ! those may be used by user
  public:: Kreal, Kreal2, Kint, Pi, Pi_d, Gdrtes, Identity, Invmx, Amirror


  ! digital precision constant
  integer, parameter:: Kreal1 = selected_real_kind (P = 6, R = 37) ! single 
  integer, parameter:: Kreal2 = selected_real_kind (P = 13, R = 200) ! double
  integer, parameter:: Kint = selected_int_kind (9) ! integer

  integer, parameter:: Kreal = Kreal1 ! or Kreal2 for double precision

  real(Kreal2), parameter:: Pi_d = acos(-1.0_Kreal2)
  real(Kreal), parameter:: Pi = acos(-1.0_Kreal)

  logical:: Amirror

  ! interfaces

  interface Abba
    module procedure Abba_i,Abba_r,Abba_rv,Abba_d, Abba_dv
  end interface

  interface Cab
    module procedure Cab_r, Cab_d
  end interface

  interface Fbsub
    module procedure Fbsub_r, Fbsub_d
  end interface

  interface Identity
    module procedure Identity_r, Identity_d
  end interface

  interface Indmax1
    module procedure Indmax1_r,Indmax1_d
  end interface

  interface Invmx
    module procedure Invmx_r, Invmx_d
  end interface

  interface Ludcp
    module procedure Ludcp_r, Ludcp_d
  end interface

contains

!-----------------------------------------------------------------------

subroutine Gdrtes(Nlyr, Nup, Ndn, Ngm, Nhm, AA, AAinv, UUinv, Odd, &
  & Ssad, Glm, PPm, Beam, Mu0, Pmu0, Thermal, ATint, B0, PPf, &
  & Hm, Bsurf, Itild)
  !use:
  implicit none

  integer(Kint),intent(in):: Nlyr, Nup, Ndn, Ngm, Nhm
  real(Kreal2),intent(in):: AA(Nup+Ndn,Nup+Ndn), AAinv(Nup+Ndn,Nup+Ndn), UUinv(Nup+Ndn,Nup+Ndn),&
  & Glm(Ngm,Nlyr), PPm(max(Ngm,Nhm),Nup+Ndn), ATint(Nup+Ndn), PPf(Nhm,Ndn), Pmu0(max(Ngm,Nhm))
  real(Kreal),intent(in):: Odd(Nlyr), Ssad(Nlyr), &
  & Beam(0:Nlyr), Mu0, B0(0:Nlyr), Hm(Nhm,Nhm), Bsurf
  logical, intent(in):: Thermal
  real(Kreal),intent(out):: Itild(Nup+Ndn,0:Nlyr)

  integer(Kint):: Nstr, Nd1

  real(Kreal2):: Mu0d
  real(Kreal2), dimension(Nup+Ndn,Nup+Ndn):: App, Dp, Gout, Gin, D, D1
  real(Kreal2), dimension(Nup+Ndn):: Dtred2, Z, S0, EE
  real(Kreal), dimension(Nup+Ndn,Nup+Ndn, Nlyr+1):: S
  real(Kreal), dimension(Nup+Ndn):: Yt, Y1, Y0, Bpm
  real(Kreal), dimension(Nup+Ndn,Nlyr+1):: Xpm
  real(Kreal), dimension(Nup,Ndn):: Rs
  real(Kreal), dimension(Nup):: Sigsup
  
  !For solving band diagonal system
  real(Kreal2):: Gt(Nup+Ndn,Nup+Ndn,Nlyr), Gb(Nup+Ndn,Nup+Ndn,Nlyr)
  real(Kreal):: Qt(Nup+Ndn,Nlyr), Qb(Nup+Ndn,Nlyr)


  integer(Kint):: l, i, j
  real(Kreal)::  Dump0 !, Dump1(Nlyr)

  Nd1 = Nup + 1
  Nstr = Nup + Ndn

  Mu0d = Mu0

  ! Loop for layers
  do l = 1, Nlyr
    
    App = matmul(transpose(PPm(1:Ngm,:)) * spread(Glm(:,l),dim=1,ncopies=Nstr), PPm(1:Ngm,:))
    App = AA - Ssad(l) * App
    D(:,:) = matmul(UUinv, App)

    Amirror = .true.
    if (Amirror) then
      ! use Asymtx for hemisphere-confined and asymmetric base
      call Deigendst(D, Dp, Dtred2, Nstr)
      !print*, 'Dtred2a', Dtred2(1:Nstr)
      !print*, Dp(1:Nstr,1:Nstr)
    else
      ! use Cholesky decomposition for general case
      call Deigen1d(App, UUinv, Dp, Dtred2)
      !print*, 'Dtred2b', Dtred2(1:Nstr)
      !print*, Dp(1:Nstr,1:Nstr)
    endif
    !stop
    
    do j = 1, Nstr
      if (Dtred2(j).gt.0.0_Kreal2) then
        Dump0 = - Odd(l)
      else
        Dump0 = 0.0_Kreal2
      endif
      do i = 1, Nup
        Gout(i,j) = Dp(i,j) * exp(Dtred2(j) * Dump0)
        Gin(i,j) = Dp(i,j) * exp(Dtred2(j) * (Odd(l) + Dump0))

        Gt(i,j,l) = Gout(i,j)
        Gb(i,j,l) = Gin(i,j)
      enddo
      do i = Nd1, Nstr
        Gout(i,j) = Dp(i,j) * exp(Dtred2(j) * (Odd(l) + Dump0))
        Gin(i,j) = Dp(i,j) * exp(Dtred2(j) * Dump0)

        Gb(i,j,l) = Gout(i,j)
        Gt(i,j,l) = Gin(i,j)
      enddo
    enddo

    S(:,:,l) = real(matmul(Gout, Invmx(Gin)), Kreal)
    D1(:,:) = Invmx(D(:,:))

    ! source from solar beam
    do i = 1, Nstr
      D(i,i) = D(i,i) + 1.0_Kreal2 / Mu0d
    enddo
    Z = matmul(transpose(PPm(1:Ngm,:))*spread(Glm(:,l),dim=1,ncopies=Nstr),Pmu0(1:Ngm)) !P0m
    Z = matmul(UUinv, Z) ! D0m
    Z = matmul(Invmx(D), Z) * Ssad(l) / (4.0_Kreal2 * Pi)

    Qt(:,l) = real(Z, Kreal) * Beam(l-1)
    Qb(:,l) = real(Z, Kreal) * Beam(l)
  
    S0(1:Nup) = Z(1:Nup) * exp(-Odd(l)/Mu0)
    S0(Nd1:Nstr) = Z(Nd1:Nstr)
    S0 = matmul(S(:,:,l), S0)
    S0(1:Nup) = Z(1:Nup) - S0(1:Nup)
    S0(Nd1:Nstr) = Z(Nd1:Nstr) * exp(-Odd(l)/Mu0) - S0(Nd1:Nstr)

    ! thermal source
    if (Thermal) then
      EE = matmul(UUinv, ATint)
      Yt = (1.0_Kreal - Ssad(l)) * real(matmul(D1,EE), Kreal) * &
      &(B0(l) - B0(l-1))
      Y1 = Yt / Odd(l)
      Y0 = real(matmul(D1,(Y1+B0(l-1)*(1.0_Kreal-Ssad(l))*EE)), Kreal)

      Qt(:,l) = Qt(:,l) + Y0
      Qb(:,l) = Qb(:,l) + Y0 + Yt

      Bpm(1:Nup) = Y0(1:Nup) + Yt(1:Nup)
      Bpm(Nd1:Nstr) = Y0(Nd1:Nstr)
      Bpm = matmul(S(:,:,l), Bpm)
      Bpm(1:Nup) = Y0(1:Nup) - Bpm(1:Nup)
      Bpm(Nd1:Nstr) = Y0(Nd1:Nstr) + Yt(Nd1:Nstr) - Bpm(Nd1:Nstr)
    else
      Bpm(1:Nstr) = 0.0_Kreal
    endif
    Xpm(:,l) = real(S0,Kreal) * Beam(l-1) + Bpm
    
  enddo

  ! Surface
  if (Nhm.gt.0) then
    call Surf(Nup, Ndn, Nhm, AAinv, PPm(1:Nhm,1:Nstr), PPf, Pmu0(1:Nhm), &
    & Hm, Thermal, ATint, Beam(Nlyr), Bsurf, Mu0, Rs, Sigsup)
  else
    Rs = 0.0_Kreal
    Sigsup = 0.0_Kreal
  endif

  S(:,:,Nlyr+1) = 0.0_Kreal
  Xpm(:,Nlyr+1) = 0.0_Kreal
  S(1:Nup,Nd1:Nstr,Nlyr+1) = Rs(1:Nup,1:Ndn)
  Xpm(1:Nup,Nlyr+1) = Sigsup

  ! Bounary conditions can be solved with Adding or Band Diagonal System
  ! Option 1: Adding
  call Adding2(Nlyr, Nup, Ndn, S, Xpm, Itild)

  ! Option 2: Band Diagonal System
  ! call BandDiag(Nlyr, Nup, Ndn, Gt, Gb, Qt, Qb, Rs, Sigsup, Itild)

endsubroutine Gdrtes

! inner subroutines ======================================================

subroutine Adding2(Nlyr, Nup, Ndn, S, Xpm, Itild)
    !use : Invmx
  implicit none
  integer(Kint), intent(in):: Nlyr, Nup, Ndn
  real(Kreal), intent(in):: S(Nup+Ndn,Nup+Ndn,Nlyr+1), Xpm(Nup+Ndn,Nlyr+1)
  real(Kreal), intent(out):: Itild(Nup+Ndn,0:Nlyr)

  real(Kreal):: Su(Nup+Ndn,Nup+Ndn,Nlyr+1), Sd(Nup+Ndn,Nup+Ndn,0:Nlyr)
  real(Kreal):: Xu(Nup+Ndn,Nlyr+1), Xd(Nup+Ndn,0:Nlyr)
  real(Kreal), dimension(Nup+Ndn):: Xustar

  integer(Kint):: l, Nd1, Nstr

  real(Kreal):: W1(Ndn,Ndn), W2(Nup,Nup)

  Nd1 = Nup + 1
  Nstr = Nup + Ndn

  Su(:,:,1) = S(:,:,1)
  Xu(:,1) = Xpm(:,1)
  do l = 2, Nlyr
  
    Su(:,:,l) = 0.0_Kreal
    Xu(:,l) = 0.0_Kreal
  
    W1 = matmul(S(Nd1:Nstr,Nd1:Nstr,l), Invmx(Identity(Ndn,1.0_Kreal) - &
    & matmul(Su(Nd1:Nstr,1:Nup,l-1), S(1:Nup,Nd1:Nstr,l))))
  
    Su(Nd1:Nstr,1:Nup,l) = S(Nd1:Nstr,1:Nup,l) + &
    & matmul(W1, matmul(Su(Nd1:Nstr,1:Nup,l-1), S(1:Nup,1:Nup,l)))
  
    Xu(Nd1:Nstr,l) = Xpm(Nd1:Nstr,l) + matmul(W1, (matmul(Su(Nd1:Nstr,1:Nup,l-1), &
    & Xpm(1:Nup,l)) + Xu(Nd1:Nstr,l-1)))
  
  enddo

  Sd(:,:,Nlyr) = S(:,:,Nlyr+1)
  Xd(:,Nlyr) = Xpm(:,Nlyr+1)
  do l = Nlyr-1, 0, -1

    Sd(:,:,l) = 0.0_Kreal
    Xd(:,l) = 0.0_Kreal
  
    W2 = matmul(S(1:Nup,1:Nup,l+1), Invmx(Identity(Nup,1.0_Kreal) - &
      matmul(Sd(1:Nup,Nd1:Nstr,l+1), S(Nd1:Nstr,1:Nup,l+1))))
  
    Sd(1:Nup,Nd1:Nstr,l) = S(1:Nup,Nd1:Nstr,l+1) + &
      matmul(W2, matmul(Sd(1:Nup,Nd1:Nstr,l+1), S(Nd1:Nstr,Nd1:Nstr,l+1)))
  
    Xd(1:Nup,l) = Xpm(1:Nup,l+1) + matmul(W2, (Xd(1:Nup,l+1) + &
      matmul(Sd(1:Nup,Nd1:Nstr,l+1), Xpm(Nd1:Nstr,l+1))))
  
  enddo

  ! Itild
  Itild(Nd1:Nstr,0) = 0.0_Kreal ! no downward diffuse at TOA
  Itild(1:Nup,0) = Xd(1:Nup,0) + matmul(Sd(1:Nup, &
  & Nd1:Nstr,0), Itild(Nd1:Nstr,0))
  do l = 1, Nlyr
    Xustar(1:Ndn) = Xu(Nd1:Nstr,l) + matmul( &
    & Su(Nd1:Nstr,Nd1:Nstr,l), Itild(Nd1:Nstr,0))
    Itild(1:Nup,l) = matmul(Invmx(Identity(Nup,1.0_Kreal) &
    & - matmul(Sd(1:Nup,Nd1:Nstr,l), Su(Nd1:Nstr, 1:Nup,l))), &
    & matmul(Sd(1:Nup,Nd1:Nstr,l), Xustar(1:Ndn)) + Xd(1:Nup,l))
    Itild(Nd1:Nstr,l) = matmul(Invmx(Identity(Ndn,1.0_Kreal) &
    & - matmul(Su(Nd1:Nstr,1:Nup,l), Sd(1:Nup,Nd1:Nstr,l))), &
    & matmul(Su(Nd1:Nstr,1:Nup,l), Xd(1:Nup,l)) + Xustar(1:Ndn))
  enddo

end subroutine Adding2

!
!subroutine Adding(Nlyr, Nup, Ndn, S, Xpm, Itild)
!  !use : Invmx
!  implicit none
!  integer(Kint), intent(in):: Nlyr, Nup, Ndn
!  real(Kreal), intent(in):: S(Nup+Ndn,Nup+Ndn,Nlyr+1), Xpm(Nup+Ndn,Nlyr+1)
!  real(Kreal), intent(out):: Itild(Nup+Ndn,0:Nlyr)
!
!  real(Kreal):: Su(Nup+Ndn,Nup+Ndn,Nlyr+1), Sd(Nup+Ndn,Nup+Ndn,0:Nlyr)
!  real(Kreal):: Xu(Nup+Ndn,Nlyr+1), Xd(Nup+Ndn,0:Nlyr)
!  real(Kreal), dimension(Nup+Ndn):: Xustar
!
!  integer(Kint):: l, Nd1, Nstr
!
!  Nd1 = Nup + 1
!  Nstr = Nup + Ndn
!
!  Su(:,:,1) = S(:,:,1)
!  Xu(:,1) = Xpm(:,1)
!  do l = 2, Nlyr
!    call Addinga(Nup, Ndn, Su(:,:,l-1), Xu(:,l-1), &
!    & S(:,:,l), Xpm(:,l), Su(:,:,l), &
!    & Xu(:,l))
!  enddo
!  Sd(:,:,Nlyr) = S(:,:,Nlyr+1)
!  Xd(:,Nlyr) = Xpm(:,Nlyr+1)
!  do l = Nlyr-1, 0, -1
!    call Addingb(Nup, Ndn, S(:,:,l+1), Xpm(:,l+1), &
!    & Sd(:,:,l+1), Xd(:,l+1), Sd(:,:,l), &
!    & Xd(:,l))
!  enddo
!
!  ! Itild
!  Itild(Nd1:Nstr,0) = 0.0_Kreal ! no downward diffuse at TOA
!  Itild(1:Nup,0) = Xd(1:Nup,0) + matmul(Sd(1:Nup, &
!  & Nd1:Nstr,0), Itild(Nd1:Nstr,0))
!  do l = 1, Nlyr
!    Xustar(1:Ndn) = Xu(Nd1:Nstr,l) + matmul( &
!    & Su(Nd1:Nstr,Nd1:Nstr,l), Itild(Nd1:Nstr,0))
!    Itild(1:Nup,l) = matmul(Invmx(Identity(Nup,1.0_Kreal) &
!    & - matmul(Sd(1:Nup,Nd1:Nstr,l), Su(Nd1:Nstr, 1:Nup,l))), &
!    & matmul(Sd(1:Nup,Nd1:Nstr,l), Xustar(1:Ndn)) + Xd(1:Nup,l))
!    Itild(Nd1:Nstr,l) = matmul(Invmx(Identity(Ndn,1.0_Kreal) &
!    & - matmul(Su(Nd1:Nstr,1:Nup,l), Sd(1:Nup,Nd1:Nstr,l))), &
!    & matmul(Su(Nd1:Nstr,1:Nup,l), Xd(1:Nup,l)) + Xustar(1:Ndn))
!  enddo
!
!end subroutine Adding
!
!subroutine Addinga(Nu0, Nd0, Sa, Xa, Sb, Xb, Sc, Xc)
!  !use : Invmx, Identity
!  implicit none
!  integer(Kint), intent(in):: Nu0, Nd0
!  real(Kreal), intent(in):: Sa(:,:), Xa(:), Sb(:,:), Xb(:)
!  real(Kreal), intent(out):: Sc(:,:), Xc(:)
!
!  real(Kreal):: W1(Nd0,Nd0)
!  integer(Kint):: N, Nd1
!
!  N = Nu0 + Nd0
!  Nd1 = Nu0 + 1
!
!  Sc = 0.0_Kreal
!  Xc = 0.0_Kreal
!
!  W1 = matmul(Sb(Nd1:N,Nd1:N), Invmx(Identity(Nd0,1._Kreal) - &
!  & matmul(Sa(Nd1:N,1:Nu0), Sb(1:Nu0,Nd1:N))))
!
!  Sc(Nd1:N,1:Nu0) = Sb(Nd1:N,1:Nu0) + &
!  & matmul(W1, matmul(Sa(Nd1:N,1:Nu0), Sb(1:Nu0,1:Nu0)))
!
!  Xc(Nd1:N) = Xb(Nd1:N) + matmul(W1, (matmul(Sa(Nd1:N,1:Nu0), &
!  & Xb(1:Nu0)) + Xa(Nd1:N)))
!
!end subroutine Addinga
!
!subroutine Addingb(Nu0, Nd0, Sa, Xa, Sb, Xb, Sc, Xc)
!  !use : Invmx, Identity
!  implicit none
!  integer(Kint), intent(in):: Nu0, Nd0
!  real(Kreal), intent(in):: Sa(:,:), Xa(:), Sb(:,:), Xb(:)
!  real(Kreal), intent(out):: Sc(:,:), Xc(:)
!
!  real(Kreal):: W2(Nu0,Nu0)
!  integer(Kint):: N, Nd1
!
!  N = Nu0 + Nd0
!  Nd1 = Nu0 + 1
!
!  Sc = 0.0_Kreal
!  Xc = 0.0_Kreal
!
!  W2 = matmul(Sa(1:Nu0,1:Nu0), Invmx(Identity(Nu0,1._Kreal) - &
!    matmul(Sb(1:Nu0,Nd1:N), Sa(Nd1:N,1:Nu0))))
!
!  Sc(1:Nu0,Nd1:N) = Sa(1:Nu0,Nd1:N) + &
!    matmul(W2, matmul(Sb(1:Nu0,Nd1:N), Sa(Nd1:N,Nd1:N)))
!
!  Xc(1:Nu0) = Xa(1:Nu0) + matmul(W2, (Xb(1:Nu0) + &
!    matmul(Sb(1:Nu0,Nd1:N), Xa(Nd1:N))))
!
!end subroutine Addingb

!-----------------------------------------------------------------------
! Band Diagonal System
subroutine BandDiag(Nlyr, Nup, Ndn, Gt, Gb, Qt, Qb, Rs, Sigsup, Itild)
  implicit none
  integer(Kint), intent(in):: Nlyr, Nup, Ndn
  real(Kreal2), intent(in):: Gt(Nup+Ndn,Nup+Ndn,Nlyr), Gb(Nup+Ndn,Nup+Ndn,Nlyr)
  real(Kreal), intent(in):: Qt(Nup+Ndn,Nlyr), Qb(Nup+Ndn,Nlyr), Rs(Nup,Ndn), Sigsup(Nup)
  real(Kreal), intent(out):: Itild(Nup+Ndn,0:Nlyr)

  real(Kreal):: A((Nup+Ndn)*Nlyr,(Nup+Ndn)*3-1), Al((Nup+Ndn)*Nlyr,Nup+Ndn+Ndn-1), &
  & B((Nup+Ndn)*Nlyr), D
  real(Kreal2):: Glp(Ndn,Nup+Ndn)
  integer(Kint):: M1, M2, Nstr, n, l, i, Indx((Nup+Ndn)*Nlyr)

  Nstr = Nup + Ndn
  M1 = Nstr + Ndn - 1
  M2 = Nstr + Nup - 1
  Glp = matmul(Rs, Gb(Nup+1:Nstr,:,Nlyr))

  do n = 1, Ndn
    A(n,:) = 0.0_Kreal
    A(n,M1+2-n:M1+Nstr+1-n) = real(Gt(Nup+n,1:Nstr,1),Kreal)
    B(n) = - Qt(Nup+n,1)
  enddo
  do l = 1, Nlyr-1
    do n = 1, Nstr
      i = Ndn + (l-1)*Nstr + n
      A(i,:) = 0.0_Kreal
      A(i,Nstr+1-n:Nstr*2-n) = real(Gb(n,1:Nstr,l),Kreal)
      A(i,Nstr*2-n+1:Nstr*3-n) =  real(- Gt(n,1:Nstr,l+1), Kreal)
      B(i) = - Qb(n,l) + Qt(n,l+1)
    enddo
  enddo
  do n = 1, Nup
    i = Ndn + (Nlyr-1)*Nstr + n
    A(i,:) = 0.0_Kreal
    A(i,Nstr+1-n:Nstr*2-n) = real(Gb(n,:,Nlyr) - Glp(n,:), Kreal)
    B(i) = - Qb(n,Nlyr) + dot_product(Rs(n,:),Qb(Nup+1:Nstr,Nlyr)) + Sigsup(n)
  enddo

  call Bandec(A, (Nup+Ndn)*Nlyr, M1, M2, Al, Indx, D)
  call Banbks(A, (Nup+Ndn)*Nlyr, M1, M2, Al, Indx, B)

  do l = 1, Nlyr
    Itild(:,l-1) = matmul(real(Gt(:,:,l),Kreal), B((Nup+Ndn)*(l-1)+1:(Nup+Ndn)*(l))) + Qt(:,l)
  enddo
  Itild(:,Nlyr) = matmul(real(Gb(:,:,Nlyr),Kreal), B((Nup+Ndn)*(Nlyr-1)+1:(Nup+Ndn)*(Nlyr))) + Qb(:,Nlyr)

endsubroutine BandDiag

SUBROUTINE Bandec(A, N, M1, M2, Al, Indx, D)
  !USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,swap,arth
  implicit none
  real(Kreal), intent(inout):: A(N,M1+M2+1)
  real(Kreal), intent(out):: Al(N,M1)
  integer(Kint), intent(in):: N, M1, M2
  integer(Kint), intent(out):: Indx(N)
  real(Kreal), intent(out) :: D
  
  real(Kreal), parameter:: Tiny = 1.0e-20_Kreal
  
  integer(Kint):: i, k, l, Mdum, Mm, Imax(1)
  real(Kreal) :: Dum

  !n=assert_eq(size(a,1),size(al,1),size(indx),'bandec: n')
  !mm=assert_eq(size(a,2),m1+m2+1,'bandec: mm')
  Mm = M1 + M2 + 1
  !mdum=assert_eq(size(al,2),m1,'bandec: mdum')
  Mdum = M1
  A(1:M1,:) = eoshift(A(1:M1,:), dim=2, shift=Arth_i(M1,-1,M1))
  D = 1.0_Kreal
  do k = 1, N
    l = min(M1+k, N)
    !i = Imaxloc(abs(A(k:l,1))) + k - 1
    Imax = maxloc(abs(A(k:l,1)))
    i = Imax(1) + k - 1
    Dum = A(i,1)
    if (Dum == 0.0_Kreal) A(k,1) = Tiny
    Indx(k) = i
    if (i /= k) then
      D = -D
      call Abba(A(k,1:Mm), A(i,1:Mm))
    endif
    do i = k+1, l
      Dum = A(i,1) / A(k,1)
      Al(k,i-k) = Dum
      A(i,1:Mm-1) = A(i,2:Mm) - Dum*A(k,2:Mm)
      A(i,Mm) = 0.0_Kreal
    enddo
  enddo
endsubroutine Bandec

function Arth_i(First, Increment, N)
  integer(Kint), intent(in) :: First, Increment, N
  integer(Kint):: Arth_i(N)
  integer(Kint):: k, k2, temp
  integer(Kint), parameter :: Npar_arth=16, Npar2_arth=8
  
  if (N > 0) Arth_i(1) = first
  if (N <= Npar_arth) then
    do k = 2, N
      Arth_i(k) = Arth_i(k-1) + Increment
    enddo
  else
    do k = 2, Npar2_arth
      Arth_i(k) = Arth_i(k-1) + Increment
    enddo
    temp = increment * Npar2_arth
    k = Npar2_arth
    do
      if (k >= N) exit
      k2 = k + k
      Arth_i(k+1:min(k2,N)) = temp + Arth_i(1:min(k,N-k))
      temp = temp + temp
      k = k2
    enddo
  endif
endfunction Arth_i

subroutine Banbks(A, N, M1, M2, Al, Indx, B)
  !USE nrtype; USE nrutil, ONLY : assert_eq,swap
  implicit none
  real(Kreal), intent(in):: A(N,M1+M2+1), Al(N,M1)
  integer(Kint), intent(in):: M1, M2, N
  integer(Kint), intent(in):: Indx(N)
  real(Kreal), intent(inout):: B(N)
  
  integer(Kint):: i, k, l, Mdum, Mm

  !n=assert_eq(size(a,1),size(al,1),size(b),size(indx),'banbks: n')
  !mm=assert_eq(size(a,2),m1+m2+1,'banbks: mm')
  MM = M1 + M2 + 1
  !mdum=assert_eq(size(al,2),m1,'banbks: mdum')
  Mdum = M1

  do k = 1, N
    l = min(N,M1+k)
    i = Indx(k)
    if (i /= k) call Abba(b(i),b(k))
    B(k+1:l) = b(k+1:l) - Al(k,1:l-k) * B(k)
  enddo
  do i = N, 1, -1
    l = min(Mm,N-i+1)
    B(i) = (B(i) - dot_product(A(i,2:l), B(1+i:i+l-1))) / A(i,1)
  enddo

endsubroutine Banbks
  

!--------------------------------------------------------------------

subroutine Surf(Nup, Ndn, Nhm, AAinv, PPh, PPf, Pmu0, Hm, Thermal, ATint, Beamsurf, Bs, Mu0, Rs, Sigsup)

  integer(Kint), intent(in):: Nup, Ndn, Nhm
  real(Kreal2), intent(in):: AAinv(Nup+Ndn,Nup+Ndn), PPh(Nhm,Nup+Ndn), PPf(Nhm,Ndn), Pmu0(Nhm), ATint(Nup+Ndn)
  real(Kreal), intent(in):: Hm(Nhm,Nhm), Beamsurf, Bs, Mu0
  logical, intent(in):: Thermal
  real(Kreal), intent(out):: Rs(Nup,Ndn), Sigsup(Nup)

  real(Kreal2):: temp(Nup,Nhm), Rh(Nup), Pint(0:Nhm), PPf1(Nhm)
  real(Kreal):: Ps0m(Nup), Epsm(Nup)
  integer(Kint):: l, i

  ! constant PPf1
  Pint(0) = 1.0_Kreal2
  Pint(1) = 0.5_Kreal2
  Pint(2:Nhm:2) = 0.0_Kreal2
  do l = 3, Nhm, 2
    Pint(l) = - (l-2) * Pint(l-2) / (l+1)
  enddo

  PPf1(1) = Pint(1)
  if (Nhm.gt.1) then
    PPf1(2) = - Pint(0) / 3.0_Kreal2
    do i = 3, Nhm, 2
      PPf1(i) = (i*Pint(i) + (i-1)*Pint(i-2)) / (2*i - 1)
    enddo
    PPf1(4:Nhm:2) = 0.0_Kreal2
  endif

  temp = matmul(AAinv(1:Nup,:), matmul(transpose(PPh), Hm))

  Rs = real(matmul(temp, PPf), Kreal) * 2 * Pi

  Ps0m = real(matmul(temp, Pmu0), Kreal)

  ! thermal source
  if (Thermal) then

    Rh = matmul(temp, PPf1) * 2 * Pi_d

    Epsm = real(matmul(AAinv(1:Nup,:), ATint) - Rh, Kreal)

  else

    Epsm = 0.0_Kreal

  endif

  Sigsup = Beamsurf * Mu0 * Ps0m + Epsm * Bs

endsubroutine Surf

! math utilities ====================================================

!-------------------------------------------------------------------- 
! Abba

subroutine Abba_i(A,B)
  integer(Kint), intent(inout):: A,B
  integer(Kint):: Tmp
  Tmp = A
  A = B
  B = Tmp
endsubroutine Abba_i
!BL
subroutine Abba_r(A,B)
  real(Kreal1), intent(inout):: A,B
  real(Kreal1):: Tmp
  Tmp = A
  A = B
  B = Tmp
endsubroutine Abba_r
!BL
subroutine Abba_d(A,B)
  real(Kreal2), intent(inout):: A,B
  real(Kreal2):: Tmp
  Tmp = A
  A = B
  B = Tmp
endsubroutine Abba_d
!BL
subroutine Abba_rv(A,B)
  real(Kreal1), dimension(:), intent(inout):: A,B
  real(Kreal1), dimension(SIZE(A)):: Tmp
  Tmp = A
  A = B
  B = Tmp
endsubroutine Abba_rv
!BL
subroutine Abba_dv(A,B)
  real(Kreal2), dimension(:), intent(inout):: A,B
  real(Kreal2), dimension(SIZE(A)):: Tmp
  Tmp = A
  A = B
  B = Tmp
endsubroutine Abba_dv

!---------------------------------------------------------------------
SUBROUTINE Asymtx( AA, EVEC, EVAL, M, IA, IEVEC, IER)
! Solves eigenfunction problem for real asymmetric matrix
!   for which it is known a priori that the eigenvalues are real.
!   Adopted from DISORT4.0.99.
!
!   I N P U T    V A R I A B L E S:
!       AA    :  input asymmetric matrix, destroyed after solved
!        M    :  order of  AA
!       IA    :  first dimension of  AA
!    IEVEC    :  first dimension of  EVEC
!
!   O U T P U T    V A R I A B L E S:
!       EVEC  :  (unnormalized) eigenvectors of  AA
!                   ( column J corresponds to EVAL(J) )
!       EVAL  :  (unordered) eigenvalues of AA ( dimension at least M )
!       IER   :  if .NE. 0, signals that EVAL(IER) failed to converge;
!                   in that case eigenvalues IER+1,IER+2,...,M  are
!                   correct but eigenvalues 1,...,IER are set to zero.
!
!   S C R A T C H   V A R I A B L E S:
!       WKD   :  work area ( dimension at least 2*M )
!       AAD   :  double precision stand-in for AA
!       EVECD :  double precision stand-in for EVEC
!       EVALD :  double precision stand-in for EVAL
!
  !   Calls- D1MACH, ERRMSG
! +-------------------------------------------------------------------+
  implicit none

! .. Scalar Arguments ..

  integer(Kint), intent(in):: IA, IEVEC, M
  integer(Kint), intent(out):: IER
! ..
! .. Array Arguments ..

  real(Kreal2), intent(inout):: AA( IA, M )
  real(Kreal2), intent(out):: EVAL( M ), EVEC( IEVEC, M )

  real(Kreal2):: AAD( IA, M ), EVALD( M ), EVECD( IA, M ), WKD( 2*M )
! ..
! .. Local Scalars ..

  LOGICAL::   NOCONV, NOTLAS
  INTEGER::   I, II, IN, J, K, KA, KKK, L, LB, LLL, N, N1, N2
  real(Kreal2):: COL, DISCRI, F, G, H, &
    &           P, Q, R, REPL, RNORM, ROW, S, SCALE, SGN, T, &
    &           TOL, UU, VV, W, X, Y, Z
! ..
! .. External Functions ..

  DOUBLE PRECISION D1MACH
  EXTERNAL  D1MACH
!     ..
!     .. External Subroutines ..

!      EXTERNAL  ERRMSG
!     ..
!     .. Intrinsic Functions ..

  INTRINSIC ABS, MIN, SIGN, SQRT
! ..
  real(Kreal2), save:: C1=0.4375D0, C2=0.5D0, C3=0.75D0, &
     &                 C4=0.95D0, C5=16.D0, C6=256.D0, ZERO=0.D0, ONE=1.D0

  P = 0d0
  Q = 0d0
  R = 0d0

  IER = 0
  TOL = real(Radix(1.D0),Kind=Kreal2) ** (-DIGITS(1.D0)) !TOL  = D1MACH( 4 )
  LB = 0

  if ( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M ) then
    print*, 'ASYMTX--bad input variable(s), M, IA, IEVEC:', M, IA, IEVEC
    stop
  endif


!                           ** Handle 1x1 and 2x2 special cases
  IF( M.EQ.1 ) THEN
    EVAL( 1 )   = AA( 1,1 )
    EVEC( 1,1 ) = 1.0_Kreal2
    RETURN
  ELSE IF( M.EQ.2 ) THEN
    DISCRI = ( AA( 1,1 ) - AA( 2,2 ) )**2 + 4.D0*AA( 1,2 )*AA( 2,1 )
    if ( DISCRI .LT. 0.0_Kreal2 ) then
      print*, 'ASYMTX--complex evals in 2x2 case'
      stop
    endif

    SGN  = ONE
    IF( AA( 1,1 ) .LT. AA( 2,2 ) ) SGN  = - ONE

    EVAL( 1 ) = 0.5_Kreal2*( AA( 1,1 ) + AA( 2,2 ) + SGN*SQRT( DISCRI ) )
    EVAL( 2 ) = 0.5_Kreal2*( AA( 1,1 ) + AA( 2,2 ) - SGN*SQRT( DISCRI ) )
    EVEC( 1,1 ) = 1.0_Kreal2
    EVEC( 2,2 ) = 1.0_Kreal2

    IF( AA( 1,1 ) .EQ. AA( 2,2 ) .AND. &
    & ( AA( 2,1 ).EQ.0.0_Kreal2 .OR. AA( 1,2 ).EQ.0.0_Kreal2 ) ) THEN

      RNORM = ABS( AA( 1,1 ) ) + ABS( AA( 1,2 ) ) + &
     &        ABS( AA( 2,1 ) ) + ABS( AA( 2,2 ) )
      W     = TOL * RNORM
      EVEC( 2,1 ) =   REAL( AA( 2,1 ) / W )
      EVEC( 1,2 ) = - REAL( AA( 1,2 ) / W )

    ELSE

      EVEC( 2,1 ) = AA( 2,1 ) / ( EVAL( 1 ) - AA( 2,2 ) )
      EVEC( 1,2 ) = AA( 1,2 ) / ( EVAL( 2 ) - AA( 1,1 ) )

    END IF

    RETURN

  END IF

!                               ** Convert single-prec. matrix to double
  AAD( 1:M,1:M ) = AA( 1:M,1:M )

!                               ** Initialize output variables
  IER  = 0
  EVALD( 1:M ) = ZERO
  EVECD( 1:M, 1:M ) = ZERO
  DO I = 1, M
    EVECD( I, I ) = ONE
  ENDDO

!                  ** Balance the input matrix and reduce its norm by
!                  ** diagonal similarity transformation stored in WK;
!                  ** then search for rows isolating an eigenvalue
!                  ** and push them down
  RNORM  = ZERO
  L  = 1
  K  = M

  50 CONTINUE
  KKK  = K

  DO J = KKK, 1, -1

    ROW  = ZERO

    DO I = 1, K
      IF( I.NE.J ) ROW  = ROW + ABS( AAD( J,I ) )
    enddo

    IF( ROW.EQ.ZERO ) THEN

      WKD( K ) = J

      IF( J.NE.K ) THEN

        DO I = 1, K
          REPL        = AAD( I, J )
          AAD( I, J ) = AAD( I, K )
          AAD( I, K ) = REPL
        enddo

        DO I = L, M
          REPL        = AAD( J, I )
          AAD( J, I ) = AAD( K, I )
          AAD( K, I ) = REPL
        enddo

      END IF

      K  = K - 1
      GO TO  50

    END IF

  enddo
!                                ** Search for columns isolating an
!                                ** eigenvalue and push them left
  100 CONTINUE
  LLL  = L

  DO J = LLL, K

    COL  = ZERO

    DO I = L, K
      IF( I.NE.J ) COL  = COL + ABS( AAD( I,J ) )
    enddo

    IF( COL.EQ.ZERO ) THEN

      WKD( L ) = J

      IF( J.NE.L ) THEN

        DO I = 1, K
          REPL        = AAD( I, J )
          AAD( I, J ) = AAD( I, L )
          AAD( I, L ) = REPL
        enddo

        DO I = L, M
          REPL        = AAD( J, I )
          AAD( J, I ) = AAD( L, I )
          AAD( L, I ) = REPL
        enddo

      END IF

      L  = L + 1
      GO TO  100

    END IF

  enddo

!                           ** Balance the submatrix in rows L through K
  WKD( L:K ) = ONE
  
  160 CONTINUE
  NOCONV = .FALSE.

  DO I = L, K

    COL  = ZERO
    ROW  = ZERO

    DO J = L, K

      IF( J.NE.I ) THEN
        COL  = COL + ABS( AAD( J,I ) )
        ROW  = ROW + ABS( AAD( I,J ) )
      END IF

    enddo

    F  = ONE
    G  = ROW / C5
    H  = COL + ROW

    180 CONTINUE
    IF( COL.LT.G ) THEN

      F    = F*C5
      COL  = COL*C6
      GO TO  180

    END IF

    G  = ROW*C5

    190    CONTINUE
    IF( COL.GE.G ) THEN

      F    = F / C5
      COL  = COL / C6
      GO TO  190

    END IF
!                                                ** Now balance
    IF( ( COL + ROW ) / F.LT.C4*H ) THEN

      WKD( I ) = WKD( I )*F
      NOCONV = .TRUE.

      AAD( I, L:M ) = AAD( I, L:M ) / F

      AAD( 1:K, I ) = AAD( 1:K, I )*F

    END IF

  enddo


  IF( NOCONV ) GO TO  160
!                                   ** Is A already in Hessenberg form?
  IF( K-1 .LT. L+1 ) GO TO  370

!                                   ** Transfer A to a Hessenberg form
  DO N = L + 1, K - 1

    H  = ZERO
    WKD( N + M ) = ZERO
    SCALE  = ZERO
!                                                 ** Scale column
    DO I = N, K
      SCALE  = SCALE + ABS( AAD( I,N - 1 ) )
    enddo

    IF( SCALE.NE.ZERO ) THEN

      DO I = K, N, -1
        WKD( I + M ) = AAD( I, N - 1 ) / SCALE
        H  = H + WKD( I + M )**2
      enddo

      G    = - SIGN( SQRT( H ), WKD( N + M ) )
      H    = H - WKD( N + M )*G
      WKD( N + M ) = WKD( N + M ) - G
!                                            ** Form (I-(U*UT)/H)*A
      DO J = N, M

        F  = ZERO

        DO I = K, N, -1
          F  = F + WKD( I + M )*AAD( I, J )
        enddo

        DO I = N, K
          AAD( I, J ) = AAD( I, J ) - WKD( I + M )*F / H
        enddo

      enddo
!                                    ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
      DO I = 1, K

        F  = ZERO

        DO J = K, N, -1
          F  = F + WKD( J + M )*AAD( I, J )
        enddo

        DO J = N, K
          AAD( I, J ) = AAD( I, J ) - WKD( J + M )*F / H
        enddo

      enddo

      WKD( N + M ) = SCALE*WKD( N + M )
      AAD( N, N - 1 ) = SCALE*G

    END IF

  enddo


  DO N = K - 2, L, -1

    N1   = N + 1
    N2   = N + 2
    F  = AAD( N + 1, N )

    IF( F.NE.ZERO ) THEN

      F  = F*WKD( N + 1 + M )

      WKD( N+2+M:K+M ) = AAD( N+2:K, N )

      IF( N + 1.LE.K ) THEN

        DO J = 1, M

          G  = ZERO

          DO I = N + 1, K
            G  = G + WKD( I + M )*EVECD( I, J )
          enddo

          G  = G / F

          DO I = N + 1, K
            EVECD( I, J ) = EVECD( I, J ) + G*WKD( I + M )
          enddo

        enddo

      END IF

    END IF

  enddo


  370 CONTINUE

  N  = 1

  DO I = 1, M

    DO J = N, M
      RNORM  = RNORM + ABS( AAD( I,J ) )
    enddo

    N  = I

    IF( I.LT.L .OR. I.GT.K ) EVALD( I ) = AAD( I, I )

  enddo

  N  = K
  T  = ZERO

!                                      ** Search for next eigenvalues
  400 CONTINUE
  IF( N.LT.L ) GO TO  550

  IN  = 0
  N1  = N - 1
  N2  = N - 2
!                          ** Look for single small sub-diagonal element
  410 CONTINUE

  DO I = L, N

    LB  = N + L - I

    IF( LB.EQ.L ) GO TO  430

    S  = ABS( AAD( LB - 1,LB - 1 ) ) + ABS( AAD( LB,LB ) )

    IF( S.EQ.ZERO ) S  = RNORM

    IF( ABS( AAD( LB, LB-1 ) ).LE. TOL*S ) GO TO  430

  enddo


  430 CONTINUE
  X  = AAD( N, N )

  IF( LB.EQ.N ) THEN
    !                                        ** One eigenvalue found
    AAD( N, N ) = X + T
    EVALD( N ) = AAD( N, N )
    N  = N1
    GO TO  400

  END IF

  Y  = AAD( N1, N1 )
  W  = AAD( N, N1 )*AAD( N1, N )

  IF( LB.EQ.N1 ) THEN
    !                                        ** Two eigenvalues found
    P  = ( Y - X )*C2
    Q  = P**2 + W
    Z  = SQRT( ABS( Q ) )
    AAD( N, N ) = X + T
    X  = AAD( N, N )
    AAD( N1, N1 ) = Y + T
    !                                        ** Real pair
    Z  = P + SIGN( Z, P )
    EVALD( N1 ) = X + Z
    EVALD( N ) = EVALD( N1 )

    IF( Z.NE.ZERO ) EVALD( N ) = X - W / Z

    X  = AAD( N, N1 )
    !                                  ** Employ scale factor in case
    !                                  ** X and Z are very small
    R  = SQRT( X*X + Z*Z )
    P  = X / R
    Q  = Z / R
    !                                             ** Row modification
    DO J = N1, M
      Z  = AAD( N1, J )
      AAD( N1, J ) = Q*Z + P*AAD( N, J )
      AAD( N, J ) = Q*AAD( N, J ) - P*Z
    enddo
    !                                             ** Column modification
    DO I = 1, N
      Z  = AAD( I, N1 )
      AAD( I, N1 ) = Q*Z + P*AAD( I, N )
      AAD( I, N ) = Q*AAD( I, N ) - P*Z
    enddo
    !                                          ** Accumulate transformations
    DO I = L, K
      Z  = EVECD( I, N1 )
      EVECD( I, N1 ) = Q*Z + P*EVECD( I, N )
      EVECD( I, N ) = Q*EVECD( I, N ) - P*Z
    enddo

    N  = N2
    GO TO  400

  END IF


  IF( IN.EQ.30 ) THEN

    !                    ** No convergence after 30 iterations; set error
    !                    ** indicator to the index of the current eigenvalue
    IER  = N
    GO TO  700

  END IF
  !                                                  ** Form shift
  IF( IN.EQ.10 .OR. IN.EQ.20 ) THEN

    T  = T + X

    DO I = L, N
      AAD( I, I ) = AAD( I, I ) - X
    enddo

    S  = ABS( AAD( N,N1 ) ) + ABS( AAD( N1,N2 ) )
    X  = C3*S
    Y  = X
    W  = -C1*S**2

  END IF


  IN  = IN + 1

  !                ** Look for two consecutive small sub-diagonal elements

  DO J = LB, N2
    I  = N2 + LB - J
    Z  = AAD( I, I )
    R  = X - Z
    S  = Y - Z
    P  = ( R*S - W ) / AAD( I + 1, I ) + AAD( I, I + 1 )
    Q  = AAD( I + 1, I + 1 ) - Z - R - S
    R  = AAD( I + 2, I + 1 )
    S  = ABS( P ) + ABS( Q ) + ABS( R )
    P  = P / S
    Q  = Q / S
    R  = R / S

    IF( I.EQ.LB ) GO TO  490

    UU   = ABS( AAD( I, I-1 ) )*( ABS( Q ) + ABS( R ) )
    VV   = ABS( P ) * ( ABS( AAD( I-1, I-1 ) ) + ABS( Z ) + &
    &                   ABS( AAD( I+1, I+1 ) ) )

    IF( UU .LE. TOL*VV ) GO TO  490

  enddo

  490 CONTINUE
  AAD( I+2, I ) = ZERO

  DO J = I + 3, N
    AAD( J, J - 2 ) = ZERO
    AAD( J, J - 3 ) = ZERO
  enddo

!             ** Double QR step involving rows K to N and columns M to N

  DO KA = I, N1

    NOTLAS = KA.NE.N1

    IF( KA.EQ.I ) THEN

      S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )

      IF( LB.NE.I ) AAD( KA, KA - 1 ) = -AAD( KA, KA - 1 )

    ELSE

      P  = AAD( KA, KA - 1 )
      Q  = AAD( KA + 1, KA - 1 )
      R  = ZERO

      IF( NOTLAS ) R  = AAD( KA + 2, KA - 1 )

      X  = ABS( P ) + ABS( Q ) + ABS( R )

      IF( X.EQ.ZERO ) cycle !GO TO  540

      P  = P / X
      Q  = Q / X
      R  = R / X
      S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
      AAD( KA, KA - 1 ) = -S*X

    END IF

    P  = P + S
    X  = P / S
    Y  = Q / S
    Z  = R / S
    Q  = Q / P
    R  = R / P
    !                                              ** Row modification
    DO J = KA, M

      P  = AAD( KA, J ) + Q*AAD( KA + 1, J )

      IF( NOTLAS ) THEN

        P  = P + R*AAD( KA + 2, J )
        AAD( KA + 2, J ) = AAD( KA + 2, J ) - P*Z

      END IF

      AAD( KA + 1, J ) = AAD( KA + 1, J ) - P*Y
      AAD( KA, J ) = AAD( KA, J ) - P*X

    enddo
    !                                                 ** Column modification
    DO II = 1, MIN( N, KA + 3 )

      P  = X*AAD( II, KA ) + Y*AAD( II, KA + 1 )

      IF( NOTLAS ) THEN

        P  = P + Z*AAD( II, KA + 2 )
        AAD( II, KA + 2 ) = AAD( II, KA + 2 ) - P*R

      END IF

      AAD( II, KA + 1 ) = AAD( II, KA + 1 ) - P*Q
      AAD( II, KA ) = AAD( II, KA ) - P

    enddo
    !                                          ** Accumulate transformations
    DO II = L, K

      P  = X*EVECD( II, KA ) + Y*EVECD( II, KA + 1 )

      IF( NOTLAS ) THEN

        P  = P + Z*EVECD( II, KA + 2 )
        EVECD( II, KA + 2 ) = EVECD( II, KA + 2 ) - P*R

      END IF

      EVECD( II, KA + 1 ) = EVECD( II, KA + 1 ) - P*Q
      EVECD( II, KA ) = EVECD( II, KA ) - P

    enddo

  enddo !540

  GO TO  410
  !                     ** All evals found, now backsubstitute real vector
  550 CONTINUE

  IF( RNORM.NE.ZERO ) THEN

    DO N = M, 1, -1

      N2   = N
      AAD( N, N ) = ONE

      DO I = N - 1, 1, -1

        W  = AAD( I, I ) - EVALD( N )

        IF( W.EQ.ZERO ) W  = TOL*RNORM

        R  = AAD( I, N )

        DO J = N2, N - 1
          R  = R + AAD( I, J )*AAD( J, N )
        enddo

        AAD( I, N ) = -R / W
        N2   = I

      enddo !570

    enddo
    !                      ** End backsubstitution vectors of isolated evals
    DO I = 1, M

      IF( I.LT.L .OR. I.GT.K ) THEN

        EVECD( I, I:M ) = AAD( I, I:M )

      END IF

    enddo !600
    !                                   ** Multiply by transformation matrix
    IF( K.NE.0 ) THEN

      DO J = M, L, -1

        DO I = L, K

          Z  = ZERO

          DO N = L, MIN( J, K )
            Z  = Z + EVECD( I, N )*AAD( N, J )
          enddo

                  EVECD( I, J ) = Z

        enddo !620          CONTINUE

      enddo !630

    END IF

  END IF


  DO I = L, K

      EVECD( I, 1:M ) = EVECD( I, 1:M ) * WKD( I )

  enddo !650 CONTINUE

  !                           ** Interchange rows if permutations occurred
  DO I = L-1, 1, -1

    J  =  INT( WKD( I ) )

    IF( I.NE.J ) THEN

      DO N = 1, M
        REPL   = EVECD( I, N )
        EVECD( I, N ) = EVECD( J, N )
        EVECD( J, N ) = REPL
      enddo !660       CONTINUE

    END IF

  enddo !670 CONTINUE


  DO I = K + 1, M

    J  = INT( WKD( I ) )

    IF( I.NE.J ) THEN

      DO N = 1, M
        REPL   = EVECD( I, N )
        EVECD( I, N ) = EVECD( J, N )
        EVECD( J, N ) = REPL
      enddo !680       CONTINUE

    END IF

  enddo !690 CONTINUE

  !                         ** Put results into output arrays
  700 CONTINUE

  EVAL( 1:M ) = EVALD( 1:M )
  EVEC( 1:M, 1:M ) = EVECD( 1:M, 1:M )

END SUBROUTINE ASYMTX

! ------------------------------------------------------------------
! Cab

function Cab_r(A, B)

  implicit none
  real(Kreal1), intent(in):: A, B
  real(Kreal1):: Cab_r
  real(Kreal1):: Absa, Absb
  Absa = abs(A)
  Absb = abs(B)
  if (Absa > Absb) then
    Cab_r = Absa * sqrt(1.0_Kreal1 + (Absb / Absa) ** 2)
  else
    if (Absb ==0.0_Kreal1) then
      Cab_r = 0.0_Kreal1
    else
      Cab_r = Absb * sqrt(1.0_Kreal1 + (Absa / Absb) ** 2)
    endif
  endif
end function Cab_r

function Cab_d(A, B)

  implicit none
  real(Kreal2), intent(in):: A, B
  real(Kreal2):: Cab_d
  real(Kreal2):: Absa, Absb
  Absa = abs(A)
  Absb = abs(B)
  if (Absa > Absb) then
    Cab_d = Absa * sqrt(1.0_Kreal2 + (Absb / Absa) ** 2)
  else
    if (Absb ==0.0_Kreal2) then
      Cab_d = 0.0_Kreal2
    else
      Cab_d = Absb * sqrt(1.0_Kreal2 + (Absa / Absb) ** 2)
    endif
  endif
end function Cab_d

! --------------------------------------------------------------------
! Cskdcp

subroutine Cskdcp(XX, P)

  implicit none
  real(Kreal2), dimension(:,:), intent(inout):: XX
  real(Kreal2), dimension(:), intent(out):: P
  integer(Kint):: i, N
  real(Kreal2):: Summ

  N = size(P)
  do i = 1, N
    Summ = XX(i,i) - dot_product(XX(i,1:i-1),XX(i,1:i-1))
    if (Summ <= 0.0_Kreal2) then
      print * , 'Error: Summ <= 0.0, Summ = ', Summ
      print*, N, i, XX(i,i), dot_product(XX(i,1:i-1),XX(i,1:i-1))
      stop 'Stopped by Cskdcp.'
    endif
    P(i) = sqrt(Summ)
    XX(i+1:N,i) = (XX(i,i+1:N) - matmul(XX(i+1:N,1:i-1),XX(i,1:i-1))) &
      / P(i)
  enddo
endsubroutine Cskdcp

! --------------------------------------------------------------------
subroutine Deigen1d(App, Uinv, Dp, D)
  !use: Invmx, Cskdcp, Evtrisym, Rdcsym

  real(Kreal2), intent(in):: App(:,:), Uinv(:,:)
  real(Kreal2), intent(out):: Dp(:,:), D(:)

  real(Kreal2):: X(size(App,1),size(App,2)), Xdiag(size(App,1)), &
                E(size(App,1)), Dump(size(App,1),size(App,1)), &
                Dp_d(size(App,1),size(App,1)), Dq(size(App,1))
  integer(Kint):: i, j, Nstr
  
  Nstr=size(App,1)

  X(1:Nstr,1:Nstr) = App(1:Nstr,1:Nstr)
  call Cskdcp(X(1:Nstr,1:Nstr), Xdiag(1:Nstr))

  do i = 1, Nstr
    do j = i+1, Nstr
      X(i,j) = 0.0_Kreal2
    enddo
    X(i,i) = Xdiag(i)
  enddo

  Dp_d(1:Nstr,1:Nstr) = matmul(transpose(X(1:Nstr,1:Nstr)), &
    Uinv(1:Nstr,1:Nstr))
  Dp_d(1:Nstr,1:Nstr) = matmul(Dp_d(1:Nstr,1:Nstr), X(1:Nstr,1:Nstr))

  call Rdcsym(Dp_d(1:Nstr,1:Nstr), Dq(1:Nstr), E(1:Nstr))
  call Evtrisym(Dq(1:Nstr), E(1:Nstr), Dp_d(1:Nstr,1:Nstr))

  Dump(1:Nstr,1:Nstr) = transpose(X(1:Nstr,1:Nstr))
  Dp(1:Nstr,1:Nstr) = matmul(Invmx(Dump(1:Nstr,1:Nstr)), &
  Dp_d(1:Nstr,1:Nstr))
  D = Dq

endsubroutine Deigen1d

!--------------------------------------------------------------------
subroutine Deigendst(D, Dp, Dtred2, Nstr)
  ! Solve eigenvalue problem using the DST method, i.e., Asymtx
  implicit none
  integer(Kint), intent(in):: Nstr
  real(Kreal2), intent(in):: D(Nstr,Nstr)
  real(Kreal2), intent(out):: Dp(Nstr,Nstr), Dtred2(Nstr)
  
  integer(Kint):: Ier, Nup, Ndn
  real(Kreal2), dimension(Nstr/2,Nstr/2):: Array, Alpha, Beta
  real(Kreal2), dimension(Nstr,Nstr/2):: Evecp, Evecm
  real(Kreal2), dimension(Nstr/2):: Eval

  Nup = Nstr / 2
  Ndn = Nstr - Nup

  Alpha = D(1:Nup,1:Nup)
  Beta = D(1:Nup,Nup+1:Nstr)
  Array = matmul(Alpha - Beta, Alpha + Beta)
  call Asymtx(Array, Evecp, Eval, Ndn, Nup, Nstr, Ier)
  if (Ier /= 0) then
    print*, 'Error: Ier /= 0, Ier = ', Ier
    stop 'Stopped by Deigendst.'
  endif

  Dtred2(1:Nstr:2) = sqrt(abs(Eval))
  !Dtred2(Nup:1:-1) = - Dtred2(Nup+1:Nstr)
  Dtred2(2:Nstr:2) = - Dtred2(1:Nstr:2)

  Evecm(1:Nup,1:Nup) = matmul(Alpha + Beta, Evecp(1:Nup,1:Nup)) / spread(Dtred2(1:Nstr:2), 1, Nup)
  Dp(1:Nup,1:Nstr:2) = (Evecp(1:Nup,1:Nup) + Evecm(1:Nup,1:Nup)) * 0.5_Kreal2
  Dp((Nup+1):Nstr,1:Nstr:2) = (Evecp(1:Nup,1:Nup) - Evecm(1:Nup,1:Nup)) * 0.5_Kreal2
  Dp(1:Nup,2:Nstr:2) = Dp((Nup+1):Nstr,1:Nstr:2)
  Dp((Nup+1):Nstr,2:Nstr:2) = Dp(1:Nup,1:Nstr:2)

endsubroutine Deigendst

! -------------------------------------------------------------------
! Evtrisym 

subroutine Evtrisym(D1, E1, ZZ)
  !USE : Cab
  implicit none
  real(Kreal2), dimension(:), intent(inout):: D1, E1
  real(Kreal2), dimension(:,:), optional, intent(inout):: ZZ
  integer(Kint):: i, Iter, l, M, N, Ndum
  real(Kreal2):: B0, C0, DD1, G0, F0, R0, P0, Q0
  real(Kreal2), dimension(size(E1)):: Ff

  N = size(D1)
  if (present(ZZ)) Ndum = N
  E1(:) = eoshift(E1(:),1)
  do l = 1, N
    Iter = 0
    Iterate: do
      do M = l, N-1
        DD1 = abs(D1(M)) + abs(D1(M+1))
        if (abs(E1(M))+DD1 == DD1) exit
      enddo
      if (M==l) exit Iterate
      if (Iter==30) then
        print*, 'Too many iterations, Iter = ', Iter
        stop 'Stopped by Evtrisym.'
      endif
      Iter = Iter + 1
      F0 = (D1(l+1) - D1(l)) / (2.0_Kreal2 * E1(l))
      P0 = Cab(F0,1.0_Kreal2)
      F0 = D1(M) - D1(l) + E1(l) / (F0 + sign(P0,F0))
      Q0 = 1.0_Kreal2
      C0 = 1.0_Kreal2
      R0 = 0.0_Kreal2
      do i = M-1, l, -1
        G0 = Q0 * E1(i)
        B0 = C0 * E1(i)
        P0 = Cab(G0,F0)
        E1(i+1) = P0
        if (P0 == 0.0_Kreal2) then
          D1(i+1) = D1(i+1) - R0
          E1(M) = 0.0_Kreal2
          cycle Iterate
        endif
        Q0 = G0 / P0
        C0 = F0 / P0
        F0 = D1(i+1) - R0
        P0 = (D1(i) - F0) * Q0 + 2.0_Kreal2 * C0 * B0
        R0 = Q0 * P0
        D1(i+1) = F0 + R0
        F0 = C0 * P0 - B0
        if (present(ZZ)) then
          Ff(1:N) = ZZ(1:N,i+1)
          ZZ(1:N,i+1) = Q0 * ZZ(1:N,i) + C0 * Ff(1:N)
          ZZ(1:N,i) = C0 * ZZ(1:N,i) - Q0 * Ff(1:N)
        endif
      enddo
      D1(l) = D1(l) - R0
      E1(l) = F0
      E1(M) = 0.0_Kreal2
    enddo Iterate
  enddo

endsubroutine Evtrisym


!-------------------------------------------------------------------- - 
! Fbsub

subroutine Fbsub_r(X, Ind, Y)
    
  implicit none
  real(Kreal1), dimension(:,:), intent(in):: X
  integer(Kint), dimension(:), intent(in):: Ind
  real(Kreal1), dimension(:), intent(inout):: Y
  integer(Kint):: i, N, ii, ll
  real(Kreal1):: Summ

  N = size(X,1)
  ii = 0
  do i = 1, N
    ll = Ind(i)
    Summ = Y(ll)
    Y(ll) = Y(i)
    if (ii /= 0) then
      Summ = Summ - dot_product(X(i,ii:i-1),Y(ii:i-1))
    elseif (Summ /= 0.0_Kreal1) then
      ii = i
    endif
    Y(i) = Summ
  enddo
  do i = N, 1, -1
    Y(i) = (Y(i) - dot_product(X(i,i+1:N),Y(i+1:N))) / X(i,i)
  enddo
endsubroutine Fbsub_r

subroutine Fbsub_d(X,Ind,Y)
  
  implicit none
  real(Kreal2), dimension(:,:), intent(in):: X
  integer(Kint), dimension(:), intent(in):: Ind
  real(Kreal2), dimension(:), intent(inout):: Y
  integer(Kint):: i, N, ii, ll
  real(Kreal2):: Summ

  N = size(X,1)
  ii = 0
  do i = 1,N
    ll = Ind(i)
    Summ = Y(ll)
    Y(ll) = Y(i)
    if (ii /= 0) then
      Summ = Summ - dot_product(X(i,ii:i-1),Y(ii:i-1))
    elseif (Summ /= 0.0_Kreal2) then
      ii = i
    endif
    Y(i) = Summ
  enddo
  do i = N, 1, -1
    Y(i) = (Y(i) - dot_product(X(i,i+1:N),Y(i+1:N))) / X(i,i)
  enddo
endsubroutine Fbsub_d

!--------------------------------------------------------------------
! Identity

function Identity_r (N, R)
  implicit none
  real(Kreal1):: R
  integer(Kint), intent(in):: N
  real(Kreal1):: Identity_r(1:N,1:N)

  integer(Kint):: i

  Identity_r = 0.0_Kreal1
  do i = 1, N
    Identity_r(i,i) = R
  enddo

endfunction Identity_r

function Identity_d (N, R)
  implicit none
  real(Kreal2):: R
  integer(Kint), intent(in):: N
  real(Kreal2):: Identity_d(1:N,1:N)

  integer(Kint):: i

  Identity_d = 0.0_Kreal2
  do i = 1, N
    Identity_d(i,i) = R
  enddo

endfunction Identity_d

! ------------------------------------------------------------------ - 
  ! Indmax1
function Indmax1_r(Arr)
  real(Kreal1), dimension(:), intent(in):: Arr
  integer(Kint):: Indmax1_r
  integer(Kint), dimension(1):: Imax
  Imax = maxloc(Arr(:))
  Indmax1_r = Imax(1)
end function Indmax1_r

function Indmax1_d(Arr)
  real(Kreal2), dimension(:), intent(in):: Arr
  integer(Kint):: Indmax1_d
  integer(Kint), dimension(1):: Imax
  Imax = maxloc(Arr(:))
  Indmax1_d = Imax(1)
end function Indmax1_d

! -------------------------------------------------------------------
! Invmx

function Invmx_r(XX)
  !use: Ludcp, Fbsub, Identity
  implicit none
  real(Kreal1), intent(in):: XX(:,:)
  real(Kreal1):: Invmx_r(size(XX,1),size(XX,2))

  real(Kreal1):: D0, XX0(size(XX,1),size(XX,2))
  integer(Kint):: Ind(size(XX,1)), j

  Invmx_r = Identity(size(XX,1),1._Kreal1)
  XX0 = XX

  call Ludcp(XX0, Ind, D0)
  do j = 1, size(Invmx_r,2)
    call Fbsub(XX0, Ind, Invmx_r(:,j))
  enddo

endfunction Invmx_r

function Invmx_d(XX)
  !use: Ludcp, Fbsub
  implicit none
  real(Kreal2), intent(in):: XX(:,:)
  real(Kreal2):: Invmx_d(size(XX,1),size(XX,2))

  real(Kreal2):: D0, XX0(size(XX,1),size(XX,2))
  integer(Kint):: Ind(size(XX,1)), j

  Invmx_d = Identity(size(XX,1),1._Kreal2)
  XX0 = XX

  call Ludcp(XX0, Ind, D0)
  do j = 1, size(Invmx_d,2)
    call Fbsub(XX0, Ind, Invmx_d(:,j))
  enddo

endfunction Invmx_d

! ------------------------------------------------------------------ - 
! Ludcp

subroutine Ludcp_r(X, Ind, D0)
  !USE: Indmax1, Abba
  implicit none
  real(Kreal1), dimension(:,:), intent(inout):: X
  integer(Kint), dimension(:), intent(out):: Ind
  real(Kreal1), intent(out):: D0
  real(Kreal1), dimension(size(X,1)):: Vv
  real(Kreal1), parameter:: Tiny = 1.0e-20_Kreal1
  integer(Kint):: j, N, Imax

  N = size(X,1)
  D0 = 1.0_Kreal1
  Vv = maxval(abs(X),dim = 2)
  if (any(Vv == 0.0_Kreal1)) then
    print * , 'Error: singular matrix, Vv = ', Vv
    stop 'Stopped by Ludcp_r.'
  endif
  Vv = 1.0_Kreal1 / Vv
  do j = 1, N
    Imax = j - 1  +  Indmax1(Vv(j:N) * abs(X(j:N,j)))
    if (j /= Imax) then
      call Abba(X(Imax,:),X(j,:))
      D0 =  - D0
      Vv(Imax) = Vv(j)
    endif
    Ind(j) = Imax
    if (X(j,j) ==0.0) X(j,j) = Tiny
    X(j+1:N,j) = X(j+1:N,j) / X(j,j)
    X(j+1:N,j+1:N) = X(j+1:N,j+1:N) - spread(X(j+1:N,j),dim=2,  &
      ncopies=N-j) * spread(X(j,j+1:N),dim=1,ncopies=N-j)
  enddo
endsubroutine Ludcp_r

subroutine Ludcp_d(X, Ind, D0)
  !USE : Indmax1, Abba
  implicit none
  real(Kreal2), dimension(:,:), intent(inout):: X
  integer(Kint), dimension(:), intent(out):: Ind
  real(Kreal2), intent(out):: D0
  real(Kreal2), dimension(size(X,1)):: Vv
  real(Kreal2), parameter:: Tiny = 1.0e-20_Kreal2
  integer(Kint):: j, N, Imax

  N = size(X,1)
  D0 = 1.0_Kreal2
  Vv = maxval(abs(X),dim = 2)
  if (any(Vv == 0.0_Kreal2)) then
    print * , 'Error: singular matrix, Vv = ', Vv
    stop 'Stopped by Ludcp_d.'
  endif
  Vv = 1.0_Kreal2 / Vv
  do j = 1, N
    Imax = j - 1  +  Indmax1(Vv(j:N) * abs(X(j:N,j)))
    if (j /= Imax) then
      call Abba(X(Imax,:),X(j,:))
      D0 = - D0
      Vv(Imax) = Vv(j)
    endif
    Ind(j) = Imax
    if (X(j,j) == 0.0_Kreal2) X(j,j) = Tiny
    X(j+1:N,j) = X(j+1:N,j) / X(j,j)
    X(j+1:N,j+1:N) = X(j+1:N,j+1:N) - spread(X(j+1:N,j),dim=2,  &   
      ncopies=N-j) * spread(X(j,j+1:N),dim=1,ncopies=N-j)
  enddo
endsubroutine Ludcp_d

!--------------------------------------------------------------------
! Rdcsym

subroutine Rdcsym(XX, D1, E1, Novec)

  implicit none
  real(Kreal2), dimension(:,:), intent(inout):: XX
  real(Kreal2), dimension(:), intent(out):: D1, E1
  logical, optional, intent(in):: Novec
  integer(Kint):: i, j, l, N
  real(Kreal2):: G0, H0, F0, Hh, Scale
  real(Kreal2), dimension(size(XX,1)):: Gg
  logical:: Yesvec

  N = size(XX,1)
  if (present(Novec)) then
    Yesvec = .not. Novec
  else
    Yesvec = .true.
  endif
  do i = N, 2, -1
    l = i - 1
    F0 = 0.0_Kreal2
    if (l > 1) then
      Scale = sum(abs(XX(i,1:l)))
      if (Scale==0.0_Kreal2) then
        E1(i) = XX(i,l)
      else
        XX(i,1:l) = XX(i,1:l) / Scale
        F0 = sum(XX(i,1:l) ** 2)
        G0 = XX(i,l)
        H0 = - sign(sqrt(F0),G0)
        E1(i) = Scale * H0
        F0 = F0 - G0 * H0
        XX(i,l) = G0 - H0
        if (Yesvec) XX(1:l,i) = XX(i,1:l) / F0
        do j = 1,l
          E1(j) = (dot_product(XX(j,1:j),XX(i,1:j)) &
            + dot_product(XX(j+1:l,j),XX(i,j+1:l))) / F0
        enddo
        G0 = dot_product(E1(1:l),XX(i,1:l))
        Hh = G0 / (F0 + F0)
        E1(1:l) = E1(1:l) - Hh * XX(i,1:l)
        do j = 1,l
          XX(j,1:j) = XX(j,1:j) - XX(i,j) * E1(1:j) - E1(j) * XX(i,1:j)
        enddo
      endif
    else
      E1(i) = XX(i,l)
    endif
    D1(i) = F0
  enddo
  if (Yesvec) D1(1) = 0.0_Kreal2
  E1(1) = 0.0_Kreal2
  do i = 1,N
    if (Yesvec) then
      l = i - 1
      if (D1(i) /= 0.0_Kreal2) then
        Gg(1:l) = matmul(XX(i,1:l),XX(1:l,1:l))
        XX(1:l,1:l) = XX(1:l,1:l) - spread(XX(1:l,i),dim=2,ncopies=l) &
         * spread(Gg(1:l),dim=1,ncopies=l)
      endif
      D1(i) = XX(i,i)
      XX(i,i) = 1.0_Kreal2
      XX(i,1:l) = 0.0_Kreal2
      XX(1:l,i) = 0.0_Kreal2
    else
      D1(i) = XX(i,i)
    endif
  enddo

endsubroutine Rdcsym

end module Gdrtesv0