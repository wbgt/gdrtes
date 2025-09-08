module Hshmv4     ! v4.4.3

  use Gdrtesv0, only: Kreal, Kreal2, Kint, Gdrtes, Identity, &
  & Invmx, Pi, Pi_d, Amirror

  implicit none

  private

  public:: Nm, Ng, Nh, Lup, Ldn, Numu, Nphi, Umu, Phi, & ! input
  & Hshminit, Hshminit1, Hshmfin, Hshm, Planckint, & !Maxm, Maxg, Maxl, &
  & Lmup, Lmdn, Itild ! Lmup, Lmdn only for diagnostic purpose, can be removed.

  ! Nm: user-specified number of azimuth expansion (m=0~Nm-1=Mm<=Maxm)
  ! Lup and Ldn: user-specified maximum number of zenith expansion
  !     >= 0, triangular or trapezoidal truncation,
  !            number of terms, Lmup(dn), decreases with m
  !     <0,   parallelogram truncation; number of terms, -Lup(dn)
  ! Ng: user-specified number of phase function terms (Ng-1<=Maxg)
  ! Nh: user-specified number terms of BDRF of surface (<=Mm assumed)
  integer(Kint):: Nm=64, Lup=0, Ldn=0, Ng=64, Nh=1, Numu, Nphi !default values can be changed

  ! Maxm: maximum order of azimuth expansion
  ! Maxl: maximum number of zenith expansion in a hemisphere
  ! Maxg: maxiumu order of phase function
  !integer(Kint), parameter:: Maxm=420, Maxl=Maxm/2, Maxg=Maxm, &  ! Maxm up to 320
  !& Maxlg=max(Maxl*2,Maxg)

  integer(Kint):: Mm, Mg, Ml
  integer(Kint), allocatable:: Lmup(:), Lmdn(:)

  !real(Kreal2):: Clm(0:Maxlg,0:Maxm), Pllm(0:Maxlg+1,0:Maxlg+1,0:Maxm), &
  !& Ullm(0:Maxl*2+Maxm,0:Maxl*2+Maxm,0:Maxm)
  real(Kreal2), allocatable:: Clmh(:,:), Pllm(:,:,:), Ullm(:,:,:), Ullmh(:,:,:)
  real(Kreal), allocatable:: Umu(:), Phi(:)

  !real(Kreal2):: &
  !& AA(Maxl*2,Maxl*2,0:Maxm), AAinv(Maxl*2,Maxl*2,0:Maxm), ATint(Maxl*2), &
  !& UUinv(Maxl*2,Maxl*2,0:Maxm), PP(Maxlg+1,Maxl*2,0:Maxm)
  real(Kreal2), allocatable:: AA(:,:), ATint(:), UUinv(:,:,:), PP(:,:,:)
  real(Kreal2), allocatable:: Legr(:,:,:), Cosphi(:,:)
  real(Kreal), allocatable:: Itild(:,:,:) !(Maxl*2,0:Nlyr,0:Mm)

  contains

!-----------------------------------------------------------

subroutine Hshminit(First)
  ! If First, the whole subroutine will run; else, only those related to
  ! Nm, Ng, Nh, Lup, Ldn will be affected. That is, if not First, the changes
  ! of Numu, Nphi, Umu, Phi will be ignored.

  logical, intent(in):: First

  integer(Kint):: i, j, l, m, Ntmp, Nup, Ndn, Nd1, Nstr, Mlg, Lmm
  !real(Kreal2), allocatable:: X1(Maxlg+Maxm), W1(Maxlg+Maxm), &
  !& P2(Maxlg+Maxm,0:Maxlg+Maxm), UU(Maxl*2,Maxl*2,0:Maxm)
  real(Kreal2), allocatable:: Clm(:,:)
  real(Kreal2), allocatable:: X1(:), W1(:), P2(:,:)
  real(Kreal2):: Temp1, Temp2

  integer(Kint), save:: Nm0, Ng0, Lup0, Ldn0
  
  Mm = Nm - 1
  Mg = Ng - 1
  
  ! To be replaced with a subroutine.
  !if (Mm > Maxm) then
  !  print*, 'Error in Hshminit: Mm > Maxm, ', Mm, Maxm
  !  stop
  !endif
  
  !if (Mg < Mm) then
  !  print*, 'Error in Hshint: Mg < Mm, ', Mg, Mm
  !  stop
  !endif
  
 !Lm(0:Maxm) = length of A^{pm} ! not maximum l for each m

  if (First) then
    allocate(Lmup(0:Mm), Lmdn(0:Mm))
  endif

  if (Lup.ge.0) then
    !Lmup(0:Mm) = (Nm + Lup - (/(i,i=0,Mm)/) + 1) /2
    Lmup(0:Mm) = (Nm - (/(i,i=0,Mm)/) + 1) /2  + Lup
  else
    Lmup(0:Mm) = - Lup
  endif
  if (Ldn.ge.0) then
    !Lmdn(0:Mm) = (Nm + Ldn - (/(i,i=0,Mm)/) + 1) /2
    Lmdn(0:Mm) = (Nm - (/(i,i=0,Mm)/) + 1) /2 + Ldn
  else
    Lmdn(0:Mm) = - Ldn
  endif

  if (Lup.eq.Ldn) then
    Amirror = .true.
  else
    Amirror = .false.
  endif
  !Amirror = .false.

  Mlg = max(2*Lmup(0)-2, 2*Lmdn(0)-2, Mg, Nh-1)
  Ml = max(Lmup(0), Lmdn(0))

  ! UUinv, PP, ATint
  if (allocated(ATint)) then
    deallocate(ATint, UUinv, PP)
  endif
  allocate(ATint(Ml*2), UUinv(Ml*2,Ml*2,0:Mm), PP(Mlg+1,Ml*2,0:Mm))

  if (First) then

    ! global allocated arrays
    allocate(Clm(0:Mlg+1,0:Mm), Pllm(0:Mlg+1,0:Mlg+1,0:Mm), &
    & Ullm(0:Ml*2+Mm,0:Ml*2+Mm,0:Mm))
    allocate(Clmh(0:Nh,0:Mm), Ullmh(0:Nh,0:Ml*2+Mm,0:Nh))

    ! local allocated arrays
    allocate(X1(Mlg+Mm+1), W1(Mlg+Mm+1), P2(Mlg+Mm+1,0:Mlg+Mm+1))

    Clm(0,0) = sqrt(2.0_Kreal2)
    do l = 1, Mlg+1
      Clm(l,0) = sqrt((2*l+1.0_Kreal2) / (2*l-1)) * Clm(l-1,0)
    enddo

    do m = 1, Mm
      Mlg = max(2*Lmup(m)+m-2, 2*Lmdn(m)+m-2, Mg, Nh-1)
      Clm(0,m) = sqrt((2*m+1.0_Kreal2)/(2*m)) / (2*m-1) * &
      & Clm(0,m-1)
      do l = m+1, Mlg+1
        Clm(l-m,m) = sqrt((2*l+1.0_Kreal2)/(2*l-1)*(l-m)/(l+m)) &
        &  * Clm(l-1-m,m)
      enddo
    enddo
    clmh(0:Nh,0:Mm) = Clm(0:Nh,0:Mm)
    
    do m = 0, Mm
      Mlg = max(2*Lmup(m)+m-2, 2*Lmdn(m)+m-2, Mg, Nh-1)
      Ntmp = Mlg + 1
      call Glquad(0.0_Kreal2, 1.0_Kreal2, X1(1:Ntmp), W1(1:Ntmp))
      P2(1:Ntmp,m:Mlg+1) = Phat(Mlg+1,m,X1(1:Ntmp))
      do l = m, Mlg+1
        Pllm(l-m,l-m,m) = 1.0_Kreal2
        do i = l+2, Mlg+1, 2
          Pllm(l-m,i-m,m) = 0.0_Kreal2
          Pllm(i-m,l-m,m) = 0.0_Kreal2
        enddo
        do i = l+1, Mlg+1, 2
          Pllm(l-m,i-m,m) = sum(P2(1:Ntmp,l) &
          & * W1(1:Ntmp) * P2(1:Ntmp,i)) * 0.5_Kreal2
          Pllm(i-m,l-m,m) = Pllm(l-m,i-m,m)
        enddo
      enddo
    enddo
  
    do m = 0, Mm
      ! when l = m
      Temp1 = sqrt((m+m+1.0_Kreal2)/((2*m+1)*(2*m+3)))
      Mlg = max(2*Lmup(m)+m-2, 2*Lmdn(m)+m-2, Nh-1)
      do i = m, Mlg, 2
        Ullm(0,i-m,m) = Pllm(1,i-m,m) * Temp1
        Ullm(i-m,0,m) = Ullm(0,i-m,m)
      enddo
      i = m + 1
      Ullm(0,i-m,m) = Pllm(1,i-m,m) * Temp1
      Ullm(i-m,0,m) = Ullm(0,i-m,m)
      do i = m+3, Mlg, 2
        Ullm(0,i-m,m) = 0.0_Kreal2
        Ullm(i-m,0,m) = Ullm(0,i-m,m)
      enddo
    ! when l > m
      do l = m+1, Mlg
        ! when l' = l + 1
        Temp1 = sqrt((l-m+1.0_Kreal2)*(l+m+1)/((2*l+1)*(2*l+3)))
        Temp2 = sqrt(real((l-m)*(l+m),Kreal2) / ((2*l-1)*(2*l+1)))
        Ullm(l-m,l+1-m,m) = Pllm(l+1-m,l+1-m,m) * Temp1
        Ullm(l+1-m,l-m,m) = Ullm(l-m,l+1-m,m)
        ! when l-l' is even
        do i = l, Mlg, 2
          Ullm(l-m,i-m,m) = Pllm(l+1-m,i-m,m)*Temp1 + Pllm(l-1-m,i-m,m)*Temp2
          Ullm(i-m,l-m,m) = Ullm(l-m,i-m,m)
        enddo
        ! when l-l' is old
        do i = l+3, Mlg, 2
          Ullm(l-m,i-m,m) = 0.0_Kreal2
          Ullm(i-m,l-m,m) = Ullm(l-m,i-m,m)
        enddo
      enddo
    enddo
    Ullmh(0:Nh,0:Ml*2+Mm,0:Nh) = Ullm(0:Nh,0:Ml*2+Mm,0:Nh)

    ! AA, AAinv
    allocate(AA(Ml*2,Ml*2))
    AA(1:Ml*2,1:Ml*2) = Identity(Ml*2, 1.0_Kreal2)
    !AAinv(1:Nstr,1:Nstr,m) = AA(1:Nstr,1:Nstr,m)

    Nm0 = Nm
    Ng0 = Ng
    Lup0 = Lup
    Ldn0 = Ldn

    deallocate(X1, W1, P2, Clm)

    allocate(Legr(Numu,0:Ml*2+Mm,0:Mm), Cosphi(Nphi,0:Mm))

    do m = 0, Mm
      Cosphi(:,m) = cos(m*Phi(:))
      Lmm = max(Lmup(m), Lmdn(m))
      l = m + 2*Lmm - 1
      Legr(:,m:l,m) = Phat(l,m,real(abs(Umu(1:Numu)), Kreal2))
    enddo

  endif

  ! Checking array sizes
  if (Nm > Nm0) then
    print*, 'Nm > Nm0 in Hshminit: ', Nm, Nm0
    stop
  endif
  if (Ng > Ng0) then
    print*, 'Ng > Ng0 in Hshminit: ', Nm, Nm0
    stop
  endif
  if (Lup0>=0) then
    if (Lup>Lup0 .or. Lup<0) then
      print*, 'Lup exceeds initial setting:', Lup, Lup0
      stop
    endif
  else
    if (Lup>=0 .or. Lup<Lup0) then 
      print*, 'Lup exceeds initial setting:', Lup, Lup0
      stop
    endif
  endif

  do m = 0, Mm
  
    Nup = Lmup(m)
    Ndn = Lmdn(m)
    Nd1 = Nup + 1
    Nstr = Nup + Ndn

    UUinv(1:Nstr,1:Nstr,m) = 0.0_Kreal2
    UUinv(1:Nup,1:Nup,m) = Invmx(Ullm(0:2*Nup-2:2,0:2*Nup-2:2,m))
    UUinv(Nd1:Nstr,Nd1:Nstr,m) = - Invmx(Ullm(0:2*Ndn-2:2,0:2*Ndn-2:2,m))

    Mlg = max(Ng-m, Nh-m)
    do i = 1, Mlg
      do j = 1, Nup
        PP(i,j,m) = Pllm(i-1,2*j-2,m)
      enddo
      do j = 1, Ndn
        PP(i,Nup+j,m) = Pllm(i-1,2*j-2,m) * (-1)**(i+1)
      enddo
    enddo

  enddo

  ! Atint for m=0. The index is indeed l+1
  ATint(1:Lmup(0)+Lmup(0)) = 0.0_Kreal2
  ATint(1) = sqrt(2.0_Kreal2) * 0.5_Kreal2
  ATint(Lmup(0)+1) = ATint(1)

endsubroutine Hshminit

subroutine Hshminit1()
  ! A light versition of Hshminit, trying to save memory, it runs only once.

  integer(Kint):: i, j, l, m, Ntmp, Nup, Ndn, Nd1, Nstr, Mlg, Lmm
  !real(Kreal2), allocatable:: X1(Maxlg+Maxm), W1(Maxlg+Maxm), &
  !& P2(Maxlg+Maxm,0:Maxlg+Maxm), UU(Maxl*2,Maxl*2,0:Maxm)
  real(Kreal2), allocatable:: Clm(:,:)
  real(Kreal2), allocatable:: X1(:), W1(:), P2(:,:)
  real(Kreal2):: Temp1, Temp2

  integer(Kint), save:: Nm0, Ng0, Lup0, Ldn0
  
  Mm = Nm - 1
  Mg = Ng - 1
  
  ! To be replaced with a subroutine.
  !if (Mm > Maxm) then
  !  print*, 'Error in Hshminit: Mm > Maxm, ', Mm, Maxm
  !  stop
  !endif
  
  !if (Mg < Mm) then
  !  print*, 'Error in Hshint: Mg < Mm, ', Mg, Mm
  !  stop
  !endif
  
 !Lm(0:Maxm) = length of A^{pm} ! not maximum l for each m

  allocate(Lmup(0:Mm), Lmdn(0:Mm))

  if (Lup.ge.0) then
    Lmup(0:Mm) = (Nm + Lup - (/(i,i=0,Mm)/) + 1) /2
  else
    Lmup(0:Mm) = - Lup
  endif
  if (Ldn.ge.0) then
    Lmdn(0:Mm) = (Nm + Ldn - (/(i,i=0,Mm)/) + 1) /2
  else
    Lmdn(0:Mm) = - Ldn
  endif

  if (Lup.eq.Ldn) then
    Amirror = .true.
  else
    Amirror = .false.
  endif
  !Amirror = .false.

  Mlg = max(2*Lmup(0)-2, 2*Lmdn(0)-2, Mg, Nh-1)
  Ml = max(Lmup(0), Lmdn(0))

  ! UUinv, PP, ATint
  !if (allocated(ATint)) then
  !  deallocate(ATint, UUinv, PP)
  !endif
  allocate(ATint(Ml*2), UUinv(Ml*2,Ml*2,0:Mm), PP(Mlg+1,Ml*2,0:Mm))

  !if (First) then

    ! global allocated arrays
    allocate(Clm(0:Mlg+1,0:Mm), Pllm(0:Mlg+1,0:Mlg+1,0:Mm), &
    & Ullm(0:Ml*2+Mm,0:Ml*2+Mm,0:Mm))
    allocate(Clmh(0:Nh,0:Nh), Ullmh(0:Nh,0:Ml*2+Mm,0:Nh))

    ! local allocated arrays
    allocate(X1(Mlg+Mm+1), W1(Mlg+Mm+1), P2(Mlg+Mm+1,0:Mlg+Mm+1))

    Clm(0,0) = sqrt(2.0_Kreal2)
    do l = 1, Mlg+1
      Clm(l,0) = sqrt((2*l+1.0_Kreal2) / (2*l-1)) * Clm(l-1,0)
    enddo

    do m = 1, Mm
      Mlg = max(2*Lmup(m)+m-2, 2*Lmdn(m)+m-2, Mg, Nh-1)
      Clm(0,m) = sqrt((2*m+1.0_Kreal2)/(2*m)) / (2*m-1) * &
      & Clm(0,m-1)
      do l = m+1, Mlg+1
        Clm(l-m,m) = sqrt((2*l+1.0_Kreal2)/(2*l-1)*(l-m)/(l+m)) &
        &  * Clm(l-1-m,m)
      enddo
    enddo
    clmh(0:Nh,0:Nh) = Clm(0:Nh,0:Nh)
    
    do m = 0, Mm
      Mlg = max(2*Lmup(m)+m-2, 2*Lmdn(m)+m-2, Mg, Nh-1)
      Ntmp = Mlg + 1
      call Glquad(0.0_Kreal2, 1.0_Kreal2, X1(1:Ntmp), W1(1:Ntmp))
      P2(1:Ntmp,m:Mlg+1) = Phat(Mlg+1,m,X1(1:Ntmp))
      do l = m, Mlg+1
        Pllm(l-m,l-m,m) = 1.0_Kreal2
        do i = l+2, Mlg+1, 2
          Pllm(l-m,i-m,m) = 0.0_Kreal2
          Pllm(i-m,l-m,m) = 0.0_Kreal2
        enddo
        do i = l+1, Mlg+1, 2
          Pllm(l-m,i-m,m) = sum(P2(1:Ntmp,l) &
          & * W1(1:Ntmp) * P2(1:Ntmp,i)) * 0.5_Kreal2
          Pllm(i-m,l-m,m) = Pllm(l-m,i-m,m)
        enddo
      enddo
    enddo
  
    do m = 0, Mm
      ! when l = m
      Temp1 = sqrt((m+m+1.0_Kreal2)/((2*m+1)*(2*m+3)))
      Mlg = max(2*Lmup(m)+m-2, 2*Lmdn(m)+m-2, Nh-1)
      do i = m, Mlg, 2
        Ullm(0,i-m,m) = Pllm(1,i-m,m) * Temp1
        Ullm(i-m,0,m) = Ullm(0,i-m,m)
      enddo
      i = m + 1
      Ullm(0,i-m,m) = Pllm(1,i-m,m) * Temp1
      Ullm(i-m,0,m) = Ullm(0,i-m,m)
      do i = m+3, Mlg, 2
        Ullm(0,i-m,m) = 0.0_Kreal2
        Ullm(i-m,0,m) = Ullm(0,i-m,m)
      enddo
    ! when l > m
      do l = m+1, Mlg
        ! when l' = l + 1
        Temp1 = sqrt((l-m+1.0_Kreal2)*(l+m+1)/((2*l+1)*(2*l+3)))
        Temp2 = sqrt(real((l-m)*(l+m),Kreal2) / ((2*l-1)*(2*l+1)))
        Ullm(l-m,l+1-m,m) = Pllm(l+1-m,l+1-m,m) * Temp1
        Ullm(l+1-m,l-m,m) = Ullm(l-m,l+1-m,m)
        ! when l-l' is even
        do i = l, Mlg, 2
          Ullm(l-m,i-m,m) = Pllm(l+1-m,i-m,m)*Temp1 + Pllm(l-1-m,i-m,m)*Temp2
          Ullm(i-m,l-m,m) = Ullm(l-m,i-m,m)
        enddo
        ! when l-l' is old
        do i = l+3, Mlg, 2
          Ullm(l-m,i-m,m) = 0.0_Kreal2
          Ullm(i-m,l-m,m) = Ullm(l-m,i-m,m)
        enddo
      enddo
    enddo
    Ullmh(0:Nh,0:Ml*2+Mm,0:Nh) = Ullm(0:Nh,0:Ml*2+Mm,0:Nh)

    ! AA, AAinv
    allocate(AA(Ml*2,Ml*2))
    AA(1:Ml*2,1:Ml*2) = Identity(Ml*2, 1.0_Kreal2)
    !AAinv(1:Nstr,1:Nstr,m) = AA(1:Nstr,1:Nstr,m)

    Nm0 = Nm
    Ng0 = Ng
    Lup0 = Lup
    Ldn0 = Ldn

    deallocate(X1, W1, P2, Clm)

    allocate(Legr(Numu,0:Ml*2+Mm,0:Mm), Cosphi(Nphi,0:Mm))

    do m = 0, Mm
      Cosphi(:,m) = cos(m*Phi(:))
      Lmm = max(Lmup(m), Lmdn(m))
      l = m + 2*Lmm - 1
      Legr(:,m:l,m) = Phat(l,m,real(abs(Umu(1:Numu)), Kreal2))
    enddo

  !endif

  do m = 0, Mm
  
    Nup = Lmup(m)
    Ndn = Lmdn(m)
    Nd1 = Nup + 1
    Nstr = Nup + Ndn

    UUinv(1:Nstr,1:Nstr,m) = 0.0_Kreal2
    UUinv(1:Nup,1:Nup,m) = Invmx(Ullm(0:2*Nup-2:2,0:2*Nup-2:2,m))
    UUinv(Nd1:Nstr,Nd1:Nstr,m) = - Invmx(Ullm(0:2*Ndn-2:2,0:2*Ndn-2:2,m))

    Mlg = max(Ng-m, Nh-m)
    do i = 1, Mlg
      do j = 1, Nup
        PP(i,j,m) = Pllm(i-1,2*j-2,m)
      enddo
      do j = 1, Ndn
        PP(i,Nup+j,m) = Pllm(i-1,2*j-2,m) * (-1)**(i+1)
      enddo
    enddo

  enddo

  ! Atint for m=0. The index is indeed l+1
  ATint(1:Lmup(0)+Lmup(0)) = 0.0_Kreal2
  ATint(1) = sqrt(2.0_Kreal2) * 0.5_Kreal2
  ATint(Lmup(0)+1) = ATint(1)

  deallocate(Pllm, Ullm)

endsubroutine Hshminit1

!--------------------------------------------------------------------

subroutine Hshmfin()
  
  if (allocated(Pllm)) then
    deallocate(Pllm, Ullm)
  endif

  deallocate(AA, ATint, UUinv, PP, Lmup, Lmdn, Clmh, Ullmh, &
  & Legr, Cosphi, Itild)

end subroutine Hshmfin

!--------------------------------------------------------------------

subroutine Hshm(Nlyr, Od0, Ssa0, Pf, I0, Mu0, T, Tsurf, GS, &
  & Thermal, Deltam, Wnumlo, Wnumhi, &
  Beam, Rad, Fup, Fdn, Act, Crt, Cpu)

  integer(Kint), intent(in):: Nlyr!, Numu, Nphi
  real(Kreal), intent(in):: Od0(Nlyr), Ssa0(Nlyr), Pf(0:Mg,Nlyr), I0, Mu0, &
  & T(0:Nlyr), Gs(0:Nh-1), Tsurf, Wnumlo, Wnumhi, Crt
  logical, intent(in):: Deltam, Thermal
  real(Kreal), intent(out):: Beam(0:Nlyr), & !Itild(Maxl*2,0:Nlyr,0:Mm), &
  & Rad(Numu,Nphi,0:Nlyr), Fup(0:Nlyr), Fdn(0:Nlyr), Act(0:Nlyr), Cpu

  integer(Kint):: Ngd, Ngm, i, j, l, m, Nup, Ndn, Nstr, Nhm, Nghm
  real(Kreal):: Od(Nlyr), Ssa(Nlyr), Ssad(Nlyr), &
  & Pfd(0:Mg,Nlyr), Odd(Nlyr), Dump1(Nlyr), B0(0:Nlyr), Hm(Nh,Nh), Bsurf
  real(Kreal2):: Glm(0:Mg,Nlyr), PPfm(Nh,Ml), Pmu0(max(Ng,Nh)), Mu0d, &
  & PPm(max(Ng,Nh),Ml*2), Phat0(1,0:max(Ng,Nh))
  logical:: Thermalm

  real(Kreal):: Icr(2,0:Nlyr), Temp, Cpu1, Cpu2
  integer(Kint):: Ncr, Mcr
  logical:: Lcr

  integer(Kint), save:: Nlyr0, Ml0, Mm0

  ! minimum significant float relative to 1 (about 0.00...1)
  real(Kreal), parameter:: Epsil10 = &
      float(RADIX(1.0_Kreal))**(1-DIGITS(1.0_Kreal)) * 10._Kreal

  call cpu_time(Cpu1)

  if (.not.allocated(Itild)) then
    allocate(Itild(0:Ml*2,0:Nlyr,0:Mm))
    Nlyr0 = Nlyr
    Ml0 = Ml
    Mm0 = Mm
  elseif (Nlyr.ne.Nlyr0.or.Ml.ne.Ml0.or.Mm.ne.Mm0) then
    deallocate(Itild)
    allocate(Itild(0:Ml*2,0:Nlyr,0:Mm))
    Nlyr0 = Nlyr
    Ml0 = Ml
    Mm0 = Mm
  endif
    
  ! Checking Od and Ssa

  Od = Od0
  if (any(Od(1:Nlyr)<0.0_Kreal)) then
    print*, 'Error in Hshmrun: Negative Od:', Od(1:Nlyr)
    stop
  endif
  where (Od < Epsil10)
    Od = Epsil10
  endwhere

  Ssa = Ssa0
  if (any(Ssa(1:Nlyr)<0.0_Kreal .or. Ssa(1:Nlyr)>1.0_Kreal)) then
    print*, 'Error in Hshmrun: Ssa out of range(0,1):', Ssa(1:Nlyr)
    stop
  endif
  where (Ssa(1:Nlyr) > 1.0_Kreal - Epsil10)
    Ssa(1:Nlyr) = 1.0_Kreal - Epsil10
  endwhere

  ! delta-M scaling or not
  ! Ngd is the number of terms of phase function
  if (Deltam) then
    Ngd = Mm + 1
    Dump1 = 1.0_Kreal - Ssa(1:Nlyr) * Pf(Ngd,1:Nlyr)
    Odd(1:Nlyr) = Dump1 * Od(1:Nlyr)
    Ssad(1:Nlyr) = Ssa(1:Nlyr) * (1.0_Kreal - Pf(Ngd,1:Nlyr)) / Dump1
    do l = 1, Nlyr
      Pfd(0:Ngd-1,l) = (Pf(0:Ngd-1,l) - Pf(Ngd,l)) / (1.0_Kreal - Pf(Ngd,l))
    enddo
  else
    Ngd = Ng
    Odd(1:Nlyr) = Od(1:Nlyr)
    Ssad(1:Nlyr) = Ssa(1:Nlyr)
    Pfd(0:Mg,1:Nlyr) = Pf(0:Mg,1:Nlyr)
  endif

  ! in principle, this should be done in Gdrtes, here is for consistency and efficiency
  Beam(0) = I0
  do l = 1, Nlyr
    Beam(l) = Beam(l-1) * exp(-Odd(l) / Mu0)
  enddo

  do m = 0, Mm

    Nup = Lmup(m)
    Ndn = Lmdn(m)
    Nstr = Nup + Ndn

    ! phase function for moment m -- scaled
    Glm(0:Ngd-1,1:Nlyr) = 0.0_Kreal2
    do l = m, Ngd-1
      Glm(l,1:Nlyr) = Pfd(l,1:Nlyr) * 0.5_Kreal2
    enddo

    Ngm = Ngd-m
    Nhm = Nh - m

    Nghm = max(Ngm, Nhm)
    ! PPm original or pre-scaled
    PPm(1:Nghm,1:Nstr) = PP(1:Nghm,1:Nstr,m)

    ! solar beam component m
    Mu0d = Mu0
    Phat0(:,m:Nghm+m-1) = Phat(Nghm+m-1,m, (/-Mu0d/))
    Pmu0(1:Nghm) = Phat0(1,m:Nghm+m-1)
    if (m .ne. 0) then
      Pmu0(1:Nghm) = Pmu0(1:Nghm) * 2
    endif

    ! Thermal source for m = 0
    if (m.eq.0) then
      Thermalm = Thermal
      if (Thermal) then
        do l = 0, Nlyr
          B0(l) = Planckint(Wnumlo, Wnumhi, T(l))
        enddo
        Bsurf = Planckint(Wnumlo, Wnumhi, Tsurf)
      else
        B0(0:Nlyr) = 0.0_Kreal
        Bsurf = 0.0_Kreal
      endif
    else
      Thermalm = .false.
    endif

    ! Surface

    if (Nhm.gt.0) then

      do j = 1, Ndn
        do i = 1, Nhm
          PPfm(i,j) = 2*(-1)**(i-1) / Clmh(i-1,m) * Ullmh(i-1,2*j-2,m)
        enddo
      enddo

      Hm(1:Nhm,1:Nhm) = 0.0_Kreal
      do i = 1, Nhm
        l = m - 1 + i
        Hm(i,i) = (2*l + 1) * Gs(l)
        do j = l-m+1, l+m
          Hm(i,i) = Hm(i,i) / j
        enddo
      enddo

    endif

    call Gdrtes(Nlyr, Nup, Ndn, Ngm, Nhm, AA(1:Nstr,1:Nstr), &
    & AA(1:Nstr,1:Nstr), UUinv(1:Nstr,1:Nstr,m), &
    & Odd, Ssad, Glm(m:Ngd-1,1:Nlyr), PPm(1:Nghm, 1:Nstr), &
    & Beam, Mu0, Pmu0(1:Nghm), Thermalm, ATint(1:Nstr), B0(0:Nlyr), &
    & PPfm(1:Nhm,1:Ndn), Hm(1:Nhm,1:Nhm), Bsurf, &
    & Itild(1:Nstr,0:Nlyr,m))

    ! Truncation of m
    if (m.eq.0) then
      Icr(1,0:Nlyr) = abs(Itild(1,0:Nlyr,0))
      Icr(2,0:Nlyr) = abs(Itild(Nup+1,0:Nlyr,0))
      Ncr = 0
      Mcr = 0
      if (Crt < 0.0_Kreal) then
        Icr(1:2,0:Nlyr) = Icr(1:2,0:Nlyr) * &
        & (Epsil10 * 0.1_Kreal * (Pi * 0.5_Kreal)**0.75_Kreal)
      else
        Icr(1:2,0:Nlyr) = Icr(1:2,0:Nlyr) * &
        & (Crt * 0.1_Kreal * (Pi * 0.5_Kreal)**0.75_Kreal)
      endif
    else

      Lcr = .true.

      do j = 1, Nup
        l = m + j*2 - 2
        Temp = sqrt(2.0_Kreal*(2*l+1)) / real(m,Kreal)**0.25_Kreal
        do l = 0, Nlyr
          if (abs(Itild(j,l,m)) * Temp > Icr(1,l)) then
            Lcr = .false.
            exit
          endif
        enddo
        if (.not.Lcr) exit
      enddo

      if (Lcr) then
        do j = 1, Ndn
          l = m + j*2 -2
          Temp = sqrt(2.0_Kreal*(2*l+1)) / m**0.25_Kreal
          do l = 0, Nlyr
            if (abs(Itild(j+Nup,l,m)) * Temp > Icr(2,l)) then
              Lcr = .false.
              exit
            endif
          enddo
          if (.not.Lcr) exit
        enddo
      endif

      if (Lcr) then

        if (Mcr == m-1) then
          Ncr = Ncr + 1
          Mcr = m
        else
          Ncr = 1
          Mcr = m
        endif
        if (Ncr==2) then
          Mcr = m
          exit
        endif

      else
        Ncr = 0
        Mcr = 0
      endif

    endif

    !print*, 'Gdrtes done for m = ', m, Ncr, Lcr

  enddo

  if (Ncr < 2) Mcr = Mm

  Itild(:,:,Mcr+1:Mm) = 0.0_Kreal

  ! Fluxes
  Nup = Lmup(0)
  Ndn = Lmdn(0)
  Nstr = Nup + Ndn
  !Fup(0:Nlyr) = real(4*Pi / Clmh(1,0) * matmul( &
  !& Pllm(0:2*Nup-2:2,1,0), Itild(1:Nup,0:Nlyr,0) ), Kreal)
  Fup(0:Nlyr) = 4*Pi / sqrt(6.0_Kreal) * real(matmul( &
  & PP(2,1:Nup,0), Itild(1:Nup,0:Nlyr,0) ), Kreal)
  !Fdn(0:Nlyr) = real(4*Pi / Clmh(1,0) * matmul( &
  !& Pllm(0:2*Ndn-2:2,1,0), Itild(Ndn+1:Nstr,0:Nlyr,0)), Kreal)
  Fdn(0:Nlyr) = 4*Pi / sqrt(6.0_Kreal) * real(matmul( &
  & - PP(2,Ndn+1:Nstr,0), Itild(Ndn+1:Nstr,0:Nlyr,0)), Kreal)
  Act(0:Nlyr) = 4*Pi * matmul(real(ATint(1:Nstr),Kreal), &
  & Itild(1:Nstr,0:Nlyr,0)) + Beam(0:Nlyr)

  call cpu_time(Cpu2)
  Cpu = Cpu2 - Cpu1

  ! Intensities for specified directions

  !print*, 'Mcr = ', Mcr
  call Intangle3(Nlyr, Rad, Mcr)


end subroutine Hshm

!--------------------------------------------------------------------
subroutine Intangle3(Nlay, Rad, Mcr)
  ! use : Phat
  implicit none

  integer(Kint), intent(in):: Nlay, Mcr!, Numu, Nphi
  !real(Kreal), intent(in):: Itild(Maxl*2,0:Nlay,0:Mm)!, Umu(Numu), Phi(Nphi)
  real(Kreal), intent(out):: Rad(Numu,Nphi,0:Nlay)

  real(Kreal2):: Im(0:Mcr,0:Nlay) !Cosphi(Nphi,0:Mm)
  integer(Kint):: Nu0, Nd0, m, j

  !do m = 0, Mm
  !  Cosphi(:,m) = cos(m*Phi(:))
  !  Lmm = max(Lmup(m), Lmdn(m))
  !  l = m + 2*Lmm - 1
  !  Legr(:,m:l,m) = Phat(l,m,real(abs(Umu(1:Numu)), Kreal2))
  !enddo

  do j = 1, Numu
    do m = 0, Mcr
      Nu0 = Lmup(m)
      Nd0 = Lmdn(m)
      if (Umu(j).lt.0.0_Kreal) then
        Im(m,:) = matmul(Legr(j,m:m+2*Nd0-2:2,m), Itild(Nu0+1:Nu0+Nd0,:,m))
      else
        Im(m,:) = matmul(Legr(j,m:m+2*Nu0-2:2,m), Itild(1:Nu0,:,m))
      endif
    enddo
    Rad(j,:,:) = real(matmul(Cosphi(:,0:Mcr), Im(0:Mcr,:)), Kreal)
  enddo  

endsubroutine Intangle3

!---------------------------------------------------------------
function Phat(L, M, Mu)
  implicit none
  integer(Kint), intent(in):: L, M
  real(Kreal2), dimension(:), intent(in):: Mu
  real(Kreal2):: Phat(size(Mu),M:L), Somx(size(Mu))

  integer(Kint):: i

  Phat(:,M) = sqrt(2.0_Kreal2)
  if (M > 0) then
    Somx = sqrt((1.0_Kreal2 - Mu) * (1.0_Kreal2 + Mu))
    do i = 2, M*2, 2
      Phat(:,M) = Phat(:,M) * sqrt((i+1.0_Kreal2)/i) * Somx
    enddo
    if (mod(M,2)==1) Phat(:,M) = - Phat(:,M)
  endif

  if (L==M) then
    return
  else
    Phat(:,M+1) = Phat(:,M) * Mu * sqrt(2.0_Kreal2*M + 3)
    if (L==M+1) then
      return
    else
      do i = M+2, L
        Phat(:,i) = sqrt((2.0_Kreal2*i+1)/(i+M)*(2*i-1)/(i-M)) * Mu * &
        & Phat(:,i-1) &
        & - sqrt((2.0_Kreal2*i+1) / (2*i-3) * (i+m-1) / (i+m) * (i-m-1) &
        & / (i-m)) * Phat(:,i-2)
      enddo
    endif
  endif

endfunction Phat

!--------------------------------------------------------------------

function Planckint( WNUMLO, WNUMHI, T )

  ! The program is adapted from DISORT package, computes Planck function
  !   integrated between two wavenumbers.
  !---------------------------------------------------------------------
  !  INPUT :  WNUMLO : Lower wavenumber (inv cm) of spectral interval
  !           WNUMHI : Upper wavenumber
  !           T      : Temperature (K)
  !
  !  OUTPUT : Planckint : Integrated Planck function ( Watts/sq m )
  !--------------------------------------------------------------------
  implicit none
  real(Kreal), intent(in):: T, WNUMHI, WNUMLO
  real(Kreal):: Planckint

  REAL(Kreal), PARAMETER:: A1 = 1. / 3., A2 = -1. / 8., A3 = 1. / 60., &
            A4 = -1. / 5040., A5 = 1. / 272160., &
            A6 = -1. / 13305600.

  INTEGER(Kint):: I, K, M, MMAX, N, SMALLV
  REAL(Kreal):: C2 = 1.438786, CONC, DEL, EPSIL, EX, EXM, HH, MV, &
    OLDVAL, PI0 = 0.0, SIGDPI, SIGMA = 5.67032E-8, VAL1, VAL0, &
    VCUT = 1.5, VMAX, VSQ

  REAL(Kreal):: D(2), P(2), V(2), &
                VCP(7) = (/10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0/)
  !     ..
  SAVE::    PI0, CONC, VMAX, EPSIL, SIGDPI

  IF( PI0 .EQ. 0.0_Kreal ) THEN

    PI0     = 2._Kreal * ASIN( 1.0_Kreal )
    Vmax = log(huge(1.0_Kreal))
    Epsil = real(radix(1.0_Kreal), kind=Kreal) ** (1 - digits(1.0_Kreal))
    SIGDPI = SIGMA / PI0
    CONC   = 15._Kreal / PI0**4

  ENDIF

  IF( T.LT.0.0_Kreal .OR. WNUMHI.LE.WNUMLO .OR. WNUMLO.LT.0._Kreal ) &
    then
    write(*,*) 'T, WNUMHI, WNUMLO = ', T, WNUMHI, WNUMLO
    stop 'Planckint--temperature or wavenums. wrong'
  endif

  IF( T .LT. 1.E-4 ) THEN

    Planckint = 0.0
    RETURN

  ENDIF


  V( 1 ) = C2*WNUMLO / T
  V( 2 ) = C2*WNUMHI / T

  IF( V( 1 ).GT.EPSIL .AND. V( 2 ).LT.VMAX .AND.            &
      ( WNUMHI - WNUMLO ) / WNUMHI .LT. 1.E-2 ) THEN

    HH     = V( 2 ) - V( 1 )
    OLDVAL = 0.0_Kreal
    VAL0   = PLKF( V( 1 ) ) + PLKF( V( 2 ) )

    DO N = 1, 10

      DEL  = HH / ( 2*N )
      VAL1  = VAL0

      DO K = 1, 2*N - 1
        VAL1  = VAL1 + 2*( 1 + MOD( K,2 ) )*PLKF( V( 1 ) + K*DEL )
      enddo

      VAL1  = DEL / 3.*VAL1
      IF( ABS( ( VAL1 - OLDVAL ) / VAL1 ).LE.1.E-6 ) then
        GO TO  30
      endif
      OLDVAL = VAL1

    enddo

    stop 'Planckint--Simpson rule didnt converge'

    30 CONTINUE

    Planckint = SIGDPI * T**4 * CONC * VAL1

    RETURN

  ENDIF

  !                          *** General case ***
  SMALLV = 0

  DO I = 1, 2

    IF( V( I ).LT.VCUT ) THEN
      !                                   ** Use power series
      SMALLV = SMALLV + 1
      VSQ    = V( I )**2
      P( I ) = CONC*VSQ*V( I )*( A1 +                            &
              V( I )*( A2 + V( I )*( A3 + VSQ*( A4 + VSQ*( A5 +  &
              VSQ*A6 ) ) ) ) )

    ELSE
      !                      ** Use exponential series
      MMAX  = 0
      !                                ** Find upper limit of series
  40  CONTINUE
      MMAX  = MMAX + 1

      IF( V(I) .LT. VCP( MMAX ) ) then
        GO TO  40
      endif

      EX     = EXP( - V(I) )
      EXM    = 1.0_Kreal
      D( I ) = 0.0_Kreal

      DO M = 1, MMAX
        MV     = M*V( I )
        EXM    = EX*EXM
        D( I ) = D( I ) + EXM*( 6._Kreal+ MV*( 6._Kreal+ MV*( &
          3._Kreal+ MV ) ) ) / M**4
      enddo

      D( I ) = CONC*D( I )

    ENDIF

  enddo

  !                              ** Handle ill-conditioning
  IF( SMALLV.EQ.2 ) THEN
    !                                 ** WNUMLO and WNUMHI both small
    Planckint = P( 2 ) - P( 1 )

  ELSEIF( SMALLV.EQ.1 ) THEN
    !                                    ** WNUMLO small, WNUMHI large
    Planckint = 1.- P( 1 ) - D( 2 )

  ELSE
    !                                 ** WNUMLO and WNUMHI both large
    Planckint = D( 1 ) - D( 2 )

  ENDIF

  Planckint = SIGDPI * T**4 * Planckint

  IF( Planckint.EQ.0.0_Kreal ) then
    stop 'Planckint--returns zero; possible underflow'
  endif

contains

  function PLKF(X)
    real(Kreal), intent(in):: X
    real(Kreal):: PLKF

    PLKF = X**3 / ( EXP( X ) - 1 )

  end function PLKF

end Function Planckint

! ------------------------------------------------------------------- 
! Glquad: Gaussian quadrature points X1 and weight W1, over [A,B]
subroutine Glquad(A, B, X1, W1)
  !USE : Ainc
  implicit none
  real(Kreal2), intent(in):: A, B
  real(Kreal2), dimension(:), intent(out):: X1, W1
  real(Kreal2), parameter:: Eps = 3.0e-14_Kreal2  !wangbiao
  integer(Kint):: i, j, M, N
  integer(Kint), parameter:: Maxit = 10   !wangbiao
  real(Kreal2):: Xl, Xm
  real(Kreal2), dimension((size(X1)+1)/2):: P1, P2, P3, PPP, Z, Z1
  logical, dimension((size(X1)+1)/2):: No

  N = size(X1)
  M = (N + 1) / 2
  Xm = 0.5_Kreal2 * (B + A)
  Xl = 0.5_Kreal2 * (B - A)
  Z = cos(Pi_d * (Ainc(1,1,M) - 0.25_Kreal2) / (N + 0.5_Kreal2))
  No = .true.
  do i = 1,Maxit
    where (No)
      P1 = 1.0_Kreal2
      P2 = 0.0_Kreal2
    endwhere
    do j = 1,N
      where (No)
        P3 = P2
        P2 = P1
        P1 = ((2.0_Kreal2 * j - 1.0_Kreal2) * Z * P2 - (j - &
              1.0_Kreal2) * P3) / j
      endwhere
    enddo
    where (No)
      PPP = N * (Z * P1 - P2) / (Z * Z - 1.0_Kreal2)
      Z1 = Z
      Z = Z1 - P1 / PPP
      No = (abs(Z - Z1) > Eps)
    endwhere
    if (.not. any(No)) exit
  enddo
  if (i == Maxit+1) then
    print*, 'Too many iterations, i > Maxit: ', i, Maxit
    stop 'Stopped by Glquad.'
  endif

  X1(1:M) = Xm - Xl * Z
  X1(N:N-M+1:-1) = Xm + Xl * Z
  W1(1:M) = 2.0_Kreal2 * Xl / ((1.0_Kreal2 - Z ** 2) * PPP ** 2)
  W1(N:N-M+1:-1) = W1(1:M)

endsubroutine Glquad

! Ainc

function Ainc(First,Incre,N)
  integer(Kint), intent(in):: First, Incre, N
  integer(Kint), dimension(N):: Ainc
  integer(Kint), parameter:: N16 = 16, N8 = 8
  integer(Kint):: k, k2, temp

  if (N > 0) Ainc(1) = First
  if (N <= N16) then
    do k = 2,N
      Ainc(k) = Ainc(k-1) + Incre
    enddo
  else
    do k = 2,N8
      Ainc(k) = Ainc(k-1) + Incre
    enddo
    temp = Incre * N8
    k = N8
    do
      if (k >= N) exit
      k2 = k + k
      Ainc(k+1:min(k2,N)) = temp + Ainc(1:min(k,N-k))
      temp = temp + temp
      k = k2
    enddo
  endif
end function Ainc


end module Hshmv4
