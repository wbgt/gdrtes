!===============================================================================
! Module for the Shifted Jacobi Hemispherical Harmonic (SJHM) method of 
! radiative transfer in plane-parallel media.
! Date: 2026-06-25
! Last revision: 2026-10-10
!-------------------------------------------------------------------------------
module Sjhm_module
  use Mod_Gdrtes, only: Kreal, Kreal2, Kint, Gdrtes, Identity, Invmx, Amirror
  implicit none
  private
  public :: Nm, Ng, Nh, Lup, Ldn, Numu, Nphi, Umu, Phi
  public :: Sjhminit, Sjhmfin, Sjhm
  public :: PI, PId

  !integer, parameter :: Kreal2 = kind(1.0d0), Kint = kind(1), Kreal = kind(1.0)
  real(Kreal2), parameter :: PId = 3.14159265358979323846d0
  real(Kreal), parameter:: PI = real(PId,Kreal)
  real(Kreal), parameter:: Epsil10 = &
      float(RADIX(1.0_Kreal))**(1-DIGITS(1.0_Kreal)) * 10._Kreal

  ! These should be set by the user before calling Sjhminit()
  integer(Kint):: Nm=64, Lup=0, Ldn=0, Ng=64, Nh=1, Numu, Nphi, Nq=2000
  ! Nm: number of azimuth expansion terms (m=0,1,...,Nm-1)
  ! Lup/Ldn: signs of zenith expansion terms for upward/downward directions
  !   if =0: l = m,m+1,..., Nm-1, i.e., Lmup/Lmdn=Nm-m, triangular
  !   if >0: l = m,m+1,..., Nm-1+Lup/Ldn, i.e,, trapezoidal truncation
  !   if <0: l = m,m+1,..., m-Lup/Ldn-1, i.e., parallelogram truncation
  ! Ng: number of phase function expansion terms (n=0,1,...,Ng-1)
  ! Nh: number of BDRF expansion terms of surface (n=0,1,...,Nh-1)
  ! Numu: number of zenith points for radiance output
  ! Nphi: number of azimuth points for radiance output

  integer(Kint):: Mm, Mg, Mh, Mlg, Ml, Lm2
  real(Kreal), allocatable:: Umu(:), Phi(:)
  integer(Kint), allocatable:: Lmup(:), Lmdn(:)
  real(Kreal2), allocatable :: AA(:,:), ATint(:), UUinv(:,:,:), PP(:,:,:), &
  & PPf(:,:,:)
  real(Kreal), allocatable:: Itild(:,:,:)

  real(Kreal2), allocatable :: gl_nodes(:), gl_weights(:)
  integer(Kint), private :: current_n = 0

  contains

  !-----------------------------------------------------------------------------
  ! 工业级最稳健高斯求积求解器：初值预测 + 本地高阶牛顿迭代修正
  ! 彻底抛弃复杂的特征值矩阵消去，在 N=3000 时提供全硬件极限精度（~10^-16）
  !-----------------------------------------------------------------------------
  subroutine Glnewton(n)
    integer, intent(in) :: n
    real(Kreal2) :: t, p1, p2, p3, dp1, delta
    integer :: i, j, l

    current_n = n
    if (allocated(gl_nodes)) deallocate(gl_nodes)
    if (allocated(gl_weights)) deallocate(gl_weights)
    allocate(gl_nodes(n), gl_weights(n))

    ! 利用对称性，只需要计算前一半的节点
    do i = 1, (n + 1) / 2
      ! 1. 经典切比雪夫初值预测公式
      t = cos(PId * (real(i, Kreal2) - 0.25d0) / (real(n, Kreal2) + 0.5d0))
            
      ! 2. 進行高階牛頓-拉弗森迭代精確修正（15次迭代確保徹底鎖死最後一位浮點數誤差）
      do j = 1, 15
        p1 = 1.0d0
        p2 = 0.0d0
        ! 标准勒让德多项式的 Clenshaw 三项递推求值
        do l = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2.0d0 * real(l, Kreal2) - 1.0d0) * t * p2 - &
          & (real(l, Kreal2) - 1.0d0) * p3) / real(l, Kreal2)
        end do
        ! 精确导数计算
        dp1 = real(n, Kreal2) * (t * p1 - p2) / (t**2 - 1.0d0)
                
        delta = p1 / dp1
        t = t - delta
        if (abs(delta) < 1.0d-16) exit
      end do
            
      ! 3. 利用镜像对称性完美赋值（保证完全正交与对称）
      gl_nodes(i) = -t
      gl_nodes(n + 1 - i) = t
            
      ! 4. 根据导数解析计算严格的高斯-勒让德权重
      gl_weights(i) = 2.0d0 / ((1.0d0 - t**2) * dp1**2)
      gl_weights(n + 1 - i) = gl_weights(i)
    end do

  end subroutine Glnewton

  !-----------------------------------------------------------------------------
  ! 函数名称: Sjf
  ! 目的: 计算兼顾极点与赤道收敛性的移位雅可比函数 (Shifted Jacobi Function)
  ! 输入: l, m   - 基函数阶数 (l >= m, m >= 0)
  !       mu     - 天顶角余弦 cos(theta), 范围 [0.0, 1.0]
  ! 返回: 该点上的基函数数值 (双精度浮点数)
  !---------------------------------------------------------------------------
  FUNCTION Sjf(l, m, mu)
    INTEGER, INTENT(IN) :: l, m
    REAL(Kreal2), INTENT(IN) :: mu(:)
    real(Kreal2) :: Sjf(size(mu))
    
    INTEGER :: n, i
    REAL(Kreal2) :: q_n(size(mu)), q_prev1(size(mu)), q_prev2(size(mu))
    REAL(Kreal2) :: den, A, B, C
    REAL(Kreal2) :: norm
    REAL(Kreal2) :: ri, rm

    integer(Kint) :: Consumed_m, Remaining_m

    Consumed_m = 0
        
    ! 边界条件检查
    IF (l < m .OR. m < 0 .OR. any(mu < 0.0_Kreal2) .OR. any(mu > 1.0_Kreal2)) THEN
      print*, " [Sjf] 错误: 输入参数不合法 (l < m 或 m < 0 或 mu 不在 [0,1] 范围内)"
      STOP
    END IF
      
    n = l - m ! 多项式的阶数
      
    ! 1. 计算三项递推多项式部分 q_n(mu)
    q_prev2 = 1.0_Kreal2 ! q_0
      
    IF (n == 0) THEN
      q_n = q_prev2
    ELSE IF (n == 1) THEN
      q_n = REAL(2 * m + 2, Kreal2) * mu - 1.0_Kreal2 ! q_1
    ELSE
      q_prev1 = REAL(2 * m + 2, Kreal2) * mu - 1.0_Kreal2
      DO i = 1, n - 1
        ! 定义双精度变量以防整数截断
        ri = REAL(i, Kreal2)
        rm = REAL(m, Kreal2)
    
        ! 计算严格化简后的 A_n, B_n, C_n 递推系数
        den = (ri + 1.0_Kreal2) * (ri + 2.0_Kreal2 * rm + 1.0_Kreal2)
    
        !A = 2.0_Kreal2 * (ri + rm + 1.0_Kreal2) * (2.0_Kreal2 * ri + 2.0_Kreal2 * rm + 1.0_Kreal2) / den
        A = 2.0_Kreal2 * (ri + rm + 1.0_Kreal2) * (2.0_Kreal2 * &
        & (ri + rm) + 1.0_Kreal2) / den
    
        ! 严格修正版 B_n
        B = -(2.0_Kreal2 * (ri + rm) + 1.0_Kreal2) * &
        & (ri**2 + 2.0_Kreal2 * rm * ri + ri + rm) / (den * (ri + rm))
    
        C = ri * (ri + 2.0_Kreal2 * rm) * (ri + rm + 1.0_Kreal2) / &
          (den * (ri + rm))
    
        ! 三项递推内核公式
        q_n = (A * mu + B) * q_prev1 - C * q_prev2

        !IF (any(ABS(q_n) > 1.0D+100) .AND. consumed_m < m) THEN
        IF (consumed_m < m) THEN
                ! 强行把一个 (1-mu) 的 1 次方乘进去，实行降维打击
                q_n     = q_n     * (1.0_Kreal2 - mu)
                q_prev1 = q_prev1 * (1.0_Kreal2 - mu)
                q_prev2 = q_prev2 * (1.0_Kreal2 - mu)
                
                ! 记录消耗：我们在内部已经帮外层消灭了 1 次幂
                consumed_m = consumed_m + 1
        END IF
    
        ! 滚动更新状态
        q_prev2 = q_prev1
        q_prev1 = q_n
      END DO
    END IF
        
    ! 2. 计算极点约束因子 (1 - mu)^m
    !pole_factor = (1.0_Kreal2 - mu)**m

    remaining_m = m - consumed_m
        
    ! 3. 计算方位角对称分量 cos(m * phi)
    !phi_factor = COS(REAL(m, Kreal2) * phi)
        
    ! 4. 计算归一化系数 N_n^m
    !IF (m == 0) THEN
    !  norm = SQRT(REAL(2 * n + 1, Kreal2) / (2.0_Kreal2 * PId))
    !ELSE
    !  norm = SQRT(REAL(2 * n + 2 * m + 1, Kreal2) / PId)
    !END IF
    norm = SQRT(REAL(2 * n + 2 * m + 1, Kreal2) * 2.0_Kreal2)

    ! 5. 组合最终的正交基函数值
    !Sjf = norm * pole_factor * q_n ! * phi_factor
    Sjf = norm * q_n * ((1.0_Kreal2 - mu)**remaining_m)
        
  END FUNCTION Sjf

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

  !-----------------------------------------------------------------------------
  subroutine Invert_sym_tridiag(n, d, e, Ainv)
    integer, intent(in)  :: n
    real(Kreal2), intent(in) :: d(:)
    real(Kreal2), intent(in) :: e(:)
    real(Kreal2), intent(out):: Ainv(:,:)

    real(Kreal2) :: f(n), g(n)  ! 前向和后向比例数组
    real(Kreal2) :: diag_inv(n) ! 逆矩阵对角线元素缓冲
    integer  :: i, j

    ! 1. 前向递推比例 f(i) = theta(i) / theta(i-1)
    f(1) = d(1)
    do i = 2, n
        if (abs(f(i-1)) < 1.0d-300) then ! 双精度下溢保底值降至 10^-300
            print *, "Error: Forward breakdown, matrix is likely singular."
            stop
        end if
        f(i) = d(i) - (e(i-1)**2) / f(i-1)
    end do

    ! 2. 后向递推比例 g(i) = phi(i) / phi(i+1)
    g(n) = d(n)
    do i = n-1, 1, -1
        if (abs(g(i+1)) < 1.0d-300) then
            print *, "Error: Backward breakdown."
            stop
        end if
        g(i) = d(i) - (e(i)**2) / g(i+1)
    end do

    ! 3. 检查最终收敛元素（【关键改动】阈值提高到双精度极限 1.0d-15）
    if (abs(f(n)) < 1.0d-15) then
        print *, "Error: Matrix is ill-conditioned or singular!"
        stop
    end if

    ! 4. 计算逆矩阵的对角线元素 Ainv(i,i)
    if (n == 1) then
        diag_inv(1) = 1.0_Kreal2 / d(1)
    else
        diag_inv(1) = 1.0_Kreal2 / g(1)
        do i = 2, n-1
            diag_inv(i) = 1.0_Kreal2 / (d(i) - (e(i-1)**2)/f(i-1) - (e(i)**2)/g(i+1))
        end do
        diag_inv(n) = 1.0_Kreal2 / f(n)
    end if

    ! 5. 填充满矩阵，并利用对称性复制
    do i = 1, n
        Ainv(i,i) = diag_inv(i)
        do j = i + 1, n
            Ainv(i,j) = -Ainv(i,j-1) * (e(j-1) / g(j))
            Ainv(j,i) = Ainv(i,j) ! 对称矩阵特性：转置相等
        end do
    end do

  end subroutine Invert_sym_tridiag

  !------------------------------------------------------------------------------
  subroutine Getuuinvm(m,UUinvm)
    integer(Kint), intent(in):: m
    real(Kreal2), intent(out):: UUinvm(:,:)

    real(Kreal2):: UU1(max(Lmup(m), Lmdn(m)),2)
    integer(Kint):: i, n, Mlm, Lm2m

    Mlm = max(Lmup(m), Lmdn(m))
    Lm2m = Lmup(m) + Lmdn(m)
    UUinvm = 0.0_Kreal2
 
    do i = 1, Mlm
      if (i==1) then
        UU1(i,1) = 0.5_Kreal2 / (real(m + 1, kind=Kreal2))
      else
        n = i - 1
        UU1(i,1) = real(n*n + 2*m*n + n + m, kind=Kreal2) / &
        & (real(n + m, kind=Kreal2) * real(n + m + 1, kind=Kreal2) * 2.0_Kreal2)
      endif

      if (Mlm>i) then
        UU1(i,2) = real(i*(i + 2*m), kind=Kreal2) / &
        & (2.0_Kreal2 * real(i+m, kind=Kreal2) * &
        & real(2*i+2*m+1, kind=Kreal2)) * &
        & sqrt(real(2*i+2*m+1,kind=Kreal2) / real(2*i+2*m-1, kind=Kreal2))
      endif
    enddo

    call Invert_sym_tridiag(Lmup(m), UU1(1:Lmup(m),1), UU1(1:Lmup(m)-1,2), &
    & UUinvm(1:Lmup(m),1:Lmup(m)))
    if (Lmdn(m) == Lmup(m)) then
      UUinvm(Lmup(m)+1:Lm2m,Lmup(m)+1:Lm2m) = - UUinvm(1:Lmdn(m),1:Lmdn(m))
    else
      call Invert_sym_tridiag(Lmdn(m), UU1(1:Lmdn(m),1), UU1(1:Lmdn(m)-1,2), &
      & UUinvm(Lmup(m)+1:Lm2m,Lmup(m)+1:Lm2m))
      UUinvm(Lmup(m)+1:Lm2m,Lmup(m)+1:Lm2m) = &
      & - UUinvm(Lmup(m)+1:Lm2m,Lmup(m)+1:Lm2m)
    endif

  end subroutine Getuuinvm

  !-----------------------------------------------------------------------------
  subroutine Getppm(m, PPm)
    integer(Kint), intent(in):: m
    real(Kreal2), intent(out):: PPm(:,:)

    real(Kreal2):: Prod(Nq), Phat1(Nq, 0:Mg), Psft(Nq), Sinp(Nq), Cosp(Nq), Sjf1(Nq)
    integer(Kint):: Mlm, i, j

    if (mod(m,2)==0) then

      Psft = 0.5d0 * (gl_nodes + 1.0d0)

      if (Mg>=m) then
        Mlm = max(Lmup(m), Lmdn(m))
        Phat1(1:Nq,m:Mg) = Phat(Mg, m, Psft)
        do j = 1, Mlm
          Sjf1(1:Nq) = Sjf(j-1+m, m, Psft)
          do i = 1, Ng-m
            Prod(1:Nq) = Phat1(1:Nq,m+i-1) * Sjf1(1:Nq)
            PPm(i,j) = sum(gl_weights * 0.5d0 * Prod) * 0.5d0
          enddo
        enddo
      endif

    else

      Psft = 0.25d0 * PId * (gl_nodes + 1.0d0)
      Sinp = sin(Psft)
      Cosp = cos(Psft)

      if (Mg>=m) then
        Mlm = max(Lmup(m), Lmdn(m))
        Phat1(1:Nq,m:Mg) = Phat(Mg, m, Sinp)
        do j = 1, Mlm
          Sjf1(1:Nq) = Sjf(j-1+m, m, Sinp) / Cosp
          do i = 1, Ng-m
            Prod(1:Nq) = Phat1(1:Nq,m+i-1) * Sjf1(1:Nq)
            PPm(i,j) = sum(gl_weights * Prod * (Cosp**2)) * (PId / 8.0d0)
          enddo
        enddo
      endif

    endif

  end subroutine Getppm

  !-----------------------------------------------------------------------------
  subroutine Getppfm(m, PPfm)
    integer(Kint), intent(in):: m
    real(Kreal2), intent(out):: PPfm(:,:)

    real(Kreal2):: Prod(Nq), Phat1(Nq, 0:Mh), Psft(Nq), Sinp(Nq), Cosp(Nq), Sjf1(Nq)
    integer(Kint):: i, j

    if (mod(m,2)==0) then
      Psft = 0.5d0 * (gl_nodes + 1.0d0)
      if (Mh>=m) then
        Phat1(1:Nq,m:Mh) = Phat(Mh,m,Psft)
        do j = 1, Lmdn(m)
          Sjf1 = Sjf(j-1+m,m,Psft) * Psft
          do i = 1, Nh-m
            Prod(1:Nq) = Phat1(1:Nq,m+i-1) * Sjf1 * (-1)**(i-1)
            PPfm(i,j) = sum(gl_weights*Prod) * 0.5_Kreal2
          enddo
        enddo
      else
        PPfm = 0.0_Kreal2
      endif
    
    else
      
      Psft = 0.25d0 * PId * (gl_nodes + 1.0d0)
      Sinp = sin(Psft)
      Cosp = cos(Psft)
      if (Mh>=m) then
        Phat1(1:Nq,m:Mh) = Phat(Mh, m, Sinp)
        do j = 1, Lmdn(m)
          Sjf1 = Sjf(j-1+m, m, Sinp) * Sinp / Cosp
          do i = 1, Nh-m
            Prod(1:Nq) = Phat1(1:Nq,m+i-1) * Sjf1 * (-1)**(i-1)
            PPfm(i,j) = sum(gl_weights*Prod*(Cosp**2)) * (PId * 0.25_Kreal2)
          enddo 
        enddo
      else
        PPfm = 0.0_Kreal2
      endif 
    endif

  end subroutine Getppfm

  !-----------------------------------------------------------------------------
  subroutine Sjhminit(Memoryless)

    logical, intent(in), optional:: Memoryless

    integer(Kint):: m

    Mm = Nm - 1
    Mg = Ng - 1
    Mh = Nh - 1

    ! Number of zenith expansion terms for each azimuthal order m
    allocate(Lmup(0:Mm), Lmdn(0:Mm))
    if (Lup >=0) then
      Lmup = Nm - (/(m,m=0,Mm)/) + Lup
    else
      Lmup = - Lup
    end if
    if (Ldn >=0) then
      Lmdn = Nm - (/(m,m=0,Mm)/) + Ldn
    else
      Lmdn = - Ldn
    end if

    Mlg = max(Lmup(0)-1, Lmdn(0)-1, Mg, Mh)
    Lm2 = Lmup(0) + Lmdn(0)

    ! AA: Identity matrix with enough size for any m
    allocate(AA(Lm2, Lm2))
    AA = Identity(Lm2, 1.0_Kreal2)

    ! ATint for m=0
    allocate(ATint(Lmup(0)+Lmdn(0)))
    ATint = 0.0_Kreal2
    ATint(1) = 1.0_Kreal2 / sqrt(2.0_Kreal2)
    ATint(Lmup(0)+1) = ATint(1)
    
    call Glnewton(Nq)

    if (.not.Memoryless) then

      ! UUinv: Inverse of the UU matrix for each m

      allocate(UUinv(Lm2, Lm2, 0:Mm))
      UUinv = 0.0_Kreal2

      do m = 0, Mm
      
        Lm2 = Lmup(m) + Lmdn(m)
        call Getuuinvm(m, UUinv(1:Lm2,1:Lm2,m))

      enddo

      ! PP

      allocate(PP(Ng, max(Lmup(0),Lmdn(0)), 0:Mm))

      do m = 0, Mm
        Ml = max(Lmup(m), Lmdn(m))
        call Getppm(m, PP(1:Ng-m, 1:Ml, m))
      enddo

      ! PPf

      allocate (PPf(Nh,Lmdn(0),0:Mm))

      !if (Nh>Ng) then ! assuming Ng is always bigger than Nh
      !  print*, 'Nh > Ng in Sjhminit, check the input or modify the code!'
      !  stop
      !endif

      do m = 0, Mm
        call Getppfm(m, PPf(1:Nh-m, 1:Lmdn(m), m))
      enddo

    endif
  
    if (.not.allocated(Umu)) allocate(Umu(Numu), Phi(Nphi))

  end subroutine Sjhminit

  !-----------------------------------------------------------------------------
  subroutine Sjhmfin()
    if (allocated(Lmup)) deallocate(Lmup)
    if (allocated(Lmdn)) deallocate(Lmdn)
    if (allocated(AA)) deallocate(AA)
    if (allocated(UUinv)) deallocate(UUinv)
    if (allocated(ATint)) deallocate(ATint)
    if (allocated(PP)) deallocate(PP)
    if (allocated(PPf)) deallocate(PPf)
    deallocate(gl_nodes, gl_weights)
    deallocate(Umu, Phi)
    deallocate(Itild)
  end subroutine Sjhmfin

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

  !-----------------------------------------------------------------------------
  ! Intensities in directions (Umu,Phi), calculated from Itild
  subroutine Iangles(Nlyr, Rad)
    ! use Sjf
    integer(Kint), intent(in):: Nlyr
    real(Kreal), intent(out):: Rad(Numu,Nphi,0:Nlyr)

    real(Kreal2):: Umu0(1), Im(0:Mm,0:Nlyr), Cosphi(Nphi,0:Mm), Sjf0(1)
    integer(Kint):: i, l, m

    do m = 0, Mm
      Cosphi(1:Nphi,m) = cos(Phi(1:Nphi)*m)
    enddo

    do i = 1, Numu
      do m = 0, Mm
        Im(m,:) = 0.0_Kreal
        if (Umu(i).lt.0.0_Kreal) then
          Umu0(1) = - Umu(i)
          do l = m, m+Lmdn(m)-1
            Sjf0 = Sjf(l,m,Umu0)
            Im(m,:) = Im(m,:) + Sjf0(1)*Itild(l-m+1+Lmup(m),:,m)
          enddo
        else
          Umu0(1) = Umu(i)
          do l = m, m+Lmup(m)-1
            Sjf0 = sjf(l,m,Umu0)
            Im(m,:) = Im(m,:) + Sjf0(1)*Itild(l-m+1,:,m)
          enddo
        endif
      enddo
      Rad(i,:,:) = real(matmul(Cosphi(:,:),Im(:,:)), Kreal)
    enddo

  endsubroutine Iangles

  !-----------------------------------------------------------------------------
  subroutine Merglevs(Nlyr0, Od0, Ssa0, Pf0, T0, Nutau, Utau, Nlyr2, Od, Mlev, &
    & Ssa, Pf, Thermal, Wnumlo, Wnumhi, B0)
    ! inputs:
    !   Nlyr0: 模式均匀介质层数
    !   Nutau: 用户指定的光学深度的个数
    !   Od0(Nlyr0): 均匀介质层的光学厚度
    !   Utau(Nutau): 用户指定的光学深度
    !   Ssa0(Nlyr0): 模式均匀层单次散射反照率
    !   Pf0(0:Mg,Nlyr0): 模式均匀层相函数
    !   T0(0:Nlyr0): 模式各层级温度
    !   Thermal: 是否计算热辐射
    !   Wnumlo, Wnumhi: 计算热辐射的波段两端波数

    ! outputs:
    !   Nlyr2: 模式分层与用户指定光学厚度融合后的总层数，Nlyr0<=Nlyr2<=Nlyr0+Nutau
    !   Od(Nlyr2): 各融合分层的光学厚度
    !   Mlev(Nutau): 用户指定光学深度对应的融合高度，0<=Mlev<=Nlyr2
    !   Ssa(Nlyr2): 融合层单次散射反照率
    !   Pf(0:Mg,Nlyr2)：融合层相函数
    !   B0(0:Nlyr2): 融合后各层级黑体辐射

    ! internal:
    !   Tau(0:Nlyr0): 模式层的累计光学厚度
    !   Tall(0:Nlyr2): 融合层的累计光学厚度
    !   Nuc: 用户指定的光学厚度的计数
  
    integer(Kint), intent(in):: Nlyr0, Nutau
    real(Kreal), intent(in):: Od0(Nlyr0), Utau(Nutau), Ssa0(Nlyr0), &
    & Pf0(0:Mg,Nlyr0), T0(0:Nlyr0), Wnumlo, Wnumhi
    integer(Kint), intent(out):: Nlyr2, Mlev(Nutau)
    logical, intent(in):: Thermal
    real(Kreal), intent(out):: Od(Nlyr0+Nutau), Ssa(Nlyr0+Nutau), &
    & Pf(0:Mg,Nlyr0+Nutau), B0(0:Nlyr0+Nutau)

    real(Kreal):: Tau(0:Nlyr0), Tall(0:Nlyr0+Nutau), B0m(0:Nlyr0)
    integer(Kint):: Nuc, Nuc0, l, i

    Tau(0) = 0.0_Kreal
    Tall(0) = 0.0_Kreal
    Nuc = 0
    Nlyr2 = 0
    if (Thermal) then
      do l = 0, Nlyr0
        B0m(l) = Planckint(Wnumlo, Wnumhi, T0(l))
      enddo
    endif

    B0(0) = B0m(0)
    if (Utau(1) < Epsil10) then
      Nuc = 1
      Mlev(1) = 0
    endif
 
    do l = 1, Nlyr0

      Tau(l) = Tau(l-1) + Od0(l)
      Nuc0 = Nuc
      do i = Nuc0+1, Nutau
        if (Utau(i).ge.(Tau(l-1)+Epsil10) .and. Utau(i).le.(Tau(l)-Epsil10)) then
          Nuc = Nuc + 1
          Nlyr2 = Nlyr2 + 1
          Mlev(Nuc) = Nlyr2
          Tall(Nlyr2) = Utau(i)
          Od(Nlyr2) = Tall(Nlyr2) - Tall(Nlyr2-1) 

          Ssa(Nlyr2) = Ssa0(l)
          Pf(0:,Nlyr2) = Pf0(0:,l)
          !T(Nlyr2) = T0(l-1) + (T0(l)-T0(l-1)) * Od(Nlyr2) / Od0(l)
          B0(Nlyr2) = B0m(l-1) + (B0m(l)-B0m(l-1)) * Od(Nlyr2) / Od0(l)
        elseif (abs(Tau(l)-Utau(i)) < Epsil10) then
          Nuc = Nuc + 1
          Nlyr2 = Nlyr2 + 1
          Mlev(Nuc) = Nlyr2
          Tall(Nlyr2) = Tau(l)
          Od(Nlyr2) = Tall(Nlyr2) - Tall(Nlyr2-1)
          Ssa(Nlyr2) = Ssa0(l)
          Pf(0:,Nlyr2) = Pf0(0:,l)
          !T(Nlyr2) = T0(l)
          B0(Nlyr2) = B0m(l)
          exit
        else
          Nlyr2 = Nlyr2 + 1
          Tall(Nlyr2) = Tau(l)
          Od(Nlyr2) = Tall(Nlyr2) - Tall(Nlyr2-1)
          Ssa(Nlyr2) = Ssa0(l)
          Pf(0:,Nlyr2) = Pf0(0:,l)
          !T(Nlyr2) = T0(l)
          B0(Nlyr2) = B0m(l)
          exit
        endif
      enddo
 
    enddo

  endsubroutine Merglevs

  !-----------------------------------------------------------------------------

  subroutine Sjhm(Nlay, Od0, Ssa0, Pf0, I0, Mu0, T0, Tsurf, GS, &
                & Thermal, Deltam, Wnumlo, Wnumhi, &
                & Beam0, Fisot, Ttemp, Temis, Usertau, Nutau, Utau, &
                & Rad, Fup, Fdn, Act, Dfdt, Crt)

    integer(Kint), intent(in):: Nlay
    real(Kreal), intent(in):: Od0(:), Ssa0(:), Pf0(0:,:), I0, Mu0, &
    & T0(0:), Gs(0:), Tsurf, Wnumlo, Wnumhi, Crt, Fisot, Ttemp, Temis
    logical, intent(in):: Deltam, Thermal

    logical, intent(in):: Usertau
    integer, intent(in):: Nutau
    real(Kreal),intent(in):: Utau(Nutau)

    real(Kreal), intent(out):: Beam0(0:), Rad(:,:,0:), &
    & Fup(0:), Fdn(0:), Act(0:), Dfdt(0:)

    integer(Kint),save:: Nlyr0, Ml0, Mm0
    
    integer(Kint):: Ngd, Nup, Ndn, Nstr, Ngm, Nhm, Nghm

    real(Kreal):: Ssa(Nlay+Nutau), Odd(Nlay+Nutau), Ssad(Nlay+Nutau), &
    & Pf(0:Mg,Nlay+Nutau), Pfd(0:Mg,Nlay+Nutau), T(0:Nlay+Nutau), &
    & Beam(0:Nlay+Nutau), Dump1(Nlay+Nutau), Idntop(Lmdn(0)), B0(0:Nlay+Nutau),&
    & Bsurf, Hm(Nh,Nh), Od(Nlay+Nutau)
    real(Kreal2):: Glm(0:Mg,Nlay+Nutau), PPm(max(Ng,Nh),Lmup(0)+Lmdn(0)), Mu0d, &
    & Pmu0(max(Ng,Nh)), Phat0(1,0:max(Ng,Nh)), PPfm(Nh,Lmdn(0)), &
    & PP2(max(Ng,Nh),max(Lmup(0),Lmdn(0))), &
    & U2(Lmup(0)+Lmdn(0),Lmup(0)+Lmdn(0))

    real(Kreal):: I0up(0:Nlay+Nutau), I0dn(0:Nlay+Nutau), Rup(Lmup(0)+Lmdn(0),0:Nlay+Nutau), &
    & Rdn(Lmup(0)+Lmdn(0),0:Nlay+Nutau), Ratiom

    ! if Usertau
    integer(Kint):: Nlyr, Mlev(Nutau)

    logical:: Thermalm

    integer(Kint):: l, m, i, Time(8), Ncrt

    !call cpu_time(Cpu)
    !call date_and_time(values=Time)
    !write(*,'(A,I2," - ",I2,":",I2,":",I2)') 'Running SJHM module at ', &
    !& Time(3), Time(5), Time(6), Time(7) !Cpu

    if (Usertau) then
      !call Usrlevs(Nlyr, Od0, Nutau, Utau, Levf)
      call Merglevs(Nlay, Od0, Ssa0, Pf0, T0, Nutau, Utau, Nlyr, Od, Mlev, Ssa, Pf, Thermal, Wnumlo, Wnumhi, B0)
    else
      Nlyr = Nlay
      Od(1:Nlyr) = Od0(1:Nlay)
      Ssa(1:Nlyr) = Ssa0(1:Nlay)
      Pf(0:,1:Nlyr) = Pf0(0:,1:Nlay)
      T(0:Nlyr) = T0(0:Nlay)
    endif

    ! allocate Itild

    Ml = Lmup(0)+Lmdn(0)
    if (.not.allocated(Itild)) then
      allocate(Itild(0:Ml,0:Nlyr,0:Mm))
    elseif (Nlyr.ne.Nlyr0 .or. Ml.ne.Ml0 .or. Mm.ne.Mm0) then
      deallocate(Itild)
      allocate(Itild(0:Ml,0:Nlyr,0:Mm))
    endif
    Nlyr0 = Nlyr
    Ml0 = Ml
    Mm0 = Mm

    Itild = 0.0_Kreal

    ! Checking Od and Ssa

    if (any(Od(1:Nlyr)<Epsil10)) then
      print*, 'Warning in Hshmrun: too small Od:', Od(1:Nlyr)
      !stop
    endif
    where (Od < Epsil10)
      Od = Epsil10
    endwhere

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
      Dump1(1:Nlyr) = 1.0_Kreal - Ssa(1:Nlyr) * Pf(Ngd,1:Nlyr)
      do l = 1, Nlyr
        Odd(l) = Od(l) * Dump1(l)
      enddo
      Ssad(1:Nlyr) = Ssa(1:Nlyr) * (1.0_Kreal - Pf(Ngd,1:Nlyr)) / Dump1(1:Nlyr)
      do l = 1, Nlyr
        Pfd(0:Ngd-1,l) = (Pf(0:Ngd-1,l) - Pf(Ngd,l)) / (1.0_Kreal - Pf(Ngd,l))
      enddo
    else
      Ngd = Ng
      Odd(1:Nlyr) = Od(1:Nlyr)
      Ssad(1:Nlyr) = Ssa(1:Nlyr)
      Pfd(0:Mg,1:Nlyr) = Pf(0:Mg,1:Nlyr)
    endif

    ! Solar beam. In principle, this should be in Gdrtes, 
    ! here is for consistency and efficiency.
    Beam(0) = I0
    do l = 1, Nlyr
      Beam(l) = Beam(l-1) * exp(-Odd(l) / Mu0)
    enddo

    !call date_and_time(values=Time)
    !write(*,'(A,I6," at ",I2," - ",I2,":",I2,":",I2)') 'Begin solving m = ', &
    !  & m, Time(3), Time(5), Time(6), Time(7)

    Ncrt = 0

    do m = 0, Mm

      Nup = Lmup(m)
      Ndn = Lmdn(m)
      Nstr = Nup + Ndn

      ! phase function for moment m -- scaled
      Glm(0:Ngd-1,1:Nlyr) = 0.0_Kreal2
      Glm(m:Ngd-1,1:Nlyr) = Pfd(m:Ngd-1,1:Nlyr) * 0.5_Kreal2

      Ngm = Ngd-m
      Nhm = Nh - m
      Nghm = max(Ngm, Nhm)

      if (allocated(PP)) then

        PPm(1:Nghm,1:Nup) = PP(1:Nghm,1:Nup,m)
        PPm(1:Nghm:2,Nup+1:Nstr) = PP(1:Nghm:2,1:Ndn,m)
        PPm(2:Nghm:2,Nup+1:Nstr) = -PP(2:Nghm:2,1:Ndn,m)

      else

        call Getppm(m, PP2(1:Nghm,1:max(Nup,Ndn)))
        PPm(1:Nghm,1:Nup) = PP2(1:Nghm,1:Nup)
        PPm(1:Nghm:2,Nup+1:Nstr) = PP2(1:Nghm:2,1:Ndn)
        PPm(2:Nghm:2,Nup+1:Nstr) = -PP2(2:Nghm:2,1:Ndn)

      endif

      ! solar beam component m
      Mu0d = Mu0
      Phat0(:,m:Nghm+m-1) = Phat(Nghm+m-1,m,(/-Mu0d/))
      Pmu0(1:Nghm) = Phat0(1,m:Nghm+m-1)
      if (m .ne. 0) then
        Pmu0(1:Nghm) = Pmu0(1:Nghm) * 2
      endif

      ! TOA diffuse source
      Idntop = 0.0_Kreal
      if (m.eq.0) then
        Idntop(1) = Fisot * sqrt(0.5_Kreal)
      endif

      ! Thermal source for m = 0
      if (m.eq.0) then
        Thermalm = Thermal
        if (Thermal) then
          Idntop(1) = Idntop(1) + Planckint(Wnumlo, Wnumhi, Ttemp) * Temis * sqrt(0.5_Kreal)
          if (.not.Usertau) then
            do l = 0, Nlyr
              B0(l) = Planckint(Wnumlo, Wnumhi, T(l))
            enddo
          endif
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

        if (allocated(PPf)) then
          PPfm(1:Nhm,1:Ndn) = PPf(1:Nhm,1:Ndn,m)
        else
          call Getppfm(m, PPfm(1:Nhm,1:Ndn))
        endif

        Hm(1:Nhm,1:Nhm) = 0.0_Kreal
        do i = 1, Nhm
          Hm(i,i) = Gs(i-1) * 0.5_Kreal
        enddo

      endif

      if (allocated(UUinv)) then
        U2(1:Nstr,1:Nstr) = UUinv(1:Nstr,1:Nstr,m)
      else
        call Getuuinvm(m, U2(1:Nstr,1:Nstr))
      endif

      ! Solve m-th GD-RTE

      call Gdrtes(Nlyr, Nup, Ndn, Ngm, Nhm, AA(1:Nstr,1:Nstr), &
      & AA(1:Nstr,1:Nstr), U2(1:Nstr,1:Nstr), &
      & Odd, Ssad, Glm(m:Ngd-1,1:Nlyr), PPm(1:Nghm, 1:Nstr), &
      & Beam(0:Nlyr), Mu0, Idntop(1:Ndn), Pmu0(1:Nghm), Thermalm, ATint(1:Nstr), &
      & B0(0:Nlyr), PPfm(1:Nhm,1:Ndn), Hm(1:Nhm,1:Nhm), Bsurf, &
      & Itild(1:Nstr,0:Nlyr,m))

      ! Truncation of m, when Crt /= 0

      if (m==0) then
        I0up(0:Nlyr) = abs(Itild(1,0:Nlyr,m)) * sqrt(2.0_Kreal)
        I0dn(0:Nlyr) = abs(Itild(Lmup(m)+1,0:Nlyr,m)) * sqrt(2.0_Kreal)
        Ratiom = 1.0_Kreal
      else
        do l = 0, Nlyr
          do i = 1, Lmup(m)
            Rup(i,l) = abs(Itild(i,l,m)) * sqrt(2.0_Kreal*(2*(i+m)-1)) / I0up(l)
          enddo
          do i = 1, Lmdn(m)
            Rdn(i,l) = abs(Itild(Lmup(m)+i,l,m)) * sqrt(2.0_Kreal*(2*(i+m)-1))/ I0dn(l)
          enddo
        enddo
        Ratiom = max(maxval(Rup(1:Lmup(m),0:Nlyr)),maxval(Rdn(1:Lmdn(m),0:Nlyr)))
      endif

      call date_and_time(values=Time)
      !write(*,'(A,I4," at ",I2," - ",I02,":",I02,":",I02, " Ratiom:", F10.7)') ' After solving m = ', &
      !& m, Time(3), Time(5), Time(6), Time(7), Ratiom

      if (Ratiom .lt. Crt) then
        Ncrt = Ncrt + 1
        if (Ncrt==2) exit
      else
        Ncrt = 0
      endif

    enddo ! m = 0, Mm

    ! Fluxes

    Act(0:Nlyr) = 2.0_Kreal * sqrt(2.0_Kreal) * PI &
    & * (Itild(1,0:Nlyr,0) + Itild(Lmup(0)+1,0:Nlyr,0)) + Beam(0:Nlyr)
    
    Dfdt(0:Nlyr) = Act(0:Nlyr) - B0(0:Nlyr) * 4.0_Kreal * PI
    if (Ssa(1)>=(1.0_Kreal-Epsil10)) then
      Dfdt(0) = 0.0_Kreal
    else
      Dfdt(0) = Dfdt(0) * (1.0_Kreal - Ssad(1))
    endif
    Dfdt(1:Nlyr) = Dfdt(1:Nlyr) * (1.0_Kreal - Ssad(1:Nlyr))
    where (Ssa(1:Nlyr)>=1.0_Kreal-Epsil10)
      Dfdt(1:Nlyr) = 0.0_Kreal
    endwhere

    if (Lmup(0).gt.1) then
      Fup(0:Nlyr) = sqrt(2.0_Kreal) * PI * (Itild(1,0:Nlyr,0) &
      & + Itild(2,0:Nlyr,0)/sqrt(3.0_Kreal))
    else
      Fup(0:Nlyr) = sqrt(2.0_Kreal) * PI * Itild(1,0:Nlyr,0)
    endif

    if (Lmdn(0).gt.1) then
      Fdn(0:Nlyr) = sqrt(2.0_Kreal) * PI * &
      & (Itild(Lmup(0)+1,0:Nlyr,0) + Itild(Lmup(0)+2,0:Nlyr,0)/sqrt(3.0_Kreal))
    else
      Fdn(0:Nlyr) = sqrt(2.0_Kreal) * PI * &
      & Itild(Lmup(0)+1,0:Nlyr,0)
    endif

    call Iangles(Nlyr, Rad)

    if (Usertau) then
      do i = 1, Nutau
        Act(i-1) = Act(Mlev(i))
        Dfdt(i-1) = Dfdt(Mlev(i))
        Fup(i-1) = Fup(Mlev(i))
        Fdn(i-1) = Fdn(Mlev(i))
        Rad(:,:,i-1) = Rad(:,:,Mlev(i))
        Beam0(i-1) = Beam(Mlev(i))
      enddo
    endif

  endsubroutine Sjhm

end module Sjhm_module

