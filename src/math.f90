module math
  use global_variables
  implicit none

  public :: Bvgd, Levmar, Diag, Svd
  private :: Write_info, Triu, Norm, Pinv, Qr

  interface Diag
     module procedure Diag_ele, Diag_mat
  end interface

contains
  
  function Bvgd(alt, mu1, mu2, s1, s2, rho) result(c)
    ! Bivariate Gaussian distribution
    implicit none
    real(DP), dimension(:), intent(in) :: alt
    real(DP), intent(in) :: mu1, mu2, s1, s2, rho
    real(DP), dimension(size(alt), size(alt)) :: c
    real(DP) :: d1, d2, apu1, apu2
    integer :: i, j, n
    n = size(alt)
    apu1 = 1.0/(2.0*PI*s1*s2*sqrt(1.0-rho**2))
    apu2 = -1.0/(2.0*(1.0-rho**2))
    do i = 1, n 
       do j = 1, n 
          d1 = alt(i) - mu1
          d2 = alt(j) - mu2
          C(i, j) = apu1*exp(apu2*((d1/s1)**2+(d2/s2)**2-(2.0*rho*d1*d2)/(s1*s2)))
       end do
    end do
  end function Bvgd
  
  function Norm(x) result(valu)
    ! L2-norm
    implicit none
    real(DP), dimension(:), intent(in) :: x
    real(DP) :: valu
    valu = sqrt(sum(x**2))
  end function Norm
  
  function Pinv(a) result(ai)
    ! Pseudo-inversion using SVD
    implicit none
    real(DP), dimension (:,:), intent(in) :: a
    real(DP), dimension(size(a,2), size(a,1)) :: ai
    real(DP), dimension(size(a,1), size(a,2)) :: u
    real(DP), dimension(size(a,2), size(a,2)) :: s, v
    integer :: r, nite
    real(DP), dimension(size(a,2)) :: d
    real(DP), parameter :: TOL = 1e-4
    real(DP), dimension(:,:), allocatable :: d2
    call Svd(a, u, s, v, nite)
    d = Diag(s)
    r = count(mask = d > TOL) ! take the singular values > tol
    allocate(d2(r,r))
    d2 = Diag(1.0 / d(1:r))
    ai = matmul(matmul(v(:,1:r), d2),transpose(u(:,1:r)))
    deallocate(d2)
  end function Pinv

  subroutine Svd(a, u, s, v, nite)
    ! Singular value decomposition
    implicit none
    real(DP), dimension(:,:), intent(in) :: a
    real(DP), dimension(size(a,1), size(a,2)), intent(out) :: u
    real(DP), dimension(size(a,2), size(a,2)), intent(out) :: s, v
    integer, intent(out) :: nite
    real(DP), dimension(size(a,1), size(a,2)) :: q1
    real(DP), dimension(size(a,1), size(a,1)) :: u1
    real(DP), dimension(size(a,2), size(a,2)) :: q, e 
    real(DP), dimension(size(a,2)) :: f
    real(DP), parameter :: EPS = 1e-4
    real(DP) :: err
    integer :: n    
    ! init u,v,u1
    u = 0.0
    v = 0.0
    u1 = 0.0
    forall (n=1:size(a,1)) u1(n,n) = 1.0
    forall (n=1:size(a,2)) v(n,n) = 1.0
    ! initial state:
    call Qr(a, q1, s)
    U = matmul(u1, q1)
    call Qr(transpose(s), q, s)
    v = matmul(v, q)
    ! iterate while converged:
    nite = 1
    do
       call Qr(transpose(s), q, s)
       u = matmul(u, q)
       call Qr(transpose(s), q, s)
       v = matmul(v, q)
       ! check the error:
       e = Triu(s)
       f = Diag(s)
       err = Norm(reshape(e, [size(e,1)*size(e,2)])) / Norm(f)
       nite = nite + 1
       if (err < EPS) exit
    end do
  end subroutine Svd

  function Diag_ele(A) result(v)
    ! Diagonal elements
    implicit none
    real(DP), dimension(:,:), intent(in) :: a
    real(DP), dimension(:), allocatable :: v
    integer :: i, n
    n = minval([size(a,1), size(a,2)]) 
    allocate(v(n)) 
    forall(i=1:n) v(i) = a(i,i)
  end function Diag_ele

  function Diag_mat(v) result(a)
    ! Diagonal matrix
    implicit none
    real(DP), dimension(:), intent(in) :: v
    real(DP), dimension(size(v), size(v)) :: a
    integer :: i
    a = 0.0
    forall(i=1:size(v)) a(i,i) = v(i)
  end function Diag_mat

  function Triu(a) result(au)
    ! Upper triangular part
    implicit none
    real(DP), dimension(:,:), intent(in) :: a
    real(DP), dimension(size(a,1), size(a,2)) :: au
    integer :: n, m, i, j
    au = 0.0
    m = size(a,1)
    n = size(a,2)
    do i = 1, m
       do j = i+1, n
          if (i+1 <= n) au(i,j) = a(i,j)
       end do
    end do
  end function Triu
  
  subroutine Qr(a,q,r)
    ! Modified Gram-Schmidt process
    implicit none
    real(DP), dimension(:,:), intent(in) :: a
    real(DP), dimension(size(a,1), size(a,2)), intent(out) :: q
    real(DP), dimension(size(a,2), size(a,2)), intent(out) :: r
    real(DP), dimension(size(a,1), size(a,2)) :: a0
    integer :: k,n
    n = size(a,2)
    q = 0.0
    r = 0.0
    a0 = a
    do k = 1, n
       r(k,k) = Norm(a0(:,k))
       q(:,k) = a0(:,k) / r(k,k)
       r(k,k+1:n) = matmul(q(:,k), a0(:,k+1:n))
       a0(:,k+1:n) = a0(:,k+1:n) - matmul(q(:,k:k), r(k:k,k+1:n))
    end do
  end subroutine Qr

  subroutine Levmar(Resfun, Jacfun, theta, cmat, r, rss, data)
    ! Levenberq-Marquardt optimization of "theta" parameters
    implicit none
    real(DP), dimension(:), intent(inout) :: theta
    real(DP), dimension(:,:), allocatable, intent(out) :: cmat
    real(DP), intent(out) :: rss
    type(obs), intent(in) :: data
    real(DP), dimension(data%nobs+size(theta)), intent(out) :: r    
    real(DP), dimension(data%nobs+size(theta), size(theta)) :: j
    real(DP), dimension(size(theta), size(theta)) :: jj, jj_mod
    real(DP), dimension(size(theta)) :: h, jr
    real(DP) :: rssold, lam
    real(DP), parameter :: ABSTOL = 1e-6 ! convergence tolerance
    integer :: ite
    ! this kind of interface is needed when 
    ! passing functions as arguments
    interface
       function Resfun(data, theta) result(r)
         use global_variables
         implicit none
         type(obs), intent(in) :: data
         real(DP), dimension(:), intent(in) :: theta
         real(DP), dimension(data%nobs) :: r
       end function Resfun
       function Jacfun(data, theta) result(j)
         use global_variables
         implicit none
         type(obs), intent(in) :: data
         real(DP), dimension(:), intent(in) :: theta
         real(DP), dimension(data%nobs, size(theta)) :: j
       end function Jacfun
    end interface
    ite = 0
    ! initial nudge factor
    lam = 1e-3
    ! initial residuals
    r = Resfun(data, theta)
    rss = sum(r**2) / (size(r) - size(theta))
    call Write_info(rss, ite, lam, .true.)
    ! iterate while converged
    do
       ite = ite + 1
       rssold = rss
       j = Jacfun(data, theta)
       jj = matmul(transpose(j), j)
       do
          jj_mod = (jj + lam*Diag(Diag(jj)))
          jr = matmul(transpose(j), r)
          h = matmul(Pinv(jj_mod), jr) 
          r = Resfun(data, theta+h)
          rss = sum(r**2) / (size(r) - size(theta))
          ! increase lam until rss decreases
          if (rss >= rssold) lam = lam*10.0
          if (abs((rss-rssold)/(rss+0.1)) < abstol) exit
          if (rss < rssold) exit
       end do
       call Write_info(rss,ite,lam)
       theta = theta+h
       if (abs((rss-rssold)/(rss+0.1)) < abstol) exit
       lam = lam/10.0
    end do
    ! estimate of the error covariance at the optimum
    cmat = Pinv(jj)*rss
  end subroutine Levmar

  subroutine Write_info(rss, ite, lam, isfirst)
    ! Write some info about LM-iteration
    implicit none
    real(DP), intent(in) :: rss, lam
    integer, intent(in) :: ite
    logical, optional, intent(in) :: isfirst
    character(len=*), parameter :: FORM = '(4x,i2,4x,f12.7,4x,es8.1)'
    if (present(isfirst) .and. isfirst) then
       write(*,*) '----------------------------------'       
       write(*,'(1x,a9,5x,a3,10x,a3)') 'Iteration', 'RSS', 'lam'
       write(*,*) '----------------------------------'
       write(*,FORM) ite, rss, lam
    else
       write(*,FORM) ite, rss, lam
    end if
  end subroutine Write_info
  
end module math
