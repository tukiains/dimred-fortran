module retrieval
  use global_variables
  implicit none
  
contains

  function Resfun(data, theta) result(r)
    ! Residual function
    use forward_model, only : Beer
    implicit none
    type(obs), intent(in) :: data
    real(DP), dimension(:), intent(in) :: theta
    real(DP), dimension(data%nobs+size(theta)) :: r
    real(DP), dimension(data%nprior) :: prof
    prof = Redu2full(data%P, data%prior, data%air, theta) 
    r(1:data%nobs) = data%t - Beer(data%cs, prof) 
    r = [r(1:data%nobs)/data%noise, theta] ! add parameters to the end
  end function Resfun

  function Jacfun(data, theta) result(j)
    ! Jacobian function
    use math, only : Diag
    use forward_model, only : Beer_jac
    implicit none
    type(obs), intent(in) :: data
    real(DP), dimension(:), intent(in) :: theta
    real(DP), dimension(data%nobs+size(theta), size(theta)) :: j
    real(DP), dimension(data%nprior) :: prof
    j = 0.0
    prof = Redu2full(data%p, data%prior, data%air, theta)
    j(1:data%nobs,:) = 0.0001 * matmul(matmul(Beer_jac(data%cs, prof), Diag(data%air)), data%p)
    j(1:data%nobs,:) = j(1:data%nobs,:) / data%noise
  end function Jacfun

  subroutine Reduce_dim(c, d, p, q)
    ! Truncation of prior covariance c using 
    ! d principal components so that p'p~=c
    use math, only : Svd, Diag
    implicit none
    real(DP), dimension(:,:), intent(in) :: c
    integer, intent(in) :: d
    real(DP), dimension(size(c,1), d), intent(out) :: p
    real(DP), dimension(size(c,1)), intent(out) :: q    
    real(DP), dimension(size(c,2), size(c,2)) :: s, v, u
    integer :: nite 
    call Svd(c, u, s, v, nite)
    q = Diag(s)
    p = matmul(u(:,1:d), Diag(sqrt(q(1:d))))
  end subroutine Reduce_dim

  function Redu2full(P, prior, air, theta) result(profile)
    ! New profile candidate from projection marix "P"
    ! and the parameters "theta"
    implicit none
    real(DP), dimension(:,:), intent(in) :: p
    real(DP), dimension(:), intent(in) :: theta, air, prior
    real(DP), dimension(size(prior)) :: profile
    profile = prior + matmul(p, theta)*air/1e9
  end function Redu2full

end module retrieval
