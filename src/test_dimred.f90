program test_dimred

  ! This program demonstrates how to:
  !
  ! 1. Simulate ground-based SWIR measurement (spectrum)
  ! 2. Retrieve ch4 profile from the measurement using 
  !    dimension reduction approach and Levenberg-Marquardt
  !    optimization
  !
  ! Simo Tukiainen, 2015

  use global_variables
  use math, only : Bvgd, Levmar
  use forward_model, only : Init_fmodel, Simulate_meas
  use retrieval
  implicit none

  type(obs) :: data
  real(DP), dimension(:), allocatable :: air, dens, alt, wn, t, q, theta, prior, profile, r
  real(DP), dimension(:,:), allocatable :: cs, p, c1, c2, c, cerr
  character(len=*), parameter :: CSFILE = "abs_coef.dat"
  character(len=*), parameter :: ATMOSFILE = "atmos.dat"
  integer :: nalt, d, nobs
  real(DP) :: noise, rss

  ! read cross sections and true/prior profiles
  call Init_fmodel(dens, prior, air, alt, cs, wn, CSFILE, ATMOSFILE)

  nalt = size(alt) ! number of altitudes
  nobs = size(wn) ! number of wavelengths

  ! simulate measurement
  write(*,*)
  write(*,*) 'Give std of noise. Signal is between 0-1 and the noise is additive ~ N(0,std^2)'
  write(*,*) 'Common value is e.g. 0.001 but you can try different values..'
  read(*,*) noise
  t = Simulate_meas(cs, dens, noise)

  ! Experimental prior covariance with two exponential components
  C1 = Bvgd(alt, 25.0_DP, 25.0_DP, 5.0_DP, 5.0_DP, 0.5_DP) * 1e7
  C2 = Bvgd(alt, 5.0_DP, 5.0_DP, 15.0_DP, 15.0_DP, 0.8_DP) * 5e6
  C = C1 + C2
  
  ! Its reduction
  write(*,*) ' Give number of principal components (~2-5 is good)'
  do 
     read(*,*) d
     if (d < 0) then 
        write(*,*) ' It can''t be negative! Try again...'
     elseif (d > nalt) then
        write(*,*) ' Too many! Try again...'
     else
        exit
     end if
  end do
  allocate(p(nalt, d), q(nalt))
  call Reduce_dim(c, d, p, q)

  ! Initial values for parameters
  allocate(theta(d))
  theta = 0.0
  
  ! Save needed variables to the data-structure
  data%t = t
  data%p = p
  data%d = d
  data%noise = noise
  data%cs = cs
  data%prior = prior
  data%air = air
  data%nobs = size(t)
  data%nprior = size(prior)

  ! Retrieve parameters
  allocate(r(nobs+d))
  call Levmar(resfun, jacfun, theta, cerr, r, rss, data)

  ! Inverted profile
  profile = Redu2full(P, prior, air, theta)

  ! Save results
  call Write_data(P,       'output/P.dat')
  call Write_data(profile, 'output/profile.dat')
  call Write_data(t,       'output/t.dat')
  call Write_data(dens,    'output/dens.dat')
  call Write_data(prior,   'output/prior.dat')
  call Write_data(air,     'output/air.dat')
  call Write_data(alt,     'output/alt.dat')
  call Write_data(r,       'output/r.dat')
  call Write_data(wn,      'output/wn.dat')
  call Write_data(d,       'output/d.dat')
  call Write_data(C,       'output/C.dat')
  call Write_data(q,       'output/q.dat')

end program test_dimred
