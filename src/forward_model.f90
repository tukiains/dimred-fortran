module forward_model
  use global_variables
  implicit none

contains

  function Beer(cs, dens) result(y)
    ! Beer-Lambert law
    implicit none
    real(DP), dimension(:,:) :: cs
    real(DP), dimension(:) :: dens
    real(DP), dimension(size(cs,1)) :: y
    y = exp(-matmul(cs,dens*1e5)) 
  end function Beer
  
  function Beer_jac(cs, dens) result(j)
    ! Jacobian of the Beer-Lambert law
    implicit none
    real(DP), dimension(:,:), intent(in) :: cs
    real(DP), dimension(:), intent(in) :: dens
    real(DP), dimension(size(cs,1)) :: y
    real(DP), dimension(size(cs,1), size(cs,2)) :: j
    integer :: n
    y = Beer(cs, dens)
    forall (n=1:size(cs,1)) j(n,:) = -cs(n,:)*y(n)
  end function Beer_jac

  function Simulate_meas(cs, dens, noise) result(t)
    ! Simulates a noisy measurement
    use lcg
    implicit none
    real(DP), dimension(:,:), intent(in) :: cs
    real(DP), dimension(:), intent(in) :: dens
    real(DP), intent(in) :: noise
    real(DP), dimension(size(cs,1)) :: t
    call Init_lcg(123123)
    t = Beer(cs, dens) + noise*(sum(Rand_lcg(size(t),12),2)-6.0)
  end function Simulate_meas

  subroutine Init_fmodel(dens, prior, air, alt, cs, wn, csfile, atmosfile)
    ! Read cross-sections, "true" density and prior density profiles    
    implicit none
    real(DP), dimension(:,:), allocatable, intent(out) :: cs
    real(DP), dimension(:), allocatable, intent(out) :: dens, air, alt, wn, prior
    character(len=*), intent(in) :: csfile, atmosfile
    integer :: nalt, err, n, nwn, m
    ! alt, dens, air, prior
    open(unit=10, file=atmosfile, status='old')
    nalt = 0
    do
       read(10, *, iostat=err)
       if (err < 0) exit
       nalt = nalt + 1
    end do
    rewind(10)
    allocate(dens(nalt), air(nalt), alt(nalt), prior(nalt))
    do n = 1, nalt
       read(10,*) alt(n), dens(n), air(n), prior(n)
    end do
    close(10)
    ! cross sections
    open(unit=10, file=csfile, status='old')
    nwn = 0
    do
       read(10, *, iostat=err)
       if (err < 0) exit
       nwn = nwn + 1
    end do
    rewind(10)
    allocate(cs(nwn,nalt), wn(nwn))
    do n = 1, nwn
       read(10, *) wn(n), (cs(n,m), m=1, nalt)
    end do
    close(10)
  end subroutine Init_fmodel

end module forward_model
