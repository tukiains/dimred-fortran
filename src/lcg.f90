module lcg
  ! Linear congruential random number generator
  use global_variables
  implicit none

  private :: Iterate_lcg 
  public :: Init_lcg, Rand_lcg

  ! parameters from "Numerical Recipies"
  integer(DP), parameter, private :: M = 4294967296_DP
  integer(DP), parameter, private :: A = 1664525_DP
  integer(DP), parameter, private :: C = 1013904223_DP
  integer(DP), private, save :: rand_val = 0_DP

  interface Rand_lcg
     module procedure Rand_lcg1, Rand_lcgN, Rand_lcgNN
  end interface Rand_lcg

contains

  subroutine Init_lcg(seed)
    implicit none
    integer, intent(in) :: seed
    if (seed < 0 .or. seed > m) then
       write(*,*) 'Error: bad lcg seed'
    else
       rand_val = seed
    end if
  end subroutine Init_lcg
  
  subroutine Iterate_lcg()
    implicit none
    rand_val = mod((A*rand_val+C), M)
  end subroutine Iterate_lcg
  
  function Rand_lcg1() result(r)
    implicit none
    real(DP) :: r
    call Iterate_lcg()
    r = abs(real(rand_val,kind=DP) / real(M,kind=DP))
  end function Rand_lcg1

  function Rand_lcgN(n) result(r)
    implicit none
    integer, intent(in) :: n
    real(DP), dimension(n) :: r
    integer :: i
    do i = 1, n
       call Iterate_lcg()
       r(n) = abs(real(rand_val, kind=DP) / real(M, kind=DP))
    end do
  end function Rand_lcgN

  function Rand_lcgNN(n,nn) result(r)
    implicit none
    integer, intent(in) :: n,nn
    real(DP), dimension(n,nn) :: r
    integer :: i,j
    do i = 1, n
       do j = 1, nn
          call Iterate_lcg()
          r(i,j) = abs(real(rand_val, kind=DP) / real(M, kind=DP))
       end do
    end do
  end function Rand_lcgNN

end module lcg
