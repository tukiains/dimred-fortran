module global_variables
  implicit none
  
  integer, parameter :: DP = 8
  real(DP), parameter :: PI = 3.14159265358979
  
  type :: obs
     real(DP), dimension(:,:), allocatable :: p, cs
     real(DP), dimension(:), allocatable :: prior, t, dens, air
     real(DP) :: noise
     integer :: d, nobs, nprior
  end type obs

  interface Write_data
     module procedure Write_1d_data, Write_2d_data, Write_int
  end interface Write_data
  
contains
  
  subroutine Write_int(table, filename)
    ! write integer point to file
    implicit none
    integer, intent(in) :: table
    character(len=*), intent(in) :: filename
    open(unit=10, file=filename, status='unknown')
    write(10,*) table
    close(10)
  end subroutine Write_int

  subroutine Write_1d_data(table, filename)
    ! write real array to file
    implicit none
    real(DP), dimension(:), intent(in) :: table
    character(len=*), intent(in) :: filename
    integer :: n
    open(unit=10, file=filename, status='unknown')
    do n=1, size(table,1)
       write(10,*) table(n)
    end do
    close(10)
  end subroutine Write_1d_data
  
  subroutine Write_2d_data(table, filename)
    ! write real matrix to file
    implicit none
    real(DP), dimension(:,:), intent(in) :: table
    character(len=*), intent(in) :: filename
    integer :: n,m
    open(unit=10, file=filename, status='unknown')
    do n=1,size(table,1)
       write(10,*) (table(n,m), m=1, size(table,2))
    end do
    close(10)
  end subroutine Write_2d_data

end module global_variables
