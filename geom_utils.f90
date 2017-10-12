module geom_utils
  use geom_types
  implicit none

  real(dp), parameter :: pi = atan(1._dp) * 4._dp
contains
  elemental integer function next_axis(axis) result(next)
    integer, intent(in) :: axis

    ! This maps 3->1, 1->2 and 2->3
    next = mod(axis, 3) + 1
  end function next_axis

  subroutine write_padding(depth, msg)
    integer, intent(in) :: depth
    character(len=*), intent(in) :: msg

    integer :: i
    do i = 1, depth
       write(*, '(1x, "|", 3x)', advance='no')
    end do
    write(*, *) msg
  end subroutine write_padding

end module geom_utils
