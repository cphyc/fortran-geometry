module geom_volumes
  use geom_types
  use geom_distances
  use geom_utils
  implicit none

  private
  logical, parameter :: debug = .false.
  integer :: MAX_DEPTH = 9
  public :: BoxContains, CapsuleContains
  public :: BoxVolume, CapsuleVolume
  public :: CapsuleBoxVolume, CapsuleBoxIntegrate

  public :: set_depth, get_depth

  interface
     real(dp) function callback_fun(x)
       use geom_types, only: dp
       real(dp), intent(in), dimension(3) :: x
     end function callback_fun
  end interface

contains
  subroutine CapsuleBoxVolume(box, capsule, V)
    ! Compute the volume of the intersection of a box with a capsule
    type(Box_t), intent(in) :: box
    type(Capsule_t), intent(in) :: capsule
    real(dp), intent(out) :: V

    real(dp) :: d2

    call CapsuleBoxDistanceSquared(box, capsule, d2)

    V = 0
    if (d2 > 0) then
       ! Their is no interaction, volume is null
       ! Nothing to do
    else
       ! Some interaction (capsule in box or the opposite)
       call divide(box, capsule, 1, 0, V, dummy, force_refine=.false.)
    end if
  contains
    real(dp) function dummy(X)
      real(dp), intent(in) :: X(3)
      dummy = 1._dp
    end function dummy
  end subroutine CapsuleBoxVolume

  subroutine CapsuleBoxIntegrate(box, capsule, V, integrand)
    ! Compute the volume of the intersection of a box with a capsule
    !
    ! Note
    ! ----
    ! The integrand function is called with arguments in the frame of
    ! the capsule. The starting point is at [0, 0, 0], the ending
    ! point at [0, 0, L]


    type(Box_t), intent(in) :: box
    type(Capsule_t), intent(in) :: capsule
    real(dp), intent(out) :: V

    procedure(callback_fun) :: integrand

    real(dp) :: d2

    call CapsuleBoxDistanceSquared(box, capsule, d2)

    V = 0
    if (d2 > 0) then
       ! Their is no interaction, volume is null
       ! Nothing to do
    else
       ! Some interaction (capsule in box or the opposite)
       call divide(box, capsule, 1, 0, V, integrand, force_refine=.true.)
    end if
  end subroutine CapsuleBoxIntegrate

  recursive subroutine divide(box, capsule, axis, depth, V, callback, force_refine)
    type(Box_t), intent(in) :: box
    type(Capsule_t), intent(in) :: capsule
    integer, intent(in) :: depth, axis
    procedure(callback_fun) :: callback
    logical, intent(in) :: force_refine

    real(dp), intent(inout) :: V

    type(Box_t) :: box1, box2
    real(dp) :: dir(3), d2
    integer :: naxis, nin, i, j, k

    real(dp) :: point(3), dV, u1(3), u2(3), u3(3)
    character(len=1) :: ikey(3)

    ! Compute distance from box to capsule
    call CapsuleBoxDistanceSquared(box, capsule, d2)
    if (d2 > 0) then
       if (debug) call write_padding(depth, 'disjoint')
    else
       ! Compute which point falls in capsule
       nin = 0
       do k = -1, 1, 2
          do j = -1, 1, 2
             do i = -1, 1, 2
                point = box%origin &
                     + i * box%u * box%extents(1) &
                     + j * box%v * box%extents(2) &
                     + k * box%w * box%extents(3)
                if (CapsuleContains(capsule, point)) then
                   nin = nin + 1
                end if
             end do
          end do
       end do

       ! Convert position in capsule's frame
       point = box%origin - capsule%segment%start

       ! Build direct basis (z-axis along capsule segment)
       u3 = capsule%segment%end - capsule%segment%start
       u3 = u3 / norm2(u3)
       u1 = [&
            -u3(2) - u3(3), &
             u3(1) - u3(3), &
             u3(1) + u3(2) &
             ]
       u1 = u1 / norm2(u1)
       u2 = vector_product(u3, u1)
       ! Project point on basis
       point = [dot_product(u1, point), dot_product(u2, point), dot_product(u3, point)]

       if (nin == 8 .and. .not. force_refine) then
          dV = BoxVolume(box) * callback(point)
          if (debug) call write_padding(depth, 'all inside')
          V = V + dV
       else if (depth == MAX_DEPTH) then
          ! Estimate volume by number of points in volume
          dV = BoxVolume(box) * nin / 8._dp * callback(point)
          V = V + dV
          if (debug) call write_padding(depth, 'max depth')
       else
          ! Split box in 2
          box1 = box
          box2 = box

          if      (axis == 1) then
             dir = box%u
          else if (axis == 2) then
             dir = box%v
          else  ! (axis == 3)
             dir = box%w
          end if

          box1%origin = box%origin - dir * box%extents(axis) / 2
          box1%extents(axis) = box%extents(axis) / 2

          box2%origin = box%origin + dir * box%extents(axis) / 2
          box2%extents(axis) = box%extents(axis) / 2

          naxis = next_axis(axis)
          ikey = (/'x', 'y', 'z'/)
          if (debug) call write_padding(depth, '> LEFT ' // ikey(axis))
          call divide(box1, capsule, naxis, depth+1, V, callback, force_refine)
          if (debug) call write_padding(depth, '> RIGHT ' // ikey(axis))
          call divide(box2, capsule, naxis, depth+1, V, callback, force_refine)
       end if
    end if

  end subroutine divide

  !----------------------------------------
  ! Containing functions
  !----------------------------------------
  logical function CapsuleContains(capsule, point) result(inside)
    type(Capsule_t), intent(in) :: capsule
    real(dp), intent(in) :: point(3)

    real(dp) :: t, d2

    call PointLineSegDistanceSquared(capsule%segment, point, t, d2)
    ! print*, t, sqrt(d2)
    inside = d2 <= capsule%r**2

  end function CapsuleContains

  logical function BoxContains(box, point) result(inside)
    type(Box_t), intent(in) :: box
    real(dp), intent(in) :: point(3)

    real(dp) :: vec(3)
    integer :: idim

    ! Project on 3 axis
    vec = (/ &
         dot_product(box%u, point) / box%extents(1), &
         dot_product(box%v, point) / box%extents(2), &
         dot_product(box%w, point) / box%extents(3)  &
         /)

    inside = .true.
    do idim = 1, 3
       if (vec(idim) < -1._dp .or. vec(idim) > 1._dp) then
          inside = .false.
          return
       end if
    end do

  end function BoxContains

  !----------------------------------------
  ! Volume functions
  !----------------------------------------
  real(dp) function CapsuleVolume(capsule) result(volume)
    type(Capsule_t), intent(in) :: capsule

    real(dp) :: vec(3)

    vec = capsule%segment%end - capsule%segment%start

    volume = 4*pi/3*capsule%r**3 + pi * capsule%r**2 * norm2(vec)
  end function CapsuleVolume

  real(dp) function BoxVolume(box) result(volume)
    type(Box_t), intent(in) :: box

    ! Each extent is the half extent
    volume = box%extents(1) * box%extents(2) * box%extents(3) * 8
  end function BoxVolume

  !----------------------------------------
  ! Setter/getter
  !----------------------------------------
  subroutine set_depth(depth)
    integer, intent(in) :: depth

    MAX_DEPTH = depth
  end subroutine set_depth

  subroutine get_depth(depth)
    integer, intent(out) :: depth
    depth = MAX_DEPTH
  end subroutine get_depth

end module geom_volumes
