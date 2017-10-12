module volumes
  use types
  use distances
  use utils
  implicit none

  private

  integer :: MAX_DEPTH = 9
  public :: CapsuleOBBVolume
  public :: BoxContains, CapsuleContains
  public :: BoxVolume, CapsuleVolume

  public :: set_depth, get_depth
contains
  subroutine CapsuleOBBVolume(box, capsule, V)
    ! Compute the volume of the intersection of a box with a capsule
    type(OBB_t), intent(in) :: box
    type(Capsule_t), intent(in) :: capsule
    real(dp), intent(out) :: V

    real(dp) :: d2

    call CapsuleOBBDistanceSquared(box, capsule, d2)

    V = 0
    if (d2 > 0) then
       ! Their is no interaction, volume is null
       ! Nothing to do
    else
       ! Some interaction (capsule in box or the opposite)
       call divide(box ,capsule, 1, 0, V)
    end if
  contains
    recursive subroutine divide(box, capsule, axis, depth, V)
      type(OBB_t), intent(in) :: box
      type(Capsule_t), intent(in) :: capsule
      integer, intent(in) :: depth, axis

      real(dp), intent(inout) :: V

      type(OBB_t) :: box1, box2
      real(dp) :: dir(3), d2
      integer :: naxis, nin, i, j, k

      real(dp) :: point(3)

      logical, parameter :: debug = .false.
      character(len=1) :: ikey(3)

      ! Compute distance from box to capsule
      call CapsuleOBBDistanceSquared(box, capsule, d2)
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

         if (nin == 8) then
            V = V + BoxVolume(box)
            if (debug) call write_padding(depth, 'all inside')
            ! if (debug) &
            ! write(*, '(3(es14.7, 2x),i3)') &
            !      box%origin, -1
         else if (depth == MAX_DEPTH) then
            ! Estimate volume by number of points in volume
            V = V + BoxVolume(box) * nin / 8._dp
            if (debug) call write_padding(depth, 'max depth')

            ! if (debug) &
            ! write(*, '(3(es14.7, 2x),i3)') &
            !      box%origin, nin
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
            call divide(box1, capsule, naxis, depth+1, V)
            if (debug) call write_padding(depth, '> RIGHT ' // ikey(axis))
            call divide(box2, capsule, naxis, depth+1, V)
         end if
      end if

    end subroutine divide

  end subroutine CapsuleOBBVolume

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
    type(OBB_t), intent(in) :: box
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

    volume = 4*pi/3*capsule%r**3 + pi * capsule%r**2 * norm(vec)
  end function CapsuleVolume

  real(dp) function BoxVolume(box) result(volume)
    type(OBB_t), intent(in) :: box

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

end module volumes
