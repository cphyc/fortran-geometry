program test_volume
  use distances
  use types
  use test_utils
  use volumes
  use utils
  implicit none

  type(Segment_t) :: segment
  type(Capsule_t) :: capsule
  type(OBB_t) :: box
  type(Line_t) :: line
  real(dp) :: t, d2, V

  !----------------------------------------
  ! Point to line
  !----------------------------------------
  call test_group('Point to line')
  call reset(box, segment, capsule)
  line%origin = (/0, 0, 0/)
  line%direction = (/1, 0, 0/)
  call PointLineDistanceSquared(line, (/0d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 0._dp, 1d-14), 'Point at origin')

  call PointLineDistanceSquared(line, (/1d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point at tip')

  line%origin = (/0, 0, 0/)
  line%direction = (/1, 1, 1/)
  call PointLineDistanceSquared(line, (/1d0, 1d0, 1d0/), t, d2)
  print*, t, d2
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point at tip (2)')

  line%origin = (/0, 0, 0/)
  line%direction = (/1, 0, 0/)
  call PointLineDistanceSquared(line, (/0d0, 1d0, 0d0/), t, d2)
  call assert(is_close(d2, 1._dp, 1d-14) .and. &
       is_close(t, 0._dp, 1d-14), 'Point on y-axis')

  call PointLineDistanceSquared(line, (/0d0, 0d0, 1d0/), t, d2)
  call assert(is_close(d2, 1._dp, 1d-14) .and. &
       is_close(t, 0._dp, 1d-14), 'Point on z-axis')

  call PointLineDistanceSquared(line, (/0d0, 0d0, 1d0/), t, d2)
  call assert(is_close(d2, 1._dp, 1d-14) .and. &
       is_close(t, 0._dp, 1d-14), 'Point on z-axis')

  call PointLineDistanceSquared(line, (/100d0, 10d0, 20d0/), t, d2)
  call assert(is_close(d2, 500._dp, 1d-14) .and. &
       is_close(t, 100._dp, 1d-14), 'Point far')

  line%origin = (/0, 0, 0/)
  line%direction = (/2, 0, 0/)
  call PointLineDistanceSquared(line, (/100d0, 10d0, 20d0/), t, d2)
  call assert(is_close(d2, 500._dp, 1d-14) .and. &
       is_close(t, 50._dp, 1d-14), 'Point far')
  call test_group_close()

  !----------------------------------------
  ! Point to segment
  !----------------------------------------
  call test_group('Point to segment')
  call reset(box, segment, capsule)
  call PointLineSegDistanceSquared(segment, (/0d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 0.5_dp, 1d-14), 'Point at origin')

  call PointLineSegDistanceSquared(segment, (/-1d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 0._dp, 1d-14), 'Point at start')

  call PointLineSegDistanceSquared(segment, (/1d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point at end')

  call PointLineSegDistanceSquared(segment, (/3d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 4._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point along axis')

  call PointLineSegDistanceSquared(segment, (/0d0, 1d0, 2d0/), t, d2)
  call assert(is_close(d2, 5._dp, 1d-14) .and. &
       is_close(t, 0.5_dp, 1d-14), 'Point perp. to axis')

  call PointLineSegDistanceSquared(segment, (/1d0, 1d0, 2d0/), t, d2)
  call assert(is_close(d2, 5._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point perp. to axis (at tip)')

  call PointLineSegDistanceSquared(segment, (/2d0, 1d0, 2d0/), t, d2)
  call assert(is_close(d2, 6._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point perp. to axis (beyond tip)')

  call test_group_close()

  !----------------------------------------
  ! CpasuleContains
  !----------------------------------------
  call test_group('Capsule contains')
  call reset(box, segment, capsule)
  call assert(CapsuleContains(capsule, (/0._dp, 0._dp, 0._dp/)), 'center of capsule')
  call assert(CapsuleContains(capsule, (/1._dp, 0._dp, 0._dp/)), 'in capsule (x-axis)')
  call assert(CapsuleContains(capsule, (/0._dp, 0.5_dp, 0._dp/)), 'in capsule (y-axis)')
  call assert(CapsuleContains(capsule, (/0._dp, 0._dp, 0.5_dp/)), 'in capsule (z-axis)')
  call assert(CapsuleContains(capsule, (/2._dp, 0._dp, 0._dp/)), 'side of capsule')
  call assert(.not. CapsuleContains(capsule, (/2.00001_dp, 0._dp, 0._dp/)), 'outside of capsule')

  call reset(box, segment, capsule)
  capsule%segment%start = (/-10, 0, 0/)
  capsule%segment%end = (/10, 0, 0/)
  call assert(CapsuleContains(capsule, (/0._dp, 0._dp, 0._dp/)), 'center of capsule')
  call assert(CapsuleContains(capsule, (/1._dp, 0._dp, 0._dp/)), 'in capsule')
  call assert(CapsuleContains(capsule, (/2._dp, 0._dp, 0._dp/)), 'side of capsule')
  call assert(.not. CapsuleContains(capsule, (/-10.00001_dp, 0._dp, 0._dp/)), 'outside of capsule')

  call assert(is_close(BoxVolume(box), 8._dp, 1d-14), 'box volume')

  ! Set spherical capsule
  capsule%segment%end = capsule%segment%start
  call assert(is_close(CapsuleVolume(capsule), 4*pi/3*capsule%r**3, 1d-14), 'capsule volume')

  call test_group_close()

  !----------------------------------------
  ! Volumes
  !----------------------------------------
  call test_group('Volumes')
  call reset(box, segment, capsule)
  call assert(is_close(CapsuleVolume(capsule), 4*pi/3*capsule%r**3 &
       + pi * capsule%r**2 * norm(capsule%segment%end - capsule%segment%start),&
       1d-14), 'capsule volume')

  ! Compute intersection volume
  call reset(box, segment, capsule)
  capsule%r = 3
  call CapsuleOBBVolume(box, capsule, V)
  call assert(is_close(V, BoxVolume(box), 1d-14), 'box in capsule')

  call reset(box, segment, capsule)
  capsule%segment%start = (/-10, 0, 0/)
  capsule%segment%end = (/10, 0, 0/)
  ! capsule%r = sqrt(2._dp)
  print*, CapsuleContains(capsule, (/0, 0, 0/)*1._dp)
  ! call CapsuleOBBVolume(box, capsule, V)
  print*, V, pi*capsule%r**2 * 2
  ! call reset(box, segment, capsule)
  ! capsule%r = 1._dp
  ! capsule%segment%start = (/0, 0, -10/)
  ! capsule%segment%end = (/0, 0, 10/)
  ! call CapsuleOBBVolume(box, capsule, V)
  ! ! call print_details()
  ! print*, V, pi*capsule%r**2 * 2
  ! call assert(is_close(&
  !      V, pi*capsule%r**2 * norm(capsule%segment%end - capsule%segment%start),&
  !      1d-14), 'capsule-cylinder in box, else outside')
contains

  subroutine reset(box, segment, capsule)
    type(OBB_t), intent(inout) :: box
    type(Segment_t), intent(inout) :: segment
    type(Capsule_t), intent(inout) :: capsule

    box%origin = (/0, 0, 0/)
    box%u = (/1, 0, 0/)
    box%v = (/0, 1, 0/)
    box%w = (/0, 0, 1/)
    box%extents = (/1, 1, 1/)

    segment%start = (/-1, 0, 0/)
    segment%end = (/1, 0, 0/)

    capsule%segment = segment
    capsule%r = 1._dp
  end subroutine reset

  subroutine print_details
    print*, 'Box:'
    print*, '----'
    print*, box%origin
    print*, box%u
    print*, box%v
    print*, box%w
    print*, box%extents
    print*, ''
    print*, 'Capsule:'
    print*, '--------'
    print*, capsule%segment%start
    print*, capsule%segment%end
    print*, capsule%r
    print*, ''
    print*, 'Intersection'
    print*, '------------'
    print*, 'd=', sqrt(d2)
  end subroutine print_details

end program
