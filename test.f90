program test
  use CylinderGeometry
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
       stdout=>output_unit, &
       stderr=>error_unit

  type(Cylinder_t) :: cyl
  type(Cube_t) :: cube

  logical :: intersect
  real(dp) :: volume

  cyl%center = (/1., 1., 1./)
  cyl%r = 1.0
  cyl%h = 1.0
  cyl%normal1 = (/1., 1., 1./)

  !----------------------------------------
  ! Testing cylinder
  !----------------------------------------
  call super_test('Cylinder')
  call assert(cyl%setup(), 'Initialization')

  ! Checking the normals are perpendicular one from another
  call assert(dot_product(cyl%normal1, cyl%normal2) == 0, "normals (1)")
  call assert(dot_product(cyl%normal1, cyl%normal3) == 0, "normals (2)")
  call assert(dot_product(cyl%normal2, cyl%normal3) == 0, "normals (3)")

  ! Check that the projection is correct
  block
    real(dp) :: tmp(3)
    type(ptCylinder_t) :: tmpOut

    tmp = cyl%center
    call cyl%project(tmp, tmpOut)
    call assert(is_close(tmpOut%r, 0._dp, 1.d-10) .and. &
         is_close(tmpOut%z, 0._dp, 1.d-10), "Project onto center")

    tmp = cyl%center + cyl%normal1
    call cyl%project(tmp, tmpOut)
    call assert(is_close(tmpOut%r, 0._dp, 1.d-10) &
         .and. is_close(tmpOut%z, cyl%h, 1.d-10), "Project onto axis 1")

    tmp = cyl%center + cyl%normal2
    call cyl%project(tmp, tmpOut)
    call assert(is_close(tmpOut%r, cyl%r, 1.d-10) &
         .and. is_close(tmpOut%z, 0._dp, 1.d-10), "Project onto axis 2")

    tmp = cyl%center + cyl%normal3
    call cyl%project(tmp, tmpOut)
    call assert(is_close(tmpOut%r, cyl%r, 1.d-10) &
         .and. is_close(tmpOut%z, 0._dp, 1.d-10), "Project onto axis 3")

    tmp = cyl%center + cyl%normal1 + cyl%normal2
    call cyl%project(tmp, tmpOut)
    call assert(is_close(tmpOut%r, cyl%r, 1.d-10) &
         .and. is_close(tmpOut%z, cyl%h, 1.d-10), "Project onto axis 1+2")

    ! Checking point in/out
    tmp = cyl%center
    call assert(cyl%contains(tmp), 'Point at center')
    tmp = cyl%center + cyl%normal1 * cyl%h * 0.9
    call assert(cyl%contains(tmp), 'Point on axis')
    tmp = cyl%center - cyl%normal1 * cyl%h * 1.1
    call assert(.not. cyl%contains(tmp), 'Point on axis (out)')
    tmp = cyl%center + cyl%normal2 * cyl%r * 0.99
    call assert(cyl%contains(tmp), 'Point on radius')
    tmp = cyl%center + cyl%normal2 * cyl%r * 1.01
    call assert(.not. cyl%contains(tmp), 'Point on radius (out)')
    tmp = cyl%center + cyl%normal2 * cyl%r * 0.99 + cyl%normal1 * cyl%h * 0.99
    call assert(cyl%contains(tmp), 'Point close to cap')
    tmp = cyl%center + cyl%normal2 * cyl%r * 1.01 + cyl%normal1 * cyl%h * 0.99
    call assert(.not. cyl%contains(tmp), 'Point close to cap (out) (1)')
    tmp = cyl%center + cyl%normal2 * cyl%r * 0.99 + cyl%normal1 * cyl%h * 1.01
    call assert(.not. cyl%contains(tmp), 'Point close to cap (out) (2)')
  end block
  ! Check volumes
  call assert(is_close(cyl%volume(), cyl%r**2 * pi * cyl%h, 1d-10), &
       "volume")

  !----------------------------------------
  ! Testing cube
  !----------------------------------------
  call super_test('Cube')
  cube%center = (/10., 1.10, 0./)
  cube%a = 0.5
  cube%normal1 = (/10, 1, 0/)
  cube%normal2 = (/1, -10, 0/)

  call assert(cube%setup(), 'Initialization')

  call assert(dot_product(cube%normal1, cube%normal2) == 0, "normals (1)")
  call assert(dot_product(cube%normal1, cube%normal3) == 0, "normals (2)")
  call assert(dot_product(cube%normal2, cube%normal3) == 0, "normals (3)")

  block
    real(dp) :: tmp(3)
    type(ptCube_t) :: tmpOut

    tmp = cube%center
    call cube%project(tmp, tmpOut)
    call assert(is_close(tmpOut%x, 0._dp, 1.d-10) .and. &
         is_close(tmpOut%y, 0._dp, 1.d-10) .and. &
         is_close(tmpOut%z, 0._dp, 1.d-10), &
         "Project onto center")

    tmp = cube%center + cube%normal1 * cube%a
    call cube%project(tmp, tmpOut)
    call assert(is_close(tmpOut%x, cube%a, 1.d-10), &
         "Project onto axis 1")

    tmp = cube%center + cube%normal2 * cube%a * 2
    call cube%project(tmp, tmpOut)
    call assert(is_close(tmpOut%y, cube%a * 2, 1.d-10), &
         "Project onto axis 2")

    tmp = cube%center + cube%normal3 * cube%a * (-3)
    call cube%project(tmp, tmpOut)
    call assert(is_close(tmpOut%z, cube%a * (-3), 1.d-10), &
         "Project onto axis 3")

    tmp = cube%center + cube%normal1 + 2*cube%normal2 - 3*cube%normal3
    call cube%project(tmp, tmpOut)
    call assert(is_close(tmpOut%x, 1._dp, 1.d-10) .and. &
         is_close(tmpOut%y, 2._dp, 1.d-10) .and. &
         is_close(tmpOut%z, -3._dp, 1.d-10), &
         "Project onto axis 1+2+3")

    tmp = cube%center
    call assert(cube%contains(tmp), 'Point at center')
    tmp = cube%center + cube%normal1 * cube%a * 0.49
    call assert(cube%contains(tmp), 'Point normal (x)')
    tmp = cube%center - cube%normal2 * cube%a * 0.49
    call assert(cube%contains(tmp), 'Point normal (y)')
    tmp = cube%center + cube%normal3 * cube%a * 0.49
    call assert(cube%contains(tmp), 'Point normal (z)')

    tmp = cube%center + cube%normal1 * cube%a * 0.51
    call assert(.not. cube%contains(tmp), 'Point normal (out, x)')
    tmp = cube%center - cube%normal2 * cube%a * 0.51
    call assert(.not. cube%contains(tmp), 'Point normal (out, y)')
    tmp = cube%center + cube%normal3 * cube%a * 0.51
    call assert(.not. cube%contains(tmp), 'Point normal (out, z)')

    tmp = cube%center + cube%normal1 * cube%a * 0.49 + &
         cube%normal2 * cube%a * 0.49 + &
         cube%normal3 * cube%a * 0.49
    call assert(cube%contains(tmp), 'Point at edge')
    tmp = cube%center + cube%normal1 * cube%a * 0.51 + &
         cube%normal2 * cube%a * 0.49 + &
         cube%normal3 * cube%a * 0.49
    call assert(.not. cube%contains(tmp), 'Point at edge (x, out)')
  end block

  ! Check volumes
  call assert(is_close(cube%volume(), cube%a**3, 1d-10), &
       "volume")
  !----------------------------------------
  ! Test intersections
  !----------------------------------------
  call super_test('Intersections')
  ! Setup a unit-cylinder along the z axis
  cyl%center = (/0., 0., 0./)
  cyl%r = 1.0
  cyl%h = 1.0
  cyl%normal1 = (/0., 0., 1./)
  if (.not. cyl%setup()) call yolo()

  ! Setup a cube and move it along
  cube%center = (/0., 0., 0./)
  cube%a = 0.5
  cube%normal1 = (/1, 5, 0/)
  cube%normal2 = (/-5, 1, 0/)
  if (.not. cube%setup()) call yolo()

  ! The cube is within the cylinder -> volume = cube volume
  block
    real(dp) :: volume
    logical :: intersect
    call cyl%intersectionVolume(cube, volume, intersect)
    call assert(is_close(volume, cube%volume(), 1d-10), "volume (cube in cylinder)")
    call assert(intersect, "Detection", .true.)

    ! Move cube (very close to boundary)
    cube%center = (/cyl%r - sqrt(3._dp) * cube%a, 0._dp, 0._dp/)
    call cyl%intersectionVolume(cube, volume, intersect)
    call assert(is_close(volume, cube%volume(), 1d-10), &
         "volume (cube in cylinder, not centered)")
    call assert(intersect, "Detection", .true.)

    ! Move cube outstide
    cube%center = (/2., 0., 0. /)
    call cyl%intersectionVolume(cube, volume, intersect)
    call assert(is_close(volume, 0._dp, 1d-10), &
         "volume (cube outside cylinder)")
    call assert(.not. intersect, "Detection", .true.)

    ! Put a tiny cube at the boundary of the cylinder. The volume is
    ! approx. 0.5 Vcube at ± 1/sqrt(N) (noise)
    cube%a = 0.001
    cube%normal1 = (/1, 0, 0/)
    cube%normal2 = (/0, 1, 0/)
    if (.not. cube%setup()) call yolo()
    ! Put cube at boundary
    cube%center = (/cyl%r, 0._dp, 0._dp/)
    call cyl%intersectionVolume(cube, volume, intersect)
    call assert(is_close(volume, cube%volume()/2, cube%volume()*2/sqrt(niter * 1._dp)), &
         'volume (tiny cube on cylinder edge)')
    call assert(intersect, "Detection", .true.)

    ! Put a giant cube and a cylinder inside
    cube%a = cyl%r * 10
    cube%center = (/0, 0, 0/)
    call cyl%intersectionVolume(cube, volume, intersect)
    call assert(is_close(volume, cyl%volume(), 1.d-10), &
         'volume (cylinder in cube)')
    call assert(intersect, "Detection", .true.)

    ! Move the cylinder to the cube's edge
    ! The cylinder is half in the cube, half out -> volume should be half
    cyl%center = (/cube%a/2, 0._dp, 0._dp/)
    call cyl%intersectionVolume(cube, volume, intersect)
    call assert(is_close(volume, cyl%volume()/2, cyl%volume()*2/sqrt(niter * 1._dp)), &
         "volume (cylinder at cube's edge)")
    call assert(intersect, "Detection", .true.)

    ! Tilt a little bit the cylinder
    cyl%normal1 = (/0, 5, 10/)
    if (.not. cyl%setup()) call yolo()
    call cyl%intersectionVolume(cube, volume, intersect)
    call assert(is_close(volume, cyl%volume()/2, cyl%volume()*2/sqrt(niter * 1._dp)), &
         "volume (cylinder at cube's edge, tilted)")
    call assert(intersect, "Detection", .true.)
  end block

contains

  subroutine super_test(name)
    character(len=*), intent(in) :: name
    write(stdout, *) 'Testing ' // name
  end subroutine super_test

  subroutine assert(bool, msg, silent)
    logical, intent(in) :: bool
    character(len=*), intent(in) :: msg
    logical, intent(in), optional :: silent

    logical :: print_msg

    if (present(silent)) then
       print_msg = .not. silent
    else
       print_msg = .true.
    end if

    if (.not. bool) then
       write(stdout, '(4x,a)') msg // '...failed!'
       stop
    else
       if (print_msg) &
            write(stdout, '(4x,a)') msg // '...ok!'
    end if
  end subroutine assert

  function is_close(A, B, aerr)
    real(dp), intent(in) :: A, B, aerr
    logical :: is_close

    if (abs(A - B) < aerr) then
       is_close = .true.
    else
       is_close = .false.
    end if
  end function is_close

  subroutine yolo
    write(stderr, *) 'An error occured. Crashing.'
    stop
  end subroutine yolo
end program test

