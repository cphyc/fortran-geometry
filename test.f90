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
  call assert(cyl%check(), 'Initialization')
  
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

  !----------------------------------------
  ! Testing cube
  !----------------------------------------
  call super_test('Cube')
  cube%center = (/10., 1.10, 0./)
  cube%a = 0.5
  cube%normal1 = (/10, 1, 0/)
  cube%normal2 = (/1, -10, 0/)

  call assert(cube%check(), 'Initialization')

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

  !----------------------------------------
  ! Test intersections
  !----------------------------------------
  
contains

  subroutine super_test(name)
    character(len=*), intent(in) :: name
    write(stdout, *) 'Testing ' // name
  end subroutine super_test

  subroutine assert(bool, msg)
    logical, intent(in) :: bool
    character(len=*), intent(in) :: msg

    if (.not. bool) then
       write(stdout, '(4x,a)') msg // '...failed!'
       stop
    else
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

end program test

