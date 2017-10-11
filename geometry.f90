module CylinderGeometry

  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
       stdout=>output_unit, &
       stderr=>error_unit

  implicit none
  private

  integer, parameter :: dp = selected_real_kind(15) !

  integer :: MAX_DEPTH = 10
  integer :: niter = 1000

  type :: Cylinder_t
     real(dp), dimension(3) :: center
     real(dp), dimension(3) :: normal1
     real(dp), dimension(3) :: normal2
     real(dp), dimension(3) :: normal3
     real(dp) :: r, h
   contains
     procedure :: setup => cylinderSetup
     procedure :: intersectionVolume
     procedure :: volume => cylinderVolume
     procedure :: draw => cylinderDraw
     procedure :: contains => cylinderContains
     procedure :: project => cylinderProject
  end type Cylinder_t

  type :: Cube_t
     real(dp), dimension(3) :: center
     real(dp), dimension(3) :: normal1
     real(dp), dimension(3) :: normal2
     real(dp), dimension(3) :: normal3
     real(dp) :: a
   contains
     procedure :: setup => cubeSetup
     procedure :: volume => cubeVolume
     procedure :: draw => cubeDraw
     procedure :: contains => cubeContains
     procedure :: project => cubeProject
     procedure :: corners => cubeCorners
  end type Cube_t

  type ptCube_t
     real(dp) :: x, y, z
  end type ptCube_t

  type ptCylinder_t
     real(dp) :: r, z
  end type ptCylinder_t

  public :: set_niter, get_niter, Cylinder_t, Cube_t, dp, ptCube_t, ptCylinder_t
  public :: set_depth, get_depth, intersect

contains

  function intersect(cube, cyl)
    logical :: intersect
    type(Cylinder_t), intent(in) :: cyl
    type(Cube_t), intent(in) :: cube

    real(dp) :: d2, corners(8, 3)
    integer :: i

    d2 = max(cyl%r, cyl%h)**2 + 3*cube%a**2

    ! If centers too far, just return
    if (sum((cyl%center - cube%center)**2) > d2) then
       intersect = .false.
       return
    end if

    ! If any point of cube in cylinder, intersection is true
    corners = cube%corners()
    do i = 1, 8
       if (cyl%contains(corners(i, :))) then
          intersect = .true.
          return
       end if
    end do

    ! If cylinder center in cube
    if (cube%contains(cyl%center)) then
       intersect = .true.
       return
    end if

    ! If the cube is bigger than the cylinder
    ! Case 1: the cube can be inside the cylinder
    if (cube%a < min(cyl%r, cyl%h)) then
       ! Case 2: the cube may contain the cylinder
       
    else if (cube%a > max(cyl%r, cyl%h)) then
       ! Case 3: none of the above
       ! the cube's size is between the radius and heigth of the cylinder
    end if

  end function intersect

  function cylinderVolume(self) result (volume)
    class(Cylinder_t) :: self
    real(dp) :: volume
    volume = self%r**2 * pi * self%h
  end function cylinderVolume

  function cubeVolume(self) result (volume)
    class(Cube_t) :: self
    real(dp) :: volume
    volume = self%a**3
  end function cubeVolume

  function cubeSetup(self) result (valid)
    class(Cube_t) :: self
    logical ::valid

    real(dp) :: tmp

    tmp = sqrt(dot_product(self%normal1, self%normal1))
    self%normal1 = self%normal1 / tmp

    tmp = sqrt(dot_product(self%normal2, self%normal2))
    self%normal2 = self%normal2 / tmp

    if (dot_product(self%normal1, self%normal2) /= 0._dp) then
       valid = .false.
       return
    end if

    self%normal3 = (/ &
         self%normal1(2)*self%normal2(3) - self%normal1(3)*self%normal2(2), &
         self%normal1(3)*self%normal2(1) - self%normal1(1)*self%normal2(3), &
         self%normal1(1)*self%normal2(2) - self%normal1(2)*self%normal2(1) &
         /)
    valid = .true.
  end function cubeSetup

  function cylinderSetup(self) result(valid)
    class(Cylinder_t) :: self
    logical :: valid

    real(dp) :: norm, n2, N(3)

    n2 = dot_product(self%normal1, self%normal1)
    norm = sqrt(n2)
    self%normal1 = self%normal1 / norm
    N = self%normal1

    ! Build a normal perpendicular to the first one for convenience
    if      (N(1) * N(2) /= 0.) then
       self%normal2 = (/N(2), -N(1), 0._dp/)
    else if (N(1) * N(3) /= 0.) then
       self%normal2 = (/N(3), 0._dp, -N(1)/)
    else
       self%normal2 = (/0._dp, N(3), N(2)/)
    end if
    self%normal2 = self%normal2 / sqrt(dot_product(self%normal2, self%normal2))

    self%normal3 = (/ &
         self%normal1(2)*self%normal2(3) - self%normal1(3)*self%normal2(2), &
         self%normal1(3)*self%normal2(1) - self%normal1(1)*self%normal2(3), &
         self%normal1(1)*self%normal2(2) - self%normal1(2)*self%normal2(1) &
         /)
    valid = .true.

  end function cylinderSetup

  subroutine cylinderProject(self, pt, ptOut)
    class(Cylinder_t) :: self
    real(dp), dimension(3), intent(in) :: pt
    type(ptCylinder_t), intent(out) :: ptOut

    real(dp) :: tmpPt(3)

    tmpPt = pt - self%center
    ! Get height
    ! h = X.N
    ptOut%z = dot_product(self%normal1, tmpPt)

    ! Project on plane
    ! r = | X - X.N * N |
    ptOut%r = sqrt(sum((tmpPt - self%normal1 * ptOut%z)**2))

  end subroutine cylinderProject

  subroutine cubeProject(self, pt, ptOut)
    ! Project the point in the local frame of the cube
    class(Cube_t) :: self
    real(dp), dimension(3), intent(in) :: pt
    type(ptCube_t), intent(out) :: ptOut

    real(dp) :: tmpPt(3)

    tmpPt = pt - self%center

    ptOut%x = dot_product(self%normal1, tmpPt)
    ptOut%y = dot_product(self%normal2, tmpPt)
    ptOut%z = dot_product(self%normal3, tmpPt)
  end subroutine cubeProject

  function cubeCorners(self) result(pt)
    class(Cube_t), intent(in) :: self
    real(dp), dimension(8, 3) :: pt

    integer :: i, j, k

    do k = 0, 1
       do j = 0, 1
          do i = 0, 1
             pt(1 + i + 2*j + 4*k, :) = self%center + &
                  self%normal1 * (i - 0.5_dp) * self%a + &
                  self%normal2 * (j - 0.5_dp) * self%a + &
                  self%normal3 * (k - 0.5_dp) * self%a
          end do
       end do
    end do

  end function cubeCorners

  function cylinderContains(self, pt) result (res)
    class(Cylinder_t) :: self
    real(dp), dimension(3) :: pt
    logical :: res

    type(ptCylinder_t) :: ptLoc

    ! Project point in local ref
    call self%project(pt, ptLoc)

    ! Setup point height
    if (abs(ptLoc%z) > self%h) then
       res = .false.
       return
    end if

    ! Setup point radius
    if (ptLoc%r > self%r) then
       res = .false.
       return
    end if

    res = .true.

  end function cylinderContains

  function cubeContains(self, pt) result (res)
    class(Cube_t) :: self
    real(dp), dimension(3) :: pt
    logical :: res

    type(ptCube_t) :: ptLoc

    real(dp) :: ao2

    ! Project point in local referential
    call self%project(pt, ptLoc)

    ao2 = self%a / 2._dp

    res = (abs(ptLoc%x) < ao2) .and. &
         (abs(ptLoc%y) < ao2) .and. &
         (abs(ptLoc%z) < ao2)

  end function cubeContains

  function cylinderDraw(self) result (point)
    class(Cylinder_t) :: self
    real(dp), dimension(3) :: point

    real(dp), dimension(3) :: tmp
    real(dp) :: r, theta


    call random_number(tmp)
    ! Compute r, theta
    r = sqrt(tmp(1)) * self%r
    theta = 2 * pi * tmp(2)

    ! Project in absolute referential
    point(:) = &
         self%normal2 * r * cos(theta) + &
         self%normal3 * r * sin(theta) + &
         (tmp(3) - 0.5_dp) * 2 * self%h * self%normal1 + &
         self%center

  end function cylinderDraw

  function cubeDraw(self) result (point)
    class(Cube_t) :: self
    real(dp), dimension(3) :: point

    real(dp), dimension(3) :: tmp

    call random_number(tmp)
    tmp = (tmp - 0.5) * 2 * self%a
    point = self%normal1 * tmp(1) + &
         self%normal2 * tmp(2) + &
         self%normal3 * tmp(3) + &
         self%center

  end function cubeDraw


  subroutine intersectionVolume(cyl, cube, volume, intersect, method)
    class(Cylinder_t), intent(in) :: cyl
    type(Cube_t), intent(in) :: cube
    real(dp), intent(out) :: volume
    logical, intent(out)  :: intersect
    character(len=*), intent(in), optional :: method

    real(dp) :: d2, Vcyl, Vcube
    real(dp), dimension(3) :: pt
    integer :: nin, ntot

    real(dp) :: corners(8, 3)

    if (present(method)) then
       if (trim(method) == 'MC') then
          call intersectionVolumeMC(cyl, cube, volume, intersect)
          return
       else if (trim(method) /= 'bisect') then
          write(stderr, *) 'Unknown method. Aborting.'
          stop
       end if
    end if

    d2 = max(cyl%r, cyl%h)**2 + 3*cube%a**2

    ! If centers too far, just return
    if (sum((cyl%center - cube%center)**2) > d2) then
       intersect = .false.
       volume = 0.0_dp
    end if

    ! Volumes may overlap
    corners = cube%corners()
    volume = 0._dp
    call divide(corners, 1, 1, volume)

    intersect = volume > 0._dp

  contains
    recursive subroutine divide(corners, depth, splitDim, volume)
      ! Corners are given as follow
      !       7 +--------+ 8
      !        /|       /|
      !     5 +------- +6|          z
      !       | |      | |          ^
      !       |3+------|-+ 4        | /y
      !       |/       |/           |/
      !     1 +--------+ 2          +----> x

      real(dp), intent(in) :: corners(8, 3) ! The location of the corners
      integer, intent(in) :: depth
      real(dp), intent(inout) :: volume ! The volume
      integer, intent(in) :: splitDim   ! The dimension for the next split

      real(dp) :: newCornersRight(8, 3), newCornersLeft(8, 3), tmpVolume

      integer :: i, nin, newDim
      nin = 0

      ! Count the number of corners in/out
      do i = 1, 8
         if (cyl%contains(corners(i, :))) then
            nin = nin + 1
         end if
      end do

      ! If all points are within cylinder, the volume is the volume of the box
      if (nin == 8) then                ! Box inside volume
         call compute_volume(corners, volume)
      else if (depth == MAX_DEPTH) then ! Estimate volume by number of corners
         call compute_volume(corners, tmpVolume)
         volume = volume + tmpVolume * nin / 8
      else
         ! Check that the box's closest point
         if (splitDim == 1) then         ! Split along x direction
            newCornersLeft(2, :) = (corners(1, :) + corners(2, :)) / 2
            newCornersLeft(4, :) = (corners(3, :) + corners(4, :)) / 2
            newCornersLeft(6, :) = (corners(5, :) + corners(6, :)) / 2
            newCornersLeft(8, :) = (corners(7, :) + corners(8, :)) / 2

            newCornersRight(1, :) = (corners(1, :) + corners(2, :)) / 2
            newCornersRight(3, :) = (corners(3, :) + corners(4, :)) / 2
            newCornersRight(5, :) = (corners(5, :) + corners(6, :)) / 2
            newCornersRight(7, :) = (corners(7, :) + corners(8, :)) / 2
         else if (splitDim == 2) then  ! Split along y direction
            newCornersLeft(3, :) = (corners(1, :) + corners(3, :)) / 2
            newCornersLeft(4, :) = (corners(2, :) + corners(4, :)) / 2
            newCornersLeft(7, :) = (corners(5, :) + corners(7, :)) / 2
            newCornersLeft(8, :) = (corners(6, :) + corners(8, :)) / 2

            newCornersRight(1, :) = (corners(1, :) + corners(3, :)) / 2
            newCornersRight(2, :) = (corners(2, :) + corners(4, :)) / 2
            newCornersRight(5, :) = (corners(5, :) + corners(7, :)) / 2
            newCornersRight(6, :) = (corners(6, :) + corners(8, :)) / 2
         else                          ! Split along z direction
            newCornersLeft(5, :) = (corners(1, :) + corners(5, :)) / 2
            newCornersLeft(6, :) = (corners(2, :) + corners(6, :)) / 2
            newCornersLeft(7, :) = (corners(3, :) + corners(7, :)) / 2
            newCornersLeft(8, :) = (corners(4, :) + corners(8, :)) / 2

            newCornersRight(1, :) = (corners(1, :) + corners(5, :)) / 2
            newCornersRight(2, :) = (corners(2, :) + corners(6, :)) / 2
            newCornersRight(3, :) = (corners(3, :) + corners(7, :)) / 2
            newCornersRight(4, :) = (corners(4, :) + corners(8, :)) / 2
         end if
         newDim = mod(splitDim, 3) + 1
         call divide(newCornersLeft, depth+1, newDim, volume)
         call divide(newCornersRight, depth+1, newDim, volume)
      end if

    end subroutine divide

    subroutine compute_volume(corners, volume)
      real(dp), intent(in) :: corners(8, 3) ! The location of the corners
      real(dp), intent(out) :: volume

      volume = sqrt(&
           dot_product(corners(2, :) - corners(1, :), corners(2, :) - corners(1, :)) * &
           dot_product(corners(3, :) - corners(1, :), corners(3, :) - corners(1, :)) * &
           dot_product(corners(5, :) - corners(1, :), corners(5, :) - corners(1, :)))
    end subroutine compute_volume

  end subroutine intersectionVolume

  subroutine intersectionVolumeMC(cyl, cube, volume, intersect)
    class(Cylinder_t), intent(in) :: cyl
    type(Cube_t), intent(in) :: cube
    real(dp), intent(out) :: volume
    logical, intent(out)  :: intersect

    real(dp) :: d2, Vcyl, Vcube
    real(dp), dimension(3) :: pt
    integer :: nin, ntot


    d2 = max(cyl%r, cyl%h)**2 + 3*cube%a**2

    ! If centers too far, just return
    if (sum((cyl%center - cube%center)**2) > d2) then
       intersect = .false.
       volume = 0.0_dp
    end if

    ! Volumes may overlap. Draw random points from the smallest volume
    nin = 0
    ntot = 0

    Vcyl = cyl%volume()
    Vcube = cube%volume()
    if (Vcyl < Vcube) then
       do while (ntot < niter)
          pt = cyl%draw()
          ntot = ntot + 1
          if (cube%contains(pt)) then
             nin = nin + 1
          end if
       end do
       volume = Vcyl * nin / ntot
    else
       do while (ntot < niter)
          pt = cube%draw()
          ntot = ntot + 1
          if (cyl%contains(pt)) then
             nin = nin + 1
          end if
       end do
       volume = Vcube * nin / ntot
    end if

    intersect = volume > 0._dp

  end subroutine intersectionVolumeMC

  subroutine set_niter(n)
    integer, intent(in) :: n
    niter = n
  end subroutine set_niter

  subroutine get_niter(n)
    integer, intent(out) :: n
    n = niter
  end subroutine get_niter


  subroutine set_depth(d)
    integer, intent(in) :: d
    MAX_DEPTH = d
  end subroutine set_depth

  subroutine get_depth(d)
    integer, intent(out) :: d
    d = MAX_DEPTH
  end subroutine get_depth

end module CylinderGeometry
