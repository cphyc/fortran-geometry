module CylinderGeometry

  implicit none
  ! private

  integer, parameter :: dp = selected_real_kind(15)
  real(dp), parameter :: pi = atan(1._dp) * 4._dp
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
  end type Cube_t

  type ptCube_t
     real(dp) :: x, y, z
  end type ptCube_t

  type ptCylinder_t
     real(dp) :: r, z
  end type ptCylinder_t

  public :: set_niter, Cylinder_t, Cube_t, dp

contains

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


  subroutine intersectionVolume(cyl, cube, volume, intersect)
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

  end subroutine intersectionVolume

  subroutine set_niter(n)
    integer, intent(in) :: n
    niter = n
  end subroutine set_niter

end module CylinderGeometry
