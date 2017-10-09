# Fortran-geometry

This is a very specific and very limited module in Fortran to easily volumetric properties of geometrical objects.
For now it supports:
* [x] cube
* [x] cylinder
* [ ] box (easy to implement, just not done)
* [ ] sphere (should be easy too)

## How to use it?

The soft is presented as a module that you can `use` (`use CylinderGeometry`). You can then create geometrical objects. Properties of these objects are accessed using the Fortran 2008 Object-Oriented syntax, e.g.:
```fortran
program test
  use CylinderGeometry

  type(Cylinder_t) :: cyl
  type(Cube_t) :: cube

  real(8) :: volume
  logical :: intersect

  ! Define a cylinder
  cyl%center = (/1., 1., 1./)
  cyl%r = 1.0  ! This is the radius
  cyl%h = 1.0  ! This is the height
  cyl%normal1 = (/1., 1., 1./)

  ! Define a cube
  cube%center = (/10., 8., 0./)
  cube%a = 1.5 ! This is the width of the cube

  ! Let's get the volume
  print*, 'Volume of the cube:', cube%volume()

  ! Let's get the volume that intersect between the cube and the cylinder
  call cyl%intersectionVolume(cube, volume, intersect)

  print*, 'Cube and volume intersect?', intersect
  print*, 'Volume of intersection   :', volume
end program test
```

## A note on the "algorithm"

The "algorithm" is extremely poor for now: the code just takes the smallest volume out of the two volumes, draws random points within it and counts the ratio of points that fall in the other volume w.r.t. to the number of points drawn. As such, the convergence is in `1/sqrt(N)`, where `N` is the number of draws.

As a side note, if the intersection volume is tiny w.r.t. the smallest volume, it is likely that you will get no point falling in the intersection and the code will erroneously report a non-existing intersection.

You should this algorithm as a first attempt to get a quick-and-dirty answer. For more precise codes, it will probably be required to develop a more complex code.

## Need some more examples?

You can read the test or send me an email to get more information.

## Licensing
The code follows the MIT license.
