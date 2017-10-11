FFLAGS=-std=f2008 --warn-all -Ofast -fimplicit-none
FCC=gfortran

all: test run_test



%.o: %.f90
	$(FCC) -c $(FFLAGS) $^ -o $@

test: test_utils.o types.o geometry.o test.o
	$(FCC) $(FFLAGS) $^ -o $@

test_obb: test_utils.o types.o obbandcapsule.o test_obb.o
	$(FCC) $(FFLAGS) $^ -o $@

run_test: test_obb
	./$^

clean:
	rm -rf *.o *.mod
