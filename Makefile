FFLAGS=-std=f2008 --warn-all -Ofast -fimplicit-none
FCC=gfortran

all: run_test


%.o: %.f90
	$(FCC) -c $(FFLAGS) $^ -o $@

test_obb: types.o  utils.o test_utils.o distances.o test_obb.o
	$(FCC) $(FFLAGS) $^ -o $@

test_volumes: types.o utils.o test_utils.o distances.o volumes.o test_volume.o
	$(FCC) $(FFLAGS) $^ -o $@

run_%: %
	./$^

run_test: test_obb test_volumes run_test_obb run_test_volumes


clean:
	rm -rf *.o *.mod
