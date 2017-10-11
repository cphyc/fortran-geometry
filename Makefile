FFLAGS=-std=f2008 --warn-all -Ofast -fimplicit-none
FCC=gfortran

all: test run_test



%.o: %.f90
	$(FCC) -c $(FFLAGS) $^ -o $@

test:  types.o test_utils.o geometry.o test.o
	$(FCC) $(FFLAGS) $^ -o $@

test_obb: utils.o types.o test_utils.o distances.o test_obb.o
	$(FCC) $(FFLAGS) $^ -o $@

test_volumes: utils.o types.o test_utils.o distances.o volumes.o test_volume.o
	$(FCC) $(FFLAGS) $^ -o $@

run_test: test_obb test_volumes
	./test_obb
	./test_volumes

clean:
	rm -rf *.o *.mod
