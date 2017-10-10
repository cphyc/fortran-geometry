FFLAGS=-std=f2008 --warn-all -Ofast -fimplicit-none
FCC=gfortran

all: test run_test



%.o: %.f90
	$(FCC) -c $(FFLAGS) $^ -o $@

test: geometry.o test.o
	$(FCC) $(FFLAGS) $^ -o $@

run_test: test
	./test

clean:
	rm -rf *.o *.mod
