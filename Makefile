FFLAGS=-std=f2008 --warn-all
FCC=gfortran

all: test run_test



%.o: %.f90
	$(FCC) -c $(FFLAGS) $^ -o $@

test: test.o geometry.o
	$(FCC) $^ -o $@

run_test: test
	./test

clean:
	rm -rf *.o *.mod
