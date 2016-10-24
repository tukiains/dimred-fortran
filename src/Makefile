FC = gfortran
FCFLAGS = -Wall -Wextra -Wconversion -pedantic -fbounds-check -frange-check -std=f2003
LIBS = 
PROGRAMS = test_dimred

all: $(PROGRAMS)

test_dimred.o: global_variables.o math.o lcg.o forward_model.o retrieval.o
test_dimred: global_variables.o math.o lcg.o forward_model.o retrieval.o

%: %.o  
	$(FC) $(FCFLAGS) -o $@ $^ ${LIBS}

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod test_dimred *~ output/*.dat	
