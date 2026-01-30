CXX = gfortran
CXXFLAGS = -O2

H2CO_F90 = pes_data.f90 pes_wrapper.f90 abinit_pot.f90 ref_pot.f90

H2CO_OBJ = pes_data.o pes_wrapper.o abinit_pot.o ref_pot.o

all: libh2co.so

$(H2CO_OBJ) :
	$(CXX) $(CXXFLAGS) -fPIC -c $(H2CO_F90)

libh2co.so: $(H2CO_OBJ)
	$(CXX) -shared -o libh2co.so -fPIC $(H2CO_OBJ)
	mv libh2co.so ../libh2co.so

clean:
	rm -f *.so *.o *.mod
