FC = gfortran
FFLAG = -llapack -lblas

symmetryED.so : mod_symmetryED.f90
	gfortran -shared -fPIC -o $@ $< $(FFLAG)

.PHONY : clean

clean :
	rm -f *.mod *.o *.out *.so

