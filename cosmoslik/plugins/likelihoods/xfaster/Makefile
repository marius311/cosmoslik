include ../../../../Makefile.inc

default: xfaster_likelihood_mod.so

xfaster_likelihood_mod.so: xfaster_likelihood_mod.f90
		f2py -c xfaster_likelihood_mod.f90 -m xfaster_likelihood_mod $(LAPACK) $(F2PYFLAGS)
		
clean:
	-rm -f xfaster_likelihood_mod.so
		
	
