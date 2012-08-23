include ../../../../Makefile.inc

default: recfast

recfast: recfast.for
	$(F90C) $(F90FLAGS) -std=legacy -ffixed-line-length-none recfast.for -o recfast

clean: 
	-rm -f recfast