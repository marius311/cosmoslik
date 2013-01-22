include ../../../../Makefile.inc

default: pywmap.so

wmap/libwmap7.a: 
	cd wmap && make libwmap7.a

pywmap.so:  pywmap.f90 wmap/libwmap7.a
	f2py --f90flags=-Iwmap -c pywmap.f90 -m pywmap -Lwmap -lwmap7 $(CFITSIO) $(LAPACK) $(F2PYFLAGS)

clean:
	-rm -f pywmap.so
	cd wmap && make clean
