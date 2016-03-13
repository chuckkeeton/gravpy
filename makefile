sief.so : sief.f90
	f2py --fcompiler=intelem --opt='-fast' -c sief.f90 -m sief

alphaf.so : hyp.f alphaf.f90
	f2py --fcompiler=intelem -lslatec -llapack --opt='-fast' -c hyp.f alphaf.f90 -m alphaf

nfwf.so : nfwf.f90
	f2py --fcompiler=intelem -lcuba --opt='-fast' -c nfwf.f90 -m nfwf

all : 
	make sief.so alphaf.so nfwf.so

clean :
	rm alphaf.so nfwf.so sief.so

re :
	make clean
	make all
