
# Makefile 
# exmaple:http://lasp.colorado.edu/cism/CISM_DX/code/CISM_DX-0.50/required_packages/octave-forge/main/sparse/Makefile
MKOCTFILE= mkoctfile

MKOCTFLAGS   = -lgsl -O3 
SRC = runc.cc 
SRC2 = runc_nomem.cc 
SRC3 = runc_nomem_nosnap.cc
SRC4 = runc_nomem_rand.cc

all: runc_nomem.oct runc_nomem_rand.oct

runc.oct: runc.cc
	$(MKOCTFILE) $(SRC) $(MKOCTFLAGS)

runc_nomem.oct: runc_nomem.cc
	$(MKOCTFILE) $(SRC2) $(MKOCTFLAGS)

runc_nomem_nosnap.oct: runc_nomem_nosnap.cc
	$(MKOCTFILE) $(SRC3) $(MKOCTFLAGS)

runc_nomem_rand.oct: runc_nomem_rand.cc
	$(MKOCTFILE) $(SRC4) $(MKOCTFLAGS)
