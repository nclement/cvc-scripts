include ../Makefile.def

all : test-existence

test-existence: ../Makefile.def
	@echo
	@echo "Testing Environment variables for hbond..."
	@test -f ${F2DOCK_REFACTORED_BIN}/testHbond || echo "  WARN:  Potential error in hbond build. Try turning on the HAVE_LIBMOL option in F2Dock"
	@test -f ${OBABEL} || echo "  WARN:  Potential error in hbond build: Could not find obabel. Please update Makefile.def"
	@test -f ${CHARMM} || echo "  WARN:  Potential error in hbond build: Could not find Charmm. Please update Makefile.def"
	@echo
