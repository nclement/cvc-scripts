#--------------------------------------
# Makefile for Nathan Clement's CVC-lab helpful scripts.
#
# Just run
#   make all
#
# and this will fix all the scripts. It will warn if the binaries in
# Makefile.def do not exist.

include Makefile.def

TARGETS=RMSDMulti RMSDBasic summaryStatsCol
OPT_FLAGS=-std=c++11 -O3
DEBUG_FLAGS=-std=c++11 -g -O0
FLAGS=${OPT_FLAGS}

SUBDIRS = docking/postprocessing hbond discrepancy/TA_algs ${TEXMOL_SAMPLES}
ifneq ($(wildcard ${TEXMOL_DIR}/samples),)
	SUBDIRS+=${TEXMOL_DIR}/samples
endif

all : ${TARGETS} subdirs test-existence
	@echo
	@echo "If there have been no errors, your build should be complete."

RMSDMulti : RMSDMulti.cpp
	g++ ${FLAGS} $< -o $@

RMSDBasic : RMSDBasic.cpp
	g++ ${FLAGS} $< -o $@

% : %.cpp
	g++ ${FLAGS} $< -o $@

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	echo making '$@'
	@${MAKE} -C $@

test-existence: Makefile.def
	@# Check and see which ones are in existence.
	@echo
	@echo "Checking Environment variables set in Makefile.def..."
	@test -f ${MolSurf} || echo "  Error: MolSurf does not exist [${MolSurf}]. Please install and update Makefile.def"
	@test -f ${MolEnergy} || echo "  Error: MolEnergy does not exist. Please install and update Makefile.def"
	@test -f ${TEXMOL_QT4} || echo "  Error: TexMol does not exist. Please install and update Makefile.def"
	@test -f ${F2DOCK_REFACTORED} || echo "  Error: F2Dock does not exist. Please install and update Makefile.def"
	@test -d ${F3DOCK_DIR} || echo "  Error: F3Dock does not exist. Please install and update Makefile.def"
	@test -f ${PDB2PQR} || echo "  Error: PDB2PQR does not exist. Please install and update Makefile.def"
	@test -d ${MGL_DIR} || echo "  Error: MGLTools does not exist. Please install and update Makefile.def"
