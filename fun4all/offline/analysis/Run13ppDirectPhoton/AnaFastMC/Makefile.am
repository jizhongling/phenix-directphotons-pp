AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror
INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libAnaFastMC.la

mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
dist_mydata_DATA = \
  dead_eff_run13pp500gev.dat

include_HEADERS =

noinst_HEADERS = \
  AnaFastMC.h \
  AnaToolsTowerID.h

libAnaFastMC_la_SOURCES = \
  AnaFastMC.C \
  AnaFastMC_Dict.C

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libAnaFastMC_la_LIBADD = \
  -lDirectPhotonPP \
  -lTOAD \
  -lPHPythia \
  -lphool \
  -lSubsysReco

BUILT_SOURCES = testexternals.C
noinst_PROGRAMS = testexternals
testexternals_LDADD = libAnaFastMC.la
testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main() {return 0;}" >> $@

AnaFastMC_Dict.C: \
  AnaFastMC.h \
  AnaFastMCLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

CLEANFILES = *Dict* testexternals*
