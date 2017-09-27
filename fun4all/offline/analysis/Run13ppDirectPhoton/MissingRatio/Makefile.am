AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror
INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = libMissingRatio.la

mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
dist_mydata_DATA = \
  dead_eff_run13pp500gev.dat

include_HEADERS =

noinst_HEADERS = \
  MissingRatio.h \
  AnaTrk.h \
  AnaToolsTowerID.h \
  AnaToolsCluster.h \
  MissingRatioLinkDef.h

libMissingRatio_la_SOURCES = \
  MissingRatio.C

nodist_libMissingRatio_la_SOURCES = \
  AnaTrk.C \
  MissingRatio_Dict.C

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libMissingRatio_la_LIBADD = \
  -lemc \
  -lEmcStatic \
  -lcteval \
  -lcompactCNT \
  -lfun4all \
  -lfun4allfuncs \
  -lphool \
  -lTOAD

#  -lEmcNode \
#  -lemc \
#  -lemcEmbed4all \
#  -lcteval \
#  -lcompactCNT \
#  -lfun4all \
#  -lfun4allfuncs \
#  -lphool \
#  -lTOAD

BUILT_SOURCES = testexternals.C
noinst_PROGRAMS = testexternals
testexternals_LDADD = libMissingRatio.la
testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main() {return 0;}" >> $@

MissingRatio_Dict.C: \
  MissingRatio.h \
  MissingRatioLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

CLEANFILES = *Dict* testexternals*
