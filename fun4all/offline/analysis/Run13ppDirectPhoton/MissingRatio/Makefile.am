AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror
INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = libMissingRatio.la

mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
dist_mydata_DATA = 

include_HEADERS =

noinst_HEADERS = \
  MissingRatio.h \
  PhotonEff.h \
  IsolationCut.h \
  AnaTrk.h \
  AnaToolsTowerID.h \
  AnaToolsCluster.h \
  MissingRatioLinkDef.h \
  IsolationCutLinkDef.h

libMissingRatio_la_SOURCES = \
  MissingRatio.C \
  PhotonEff.C \
  IsolationCut.C

nodist_libMissingRatio_la_SOURCES = \
  AnaTrk.C \
  MissingRatio_Dict.C \
  IsolationCut_Dict.C

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
  PhotonEff.h \
  MissingRatioLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

IsolationCut_Dict.C: \
  IsolationCut.h \
  IsolationCutLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

CLEANFILES = *Dict* testexternals*
