AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror
INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libMissingRatio.la

mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
dist_mydata_DATA = \
  ERTEff-SM.root

include_HEADERS =

noinst_HEADERS = \
  AnaToolsTrack.h \
  AnaTrk.h \
  ERTSimTrigger.h \
  MissingRatio.h \
  PhotonEff.h \
  IsolationCut.h \
  Isolation.h \
  HadronResponse.h \
  MissingRatioLinkDef.h

libMissingRatio_la_SOURCES = \
  AnaTrk.C \
  ERTSimTrigger.C \
  MissingRatio.C \
  PhotonEff.C \
  IsolationCut.C \
  Isolation.C \
  HadronResponse.C \
  MissingRatio_Dict.C

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libMissingRatio_la_LIBADD = \
  -lDirectPhotonPP \
  -lTOAD \
  -lemc \
  -lEmcStatic \
  -lcteval \
  -lcompactCNT \
  -lfun4all \
  -lfun4allfuncs \
  -lphool

BUILT_SOURCES = testexternals.C
noinst_PROGRAMS = testexternals
testexternals_LDADD = libMissingRatio.la
testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main() {return 0;}" >> $@

MissingRatio_Dict.C: \
  MissingRatio.h \
  PhotonEff.h \
  IsolationCut.h \
  Isolation.h \
  HadronResponse.h \
  MissingRatioLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

CLEANFILES = *Dict* testexternals*
