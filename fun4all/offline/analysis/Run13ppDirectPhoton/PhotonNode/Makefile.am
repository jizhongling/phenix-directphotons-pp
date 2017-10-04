AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror
INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = libPhotonNode.la

mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
dist_mydata_DATA = \
  clock-counts.root

include_HEADERS =

noinst_HEADERS = \
  FillHisto.h \
  EmcLocalRecalibrator.h \
  EmcLocalRecalibratorSasha.h \
  PhotonContainerClone.h \
  FillHistoMB.h \
  EmcLocalRecalibratorMB.h \
  PhotonContainerCloneMB.h \
  AnaToolsTowerID.h \
  AnaToolsPhoton.h \
  PhotonNodeLinkDef.h

libPhotonNode_la_SOURCES = \
  FillHisto.C \
  EmcLocalRecalibrator.C \
  EmcLocalRecalibratorSasha.C \
  PhotonContainerClone.C \
  FillHistoMB.C \
  EmcLocalRecalibratorMB.C \
  PhotonContainerCloneMB.C

nodist_libPhotonNode_la_SOURCES = \
  PhotonNode_Dict.C

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libPhotonNode_la_LIBADD = \
  -lDirectPhotonPP \
  -lfun4all \
  -lfun4allfuncs \
  -lphool \
  -lTOAD

BUILT_SOURCES = testexternals.C
noinst_PROGRAMS = testexternals
testexternals_LDADD = libPhotonNode.la
testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main() {return 0;}" >> $@

PhotonNode_Dict.C: \
  FillHisto.h \
  EmcLocalRecalibrator.h \
  EmcLocalRecalibratorSasha.h \
  PhotonContainerClone.h \
  FillHistoMB.h \
  EmcLocalRecalibratorMB.h \
  PhotonContainerCloneMB.h \
  PhotonNodeLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

CLEANFILES = *Dict* testexternals*
