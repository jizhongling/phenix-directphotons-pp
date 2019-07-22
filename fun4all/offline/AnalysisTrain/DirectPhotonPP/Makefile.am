AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror
INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libDirectPhotonPP.la

mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
dist_mydata_DATA = \
  Warnmap_Run13pp510.txt \
  Run13pp510_EMC_TOF_Correction.root \
  Run13pp_RunbyRun_Calib.dat \
  warn_all_run13pp500gev.dat \
  ecorr_run13pp500gev.txt \
  ecorr_run_run13pp500gev.txt \
  tcorr_run13pp500gev.txt

include_HEADERS = \
  DCDeadmapChecker.h \
  PhotonContainer.h \
  Photon.h \
  PhotonERT.h \
  SpinPattern.h

noinst_HEADERS = \
  AnaToolsTowerID.h \
  AnaToolsCluster.h \
  AnaToolsTrigger.h \
  HistogramBooker.h \
  EmcLocalRecalibrator.h \
  EmcLocalRecalibratorSasha.h \
  DirectPhotonPP.h \
  PhotonNode.h \
  PhotonHistos.h \
  DirectPhotonPPLinkDef.h

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libDirectPhotonPP_la_LIBADD = \
  -lTOAD \
  -luspin \
  -lodbc \
  -lemc \
  -lcompactCNT \
  -lphool \
  -lfun4all \
  -lfun4allfuncs \
  -lSubsysReco

libDirectPhotonPP_la_SOURCES = \
  HistogramBooker.cc \
  EmcLocalRecalibrator.cc \
  EmcLocalRecalibratorSasha.cc \
  DCDeadmapChecker.cc \
  DirectPhotonPP.cc \
  PhotonContainer.cc \
  Photon.cc \
  PhotonERT.cc \
  PhotonNode.cc \
  PhotonHistos.cc \
  SpinPattern.cc \
  DirectPhotonPP_Dict.C

# Rule for generating CINT dictionaries from class headers.
DirectPhotonPP_Dict.C: \
  HistogramBooker.h \
  EmcLocalRecalibrator.h \
  EmcLocalRecalibratorSasha.h \
  DCDeadmapChecker.h \
  DirectPhotonPP.h \
  PhotonContainer.h \
  Photon.h \
  PhotonERT.h \
  PhotonNode.h \
  PhotonHistos.h \
  SpinPattern.h \
  DirectPhotonPPLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

CLEANFILES = *Dict* testexternals*
