AUTOMAKE_OPTIONS = foreign

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

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib

include_HEADERS = \
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
  DirectPhotonPPLinkDef.h

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
  DirectPhotonPP.cc \
  PhotonContainer.cc \
  Photon.cc \
  PhotonERT.cc \
  PhotonNode.cc \
  SpinPattern.cc \
  DirectPhotonPP_Dict.C

# Rule for generating CINT dictionaries from class headers.

DirectPhotonPP_Dict.C: \
  HistogramBooker.h \
  EmcLocalRecalibrator.h \
  EmcLocalRecalibratorSasha.h \
  DirectPhotonPP.h \
  PhotonContainer.h \
  Photon.h \
  PhotonERT.h \
  PhotonNode.h \
  SpinPattern.h \
  DirectPhotonPPLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
