AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror
INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libAnaPHPythiaDirectPhoton.la

#mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
#dist_mydata_DATA =

include_HEADERS =

noinst_HEADERS = \
  AnaPHPythiaDirectPhoton.h \
  AnaPHPythiaHistos.h

libAnaPHPythiaDirectPhoton_la_SOURCES = \
  AnaPHPythiaDirectPhoton.C \
  AnaPHPythiaHistos.C \
  AnaPHPythiaDirectPhoton_Dict.C

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libAnaPHPythiaDirectPhoton_la_LIBADD = \
  -lTOAD \
  -lPHPythia \
  -lphool \
  -lSubsysReco

BUILT_SOURCES = testexternals.C
noinst_PROGRAMS = testexternals
testexternals_LDADD = libAnaPHPythiaDirectPhoton.la
testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main() {return 0;}" >> $@

AnaPHPythiaDirectPhoton_Dict.C: \
  AnaPHPythiaDirectPhoton.h \
  AnaPHPythiaHistos.h \
  AnaPHPythiaDirectPhotonLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

CLEANFILES = *Dict* testexternals*
