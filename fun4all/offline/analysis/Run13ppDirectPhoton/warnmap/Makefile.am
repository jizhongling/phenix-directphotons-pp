AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror
INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libdgpp_warnmap.la

#mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
#  dist_mydata_DATA =

include_HEADERS = \
  GenerateWarnmap.h

noinst_HEADERS = \
  GenerateWarnmapLinkDef.h

libdgpp_warnmap_la_SOURCES = \
  GenerateWarnmap.cc \
  GenerateWarnmap_Dict.C

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libdgpp_warnmap_la_LIBADD = \
  -lSubsysReco

#BUILT_SOURCES = testexternals.C
#noinst_PROGRAMS = testexternals
#testexternals_LDADD = libdgpp_warnmap.la
#testexternals.C:
#	echo "//*** this is a generated file. Do not commit, do not edit" > $@
#	echo "int main() {return 0;}" >> $@

# Rule for generating CINT dictionaries from class headers.
GenerateWarnmap_Dict.C: \
  GenerateWarnmap.h \
  GenerateWarnmapLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

CLEANFILES = *Dict* testexternals*
