AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libdgpp_warnmap.la

noinst_HEADERS = \
  GenerateWarnmapLinkDef.h

include_HEADERS = \
  GenerateWarnmap.h

libdgpp_warnmap_la_SOURCES = \
  GenerateWarnmap_Dict.C \
  GenerateWarnmap.cc

libdgpp_warnmap_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco


# Rule for generating CINT dictionaries from class headers.
GenerateWarnmap_Dict.C: \
  GenerateWarnmap.h \
  GenerateWarnmapLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

