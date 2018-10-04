AUTOMAKE_OPTIONS = foreign 

bin_PROGRAMS = getBBCCount getBBCCount-allruns

getBBCCount_SOURCES =  getBBCCount.C

getBBCCount_allruns_SOURCES =  getBBCCount-allruns.C

INCLUDES  =  -I$(includedir) -I/opt/phenix/core/include -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

EXLIBS = -L/opt/phenix/core/lib -lodbc++ -L$(OFFLINE_MAIN)/lib -lphool $(ROOTLIBS)

LDADD =  $(EXLIBS)
