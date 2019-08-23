#ifndef __HISTOGRAMBOOKER_H__
#define __HISTOGRAMBOOKER_H__

#include <string>

class Fun4AllHistoManager;

/* Root classes */
class TH1;
class TH2;
class THnSparse;

class HistogramBooker
{
  public:

    /**
     * Default constructor
     */
    HistogramBooker();

    /**
     * Default destructor
     */
    virtual ~HistogramBooker();


    /**
     *  create histograms
     */
    static Fun4AllHistoManager * GetHistoManager( std::string managername );


  protected:

};
#endif /* __HISTOGRAMBOOKER_H__ */
