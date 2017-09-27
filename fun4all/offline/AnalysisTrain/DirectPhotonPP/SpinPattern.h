#ifndef __SPINPATTERN_H__
#define __SPINPATTERN_H__

#include <PHObject.h>

#include <iostream>

#define CHECK_BUNCH_GET(name) \
  do { \
    if(bunch<0 || bunch>119) \
    { \
      std::cerr << "Wrong bunch number: " << bunch << std::endl; \
      return -9999; \
    } \
    return name[bunch]; \
  } while(0)

#define CHECK_BUNCH_SET(name) \
  do { \
    if(bunch<0 || bunch>119) \
    { \
      std::cerr << "Wrong bunch number: " << bunch << std::endl; \
      return; \
    } \
    name[bunch] = a_##name; \
  } while(0)

class SpinPattern: public PHObject
{
  public:
    SpinPattern() { Reset(); }
    virtual ~SpinPattern() {}

    void Reset();

    int get_runnumber() const { return runnumber; }
    int get_qa_level() const { return qa_level; }
    int get_fillnumber() const { return fillnumber; }
    int get_badrunqa() const { return badrunqa; }
    int get_crossing_shift() const { return crossing_shift; }

    int get_badbunch(int bunch) const { CHECK_BUNCH_GET(badbunch); }
    int get_spinpattern_blue(int bunch) const { CHECK_BUNCH_GET(spinpattern_blue); }
    int get_spinpattern_yellow(int bunch) const { CHECK_BUNCH_GET(spinpattern_yellow); }

    float get_pb() const { return pb; }
    float get_pbstat() const { return pbstat; }
    float get_pbsyst() const { return pbsyst; }
    float get_py() const { return py; }
    float get_pystat() const { return pystat; }
    float get_pysyst() const { return pysyst; }

    long long get_bbc_narrow(int bunch) const { CHECK_BUNCH_GET(bbc_narrow); }
    long long get_bbc_wide(int bunch) const { CHECK_BUNCH_GET(bbc_wide); }
    long long get_zdc_narrow(int bunch) const { CHECK_BUNCH_GET(zdc_narrow); }
    long long get_zdc_wide(int bunch) const  { CHECK_BUNCH_GET(zdc_wide); }

    void set_runnumber(int a_runnumber) { runnumber = a_runnumber; }
    void set_qa_level(int a_qa_level) { qa_level = a_qa_level; }
    void set_fillnumber(int a_fillnumber) { fillnumber = a_fillnumber; }
    void set_badrunqa(int a_badrunqa) { badrunqa = a_badrunqa; }
    void set_crossing_shift(int a_crossing_shift) { crossing_shift = a_crossing_shift; }

    void set_badbunch(int bunch, int a_badbunch) { CHECK_BUNCH_SET(badbunch); }
    void set_spinpattern_blue(int bunch, int a_spinpattern_blue) { CHECK_BUNCH_SET(spinpattern_blue); }
    void set_spinpattern_yellow(int bunch, int a_spinpattern_yellow) { CHECK_BUNCH_SET(spinpattern_yellow); }

    void set_pb(float a_pb) { pb = a_pb; }
    void set_pbstat(float a_pbstat) { pbstat = a_pbstat; }
    void set_pbsyst(float a_pbsyst) { pbsyst = a_pbsyst; }
    void set_py(float a_py) { py = a_py; }
    void set_pystat(float a_pystat) { pystat = a_pystat; }
    void set_pysyst(float a_pysyst) { pysyst = a_pysyst; }

    void set_bbc_narrow(int bunch, long long a_bbc_narrow) { CHECK_BUNCH_SET(bbc_narrow); }
    void set_bbc_wide(int bunch, long long a_bbc_wide) { CHECK_BUNCH_SET(bbc_wide); }
    void set_zdc_narrow(int bunch, long long a_zdc_narrow) { CHECK_BUNCH_SET(zdc_narrow); }
    void set_zdc_wide(int bunch, long long a_zdc_wide) { CHECK_BUNCH_SET(zdc_wide); }

  protected:
    int runnumber;
    int qa_level;
    int fillnumber;
    int badrunqa;
    int crossing_shift;

    int badbunch[120];
    int spinpattern_blue[120];
    int spinpattern_yellow[120];

    float pb;
    float pbstat;
    float pbsyst;
    float py;
    float pystat;
    float pysyst;

    long long bbc_narrow[120];
    long long bbc_wide[120];
    long long zdc_narrow[120];
    long long zdc_wide[120];

    ClassDef(SpinPattern, 1)
};

#endif /* __SPINPATTERN_H__ */
