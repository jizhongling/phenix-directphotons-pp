#include "SpinPattern.h"

ClassImp(SpinPattern)

void SpinPattern::Reset()
{
  runnumber = -9999;
  qa_level = -9999;
  fillnumber = -9999;
  badrunqa = -9999;
  crossing_shift = -9999;

  pb = -9999.;
  pbstat = -9999.;
  pbsyst = -9999.;
  py = -9999.;
  pystat = -9999.;
  pysyst = -9999.;

  for(int i=0; i<120; i++)
  {
    badbunch[i] = -9999;
    spinpattern_blue[i] = -9999;
    spinpattern_yellow[i] = -9999;

    bbc_narrow[i] = -9999;
    bbc_wide[i] = -9999;
    zdc_narrow[i] = -9999;
    zdc_wide[i] = -9999;
  }
}


