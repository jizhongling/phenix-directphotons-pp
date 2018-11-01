#include "BBCCounts.h"

void draw_ALL()
{
  DataBase *db = new DataBase();
  char pattern[10];
  int runnumber = 387550;
  db->GetSpinPattern(runnumber, pattern);
  cout << pattern << endl;
}
