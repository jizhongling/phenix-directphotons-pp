#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main()
{
  ifstream ftot("/phenix/plhf/zji/taxi/Run13pp510ERT/runnumber.txt");
  ifstream fexl("/phenix/plhf/zji/taxi/Run13pp510ERT/runexclude.txt");
  ofstream flist("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist.txt");

  string runno;
  string exlno;

  while( ftot >> runno )
  {
    bool islist = true;
    fexl.clear();
    fexl.seekg(0);

    while( fexl >> exlno )
    {
      size_t found = exlno.find("-");
      if( found == string::npos )
      {
        if( runno.compare(exlno) == 0 )
        {
          islist = false;
          break;
        }
        else
        {
          continue;
        }
      }
      else
      {
        if( runno.compare(0,6, exlno,0,6) >= 0 && runno.compare(0,6, exlno,found+1,6) <= 0 )
        {
          islist = false;
          break;
        }
        else
        {
          continue;
        }
      }
    }

    if(islist)
      flist << runno << " ";
  }

  ftot.close();
  fexl.close();
  flist.close();

  return 0;
}
