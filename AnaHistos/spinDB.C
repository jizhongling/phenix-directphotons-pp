#include <SpinDBOutput.hh>
#include <SpinDBContent.hh>

#include <iostream>
#include <fstream>

using namespace std;

int main()
{
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist.txt");
  int nrun = 0;
  int runnumber[1024];
  while(fin >> runnumber[nrun]) nrun++;
  fin.close();

  SpinDBOutput spin_out;
  SpinDBContent spin_cont;

  // initialize opbject to access spin DB
  spin_out.Initialize();
  spin_out.SetUserName("phnxrc");
  spin_out.SetTableName("spin");

  long long bbc15_total = 0;
  long long bbc30_total = 0;

  for(int ir=0; ir<nrun; ir++)
  {
    // retrieve entry from Spin DB
    int qa_level = spin_out.GetDefaultQA(runnumber[ir]);
    spin_out.StoreDBContent(runnumber[ir], runnumber[ir], qa_level);
    spin_out.GetDBContentStore(spin_cont, runnumber[ir]);

    long long bbc15 = 0; 
    long long bbc30 = 0;
    for(int i=0; i<120; i++)
    {
      bbc15 += spin_cont.GetScalerBbcVertexCut(i);
      bbc30 += spin_cont.GetScalerBbcNoCut(i);
    }

    bbc15_total += bbc15;
    bbc30_total += bbc30;
  }

  cout << "BBC15cm: " << bbc15_total << endl;
  cout << "BBC30cm: " << bbc30_total << endl;

  return 0;
}
