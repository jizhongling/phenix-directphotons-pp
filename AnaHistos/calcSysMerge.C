#include "GlobalVars.h"
#include "QueryTree.h"

void calcSysMerge()
{
  const double n = 7.22;

  QueryTree *qt_syserr = new QueryTree("data/syserr-merge.root", "RECREATE");

  double xpt[9], syspi0[2][9];
  int i = 0;
  ifstream fin("data/syserr-merge.txt");
  while(fin >> xpt[i] >> syspi0[0][i] >> syspi0[1][i] && ++i < 9);

  for(int i=0; i<9; i++)
  {
    int ipt = i + 21;
    double sysph[2] = {};
    for(int j=i; j<9; j++)
    {
      int jpt = j + 21;
      double weight;
      double factor = (n - 1)*(pow(pTbin[jpt],-n) - pow(pTbin[jpt+1],-n))/(pow(pTbin[ipt],-n+1) - pow(pTbin[ipt+1],-n+1));
      if(j == i)
        weight = n - pTbin[ipt]*factor;
      else
        weight = (pTbin[ipt+1] - pTbin[ipt])*factor;
      for(int part=0; part<2; part++)
        sysph[part] += weight*syspi0[part][j];
    } // j
    for(int part=0; part<2; part++)
    {
      qt_syserr->Fill(ipt, part, xpt[i], syspi0[part][i], 0.);
      qt_syserr->Fill(ipt, part+2, xpt[i], sysph[part], 0.);
    }
  } // i

  qt_syserr->Save();
}
