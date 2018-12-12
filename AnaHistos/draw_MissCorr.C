#include "GlobalVars.h"
#include "QueryTree.h"

void draw_MissCorr()
{
  const double Conv[3] = {0.849, 0.959, 0.959};
  const double eConv[3] = {0.027, 0.023, 0.023};

  QueryTree *qt_mergecorr1 = new QueryTree("data/MergeCorr-1photon.root", "RECREATE");
  QueryTree *qt_mergecorr2 = new QueryTree("data/MergeCorr-2photon.root", "RECREATE");

  QueryTree *qt_merge1 = new QueryTree("data/Merge-1photon.root");
  QueryTree *qt_merge2 = new QueryTree("data/Merge-2photon.root");
  QueryTree *qt_badpass = new QueryTree("data/MergePassRate.root");

  for(int ipt=0; ipt<npT; ipt++)
    for(int part=0; part<3; part++)
    {
      double xpt, Merge2, eMerge2, Merge1, eMerge1, BadPass, eBadPass;
      qt_merge2->Query(ipt, part, xpt, Merge2, eMerge2);
      qt_merge1->Query(ipt, part, xpt, Merge1, eMerge1);
      qt_badpass->Query(ipt, part/2, xpt, BadPass, eBadPass);

      double MergeCorr1 = Merge1 * (1.-Conv[part]) / Conv[part];
      double eMergeCorr1 = MergeCorr1 * sqrt( pow(eMerge1/Merge1,2) + pow(eConv[part]/Conv[part]/(1.-Conv[part]),2) );
      double MergeCorr2 = Merge2 * BadPass;
      double eMergeCorr2 = MergeCorr2 * sqrt( pow(eMerge2/Merge2,2) + pow(eBadPass/BadPass,2) );

      qt_mergecorr1->Fill(ipt, part, xpt, MergeCorr1, eMergeCorr1);
      qt_mergecorr2->Fill(ipt, part, xpt, MergeCorr2, eMergeCorr2);
    }

  qt_mergecorr1->Write();
  qt_mergecorr2->Write();
}
