#include "GlobalVars.h"
#include "QueryTree.h"
#include "DataBase.h"

void anaPi0ALL(const int process = 0)
{
  const int ngroup = 2;
  const int nThread = 50;
  int thread = -1;
  int runnumber;
  int spin_pattern = 0;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist-ALL.txt");

  QueryTree *qt_asym = new QueryTree(Form("histos/pion-asym-%d.root",process), "RECREATE");

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    if(runnumber == 0) { spin_pattern++; continue; }
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;
    if(spin_pattern > 3) continue;

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/15354/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    // h[evtype][icr][ipol]
    TH2 *h2_pion[3][2][3];

    int bbc10cm = 0;
    int tof = 1;
    int prob = 1;
    int ival = 1;

    TH2 *h2_pion_t = (TH2*)f->Get("h2_pion_0");
    h2_pion_t = (TH2*)h2_pion_t->Clone();
    h2_pion_t->Reset();
    for(int evtype=1; evtype<3; evtype++)
      for(int icr=0; icr<2; icr++)
        for(int ipol=0; ipol<3; ipol+=2)
        {
          h2_pion[evtype][icr][ipol] = (TH2*)h2_pion_t->Clone(Form("h2_pion_type%d_cross%d_pol%d",evtype,icr,ipol));
          for(int sector=0; sector<8; sector++)
            for(int isolated=0; isolated<2; isolated++)
            {
              int ih = sector + 8*icr + 8*2*ipol + 8*2*3*isolated + 8*2*3*2*tof + 8*2*3*2*3*prob + 8*2*3*2*3*2*evtype + 8*2*3*2*3*2*3*bbc10cm + 8*2*3*2*3*2*3*2*ival;
              TH2 *h2_tmp = (TH2*)f->Get(Form("h2_pion_%d",ih));
              h2_pion[evtype][icr][ipol]->Add(h2_tmp);
              delete h2_tmp;
            }
        }
    for(int evtype=1; evtype<3; evtype++)
      for(int icr=0; icr<2; icr++)
        for(int ipol=0; ipol<3; ipol+=2)
          h2_pion[evtype][icr][ipol]->RebinX(2);

    double pb[2], py[2];  // pb, pbstat
    double rlum[2], erlum[2];  // even/odd crossings
    db->GetPol(runnumber, pb, py);
    db->GetRelLum(runnumber, rlum, erlum);

    for(int icr=0; icr<2; icr++)
    {
      double r = rlum[icr];
      double er = erlum[icr];

      for(int ipt=0; ipt<npT/ngroup; ipt++)
      {
        int evtype = 2;
        if(ipt < 22/ngroup)  // <14GeV use ERT_4x4c
          evtype = 2;
        else  // >14GeV use ERT_4x4b
          evtype = 1;

        for(int ibg=0; ibg<2; ibg++)
        {
          double npp, npm;
          if(ibg == 0)
          {
            npp = h2_pion[evtype][icr][2]->Integral(ipt+1,ipt+1, 113,162);
            npm = h2_pion[evtype][icr][0]->Integral(ipt+1,ipt+1, 113,162);
          }
          else
          {
            npp = h2_pion[evtype][icr][2]->Integral(ipt+1,ipt+1, 48,97) +
              h2_pion[evtype][icr][2]->Integral(ipt+1,ipt+1, 178,227);
            npm = h2_pion[evtype][icr][0]->Integral(ipt+1,ipt+1, 48,97) +
              h2_pion[evtype][icr][0]->Integral(ipt+1,ipt+1, 178,227);
          }
          if( npp+npm < 10. || npp <= 1. || npm <= 1. ) continue;
          double ALL = 1./(pb[0]*py[0]) * (npp - r*npm) / (npp + r*npm);
          double eALL = sqrt( pow(2.*r*npp*npm/pb[0]/py[0],2) / pow(npp+r*npm,4) * (1./npp+1./npm+er*er/r/r)
              + (pow(pb[1]/pb[0],2)+pow(py[1]/py[0],2)) * ALL*ALL );

          if( TMath::Finite(ALL+eALL) && eALL > 0. )
          {
            int ig = ipt + npT/ngroup*icr + 2*npT/ngroup*spin_pattern + 4*2*npT/ngroup*ibg;
            qt_asym->Fill(runnumber, ig, runnumber, ALL, eALL);
          }
        } // ibg
      } // ipt
    } // icr

    delete f;
  }

  qt_asym->Save();
}
