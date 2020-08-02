#include "GlobalVars.h"
#include "QueryTree.h"
#include "DataBase.h"

void anaRawALL(const int process = 0)
{
  const int ngroup = 2;
  const int nThread = 50;
  int thread = -1;
  int runnumber;
  int spin_pattern = 0;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist-ALL.txt");

  QueryTree *qt_asym = new QueryTree(Form("histos/raw-asym-%d.root",process), "RECREATE");

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    if(runnumber == 0) { spin_pattern++; continue; }
    thread++;
    if( thread < process*nThread )
      continue;
    else if( thread >= (process+1)*nThread )
      break;
    if(spin_pattern > 3) continue;

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/15811/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    // h[icr][ipol]
    TH1 *h_1photon[2][3];

    int evtype = 2;
    int bbc10cm = 0;
    int ival = 1;

    TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
    h_1photon_t = (TH1*)h_1photon_t->Clone();
    h_1photon_t->Reset();
    for(int icr=0; icr<2; icr++)
      for(int ipol=0; ipol<3; ipol+=2)
      {
        h_1photon[icr][ipol] = (TH1*)h_1photon_t->Clone(Form("h_1photon_cross%d_pol%d",icr,ipol));
        for(int sector=0; sector<8; sector++)
          for(int isolated=0; isolated<2; isolated++)
          {
            int ih = sector + 8*icr + 8*2*ipol + 8*2*3*isolated + 8*2*3*2*evtype + 8*2*3*2*3*bbc10cm + 8*2*3*2*3*2*ival + 8*2*3*2*3*2*4*(tof-1);
            TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
            h_1photon[icr][ipol]->Add(h_tmp);
            delete h_tmp;
          }
      }
    for(int icr=0; icr<2; icr++)
      for(int ipol=0; ipol<3; ipol+=2)
        h_1photon[icr][ipol]->Rebin(2);

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
        double npp = h_1photon[icr][2]->GetBinContent(ipt+1);
        double npm = h_1photon[icr][0]->GetBinContent(ipt+1);
        if( npp+npm < 10. || npp <= 1. || npm <= 1. ) continue;
        double ALL = 1./(pb[0]*py[0]) * (npp - r*npm) / (npp + r*npm);
        double eALL = sqrt( pow(2.*r*npp*npm/pb[0]/py[0],2) / pow(npp+r*npm,4) * (1./npp+1./npm+er*er/r/r)
            + (pow(pb[1]/pb[0],2)+pow(py[1]/py[0],2)) * ALL*ALL );

        if( TMath::Finite(ALL+eALL) && eALL > 0. )
        {
          int ig = ipt + npT/ngroup*icr + 2*npT/ngroup*spin_pattern;
          qt_asym->Fill(runnumber, ig, runnumber, ALL, eALL);
        }
      } // ipt
    } // icr

    delete f;
  }

  qt_asym->Save();
}
