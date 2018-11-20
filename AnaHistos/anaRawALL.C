#include "GlobalVars.h"
#include "DataBase.h"

void anaRawALL(const int process = 0)
{
  const int nThread = 50;
  int thread = -1;
  int runnumber;
  int spin_pattern = 0;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist-ALL.txt");

  TGraphErrors *gr[4*2*npT];
  int igp[4*2*npT] = {};
  for(int pattern=0; pattern<4; pattern++)
    for(int icr=0; icr<2; icr++)
      for(int ipt=0; ipt<npT; ipt++)
      {
        int ig = ipt + npT*icr + 2*npT*pattern;
        gr[ig] = new TGraphErrors(nThread);
        gr[ig]->SetName(Form("gr_%d",ig));
      }

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    if(runnumber == 0) { spin_pattern++; continue; }
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/14194/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    // h[icr][ipol]
    TH1 *h_1photon[2][3];

    int evtype = 2;
    int bbc10cm = 0;
    int tof = 1;
    int prob = 1;
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
            int ih = sector + 8*icr + 2*8*ipol + 3*2*8*isolated + 2*3*2*8*tof + 2*2*3*2*8*prob + 2*2*2*3*2*8*evtype + 3*2*2*2*3*2*8*bbc10cm + 2*3*2*2*2*3*2*8*ival;
            TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
            h_1photon[icr][ipol]->Add(h_tmp);
            delete h_tmp;
          }
      }

    double pb[2], py[2];  // pb, pbstat
    double rlum[2], erlum[2];  // even/odd crossings
    db->GetPol(runnumber, pb, py);
    db->GetRelLum(runnumber, rlum, erlum);

    for(int icr=0; icr<2; icr++)
    {
      double r = rlum[icr];
      double er = erlum[icr];

      for(int ipt=0; ipt<npT; ipt++)
      {
        double npp = h_1photon[icr][2]->GetBinContent(ipt+1);
        double npm = h_1photon[icr][0]->GetBinContent(ipt+1);
        if( npp+npm < 10. || npp <= 1. || npm <= 1. ) continue;
        double ALL = 1./(pb[0]*py[0]) * (npp - r*npm) / (npp + r*npm);
        double eALL = sqrt( pow(2.*r*npp*npm/pb[0]/py[0],2) / pow(npp+r*npm,4) * (1./npp+1./npm+er*er/r/r)
            + (pow(pb[1]/pb[0],2)+pow(py[1]/py[0],2)) * ALL*ALL );

        if( TMath::Finite(ALL+eALL) && eALL > 0. )
        {
          int ig = ipt + npT*icr + 2*npT*spin_pattern;
          gr[ig]->SetPoint(igp[ig], (double)runnumber, ALL);
          gr[ig]->SetPointError(igp[ig], 0., eALL);
          igp[ig]++;
        }
      }
    }

    delete f;
  }

  TFile *f_out = new TFile(Form("histos/raw-asym-%d.root",process), "RECREATE");
  for(pattern=0; pattern<4; pattern++)
    for(int icr=0; icr<2; icr++)
      for(int ipt=0; ipt<npT; ipt++)
      {
        int ig = ipt + npT*icr + 2*npT*pattern;
        gr[ig]->Set(igp[ig]);
        gr[ig]->Write();
      }
  f_out->Close();
}
