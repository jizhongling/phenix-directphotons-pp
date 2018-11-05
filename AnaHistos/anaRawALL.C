#include "GlobalVars.h"
#include "DataBase.h"

void anaRawALL(const int process = 0)
{
  const int nThread = 10;
  int thread = -1;
  int runnumber;
  ifstream fin(Form("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist-%d.txt",process));

  TGraphErrors *gr[2];
  int igr[2] = {};
  for(int icr=0; icr<2; icr++)
    gr[icr] = new TGraphErrors(nThread);

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/13912/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    // h[pattern]
    TH1 *h_1photon[3];

    int evtype = 2;
    int bbc10cm = 0;
    int tof = 1;
    int prob = 1;
    int ival = 1;

    TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
    h_1photon_t = (TH1*)h_1photon_t->Clone();
    h_1photon_t->Reset();
    for(int pattern=0; pattern<3; pattern+=2)
    {
      h_1photon[pattern] = (TH1*)h_1photon_t->Clone(Form("h_1photon_pattern%d",pattern));
      for(int sector=0; sector<8; sector++)
        for(int isolated=0; isolated<2; isolated++)
        {
          int ih = sector + 8*pattern + 3*8*isolated + 2*3*8*tof + 2*2*3*8*prob + 2*2*2*3*8*evtype + 3*2*2*2*3*8*bbc10cm + 2*3*2*2*2*3*8*ival;
          TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
          h_1photon[pattern]->Add(h_tmp);
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
        double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
        double npp = h_1photon[2]->GetBinContent(ipt+1);
        double npm = h_1photon[0]->GetBinContent(ipt+1);
        double ALL = 1./(pb[0]*py[0]) * (npp - r*npm) / (npp + r*npm);
        double eALL = sqrt( pow(2.*r*npp*npm/pb[0]/py[0],2) / pow(npp+r*npm,4) * (1./npp+1./npm+er*er/r/r)
            + (pow(pb[1]/pb[0],2)+pow(py[1]/py[0],2)) * ALL*ALL );

        gr[icr]->SetPoint(igr[icr], xpt, ALL);
        gr[icr]->SetPointError(igr[icr], 0., eALL);
        igr[icr]++;
      }
    }

    delete f;
  }
}
