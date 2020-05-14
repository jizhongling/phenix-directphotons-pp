#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "SumSysErr.h"

void draw_SysErrEn()
{
  const char *sysname[4] = {"Sum", "Global scale", "Non-lin", "Geom"};

  QueryTree *qt_sys = new QueryTree("data/syserr-en-fast.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo-syserr.root");
  THnSparse *hn_photon = (THnSparse*)f->Get("hn_photon");
  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");

  const int s1 = 4;
  const int s2 = 1;
  for(int ipt=0; ipt<npT;ipt+=(ipt<20?s1:s2))
  {
    double ndir[2][4], endir[2][4];  // photon|pion, isys
    double rsys[2][3], ersys[2][3];  // photon|pion, isys-1

    for(int isys=0; isys<4; isys++)
    {
      int sysCond = isys==2 ? isys+2 : isys;
      hn_photon->GetAxis(3)->SetRange(sysCond+1,sysCond+1);
      //hn_photon->GetAxis(2)->SetRange(1,6);
      TH1 *h_pt = hn_photon->Projection(1);  // pt_reco
      ndir[0][isys] = h_pt->IntegralAndError(ipt+1,ipt+(ipt<20?s1:s2), endir[0][isys]);
      delete h_pt;

      TH1 *h_minv;
      hn_pion->GetAxis(5)->SetRange(sysCond+1,sysCond+1);
      //hn_pion->GetAxis(3)->SetRange(1,6);
      hn_pion->GetAxis(1)->SetRange(ipt+1,ipt+(ipt<20?s1:s2));  // pt_reco
      h_minv = (TH1*)hn_pion->Projection(2);  // minv
      double hpt = ipt<20 ? 0. : 0.01;
      FitMinv(h_minv, ndir[1][isys], endir[1][isys], false, 0.11-hpt,0.16+hpt);
      delete h_minv;

      if(isys)
        for(int itype=0; itype<2; itype++)
        {
          int index = itype + 2*isys;
          double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
          rsys[itype][isys-1] = ndir[itype][isys]/ndir[itype][0];
          ersys[itype][isys-1] = rsys[itype][isys-1]*sqrt(pow(endir[itype][isys]/ndir[itype][isys],2) + pow(endir[itype][0]/ndir[itype][0],2));
          if( TMath::Finite(rsys[itype][isys-1]+ersys[itype][isys-1]) )
            qt_sys->Fill(ipt, index, xpt, fabs(rsys[itype][isys-1]-1), ersys[itype][isys-1]);
        } // itype
    } // isys

    for(int itype=0; itype<2; itype++)
    {
      double sum, esum;
      SumSysErr(3, rsys[itype], ersys[itype], sum, esum);
      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      if( TMath::Finite(sum+esum) )
        qt_sys->Fill(ipt, itype, xpt, sum, esum);
    } // itype
  } // ipt

  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(4);
  qt_sys->Write();
  for(int itype=0; itype<2; itype++)
  {
    mc();
    mcd();
    for(int isys=0; isys<4; isys++)
    {
      int index = itype + 2*isys;
      TGraphErrors *gr = qt_sys->Graph(index);
      aset(gr, "p_{T} [GeV]","SysErr", 0.,30., 0.,0.15);
      style(gr, 20+isys, 1+isys);
      char *opt = isys ? "P" : "AL";
      gr->Draw(opt);
      if(itype == 0)
        leg0->AddEntry(gr, sysname[isys], "P");
    }
    leg0->Draw();
    mcw( 0, Form("sysen-pion%d",itype) );
    char *type = itype ? "pion" : "photon";
    c0->Print(Form("plots/SysErrEn-%s.pdf",type));
  }
  qt_sys->Close();
}
