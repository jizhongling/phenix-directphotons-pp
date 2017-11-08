#include "Pileup.h"
#include "MultiGraph.h"

void draw_Pileup()
{
  const char *dname[2] = {"ERT", "MB"};
  const char *cname[4] = {"PbSc without ToF", "PbGl without ToF", "PbSc with ToF", "PbGl with ToF"};
  TMultiGraph *mg[npT*8];
  for(Int_t ipt=0; ipt<npT; ipt++)
    for(Int_t id=0; id<2; id++)
      for(Int_t ic=0; ic<2; ic++)
        for(Int_t is=0; is<2; is++)
        {
          Int_t ig = ipt*8+id*4+ic*2+is;
          mg[ig] = new TMultiGraph();
        }

  for(Int_t i=0; i<44; i++)
  {
    TFile *f = new TFile(Form("pileup/Pileup-%d.root",i));
    if( f->IsZombie() ) continue;

    for(Int_t ipt=0; ipt<npT; ipt++)
      for(Int_t id=0; id<2; id++)
        for(Int_t ic=0; ic<2; ic++)
          for(Int_t is=0; is<2; is++)
          {
            Int_t ig = ipt*8+id*4+ic*2+is;
            TGraphErrors *gr = (TGraphErrors*)f->Get(Form("gr_%d",ig));
            if( gr->GetN() > 0)
              mg[ig]->Add(gr);
          }
  }

  TGraphErrors *gr_ratio[8];
  Int_t igp[8] = {};
  for(Int_t id=0; id<2; id++)
    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<2; is++)
      {
        Int_t igr = id*4+ic*2+is;
        gr_ratio[igr] = new TGraphErrors(npT);
      }

  TF1 *fn_fit = new TF1("fn_fit", "pol1");
  //TF1 *fn_fit = new TF1("fn_fit", "-[0]*log(1-[1]*x)/([1]*x)");
  mc(0, 2,2);

  for(Int_t ipt=0; ipt<npT; ipt++)
    for(Int_t id=0; id<2; id++)
    {
      for(Int_t ic=0; ic<2; ic++)
        for(Int_t is=0; is<2; is++)
        {
          Int_t cond = ic*2+is;
          Int_t igr = id*4+cond;
          Int_t ig = ipt*8+igr;

          mcd(0, cond+1);
          mg[ig]->Draw("AP");  // must before GetXaxis()
          mg[ig]->SetTitle(cname[cond]);
          mg[ig]->GetXaxis()->SetTitle("Nmb/Nclock");
          mg[ig]->GetYaxis()->SetTitle("Npi0/Nevent");
          mg[ig]->GetXaxis()->SetLimits(0., 0.2);  // Do not use SetRangeUser()
          //mg[ig]->GetYaxis()->SetRangeUser(0., 3.2e-3);  // Do not use SetLimits()

          //fn_fit->SetParameters(1e-4, 1.);
          //for(Int_t ifit=0; ifit<5; ifit++)
          //{
          //  mg[ig]->Fit(fn_fit, "Q0");
          //  fn_fit->SetParameters( fn_fit->GetParameters() );
          //}
          mg[ig]->Fit(fn_fit, "Q");

          if( ipt > 0 )
          {
            Double_t scale = sqrt( fn_fit->GetChisquare() / fn_fit->GetNDF() );
            Double_t p0 = fn_fit->GetParameter(0);
            Double_t ep0 = fn_fit->GetParError(0) * scale;

            Double_t mean, emean;
            GetMeanError<TGraphErrors>(mg[ig], mean, emean);

            Double_t xx = ( pTbin[id][ipt-1] + pTbin[id][ipt] ) / 2.;
            Double_t yy = p0 / mean;
            Double_t eyy = yy * sqrt( pow(ep0/p0,2.) + pow(emean/mean,2.) );
            if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
            {
              gr_ratio[igr]->SetPoint(igp[igr], xx, yy);
              gr_ratio[igr]->SetPointError(igp[igr], 0., eyy);
              igp[igr]++;
            }
          }
        }

      c0->Print( Form("pileup/Pileup-data%d-pt%d-%d-pol1.pdf", id, pTlow[id][ipt], pThigh[id][ipt]) );
      c0->Clear("D");
    }

  mc(1, 2,2);

  for(Int_t id=0; id<2; id++)
    for(Int_t is=0; is<2; is++)
    {
      Int_t igr = id*4+2+is;
      gr_ratio[igr]->Set(igp[igr]);

      mcd(1, id*2+is+1);
      gr_ratio[igr]->SetTitle( Form("%s %s", dname[id], cname[2+is]) );
      aset(gr_ratio[igr], "pT [GeV]", "#frac{p0}{mean}", 1.,4.25);
      style(gr_ratio[igr], 20, kRed);
      gr_ratio[igr]->Draw("AP");
      gr_ratio[igr]->Fit("pol0", "Q","", 1.5,30.);
    }

  c1->Print("pileup/Pileup-ratio-pol1.pdf");
}
