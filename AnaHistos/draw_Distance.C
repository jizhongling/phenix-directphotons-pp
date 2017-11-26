void Draw(Int_t fn)
{
  if(fn==0)
    TFile *f = new TFile("data/MissingRatio-histo.root");
  else if(fn==1)
    TFile *f = new TFile("data/AnaPHPythia-histo.root");

  THnSparse *hn_distance = (THnSparse*)f->Get("hn_distance");

  TCanvas *c = new TCanvas("c", "Canvas", 3600, 3000);
  gStyle->SetOptStat(0);
  c->Divide(6,5);

  Int_t secl[3] = {1, 5, 7};
  Int_t sech[3] = {4, 6, 8};

  Double_t scale = hn_distance->Projection(0)->GetSumOfWeights();

  for(Int_t part=0; part<3; part++)
  {
    hn_distance->GetAxis(1)->SetRange(secl[part],sech[part]);
    Int_t ipad = 1;
    for(Int_t ipt=0; ipt<30; ipt++)
    {
      c->cd(ipad++);
      hn_distance->GetAxis(0)->SetRange(ipt+1,ipt+1);
      TH1 *h_distance = (TH1*)hn_distance->Projection(2);
      h_distance->Scale(1./scale);
      Float_t pTlow = hn_distance->GetAxis(0)->GetBinLowEdge(ipt+1);
      Float_t pThigh = pTlow + hn_distance->GetAxis(0)->GetBinWidth(ipt+1);
      h_distance->SetTitle(Form("p_{T}: %3.1f-%3.1f GeV",pTlow,pThigh));
      h_distance->GetXaxis()->SetRangeUser(0.,15.);
      h_distance->SetMarkerStyle(20+part);
      h_distance->SetMarkerColor(1+part);
      if(part == 0)
        h_distance->DrawCopy("P");
      else
        h_distance->DrawCopy("PSAME");
    }
  }

  if(fn==0)
    c->Print("plots/Distance-PISA.pdf");
  else if(fn==1)
    c->Print("plots/Distance-FastMC.pdf");

  return;
}

void draw_Distance()
{
  for(Int_t fn=0; fn<2; fn++)
    Draw(fn);
}
