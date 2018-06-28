void draw_Disp()
{
  TFile *f = new TFile("data/MissingRatio-histo.root");
  THnSparse *hn_dispProb = (THnSparse*)f->Get("hn_dispProb_photon");

  TCanvas *c = new TCanvas("c", "Canvas", 3600, 3000);
  gStyle->SetOptStat(0);
  c->Divide(6,5);

  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};
  
  TAxis *axis_i = hn_dispProb->GetAxis(3);
  const char *title_i = axis_i->GetTitle();

  for(int part=0; part<3; part++)
  {
    hn_dispProb->GetAxis(1)->SetRange(secl[part],sech[part]);
    int ipad = 1;
    for(int i=0; i<60; i+=2)
    {
      c->cd(ipad++);
      int bin_i = i + 1;
      axis_i->SetRange(bin_i,bin_i);
      TH1 *h_disp = (TH1*)hn_dispProb->Projection(4)->Clone();
      float low_i = axis_i->GetBinLowEdge(bin_i);
      float high_i = low_i + axis_i->GetBinWidth(bin_i);
      h_disp->SetTitle(Form("%s: %3.1f-%3.1f GeV",title_i,low_i,high_i));
      h_disp->GetXaxis()->SetRangeUser(0.1, 1.);
      h_disp->DrawCopy();
    }
    c->Print(Form("plots/Disp-disp-photon-%d.pdf",part));
    c->Clear("D");
  }
}
