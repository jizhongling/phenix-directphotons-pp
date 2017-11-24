void draw_SimCmpAsym()
{
  TFile *f_pisa = new TFile("data/MissingRatio-histo.root");
  TFile *f_fastmc = new TFile("data/AnaPHPythia-histo.root");

  THnSparse *hn_total_pisa = (THnSparse*)f_pisa->Get("hn_total");
  THnSparse *hn_total_fastmc = (THnSparse*)f_fastmc->Get("hn_total");

  mc(0, 6,5);

  Int_t ipad = 1;
  for(Int_t ipt=0; ipt<30; ipt++)
  {
    mcd(0, ipad++);

    TAxis *axis_pt = hn_total_fastmc->GetAxis(0);
    Double_t pt_low = axis_pt->GetBinLowEdge(ipt+1);
    Double_t pt_up = axis_pt->GetBinUpEdge(ipt+1);

    hn_total_fastmc->GetAxis(0)->SetRange(ipt+1,ipt+1);
    TH1 *h_asym_fastmc = hn_total_fastmc->Projection(3);
    Double_t ss = h_asym_fastmc->GetEntries();
    aset(h_asym_fastmc, "Asym", "Yield");
    h_asym_fastmc->SetTitle(Form("p_{T}: %3.1f-%3.1f GeV",pt_low,pt_up));
    h_asym_fastmc->GetXaxis()->SetRangeUser(0.,1.);
    h_asym_fastmc->GetYaxis()->SetRangeUser(0.,ss/70.);
    h_asym_fastmc->SetLineColor(kRed);
    h_asym_fastmc->DrawCopy();
    delete h_asym_fastmc;

    hn_total_pisa->GetAxis(0)->SetRange(ipt+1,ipt+1);
    TH1 *h_asym_pisa = hn_total_pisa->Projection(3);
    h_asym_pisa->SetLineColor(kBlue);
    h_asym_pisa->DrawCopy("SAME");
    delete h_asym_pisa;
  }

  c0->Print("SimCmpAsym.pdf");
}
