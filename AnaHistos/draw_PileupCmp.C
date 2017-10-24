void draw_PileupCmp()
{
  TFile *f_mine = new TFile("Pileup.root");
  TFile *f_sasha = new TFile("Pileup-Sasha.root");

  mc();
  mcd();

  TGraph *gr_ratio[2];
  TGraphErrors *gr_mine[2];
  TGraphErrors *gr_sasha[2];
  for(Int_t igr=0; igr<2; igr++)
  {
    gr_ratio[igr] = new TGraph(1000);
    gr_mine[igr] = (TGraphErrors*)f_mine->Get(Form("gr_run_%d",igr));
    gr_sasha[igr] = (TGraphErrors*)f_sasha->Get(Form("gr_%d",igr));

    Int_t gn_mine = gr_mine[igr]->GetN();
    Int_t gn_sasha = gr_sasha[igr]->GetN();
    Double_t *runno_mine = gr_mine[igr]->GetX();
    Double_t *npi0_mine = gr_mine[igr]->GetY();
    Double_t *runno_sasha = gr_sasha[igr]->GetX();
    Double_t *npi0_sasha = gr_sasha[igr]->GetY();

    Int_t ip = 0;
    for(Int_t i=0; i<gn_mine; i++)
      for(Int_t j=0; j<gn_sasha; j++)
        if( abs(runno_mine[i]-runno_sasha[j]) < 0.1 )
        {
          Double_t xx = runno_mine[i];
          Double_t yy = npi0_mine[i] / npi0_sasha[j];
          gr_ratio[igr]->SetPoint(ip, xx, yy);
          ip++;
        }

    aset(gr_ratio[igr], "Runnumber","Mine/Sasha", 386000.,400000., 0.5, 1.5);
    style(gr_ratio[igr], igr+20, igr+1);
  }
  gr_ratio[0]->Draw("AP");
  gr_ratio[1]->Draw("P");

  c0->Print("PileupCmp.pdf");
}
