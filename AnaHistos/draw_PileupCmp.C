void draw_PileupCmp()
{
  Int_t gn_mine[2] = {}, gn_sasha[2] = {};
  Double_t runno_mine[2][1000] = {}, runno_sasha[2][1000] = {};
  Double_t npi0_mine[2][1000] = {}, npi0_sasha[2][1000] = {};

  for(Int_t i=0; i<44; i++)
  {
    TFile *f_mine = new TFile(Form("pileup/Pileup-%d.root",i));
    if( f_mine->IsZombie() ) continue;

    for(Int_t igr=0; igr<2; igr++)
    {
      TGraphErrors *gr = (TGraphErrors*)f_mine->Get(Form("gr_run_%d",4+igr));
      for(Int_t ip=0; ip<gr->GetN(); ip++)
      {
        gr->GetPoint(ip, runno_mine[igr][gn_mine[igr]], npi0_mine[igr][gn_mine[igr]]);
        gn_mine[igr]++;
      }
    }

    delete f_mine;
  }

  TFile *f_sasha = new TFile("Pileup-Sasha-CVS.root");
  TFile *f_sasha0 = new TFile("Pileup-Sasha-TAXI.root");

  for(Int_t igr=0; igr<2; igr++)
  {
    TGraphErrors *gr = (TGraphErrors*)f_sasha->Get(Form("gr_run_%d",2+igr));
    TGraphErrors *gr0 = (TGraphErrors*)f_sasha0->Get(Form("gr_run_%d",2+igr));
    for(Int_t ip=0; ip<gr->GetN(); ip++)
    {
      gr->GetPoint(ip, runno_sasha[igr][gn_sasha[igr]], npi0_sasha[igr][gn_sasha[igr]]);
      gn_sasha[igr]++;
      //gr0->GetPoint(ip, runno_mine[igr][gn_mine[igr]], npi0_mine[igr][gn_mine[igr]]);
      //gn_mine[igr]++;
    }
  }

  delete f_sasha;
  delete f_sasha0;

  for(Int_t igr=0; igr<2; igr++)
  {
    for(Int_t j=0; j<gn_sasha[igr]; j++)
    {
      bool matched = false;
      for(Int_t i=0; i<gn_mine[igr]; i++)
        if( abs(runno_mine[igr][i]-runno_sasha[igr][j]) < 0.1 )
          matched = true;
      if(!matched)
        cout << "Part " << igr << " Runnumber = " << runno_sasha[igr][j] << " not matched!!!" << endl;
    }
  }

  mc();
  mcd();

  TGraph *gr_ratio[2];
  for(Int_t igr=0; igr<2; igr++)
  {
    gr_ratio[igr] = new TGraph(1000);
    Int_t igp1 = 0;
    for(Int_t i=0; i<gn_mine[igr]; i++)
      for(Int_t j=0; j<gn_sasha[igr]; j++)
        if( abs(runno_mine[igr][i]-runno_sasha[igr][j]) < 0.1 )
        {
          Double_t xx = runno_mine[igr][i];
          Double_t yy = npi0_mine[igr][i] / npi0_sasha[igr][j];
          gr_ratio[igr]->SetPoint(igp1, xx, yy);
          igp1++;
        }
    gr_ratio[igr]->Set(igp1);
    aset(gr_ratio[igr], "Runnumber","Mine/Sasha", 386000.,400000., 0.9, 1.1);
    style(gr_ratio[igr], igr+24, igr+1);
    cout << "gn_mine[" << igr << "] = " << gn_mine[igr] << endl
      << "gn_sasha[" << igr << "] = " << gn_sasha[igr] << endl
      << "gn_ratio[" << igr << "] = " << gr_ratio[igr]->GetN() << endl;
  }
  gr_ratio[0]->Draw("AP");
  gr_ratio[1]->Draw("P");

  c0->Print("PileupCmp.pdf");
}