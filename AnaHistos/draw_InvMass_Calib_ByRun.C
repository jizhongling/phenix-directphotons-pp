void draw_InvMass_Calib_ByRun() {

  const int sectors = 8;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist.txt");
  int nrun = 0;
  int runnumber[1024];
  while(fin >> runnumber[nrun]) nrun++;
  fin.close();

  TGraphErrors* g_minv_runnumber[sectors];
  for(int is=0; is<sectors; is++)
    g_minv_runnumber[is] = new TGraphErrors(nrun);

  for(int ir=0; ir<nrun; ir++) {

    char buf[100];
    sprintf(buf, "/phenix/plhf/zji/taxi/Run13pp510ERT/10853/data/DirectPhotonPP-%d.root", runnumber[ir]);
    TFile* f = new TFile(buf);
    TH3* h3_minv = (TH3*)f->Get("h3_inv_mass_pi0calib_raw");
    //TH3* h3_minv = (TH3*)f->Get("h3_inv_mass_pi0calib");

    for(int is=0; is<sectors; is++) {
      TH1* h_minv = (TH1*)h3_minv->ProjectionZ("h_minv",is+1,is+1);
      h_minv->GetXaxis()->SetRangeUser(0.1,0.2);
      double max = h_minv->GetMaximum();
      int bin1 = h_minv->FindFirstBinAbove(max/2);
      int bin2 = h_minv->FindLastBinAbove(max/2);
      double minv_peak = h_minv->GetBinCenter((bin1+bin2)/2);
      double fwhm = h_minv->GetBinCenter(bin2) - h_minv->GetBinCenter(bin1);
      //cout << "Max=" << max << " Left bin=" << bin1 << " Right bin=" << bin2 << " Tof peak=" << minv_peak << " FWHM=" << fwhm << endl;
      g_minv_runnumber[is]->SetPoint(ir, (double)runnumber[ir], minv_peak);
      g_minv_runnumber[is]->SetPointError(ir, 0., fwhm/2);
      delete h_minv;
    }

    delete f;

  }

  mc(0, 2,4);
  gStyle->SetOptStat(1);

  for(int is=0; is<sectors; is++) {
    mcd(0, is+1);
    char buf[100];
    sprintf(buf, "Sector %d", is);
    aset(g_minv_runnumber[is]);
    g_minv_runnumber[is]->SetTitle(buf);
    g_minv_runnumber[is]->GetXaxis()->SetTitle("Runnumber");
    g_minv_runnumber[is]->GetXaxis()->SetLimits(387000.,399000.);
    g_minv_runnumber[is]->GetYaxis()->SetTitle("ToF peak [ns]");
    g_minv_runnumber[is]->GetYaxis()->SetRangeUser(0.1,0.2);
    g_minv_runnumber[is]->SetMarkerColor(1);
    g_minv_runnumber[is]->SetMarkerSize(1);
    g_minv_runnumber[is]->SetMarkerStyle(21);
    g_minv_runnumber[is]->Draw("AP");
  }

  c0->Print("plots/InvMass_NoCalib_ByRun.pdf");
  //c0->Print("plots/InvMass_Calib_ByRun.pdf");

}
