plot_tof()
{
  gStyle->SetOptStat(0);

  //TFile *fin = new TFile("taxi_test/keep/WarnmapData-CNT_MB_run13pp_510GeV_pro97-0000393483-9000.root","OPEN");
  //TFile *fin = new TFile("taxi_test/keep/WarnmapData-CNT_ERT_run13pp_510GeV_pro97-0000393483-9000.root","OPEN");

  TFile *fin = new TFile("/gpfs/mnt/gpfs02/phenix/spin3/nfeege/taxi/Run13pp510ERT/8058/data/DirectPhotonPP-387565.root","OPEN");

  TH2D* h_tof = (TH2D*) fin->Get("hitmap_tof");
  TH2D* h_tof_sector = (TH2D*) fin->Get("hitmap_tof_sector");

  h_tof->SetTitle("");
  h_tof->GetXaxis()->SetTitle("tower ID");
  h_tof->GetYaxis()->SetTitle("TOF [ns]");

  h_tof_sector->SetTitle("");
  h_tof_sector->GetXaxis()->SetTitle("sector");
  h_tof_sector->GetYaxis()->SetTitle("TOF [ns]");

  TCanvas *c1 = new TCanvas();
  c1->SetLogz(1);
  c1->SetRightMargin(0.15);
  h_tof->Draw("colz");
  c1->Print("plots_tof/tof_vs_towerid_example.png");

  TCanvas *c2 = new TCanvas();
  c2->SetLogz(1);
  c2->SetRightMargin(0.15);
  h_tof_sector->Draw("colz");
  c2->Print("plots_tof/tof_vs_sector_example.png");

  /* Project TOF distribution: single tower */
  int bin_tower = 100;
  TH1D* proj_tof_1tower = h_tof->ProjectionY("proj_tof_1tower",bin_tower,bin_tower);
  proj_tof_1tower->GetYaxis()->SetTitle("# entries");
  proj_tof_1tower->Rebin(1);

  TCanvas *c3 = new TCanvas();
  proj_tof_1tower->Draw("");
  c3->Print("plots_tof/tof_1tower_example.png");

  /* Project TOF distribution: single sector */
  int bin_sector = 1;
  TH1D* proj_tof_1sector = h_tof_sector->ProjectionY("proj_tof_1sector",bin_sector,bin_sector);
  proj_tof_1sector->GetYaxis()->SetTitle("# entries");
  proj_tof_1sector->Rebin(1);

  TCanvas *c4 = new TCanvas();
  proj_tof_1sector->Draw("");
  c4->Print("plots_tof/tof_1sector_example.png");

  cout << proj_tof_1tower->GetEntries() << " - " << proj_tof_1sector->GetEntries() << endl;
}
