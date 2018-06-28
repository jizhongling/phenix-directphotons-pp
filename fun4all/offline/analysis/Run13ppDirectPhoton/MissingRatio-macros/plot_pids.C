int
plot_pids()
{
  // Settings
  gStyle->SetOptStat(0);

  // Open file
  TFile *fpisa = new TFile("histos/IsolationCut-MB.root", "OPEN");
  //TFile *fpisa = new TFile("histos/IsolationCut-test.root", "OPEN");

  // Get tree
  TTree *ttruth = (TTree*)fpisa->Get("mcparticles");

  TH1F* h1_pid = new TH1F("h_pid","",2001,-0.5,2000.5);
  h1_pid->GetXaxis()->SetTitle("PISA particle ID");
  h1_pid->GetYaxis()->SetTitle("# entries / #Sigma entries");
  ttruth->Draw("t_pid >> h_pid", "t_pt > 2");
  h1_pid->Scale(1./h1_pid->Integral());

  TCanvas *c1 = new TCanvas();
  h1_pid->GetXaxis()->SetRangeUser(0,35);
  h1_pid->GetYaxis()->SetRangeUser(0,0.17);
  h1_pid->DrawClone();
  c1->Print("pisa_pid_1.eps");

  TCanvas *c2 = new TCanvas();
  h1_pid->GetXaxis()->SetRangeUser(1060,1095);
  h1_pid->GetYaxis()->SetRangeUser(0,0.17);
  h1_pid->DrawClone();
  c2->Print("pisa_pid_2.eps");

  // return
  return 0;
}
