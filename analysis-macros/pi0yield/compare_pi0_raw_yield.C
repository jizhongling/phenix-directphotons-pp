int
compare_pi0_raw_yield()
{
  /* get number of pi0 from runs extracted with PhotonNode + FillHistos code */
  TFile *f1 = new TFile("data/RawYield-Sasha.root", "OPEN");
  TTree *t1 = (TTree*)f1->Get("t_ert");
  t1->SetName("t1");
  t1->BuildIndex("runnumber");

  /* get number of pi0 from runs extracted with DirectPhotonPP code */
  TTree *t2 = new TTree();
  t2->ReadFile("data/raw_pi0_yield_Run13pp510ERT.txt", "runnumber/F:npi0/F");
  t2->SetName("t2");
  t2->BuildIndex("runnumber");

  /* add first tree as friend to second */
  t2->AddFriend("t1");

  TCanvas *c1 = new TCanvas();
  t2->Draw("t2.npi0 / ( t1.npion_tof[0] + t1.npion_tof[1] + t1.npion_tof[2]) : runnumber");

  TCanvas *c2 = new TCanvas();
  t2->Draw("t2.npi0 / ( t1.npion_tof[0] + t1.npion_tof[1] + t1.npion_tof[2]) : t2.npi0");

  return 0;
}
