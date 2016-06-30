void draw_Warnmap()
{
  TFile *f = new TFile("Warnmap.root", "recreate");
  TTree *T1 = new TTree("T1", "Our warnmap");
  TTree *T2 = new TTree("T2", "Inseok warnmap");

  T1->ReadFile("/phenix/plhf/zji/install/share/DirectPhotonPP/Warnmap_Run13pp510.txt", "sector/I:y:z:status");
  T2->ReadFile("/phenix/plhf/zji/install/share/DirectPhotonPP/Warnmap_Run13pp510_ISYOON.txt", "sector/I:y:z:status");

  T1->Write();
  T1->Delete();
  T2->Write();
  T2->Delete();
  f->Close();
}
