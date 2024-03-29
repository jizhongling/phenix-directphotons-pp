void draw_EventSasha()
{
  TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510MinBias/12233/data/Pi0PP-386881.root");

  TTree *treePi0 = (TTree*)f->Get("treePi0");

  int runno;
  int evtype;
  float bbc_z, bbc_t0;

  float pi0_eg1, pi0_eg2, pi0_tof1, pi0_tof2, pi0_prob1, pi0_prob2;
  float pi0_pt;
  float pi0_mc;
  int pi0_id1, pi0_id2;

  treePi0->SetBranchAddress("runno",&runno);
  treePi0->SetBranchAddress("evtype",&evtype);
  treePi0->SetBranchAddress("bbc_z",&bbc_z);
  treePi0->SetBranchAddress("bbc_t0",&bbc_t0);

  treePi0->SetBranchAddress("id1",&pi0_id1);
  treePi0->SetBranchAddress("id2",&pi0_id2);
  treePi0->SetBranchAddress("e1",&pi0_eg1);
  treePi0->SetBranchAddress("e2",&pi0_eg2);
  treePi0->SetBranchAddress("tof1",&pi0_tof1);
  treePi0->SetBranchAddress("tof2",&pi0_tof2);
  treePi0->SetBranchAddress("prob1",&pi0_prob1);
  treePi0->SetBranchAddress("prob2",&pi0_prob2);
  treePi0->SetBranchAddress("pt",&pi0_pt);
  treePi0->SetBranchAddress("mc",&pi0_mc);

  double ID1list[10];
  double ID2list[10];
  int count = 0;
  int nentries = treePi0->GetEntries();
  for(int i=0; i<nentries; i++)
  {
    treePi0->GetEntry(i);
    if( count < 10 &&
        (evtype & 0x1) && 
        TMath::Abs(bbc_z) < 10. &&
        pi0_prob1 > 0.02 &&
        pi0_prob2 > 0.02 && 
        TMath::Abs(pi0_tof1) < 10. &&
        TMath::Abs(pi0_tof2) < 10. &&
        TMath::Abs( (pi0_eg1 - pi0_eg2) / (pi0_eg1 + pi0_eg2) ) < 0.8 )
    {
      ID1list[count] = pi0_id1;
      ID2list[count] = pi0_id2;
      count++;
      cout << "Count " << count << endl
        << "Entry " << i << endl
        << "TowerID1 = " << pi0_id1 << endl
        << "TowerID2 = " << pi0_id2 << endl
        << "bbc_t0 = " << bbc_t0 << endl
        << "ToF1 = " << pi0_tof1 + bbc_t0 << endl
        << "ToF2 = " << pi0_tof2 + bbc_t0 << endl
        << "E1 = " << pi0_eg1 << endl
        << "E2 = " << pi0_eg2 << endl
        << "pi0_pT = " << pi0_pt << endl
        << "minv = " << pi0_mc << endl
        << endl;
    }
  }

  for(int i=0; i<10; i++)
    cout << ID1list[i] << ", ";
  cout << endl;
  for(int i=0; i<10; i++)
    cout << ID2list[i] << ", ";
  cout << endl;
}
