// .x RunMyMacro.C("Run_DirectPhotonPP_PhotonNode_MB.C","num.root",1000,"Run13pp510MB_Fast")

void Run_DirectPhotonPP_PhotonNode_MB(const char *filename = "num.root")
{
  gSystem->Load("libDirectPhotonPP.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  PhotonNode *my1 = new PhotonNode("PhotonNode");
  my1->SelectMB();
  se->registerSubsystem(my1);

  string outFile = "PhotonNode-";
  outFile.append(filename);
  Fun4AllOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", outFile.c_str());
  out->AddEventSelector("PHOTONNODE");
  out->AddNode("PhotonContainer");
  se->registerOutputManager(out);
}

void InputData(vector<string> &indata)
{
  indata.push_back("CNT");
  return;
}
