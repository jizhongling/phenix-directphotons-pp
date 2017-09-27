void Run_DirectPhotonPP_PhotonNodeMB(const char *filename = "TREE.root")
{
  gSystem->Load("libDirectPhotonPP.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  PhotonNodeMB *my1 = new PhotonNodeMB("PHOTONNODEMB");
  se->registerSubsystem(my1);

  string outFile = "DirectPhotonPP_PhotonNode-";
  outFile.append(filename);
  Fun4AllOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", outFile.c_str());
  out->AddEventSelector("PHOTONNODEMB");
  out->AddNode("PhotonContainerMB");
  se->registerOutputManager(out);
}

void InputData(vector<string> &indata)
{
  indata.push_back("CNT");
  return;
}
