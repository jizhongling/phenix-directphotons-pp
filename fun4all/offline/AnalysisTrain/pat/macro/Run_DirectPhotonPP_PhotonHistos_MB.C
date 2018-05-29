// .x RunMyMacro.C("Run_DirectPhotonPP_PhotonHistos_MB.C","num.root",1000,"Run13pp510MB_Fast")

void Run_DirectPhotonPP_PhotonHistos_MB(const char *filename = "num.root")
{
  gSystem->Load("libDirectPhotonPP.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  PhotonHistos *my1 = new PhotonHistos("PHOTONHISTOS", filename);
  my1->SelectMB();
  se->registerSubsystem(my1);
}

void InputData(vector<string> &indata)
{
  indata.push_back("CNT");
  return;
}
