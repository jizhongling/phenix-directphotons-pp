void anaFillHisto_Sample()
{
  // Set up Fun4All libraries
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("libPhotonNode.so");

  // Server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Reconstruction Module
  FillHisto *my1 = new FillHisto("FillHisto_Sample");
  my1->SelectERT();
  se->registerSubsystem(my1);

  // Input Manager
  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // Input DST files
  const char *dstFileName = "/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/AnalysisTrain/pat/macro/PhotonNode-num.root";

  cout << "\nfileopen for " << dstFileName << endl; 
  int openReturn = se->fileopen("DSTin1", dstFileName);
  if(openReturn)
  {
    cout << "\nAbnormal return: openReturn from Fun4All fileopen method = " << openReturn << endl;
    continue;
  }

  // Do the analysis for this DST file
  se->run(0);

  cout << "\nClosing input file, and a No Input file open message from Fun4All should appear" << endl;
  int closeReturn = se->fileclose("DSTin1");
  if(closeReturn)
    cout << "\nAbnormal return: closeReturn from Fun4All fileclose = " << closeReturn << endl;

  // Write out the histogram file
  se->End();
}
