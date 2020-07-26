void anaFillHisto_Sample()
{
  gSystem->Load("libfun4all.so");	// framework + reco modules
  gSystem->Load("libPhotonNode.so");

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",0);

  // Server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // My Reconstruction Module
  FillHisto *my1 = new FillHisto("FillHisto_Sample");
  my1->SelectERT();
  se->registerSubsystem(my1);

  // Real input from DST files
  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // Input DST file
  char dstFileName[1000];
  sprintf(dstFileName, "/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/AnalysisTrain/pat/macro/PhotonNode-num.root");

  cout << "\nfileopen for " << dstFileName << endl; 
  int openReturn = se->fileopen("DSTin1", dstFileName);
  if(openReturn)
    cout << "\nAbnormal return: openReturn from Fun4All fileopen method = " << openReturn << endl;

  // Do the analysis for this DST file
  se->run(0);

  cout << "\nClosing input file, and a No Input file open message from Fun4All should appear" << endl;
  int closeReturn = se->fileclose("DSTin1");
  if(closeReturn)
    cout << "\nAbnormal return: closeReturn from Fun4All fileclose = " << closeReturn << endl;

  // Write histograms
  se->End();

  // Delete Fun4All server
  delete se;
}
