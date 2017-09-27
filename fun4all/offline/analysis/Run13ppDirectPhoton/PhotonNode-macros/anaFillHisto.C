void anaFillHisto(const char *filelist="taxifiles.txt", const char *outfile = "histo.root", const int nEvents=0)
{
  // Set up input file location
  gSystem->Load("libFROG.so");
  FROG fr;

  // Set up Fun4All libraries
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("libPhotonNode.so");

  char dstString[5000];
  char *dstFileName;

  ifstream inFiles(filelist);
  if(!inFiles)
  {
    cerr << "\nUnable to open input file list " << filelist << endl;
    return;
  }
  cout << "\nUsing input files list " << filelist << endl << endl;

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Reconstruction Module
  //SubsysReco *my1 = new FillHisto("FILLHISTO", outfile);
  SubsysReco *my1 = new FillHistoMB("FILLHISTOMB", outfile);
  se->registerSubsystem(my1);

  // Input Manager
  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // Loop over input DST files
  while( inFiles >> dstString )
  {
    dstFileName = fr.location(dstString);  // get FROG location
    if(!dstFileName)
    {
      cerr << "\nBad return from FROG, null dstFileName" << endl;
      return;
    }

    cout << "\nfileopen for " << dstFileName << endl; 
    int openReturn = se->fileopen("DSTin1", dstFileName);
    if(openReturn)
      cout << "\nAbnormal return: openReturn from Fun4All fileopen method = " << openReturn << endl; 

    // Do the analysis for this DST file
    se->run(nEvents);

    cout << "\nClosing input file, and a No Input file open message from Fun4All should appear" << endl;
    int closeReturn = se->fileclose("DSTin1");
    if(closeReturn)
      cout << "\nAbnormal return: closeReturn from Fun4All fileclose = " << closeReturn << endl;
  }

  // Write out the histogram file
  se->End();

  delete my1;
  delete se;
  inFiles.close();
}
