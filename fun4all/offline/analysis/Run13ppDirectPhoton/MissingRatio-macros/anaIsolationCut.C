void anaIsolationCut(const int process=0)
{
  // Set up Fun4All libraries
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("librecal.so");
  gSystem->Load("libcteval.so");
  gSystem->Load("libcompactCNT.so");
  gSystem->Load("libemc.so");
  gSystem->Load("libemcEmbed4all.so");
  //gSystem->Load("libIsolationCut.so");
  gSystem->Load("libMissingRatio.so");

  const int nThread = 10;
  char dstFileName[1000];

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Setup recoConsts
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", 390039);
  rc->set_IntFlag("EMCNEW_DEBUG", 0);  // set debugging verbosity level
  rc->set_IntFlag("EMCNEW_PI0VERBOUT", 2);  // 2 - write only clean pi->2g cases

  // Register the master recal
  MasterRecalibratorManager *mr = new MasterRecalibratorManager("MASTERRECALIBRATORMANAGER");
  se->registerSubsystem(mr);

  // Import simulated data
  SubsysRecoStack *simimp = new EmcGeaContainerImporter();
  simimp->x_push_back( new EmcUnclusterizer() );
  EmcTowerScalerSmearer *emcsm = new EmcTowerScalerSmearer(1., 0.);
  emcsm->SetScale(0.98, 1.03);
  emcsm->SetSmear(0.08, 0.12);
  simimp->x_push_back(emcsm);
  se->registerSubsystem(simimp);

  // Clusterize data
  se->registerSubsystem( new EmcEmbedReclusterizer("TOP", "TOP", "TOP", "") );

  // Reconstruction Module
  //se->registerSubsystem( new EmcGeaContainerImporter() );
  SubsysReco *my1 = new IsolationCut( "test.root" );
  se->registerSubsystem(my1);

  // Input Manager
  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

//  // Process single file
//  //int openReturn = se->fileopen("DSTin1", "data_pisa/simDST_phpythia_pp510GeV_all_minpt1_seed1.root");
//  //int openReturn = se->fileopen("DSTin1", "data_pisa/simDST_phpythia_pp510GeV_dirphoton_minpt1_seed1.root");
//
//  int openReturn = se->fileopen("DSTin1", "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/zji/data/pisaRun13/simDST/simDST1.root");
//
//  if(openReturn)
//    {
//      cout << "\nAbnormal return: openReturn from Fun4All fileopen method = " << openReturn << endl;
//      continue;
//    }
//  se->run(0);
//  cout << "\nClosing input file, and a No Input file open message from Fun4All should appear" << endl;
//  int closeReturn = se->fileclose("DSTin1");
//  if(closeReturn)
//    cout << "\nAbnormal return: closeReturn from Fun4All fileclose = " << closeReturn << endl;


  // Loop over input DST files
  for(int fileid=0; fileid<=999; fileid++)
  {
    sprintf(dstFileName, "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/zji/data/pisaRun13/simDST-isolation/simDST%d.root", fileid);

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
  }

  // Write out the histogram file
  se->End();

  delete my1;
  delete se;
}
