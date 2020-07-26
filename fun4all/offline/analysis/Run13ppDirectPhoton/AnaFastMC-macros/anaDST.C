void anaDST(const int process = 0,
    const int scale = 40,
    const char *dstFileName = "simDST.root",
    const char *histoname = "histo.root")
{
  gSystem->Load("libfun4all.so");	// framework + reco modules
  gSystem->Load("librecal.so");
  gSystem->Load("libcteval");
  gSystem->Load("libemcEmbed4all.so");
  gSystem->Load("libAnaFastMC.so");
  gSystem->Load("libMissingRatio.so");

  // Setup recoConsts
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", 390039);
  rc->set_IntFlag("EMCNEW_DEBUG", 0);  // set debugging verbosity level
  rc->set_IntFlag("EMCNEW_PI0VERBOUT", 2);  // 2 - write only clean pi->2g cases

  // Server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Register Master Recalibrator
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

  // Reclusterize data
  se->registerSubsystem( new EmcEmbedReclusterizer("TOP", "TOP", "TOP", "") );

  // My Reconstruction Module
  HadronResponse *my1 = new HadronResponse("HadronResponse");
  //MissingRatio *my1 = new MissingRatio("MissingRatio");
  my1->set_outfile(histoname);
  my1->use_xsec_weight();
  se->registerSubsystem(my1);

  // Real input from DST files
  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // Loop over input DST file
  cout << "\nfileopen for " << dstFileName << endl; 
  int openReturn = se->fileopen("DSTin1", dstFileName);
  if(openReturn)
  {
    cout << "\nAbnormal return: openReturn from Fun4All fileopen method = " << openReturn << endl;
    return;
  }

  // Do the analysis for this DST file
  //double pt_start = 3. + process/scale * 0.1;
  //PtWeights *ptweights = new PtWeights();
  //double weight_pythia = ptweights->Integral(pt_start, pt_start+1., "MinBias");
  //my1->SetWeightPythia(weight_pythia);
  se->run(0);

  cout << "\nClosing input file, and a No Input file open message from Fun4All should appear" << endl;
  int closeReturn = se->fileclose("DSTin1");
  if(closeReturn)
  {
    cout << "\nAbnormal return: closeReturn from Fun4All fileclose = " << closeReturn << endl;
    return;
  }

  // Write histograms
  se->End();

  // Delete Fun4All server
  delete se;
}
