//
// This macro is an example for generating single particles written in the 
// PHPythia format (which is fun4all compliant)
//
void
phparticlegen(const int process = 0, const int nevents = 10, const char *outputname = "phparticlegen.root")
{
  gSystem->Load("libPythia6.so");
  gSystem->Load("libfun4all.so");	// framework + reco modules
  gSystem->Load("libpdbcalBase.so");
  gSystem->Load("libPHPythia.so");
  //gSystem->Load("libPHPyTrigger.so");	// only needed if triggering
  gSystem->Load("libPHParticleGen.so");
  //gSystem->Load("libPHParticleCustomGen.so");	// only needed for custom generator

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",0);

  /////////////////////////////////////////////////////////////////
  //  Server...
  Fun4AllServer *se = Fun4AllServer::instance();

  /////////////////////////////////////////////////////////////////
  //  Reconstruction Modules...

  PHParticleGen* ss = new PHParticleGen();

  //*** Put your own custom generator here, in case what's available
  //*** won't work for you.  If you have some distribution which would
  //*** be generally useful, mail winter@alum.mit.edu or chiu@bnl.gov
  //TSingleParticleCustomGenerator *spg = new TSingleDYGenerator();
  //TSingleDYGenerator *spg = new TSingleDYGenerator();
  //ss->SetGenerator(spg);

  //ss->SetSeed(1999);			// Set your own seed, otherwise, defaults to /dev/random

  //*** Set Particle Parameters. Overrides config file.

  // Particle Type
  //ss->GetGenerator()->SetParameter("pid",13);		// mu+
  //ss->GetGenerator()->SetParameter("pid",111);	// pi0
  ss->GetGenerator()->SetParameter("pid",process%2?-211,211);	// even proc: pi+; odd proc: pi-
  //ss->GetGenerator()->SetParameter("pid",221);	// eta
  //ss->GetGenerator()->SetParameter("pid",22);		// gamma
  //ss->GetGenerator()->SetParameter("pid",23);		// gamma*/Z0
  //  ss->GetGenerator()->SetParameter("pid",11);	// electron

  // Momentum/Transverse Momentum/Energy
  //ss->GetGenerator()->SetParameter("momflag",TSingleParticleGenerator::STEPS_ALL);
  ss->GetGenerator()->SetParameter("momflag",TSingleParticleGenerator::FLAT_PT);
  //ss->GetGenerator()->SetParameter("momflag",TSingleParticleGenerator::FLAT_E);
  //ss->GetGenerator()->SetParameter("momflag",TSingleParticleGenerator::FLAT_P);
  ss->GetGenerator()->SetParameter("pmin",0.);
  ss->GetGenerator()->SetParameter("pmax",40.);

  // Eta distribution
  ss->GetGenerator()->SetParameter("etamin",-0.5);
  ss->GetGenerator()->SetParameter("etamax",0.5);

  // steps and phi distribution
  ss->GetGenerator()->SetParameter("nstepse",20.);
  ss->GetGenerator()->SetParameter("nstepstheta",50.);
  ss->GetGenerator()->SetParameter("nstepsphi",50);
  ss->GetGenerator()->SetParameter("phitot",0.5);

  // Vertex Distribution (sometimes better to leave at 0)
  //ss->GetGenerator()->SetParameter("vtxflag",TSingleParticleGenerator::FLAT);
  //ss->GetGenerator()->SetParameter("vtxflag",TSingleParticleGenerator::GAUSSIAN);
  //ss->GetGenerator()->SetParameter("zwidth",10);

  //** Force Your Particles To Decay using the pythia6 decayer
  //TPythia6Decayer *tpythia6decayer = TPythia6Decayer::Instance();
  //tpythia6decayer->ForceParticleDecay(221,22,2);  // Force eta->gamgam mode
  //ss->GetGenerator()->setDecayer( tpythia6decayer );

  //** (Uncomment only one of the two lines below)
  //ss->GetGenerator()->SetParameter("decayfile","default");
  //ss->GetGenerator()->SetParameter("decayfile","decay.table");

  se->registerSubsystem(ss);

  //** Here you could put a module to trigger on selected events
  //** For more info, see
  //**  https://www.phenix.bnl.gov/WWW/offline/wikioffline/index.php/PHPythia_Trigger
  //PHPyTrigger* trigger = new PHFvtxTrigger();
  //PHPyTrigger* trigger = new PHPyGammaNccTrigger();
  //trigger->Verbosity(1);
  //se->registerSubsystem(trigger);

  // This is to read in null events (necessary for Fun4All Framework)
  Fun4AllDummyInputManager *in1 = new Fun4AllDummyInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  TString OUTPUT = outputname;
  Fun4AllDstOutputManager *dst_output_mgr  = new Fun4AllDstOutputManager("PHPYTHIA",OUTPUT.Data());
  dst_output_mgr->AddNode("Sync");
  dst_output_mgr->AddNode("PHPythiaHeader");
  dst_output_mgr->AddNode("PHPythia");
  se->registerOutputManager(dst_output_mgr);

  cout<<" PHParticleGen will generate "<<nevents<<" events"<<endl;


  se->run(nevents);  // run over all events
  se->End();

}

