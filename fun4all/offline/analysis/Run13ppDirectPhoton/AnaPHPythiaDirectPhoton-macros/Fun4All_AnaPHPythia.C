//#include <libgen.h>

void Fun4All_AnaPHPythia(const int nevents = 100000,
                 const int runnumber = 0,
                 //const char *pyfname = "keep/phpythia_directphoton.root")
		 const char *pyfname = "data_phpythia/phpythia_pp_510_all_seed1.root")
//const char *pyfname = "data_phpythia/phpythia_pp_510_dirphoton_seed1.root")
{
  gSystem->Load("libfun4allfuncs.so");	// framework + reco modules
  gSystem->Load("libAnaPHPythiaDirectPhoton.so");

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",runnumber);

  /////////////////////////////////////////////////////////////////
  //  Server...
  Fun4AllServer *se = Fun4AllServer::instance();

  /////////////////////////////////////////////////////////////////
  //  Reconstruction Modules...
 
  SubsysReco *anaphpythia = new AnaPHPythiaDirectPhoton("anaphpythia.root");
  se->registerSubsystem(anaphpythia);

  Fun4AllDstInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  in1->AddFile( pyfname );
  se->registerInputManager(in1);

  se->run(nevents);  // run over all events
  se->End();

  cout << "All done." << endl;
}

