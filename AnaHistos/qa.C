qa(char* flist="data/run_good-mine.list")
// Show pi0 mass and width vs run
{
  const bool bgamma = false; // 
  const bool read_db = true; // Normalize to MB counts from in DB

  //  char* fdb = "../../ntrig/runsummary_500gev.log";
  char* fdb = "/gpfs/mnt/gpfs02/phenix/spin/spin1/data77/phnxsp01/shura/trig/run13/summary_run13pp510gev.log";

  // ../data_gamma/pi0_nt_gamma_280187.root
  //  const int nchar = 27;
  // ../data_gamma_final/pi0_nt_gamma_280185.root
  //  const int nchar = 33;
  // ../data_mb_final/pi0_nt_mb_280185.root
  // /phenix/spin/phnxsp01/shura/taxi/Run13pp510MinBias/5096/data/387027.root

  //int nchar = 61;  // for run_good.list
  int nchar = 57;  // for run_good-mine.list
  if( bgamma ) nchar = 33;
  //  if( bgamma ) nchar = 27; // for gamma_raw

  const float ClockRate = 1.0e7; // 10 MHz

  float ptmin = 2.0; // MB
  if( bgamma )  ptmin = 6.0; // G3

  const float ptmax = 10.0;

  const int imin = ptmin*2;
  const int imax = ptmax*2;
  const int IB = 0;
  const float xmin = 0.06; // for fit
  const float xmax = 0.25; // for fit
  const int ibeg = 12;
  const int iend = 16;

  const int NMAX = 10000;
  const int NSEC = 3; // PbSc-West, PbSc-East and PbGl
  char* sname[NSEC] = {"PbSc-West","PbSc-East","PbGl"};
  float mass[NSEC][NMAX];
  float emass[NSEC][NMAX];
  float width[NSEC][NMAX];
  float ewidth[NSEC][NMAX];
  float npi0[NSEC+1][NMAX];
  float enpi0[NSEC+1][NMAX];
  float npi0_corr[NSEC+1][NMAX];
  float enpi0_corr[NSEC+1][NMAX];
  float rbg[NSEC+1][NMAX];
  float erbg[NSEC+1][NMAX];
  float rrun[NMAX];
  float rate[NMAX];
  float rejection[NMAX];

  printf("Min/Max pT= %f/%f (imin/imax=%d/%d)\n",ptmin,ptmax,imin,imax);

  // Read DB data

  int runlist_db[NMAX];
  ULong64_t ntr_mb_sc[NMAX];
  ULong64_t ntr_mb_lv[NMAX];
  ULong64_t ntr_cl_sc[NMAX];
  ULong64_t ntr_cl_lv[NMAX];
  int sdown_mb[NMAX];
  float evrate[NMAX];
  float ltime[NMAX];
  int ndb = 0;
  
  if( read_db ) {

    ifstream inFile(fdb);
    
    char tmp[256];
    if(!inFile){
      printf("Error: failed to open file %s\n",fdb);
      return;
    }
    inFile.getline(tmp, 255);
    inFile.getline(tmp, 255);
    
    int id;
    int fill, run, nraw, nlive, nscaled, scaledown, tbeg, tend;
    int idum;
    char tname[50], sdum[10];
    while(inFile){
      inFile.getline(tmp, 255);
      sscanf( tmp, "%d %s %d %s %d %s %s\n", &fill, sdum, &run, sdum,  &idum, sdum, tname );
      // Check the format below !!!!!
      sscanf( &(tmp[85]), "%d %s %d %s %d %s %d\n", &nraw, sdum, &nlive, sdum, &nscaled, sdum, &scaledown );
      sscanf( &(tmp[141]), "%d %s %d\n", &tbeg, sdum, &tend );
      //    printf("%d %s %d %d %d\n",run, tname, nraw, nlive, nscaled );
      
      id = FindRun(ndb,runlist_db,run);
      if( ndb<=0 || (ndb>0 && id<0) )  {
	runlist_db[ndb]=run;
	ndb++;
	id = ndb-1;
      }
      
      //      if( strcmp(tname,"BBCLL1(>0 tubes) novertex")==0 ) {
      if( strcmp(tname,"BBCLL1(>0")==0 ) {
	//	printf("Found!!!\n");
	ntr_mb_sc[id] = nscaled;
	//      ntr_mb_lv[id] = nlive;
	ntr_mb_lv[id] = ntr_mb_sc[id]*(scaledown+1);
	sdown_mb[id] = scaledown;
	//	evrate[id] = (float)ntr_mb_lv[id]/(tend-tbeg);
      }

      if( strcmp(tname,"CLOCK")==0 ) {
	ntr_cl_sc[id] = nscaled;
	ntr_cl_lv[id] = ntr_cl_sc[id]*(scaledown+1);
	ltime[id] = float(ntr_cl_lv[id])/ClockRate/(tend-tbeg);
      }
      
    } // while(inFile

    printf("nRuns=%d read from DB file %s\n",ndb,fdb);

    for( int id=0; id<ndb; id++ ) {
      evrate[id] = 0;
      if( ltime[id]>0.5 && ltime[id]<1.0 )
	evrate[id] = (float)ntr_mb_lv[id]/float(ntr_cl_lv[id]);
      else 
	printf("Bad Run %d: ltime %f\n", runlist_db[id], ltime[id]);
      //      printf("EvRate %e\n",evrate[id]);
    }

  }

  //  return;

  // Read data from nt files
  
  ifstream inFile(flist);

  if(!inFile){
    printf("Error: failed to open file %s\n",flist);
    return false;
  }

  TFile* fnt;
  TH1F* hm0;
  TH1F* hm;
  TH1F* htmp;
  TF1* fit = new TF1("fit","gaus+pol2(3)",xmin,xmax);
  TF1* fsub = new TF1("fsub","pol2",xmin,xmax);

  char fin[256];
  char hname[30];
  float sum;
  float sumb;
  float nev, rej;
  int irun, idb;

  int nn0 = 0;

  while(inFile){
    
    if( nn0%10 == 0 && nn0!=0 ) printf("Nfiles = %d\n",nn0);
    
    inFile.getline(fin, 255);
    if(fin[0]=='#' || fin[0] == 0) continue; // strip comment or empty lines

    sscanf(&(fin[nchar]),"%6d",&irun);
    //    printf("%d\n",irun);

    rrun[nn0] = irun;

    rate[nn0] = 0;
    if( read_db ) {
      idb = FindRun(ndb,runlist_db,irun);
      if( idb<0 ) {
	printf("No run %d in DB data: %s. May need to tune nchar parameter\n",irun,fin);
	continue;
      }
      rate[nn0] = evrate[idb];
      //      if( evrate[idb]>0.12 ) printf("*** Run %d: rate %f\n",irun,evrate[idb]);
    }

    fnt = new TFile(fin);
    if( !fnt->IsOpen() ) {
      printf("Failed to open file %s\n",fin);
      //    return false;
      continue;
    }

    htmp = (TH1F*)fnt->Get("hevtype");

    if( bgamma ) { // Gamma3
      rej  = 1;
      if( read_db )
	nev = ntr_mb_lv[idb];
	//	nev = htmp->GetBinContent(1);
      else 
	nev = htmp->GetBinContent(IB+1);
    }
    else { // MB
      rej  = htmp->GetBinContent(1)/htmp->GetBinContent(2);
      if( read_db )
	//	nev = ntr_mb_sc[idb];
	nev = htmp->GetBinContent(IB+1);
      else 
	nev = htmp->GetBinContent(IB+1);
    }

    if( nev<=0 ) {
      printf("Nmb = %d in fill %d\n",nev,irun);
      continue;
    }

    rejection[nn0] = rej;

    for( int is=0; is<3; is++ ) {

      htmp = (TH1F*)fnt->Get("mc_s0_bcc0_pt_020_p");
      hm = (TH1F*)htmp->Clone();
      hm->Reset();
      for( int ip=imin; ip<imax; ip++ ) {
	sprintf(hname,"mc_s%d_bcc%d_pt_%03d_p",is,IB,5*ip); // Pof cut
	//	sprintf(hname,"mc_s%d_bcc%d_pt_%03d_tp",is,IB,5*ip); // Prof&ToF cut
	hm0 = (TH1F*)fnt->Get(hname);
	hm->Add(hm0);
      }
      hm->Rebin(10);
      hm->Scale(0.5);
      fit->SetParameters(hm->GetMaximum(),0.140,0.010,0,0,0);
      hm->Fit(fit,"QN","",xmin,xmax);

      mass[is][nn0]  = fit->GetParameter(1);
      emass[is][nn0]  = fit->GetParError(1);
      width[is][nn0] = fit->GetParameter(2);
      ewidth[is][nn0] = fit->GetParError(2);
      if( emass[is][nn0] > 0.1 || ewidth[is][nn0] > 0.1 ) {
	mass[is][nn0]  = 0;
	emass[is][nn0]  = 0;
	width[is][nn0] = 0;
	ewidth[is][nn0] = 0;	
      }

      sum = 0;
      for( int ib=ibeg; ib<=iend; ib++ ) {
	sum += hm->GetBinContent( ib );
      }
      if( sum<=0 ) sum=1; // !!!!! To avoid zero errors
      npi0[is][nn0] = sum/nev;
      enpi0[is][nn0] = sqrt(sum)/nev;
      if( npi0[is][nn0] <= 0 ) printf("Bad Run %d: pi0 rate [%d] = %e\n",irun,is,npi0[is][nn0]);
      //      if( is==0 && npi0[is][nn0] <= 1.9e-6 ) printf("Bad Run %d (%d): pi0 rate [%d] = %e\n",irun,nn0,is,npi0[is][nn0]);
      //      printf("%f %f %f\n",sum,nev,npi0[is][nn0]);

      fsub->SetParameter(0, fit->GetParameter(3) );
      fsub->SetParameter(1, fit->GetParameter(4) );
      fsub->SetParameter(2, fit->GetParameter(5) );
      sumb = 0;
      for( int ib=ibeg; ib<=iend; ib++ ) {
	sumb += fsub->Eval( hm->GetBinCenter(ib) );
      }
      if( sumb<=0 ) sumb = 1; // !!!!! To avoid zero errors
      rbg[is][nn0] = 0;
      erbg[is][nn0] = 0;
      npi0_corr[is][nn0] = 0;
      enpi0_corr[is][nn0] = 0;
      if( sum>0 ) {	
	rbg[is][nn0] = sumb/sum;
	erbg[is][nn0] = sqrt(sumb)/sum;
	npi0_corr[is][nn0] = npi0[is][nn0]*(1-rbg[is][nn0]);
	enpi0_corr[is][nn0] = sqrt((1-rbg[is][nn0])*(1-rbg[is][nn0])*enpi0[is][nn0]*enpi0[is][nn0] + npi0[is][nn0]*npi0[is][nn0]*erbg[is][nn0]*erbg[is][nn0]);
      }

    } // for( int isec

    nn0++;
    fnt->Close();

  }  // while(inFile

  printf("Number of files read: %d\n",nn0);

  float order[NMAX];
  float zero[NMAX];
  for( int i=0; i<NMAX; i++ ) {
    order[i] = i+1;
    zero[i] = 0;
  }

  // PbSc
  for( int i=0; i<nn0; i++ ) { 
    npi0[3][i] = npi0[0][i]+npi0[1][i];
    enpi0[3][i] = sqrt(enpi0[0][i]*enpi0[0][i]+enpi0[1][i]*enpi0[1][i]);
    npi0_corr[3][i] = npi0_corr[0][i]+npi0_corr[1][i];
    enpi0_corr[3][i] = sqrt(enpi0_corr[0][i]*enpi0_corr[0][i]+enpi0_corr[1][i]*enpi0_corr[1][i]);
    rbg[3][i] = (rbg[0][i]+rbg[1][i])/2.;
    erbg[3][i] = sqrt(erbg[0][i]*erbg[0][i]+erbg[1][i]*erbg[1][i])/2.;
    //    if( npi0[3][i] < 2.7e-6  ) {npi0[3][i] *= 2.; enpi0[3][i] *= 2.;} // !!!!!
    //    if( npi0[2][i] < 0.36e-6 ) {npi0[2][i] *= 2.; enpi0[2][i] *= 2.;} // !!!!!
  }

  TGraphErrors* gm[NSEC];
  TGraphErrors* gw[NSEC];
  TGraphErrors* gnpi0[NSEC];
  TGraphErrors* gbg[NSEC];
  for( int is=0; is<NSEC; is++ ) {
    gm[is] = new TGraphErrors(nn0,order,mass[is],zero,emass[is]);
    gw[is] = new TGraphErrors(nn0,order,width[is],zero,ewidth[is]);
    gnpi0[is] = new TGraphErrors(nn0,order,npi0[is],zero,enpi0[is]);
    gbg[is] = new TGraphErrors(nn0,order,rbg[is],zero,erbg[is]);
  }

  TGraphErrors* grate[2];
  grate[0] = new TGraphErrors(nn0,rate,npi0_corr[3],zero,enpi0_corr[3]);
  grate[1] = new TGraphErrors(nn0,rate,npi0_corr[2],zero,enpi0_corr[2]);
  //  grate[0] = new TGraphErrors(nn0,rate,mass[0],zero,emass[0]);
  //  grate[1] = new TGraphErrors(nn0,rate,mass[2],zero,emass[2]);
  //  grate[0] = new TGraphErrors(nn0,rate,rbg[3],zero,erbg[3]);
  //  grate[1] = new TGraphErrors(nn0,rate,rbg[2],zero,erbg[2]);


  //  TGraphErrors* grate = new TGraphErrors(nn0,order,rate,zero,zero);
  //  TGraphErrors* grate = new TGraphErrors(nn0,rate,npi0_corr[0],zero,enpi0_corr[0]);
  //  TGraphErrors* grate = new TGraphErrors(nn0,rate,rbg[0],zero,erbg[0]);

  //  for( int i=0; i<nn0; i++ ) if( rate[i]<0.041 || rate[i]>0.12 ) printf("*** Run %d: Rate %f  Npi0=%e +/- %e\n",rrun[i],rate[i],npi0[0][i],enpi0[0][i]);

  TGraphErrors* grej = new TGraphErrors(nn0,rate,rejection,zero,zero);

  TF1* fpol1 = new TF1("fpol1","pol1");
  grej->Fit(fpol1,"Q");
  grej->SetMarkerStyle(20);
  //  grej->Draw("AP");
  //  return;

  gStyle->SetOptStat(0);
  c1 = new TCanvas("c1","The Ntuple canvas",10,10,350,600);
  c1->Range(0,0,1,1);
  TPad* pad1[3];
  pad1[0] = new TPad("pad10","This is pad1",0.0,0.66,1.0,0.99);
  pad1[1] = new TPad("pad11","This is pad2",0.0,0.33,1.0,0.66);
  pad1[2] = new TPad("pad12","This is pad2",0.0,0.00,1.0,0.33);
  for( int i=0; i<3; i++ ) pad1[i]->Draw();

  for( int is=0; is<3; is++ ) {
    pad1[is]->cd();
    sprintf(hname, "Mass: %s",sname[is]);
    gm[is]->SetTitle(hname);
    gm[is]->SetMarkerStyle(20);
    gm[is]->Draw("AP");
  }

  gStyle->SetOptStat(0);
  c2 = new TCanvas("c2","The Ntuple canvas",100,100,350,600);
  c2->Range(0,0,1,1);
  TPad* pad2[3];
  pad2[0] = new TPad("pad20","This is pad1",0.0,0.66,1.0,0.99);
  pad2[1] = new TPad("pad21","This is pad2",0.0,0.33,1.0,0.66);
  pad2[2] = new TPad("pad22","This is pad2",0.0,0.00,1.0,0.33);
  for( int i=0; i<3; i++ ) pad2[i]->Draw();

  for( int is=0; is<3; is++ ) {
    pad2[is]->cd();
    sprintf(hname, "Width: %s",sname[is]);
    gw[is]->SetTitle(hname);
    gw[is]->SetMarkerStyle(20);
    gw[is]->Draw("AP");
  }

  gStyle->SetOptStat(0);
  c3 = new TCanvas("c3","The Ntuple canvas",200,200,350,600);
  c3->Range(0,0,1,1);
  TPad* pad3[3];
  pad3[0] = new TPad("pad30","This is pad1",0.0,0.66,1.0,0.99);
  pad3[1] = new TPad("pad31","This is pad2",0.0,0.33,1.0,0.66);
  pad3[2] = new TPad("pad32","This is pad2",0.0,0.00,1.0,0.33);
  for( int i=0; i<3; i++ ) pad3[i]->Draw();

  for( int is=0; is<3; is++ ) {
    pad3[is]->cd();
    sprintf(hname, "Yield: %s",sname[is]);
    gnpi0[is]->SetTitle(hname);
    gnpi0[is]->SetMarkerStyle(20);
    gnpi0[is]->Draw("AP");
  }

  gStyle->SetOptStat(0);
  c4 = new TCanvas("c4","The Ntuple canvas",300,300,350,600);
  c4->Range(0,0,1,1);
  TPad* pad4[3];
  pad4[0] = new TPad("pad40","This is pad1",0.0,0.66,1.0,0.99);
  pad4[1] = new TPad("pad41","This is pad2",0.0,0.33,1.0,0.66);
  pad4[2] = new TPad("pad42","This is pad2",0.0,0.00,1.0,0.33);
  for( int i=0; i<3; i++ ) pad4[i]->Draw();

  for( int is=0; is<3; is++ ) {
    pad4[is]->cd();
    sprintf(hname, "Bgr: %s",sname[is]);
    gbg[is]->SetTitle(hname);
    gbg[is]->SetMarkerStyle(20);
    gbg[is]->Draw("AP");
  }

  c5 = new TCanvas("c5","The Ntuple canvas",400,400,350,700);
  c5->Range(0,0,1,1);
  TPad* pad5[2];
  pad5[0] = new TPad("pad50","This is pad1",0.,0.5,1.0,1.0);
  pad5[1] = new TPad("pad51","This is pad1",0.,0.0,1.0,0.5);
  pad5[0]->Draw();
  pad5[1]->Draw();

  TF1* fpol0 = new TF1("fpol0","pol0");
  grate[0]->Fit(fpol0,"QN");
  printf("<N>  PbSc: %e +/- %e\n",fpol0->GetParameter(0),fpol0->GetParError(0));
  grate[1]->Fit(fpol0,"QN");
  printf("<N>  PbGl: %e +/- %e\n",fpol0->GetParameter(0),fpol0->GetParError(0));

  TF1* fpol1 = new TF1("fpol1","pol1");
  grate[0]->Fit(fpol1,"Q");
  printf("<N0> PbSc: %e +/- %e\n",fpol1->GetParameter(0),fpol1->GetParError(0));
  grate[1]->Fit(fpol1,"Q");
  printf("<N0> PbSc: %e +/- %e\n",fpol1->GetParameter(0),fpol1->GetParError(0));

  pad5[0]->cd();
  grate[0]->SetTitle("PbSc: Npi0/Nmb vs Nmb");
  grate[0]->Draw("AP");

  pad5[1]->cd();
  grate[1]->SetTitle("PbGl: Npi0/Nmb vs Nmb");
  grate[1]->Draw("AP");

  TFile *f_out = new TFile("data/Pileup-Sasha.root", "RECREATE");
  TGraphErrors *gr_run[4];
  gr_run[0] = new TGraphErrors(nn0,rrun,npi0[3],zero,enpi0[3]);
  gr_run[1] = new TGraphErrors(nn0,rrun,npi0[2],zero,enpi0[2]);
  gr_run[2] = new TGraphErrors(nn0,rrun,npi0_corr[3],zero,enpi0_corr[3]);
  gr_run[3] = new TGraphErrors(nn0,rrun,npi0_corr[2],zero,enpi0_corr[2]);
  for(int igr=0; igr<4; igr++)
  {
    gr_run[igr]->SetName(Form("gr_run_%d",igr));
    gr_run[igr]->Write();
  }
  f_out->Close();
}


int FindRun(int nrun, int* runlist, int run)
{
  for( int i=0; i< nrun; i++ ) {
    if( runlist[i] == run ) return i;
  }
  return -1;
}
