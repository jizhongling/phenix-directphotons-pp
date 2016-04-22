Double_t GetRho(const TMatrixD &b, const vector<Double_t> &y,
    const Int_t nR, const Int_t nG);

Double_t BgGPR(vector<Double_t> &x1, vector<Double_t> &y1, vector<Double_t> &sigma_y1,
    vector<Double_t> &x2, vector<Double_t> &y2, vector<Double_t> &sigma_y2)
{
  const Int_t verbosity = 0;
  const char *outfile = "BgGPR.root";

  const Double_t xmin = 0.047;
  const Double_t xmax = 0.227;
  const Int_t nPredictions = 120;
  const Int_t nRange = 34;
  const Int_t nGap = 10;

  GausProc a(x1, y1, sigma_y1, xmin, xmax, nPredictions, outfile);
  //gSystem->Exec(Form("rm -f %s",outfile));

  a.SetVerbosity(verbosity);
  a.SetKernel(GausProc::RBF);
  a.warp(0);  //should be used only if you have data spanning several orders of magnitude

  GPOptimizer c(&a, 10., 1.);
  c.GPoptimize(2,0);

  GausProc d(x2, y2, sigma_y2, xmin, xmax, nPredictions, outfile);
  d.warp(0);  //see comment above
  d.SetPar(0, c.getPar(0));
  d.SetPar(1, c.getPar(1));
  d.process();
  //d.Write(-1);
  d.unwarp(0);  //see comment above
  //d.Write(-1, "_unwarp");

  TMatrixDSym *covMatrix = d.GetCovarianceMatrix();
  Double_t rho = GetRho(*covMatrix, y2, nRange, nGap);
  delete covMatrix;

  return rho;
}

Double_t GetRho(const TMatrixD &b, const vector<Double_t> &y,
    const Int_t nR, const Int_t nG)
{
  Double_t rhobar = 0.;
  Double_t count = 0.;
  for(Int_t i=0; i<b.GetNrows(); i++)
    for(Int_t j=0; j<b.GetNcols(); j++)
      if( (i>=nR+nG && i<nR*2+nG) && (j<nR || j>=nR*2+nG*2) )
      {
        Double_t rho = b(i,j) / sqrt( b(i,i) * b(j,j) );
        Double_t weight = y.at(i*3/2);
        rhobar += fabs(rho) * weight;
        count += weight;
      }
  rhobar /= count;
  return rhobar;
}
