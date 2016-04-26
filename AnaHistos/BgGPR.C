void BgGPR(vector<Double_t> &x, vector<Double_t> &y, vector<Double_t> &sigma_y,
    Double_t &Integral, Double_t &dIntegral)
{
  const Int_t verbosity = 0;
  const char *outfile = "BgGPR.root";

  const Double_t xmin = 0.047;
  const Double_t xmax = 0.227;
  const Int_t nPredictions = 90;

  GausProc a(x, y, sigma_y, xmin, xmax, nPredictions, outfile);
  //gSystem->Exec(Form("rm -f %s",outfile));

  a.SetVerbosity(verbosity);
  a.SetKernel(GausProc::RBF);
  //a.warp(0);  //should be used only if you have data spanning several orders of magnitude

  GPOptimizer c(&a, 2., 10.);
  c.GPoptimize(2,0);

  GausProc d(x, y, sigma_y, xmin, xmax, nPredictions, outfile);
  //d.warp(0);  //see comment above
  d.SetPar(0, c.getPar(0));
  d.SetPar(1, c.getPar(1));
  d.process();
  //d.Write(-1);
  //d.unwarp(0);  //see comment above
  //d.Write(-1, "_unwarp");

  d.Integral(0.112, 0.162, Integral, dIntegral);

  return;
}
