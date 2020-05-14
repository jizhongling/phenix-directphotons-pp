void BgGPR(vector<double> &x, vector<double> &y, vector<double> &sigma_y,
    double &Integral, double &dIntegral,
    const double xmin, const double xmax, const int nPredictions,
    const double xint1, const double xint2,
    const char *outfile = "default.root", const int index = 0)
{
  const int verbosity = 0;

  GausProc a(x, y, sigma_y, xmin, xmax, nPredictions, outfile);

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
  if(strcmp(outfile,"default.root") != 0)
  {
    d.Write(-1, Form("_%d",index));
    //d.unwarp(0);  //see comment above
    //d.Write(-1, Form("_unwarp_%d",index));
  }

  d.Integral(xint1, xint2, Integral, dIntegral);

  return;
}
