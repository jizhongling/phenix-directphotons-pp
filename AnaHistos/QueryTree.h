class QueryTree:
{
  public:
    QueryTree(const char *fname, const char *type = "READ");

    bool Query(const int ipt, const int part,
        double &xpt, double &value, double &error);

    void Fill(const int ipt, const int part,
        const double xpt, const double value, const double error);
    void Fill(const TH1 *h1, const TH1 *h2, const int part,
        const double kh1 = 1., const double kh2 = 1.);
    void Write(); 

  protected:
    int _ipt;
    int _part;
    double _xpt;
    double _value;
    double _error;

    TString _name;
    TFile *_file;
    TTree *_tree;
};

QueryTree::QueryTree(const char *fname, const char *type):
  _name(fname)
{
  if(strcmp(type,"READ") == 0)
  {
    _file = new TFile(fname, type);
    _tree = (TTree*)_file->Get("t1");
  }

  else if(strcmp(type,"RECREATE") == 0)
  {
    _file = new TFile(fname, type);
    _tree = new TTree("t1", "t1 tree");
    _tree->Branch("ipt", &_ipt, "ipt/I");
    _tree->Branch("part", &_part, "part/I");
    _tree->Branch("xpt", &_xpt, "xpt/D");
    _tree->Branch("value", &_value, "value/D");
    _tree->Branch("error", &_error, "error/D");
  }

  else
  {
    cout << _name << ": wrong file open type" << endl;
  }
}

bool QueryTree::Query(const int ipt, const int part,
    double &xpt, double &value, double &error)
{
  xpt = 0.;
  value = 0.;
  error = 0.;

  TSQLResult *res = _tree->Query("xpt:value:error", Form("ipt==%d&&part==%d",ipt,part));
  TSQLRow *row;
  if( row = res->Next() )
  {
    TString field0 = row->GetField(0);
    TString field1 = row->GetField(1);
    TString field2 = row->GetField(2);
    xpt = field0.Atof();
    value = field1.Atof();
    error = field2.Atof();
    delete row;
    delete res;
    return true;
  }

  cout << _name << ": no info for ipt " << ipt << " and part " << part << endl;
  delete res;
  return false;
}

void QueryTree::Fill(const int ipt, const int part,
    const double xpt, const double value, const double error)
{
  _ipt = ipt;
  _part = part;
  _xpt = xpt;
  _value = value;
  _error = error;
  _tree->Fill();
}

void QueryTree::Fill(const TH1 *h1, const TH1 *h2, const int part,
    const double kh1, const double kh2)
{
  int nbin = h1->GetXaxis()->GetNbins();
  if( nbin != h2->GetXaxis()->GetNbins() )
  {
    cout << _name << ": histograms have different bin numbers" << endl;
    return;
  }

  for(int ipt=0; ipt<nbin; ipt++)
  {
    double h1y = h1->GetBinContent(ipt+1) * kh1; 
    double h2y = h2->GetBinContent(ipt+1) * kh2; 
    double eh1y = h1->GetBinError(ipt+1) * kh1;
    double eh2y = h2->GetBinError(ipt+1) * kh2;
    if(h1y != 0. && h2y != 0.)
    {
      _ipt = ipt;
      _part = part;
      _xpt = h1->GetXaxis()->GetBinCenter(ipt+1);
      _value = h1y / h2y;
      _error = _value * sqrt( pow(eh1y/h1y,2) + pow(eh2y/h2y,2) );
      _tree->Fill();
    }
  }
}

void QueryTree::Write()
{
  _file->cd();
  _tree->Write();
  _file->Close();
}
