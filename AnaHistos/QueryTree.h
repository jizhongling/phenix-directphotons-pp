#include "GetEfficiency.h"

class QueryTree
{
  public:
    QueryTree(const char *fname, const char *option = "READ");
    void SetQuiet(bool quiet = true) { _quiet = quiet; }

    TSQLResult *Query(int part);
    bool Query(int ipt, int part,
        double &xpt, double &value, double &error);
    bool Query(int ipt, int part,
        double &xpt, double &value, double &errorlow, double &errorhigh);

    TGraphErrors *Graph(int part);
    TGraphAsymmErrors *GraphAsymm(int part);

    void Fill(int ipt, int part,
        double xpt, double value, double error);
    void Fill(int ipt, int part,
        double xpt, double value, double errorlow, double errorhigh);
    void Fill(const TH1 *h1, const TH1 *h2, int part,
        double kh1 = 1., double kh2 = 1.);

    TObject *Get(const char *namecycle);
    void cd(); 
    void Write(); 
    void Close(); 
    void Save(); 

  protected:
    bool _quiet;

    int _ipt;
    int _part;
    double _xpt;
    double _value;
    double _error;
    double _errorlow;
    double _errorhigh;

    TString _name;
    TFile *_file;
    TTree *_tree;
};

QueryTree::QueryTree(const char *fname, const char *option):
  _quiet(false),
  _name(fname)
{
  TString soption = option;
  soption.ToUpper();

  if(soption == "READ")
  {
    _file = new TFile(fname, option);
    _tree = (TTree*)_file->Get("t1");
  }

  else if(soption == "RECREATE")
  {
    _file = new TFile(fname, option);
    _tree = new TTree("t1", "t1 tree");
    _tree->Branch("ipt", &_ipt, "ipt/I");
    _tree->Branch("part", &_part, "part/I");
    _tree->Branch("xpt", &_xpt, "xpt/D");
    _tree->Branch("value", &_value, "value/D");
    _tree->Branch("error", &_error, "error/D");
    _tree->Branch("errorlow", &_errorlow, "errorlow/D");
    _tree->Branch("errorhigh", &_errorhigh, "errorhigh/D");
  }

  else
  {
    cout << _name << ": wrong file open mode" << endl;
  }
}

TSQLResult *QueryTree::Query(int part)
{
  TSQLResult *res = _tree->Query("ipt:xpt:value:error:errorlow:errorhigh", Form("part==%d",part));
  return res;
}

bool QueryTree::Query(int ipt, int part,
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

    if( !TMath::Finite(value + error) )
    {
      if(!_quiet)
        cout << _name << " at pt " << xpt << " in part " << part
          << ": not finite value = " << value << " and/or error = " << error << endl; 
      value = 0.;
      error = 0.;
      return false;
    }
    return true;
  }

  if(!_quiet)
    cout << _name << ": no info for ipt " << ipt << " and part " << part << endl;
  delete res;
  return false;
}

bool QueryTree::Query(int ipt, int part,
    double &xpt, double &value, double &errorlow, double &errorhigh)
{
  xpt = 0.;
  value = 0.;
  errorlow = 0.;
  errorhigh = 0.;

  TSQLResult *res = _tree->Query("xpt:value:errorlow:errorhigh", Form("ipt==%d&&part==%d",ipt,part));
  TSQLRow *row;
  if( row = res->Next() )
  {
    TString field0 = row->GetField(0);
    TString field1 = row->GetField(1);
    TString field2 = row->GetField(2);
    TString field3 = row->GetField(3);
    xpt = field0.Atof();
    value = field1.Atof();
    errorlow = field2.Atof();
    errorhigh = field3.Atof();
    delete row;
    delete res;

    if( !TMath::Finite(value + errorlow + errorhigh) )
    {
      if(!_quiet)
        cout << _name << " at pt " << xpt << " in part " << part
          << ": not finite value = " << value << " and/or errorlow = " << errorlow
          << " and/or errorhigh = " << errorhigh << endl; 
      value = 0.;
      errorlow = 0.;
      errorhigh = 0.;
      return false;
    }
    return true;
  }

  if(!_quiet)
    cout << _name << ": no info for ipt " << ipt << " and part " << part << endl;
  delete res;
  return false;
}

TGraphErrors *QueryTree::Graph(int part)
{
  TSQLResult *res = _tree->Query("xpt:value:error", Form("part==%d",part));
  int nrow = res->GetRowCount();
  if(nrow == 0 && !_quiet)
    cout << _name << ": no graph for part " << part << endl;

  TGraphErrors *graph = new TGraphErrors(nrow);
  int igp = 0;

  TSQLRow *row;
  while( row = res->Next() )
  {
    TString field0 = row->GetField(0);
    TString field1 = row->GetField(1);
    TString field2 = row->GetField(2);
    double xpt = field0.Atof();
    double value = field1.Atof();
    double error = field2.Atof();
    if( TMath::Finite(value + error) &&
        error > 0. )
    {
      graph->SetPoint(igp, xpt, value);
      graph->SetPointError(igp, 0., error);
      igp++;
    }
    delete row;
  }
  delete res;

  graph->Set(igp);
  return graph;
}

TGraphAsymmErrors *QueryTree::GraphAsymm(int part)
{
  TSQLResult *res = _tree->Query("xpt:value:errorlow:errorhigh", Form("part==%d",part));
  int nrow = res->GetRowCount();
  if(nrow == 0 && !_quiet)
    cout << _name << ": no graph for part " << part << endl;

  TGraphAsymmErrors *graph = new TGraphAsymmErrors(nrow);
  int igp = 0;

  TSQLRow *row;
  while( row = res->Next() )
  {
    TString field0 = row->GetField(0);
    TString field1 = row->GetField(1);
    TString field2 = row->GetField(2);
    TString field3 = row->GetField(3);
    double xpt = field0.Atof();
    double value = field1.Atof();
    double errorlow = field2.Atof();
    double errorhigh = field3.Atof();
    if( TMath::Finite(value + errorlow + errorhigh) &&
        errorlow > 0. && errorhigh > 0. )
    {
      graph->SetPoint(igp, xpt, value);
      graph->SetPointError(igp, 0., 0., errorlow, errorhigh);
      igp++;
    }
    delete row;
  }
  delete res;

  graph->Set(igp);
  return graph;
}

void QueryTree::Fill(int ipt, int part,
    double xpt, double value, double error)
{
  _ipt = ipt;
  _part = part;
  _xpt = xpt;
  _value = value;
  _error = error;
  _errorlow = error;
  _errorhigh = error;
  _tree->Fill();
}

void QueryTree::Fill(int ipt, int part,
    double xpt, double value, double errorlow, double errorhigh)
{
  _ipt = ipt;
  _part = part;
  _xpt = xpt;
  _value = value;
  _error = sqrt((errorlow*errorlow + errorhigh*errorhigh)/2.);
  _errorlow = errorlow;
  _errorhigh = errorhigh;
  _tree->Fill();
}

void QueryTree::Fill(const TH1 *h1, const TH1 *h2, int part,
    double kh1, double kh2)
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
    double Eff;

    _ipt = ipt;
    _part = part;
    _xpt = h1->GetXaxis()->GetBinCenter(ipt+1);
    _value = h1y / h2y;
    _error = _value * sqrt( pow(eh1y/h1y,2) + pow(eh2y/h2y,2) );
    if( !GetEfficiency(h2y, h1y, Eff, _errorlow, _errorhigh) )
    {
      _errorlow = _error;
      _errorhigh = _error;
    }

    if( TMath::Finite(_value + _error) ||
        TMath::Finite(_value + _errorlow + _errorhigh) )
      _tree->Fill(); 
  }
}

TObject *QueryTree::Get(const char *namecycle)
{
  return _file->Get(namecycle);
}

void QueryTree::cd()
{
  _file->cd();
}

void QueryTree::Write()
{
  _file->cd();
  _tree->Write();
}

void QueryTree::Close()
{
  _file->Close();
}

void QueryTree::Save()
{
  Write();
  Close();
}
