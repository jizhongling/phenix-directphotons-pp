template <class GraphType>
void ReadGraph(const char *name, Int_t igr, Double_t *gx, Double_t *gy, Double_t *egy)
{
  TFile *f = new TFile(name);
  GraphType *gr = (GraphType*)f->Get(Form("gr_%d",igr));
  gx = gr->GetX();
  gy = gr->GetY();
  egy = gr->GetEY();
  f->Close();

  return;
}
