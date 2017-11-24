template <class GraphType>
void ReadGraph(const char *name, Int_t igr, Double_t *gx, Double_t *gy, Double_t *egy)
{
  TFile *f = new TFile(name);
  GraphType *gr = (GraphType*)f->Get(Form("gr_%d",igr));
  for(int i=0; i<gr->GetN(); i++)
  {
    gr->GetPoint(i, gx[i], gy[i]);
    egy[i] = gr->GetErrorY(i);
  }
  f->Close();

  return;
}
