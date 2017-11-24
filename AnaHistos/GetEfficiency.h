Bool_t GetEfficiency(Double_t total, Double_t passed, Double_t &Eff, Double_t &eLow, Double_t &eHigh)
{
  const Double_t level = 0.682689492137086;

  Eff = passed / total;
  if( total >= passed && passed >= 0. )
  {
    Double_t Low = TEfficiency::ClopperPearson(total, passed, level, kFALSE);
    Double_t Up = TEfficiency::ClopperPearson(total, passed, level, kTRUE);
    eLow = TMath::Abs( Low - Eff );
    eHigh = TMath::Abs( Up - Eff );
    return kTRUE;
  }

  return kFALSE;
}
