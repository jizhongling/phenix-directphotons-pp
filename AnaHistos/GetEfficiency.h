bool GetEfficiency(double total, double passed, double &Eff, double &eLow, double &eHigh)
{
  const double level = 0.682689492137086;

  Eff = passed / total;
  if( total >= passed && passed >= 0. )
  {
    double Low = TEfficiency::ClopperPearson(total, passed, level, false);
    double Up = TEfficiency::ClopperPearson(total, passed, level, true);
    eLow = TMath::Abs( Eff - Low );
    eHigh = TMath::Abs( Up - Eff );
    return true;
  }

  return false;
}
