int ULP_distance(const double ref, const double val, const int MAX)
{
  
  double tmp=val;
  int ULP=0;

  while( (ref!=tmp) && abs(ULP)<MAX ){ ULP++; tmp=nextafter(tmp,ref);}

  return ULP;
}