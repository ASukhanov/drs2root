// Make shape from a profile histogram
{
Double_t xlo = -1.;
Double_t xhi = 10.;

Int_t nbins = d4r->fprofile[0]->GetNbinsX();
Double_t binw = d4r->fprofile[0]->GetBinWidth(0);
Double_t shape_step = 0.1;
Int_t xlo_bin = (xlo - d4r->fprofile[0]->GetBinCenter(0))/binw;
Int_t xhi_bin = (xhi - d4r->fprofile[0]->GetBinCenter(0))/binw;
Int_t nsteps = (xhi - xlo) / shape_step + 1;

Int_t ii,xi;
Double_t xx, x0, yy;
TArrayD tax(nsteps),tay(nsteps);

for(ii=0;ii<nsteps;ii++) {
  xx = xlo + shape_step*ii;
  tax.SetAt(xx,ii);
  x0 = d4r->fprofile[0]->GetBinCenter(0);
  xi = (Int_t)((xx-x0)/binw);
  yy = d4r->fprofile[0]->GetBinContent(xi);
  tay.SetAt(yy,ii);
  //cout<<xx<<", "<<d4r->fprofile[0]->GetBinContent(xi)<<endl;
}
}
