// Make shape from a profile histogram

void get_shape_from_saved_profile()
{
  Double_t xlo = -1.;
  Double_t xhi = 10.;
  Double_t kCellWidth = 200./1024;
  Double_t shape_step = kCellWidth;

  TFile tf("d4r_profile_0.root");
  TProfile *tp;  
  tp = (TProfile*) tf.Get("d4r_profile_0");
  cout<<"got profile:"<<endl;
  tp->Print();  
  int imax = tp->GetMaximumBin();
  double ymx = tp->GetMaximum();
  cout<<"max="<<ymx<<" @ "<<imax<<endl;
  
  Int_t nbins = tp->GetNbinsX();
  Double_t binw = tp->GetBinWidth(0);
  Double_t x0 = tp->GetBinCenter(0);
  Int_t xlo_bin = (xlo - x0)/binw;
  Int_t xhi_bin = (xhi - x0)/binw;
  Int_t nsteps = (xhi - xlo) / shape_step + 1;

  Int_t ii,xi;
  Double_t xx[1000], x0, yy[1000];

  //cout<<"loop from bin "<<xlo_bin<<" to "<<xhi_bin<<endl;
  
  //sample such way that that max point is sampled
  double xmx = imax*binw - (xlo - x0);
  double offset = (xmx/shape_step - floor(xmx/shape_step))*shape_step;
  //cout<<"xmx="<<xmx<<" @ "<<imax<<" ofs="<<offset<<" bw="<<binw<<" st"<<shape_step<<" "<<shape_step/binw<<endl;
  for(ii=0; ii<nsteps && ii<1000; ii++) {
    xx[ii] = shape_step*ii + offset;
    xi = (Int_t)((xx[ii] + xlo - x0)/binw);
    yy[ii] = tp->GetBinContent(xi);
    //cout<<"y="<<tp->GetBinContent(xi)<<" #"<<xi<<" @ "<<xx[ii]<<endl;
  }
  d4r->SetFilter(4,100);
  cout<<"Embedding shape type "<<d4r->mf_shape
    <<" max="<<yy[TMath::LocMax(nsteps,yy)]<<endl;
  d4r->Set_shape(xx, yy, nsteps);
  tf.Close(); cout<<"file closed"<<endl;
}
