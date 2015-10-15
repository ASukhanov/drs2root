void write_profile()
{
  Int_t ch = 0;
  TString fname("d4r_profile_");
  fname += ch;
  fname += ".root";
  TFile tf(fname,"RECREATE");
  
  d4r->GetProfile(0)->GetXaxis()->UnZoom();
  d4r->GetProfile(0)->GetYaxis()->UnZoom();
  d4r->GetProfile(0)->Write();
  tf.Close();
  cout<<"fprofile[0] is written to "<<fname<<endl;
}