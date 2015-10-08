{
  cout<<"gSystem.Load(""drs4root"")"<<endl;
  gSystem->Load("drs4root");
  drs4root* d4r = new drs4root();
  
  //assign options here
  d4r->gverb = 0;//0x2 | 0x4 | 0x8 ;// verbosity
  //d4r->baseline_npoints = 1;// <2 to suppress, >0 to calculate single cell noise, default 10
  d4r->invert[0] = 1.; // negative pulse processing of ch0
  d4r->invert[1] = 1.; // negative pulse processing of ch1

  //if FILTERING defined
  d4r->mf_shape = 0; // 0: no filtering, 1: from first event, 2: rectangle, 3: triangle
  d4r->mf_size = 10;
  
  d4r = new drs4root("~/data/1509301458.drs");
  cout<<"Functions: ";
  cout<<"d4r->Skip_events(i); ";
  cout<<"d4r->Next_event(); ";
  cout<<"d4r->Print_header(); ";
  cout<<"d4r->Print_stat(); ";
  cout<<endl;
}
