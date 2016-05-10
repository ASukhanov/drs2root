{
  TH2S *anyh2s; // need this for gSystem->Load to work
  TTree *anytree; // need this for gSystem->Load to work
  
  cout<<"gSystem.Load(""drs4root"")"<<endl;
  gSystem->Load("drs4root");
  
  /**/
  drs4root* d4r = new drs4root(); 
  //assign options here
  d4r->gverb = 1; //0x2 | 0x4 | 0x8 ;// verbosity
  //d4r->gverb = 0xffff;
  //d4r->baseline_npoints = 1;// <2 to suppress, >0 to calculate single cell noise, default 10
  d4r->invert[0] = 1.; // negative pulse processing of ch0
  d4r->invert[1] = 1.; // negative pulse processing of ch1
  d4r->ghist = 0x1; //enable histogramming
  d4r->gthreshold_relative = 0.5;// use relative threshold level, if 0. then use absolute
  d4r->threshold = 0.2; // absolute threshold

  //if FILTERING defined
  //d4r->mf_shape = 1; // 0: no filtering, 1: from first event, 2: rectangle (moving average), 3: triangle
  //d4r->mf_size = 100;
  ////d4r->gfilter_roi_length = 0; // 0 to filter all cells, it is time consuming
  
  //d4r = new drs4root("r_0429_2040.dat");
  cout<<"Open data file using: d4r = new drs4root(filename)\n";
  cout<<"Functions:\n";
  cout<<"  d4r->Skip_events(i);\n";
  cout<<"  d4r->Next_event();\n";
  cout<<"  d4r->Print_header();\n";
  cout<<"  d4r->Print_stat();\n\n";
  cout<<" Analyze 100 events: for(ii=0;ii<100;ii++) d4r->Next_event();\n";
}
