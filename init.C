{
    cout<<"gSystem.Load(""drs4root"")"<<endl;
    gSystem->Load("drs4root");
    drs4root* d4r = new drs4root();
    
    //assign options here
    d4r->gverb = 0; //0x2 | 0x4 | 0x8 ;
    d4r->invert[0] = 1.;
    d4r->invert[1] = 1.;
    
    d4r = new drs4root("~/data/1509301458.drs");
    cout<<"Functions: ";
    cout<<"d4r->Find_event(i); ";
    cout<<"d4r->Next_event(); ";
    cout<<"d4r->Print_header(); ";
    cout<<endl;
}
