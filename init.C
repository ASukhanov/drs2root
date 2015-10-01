{
    cout<<"gSystem.Load(""drs2root"")"<<endl;
    gSystem->Load("drs2root");
    //drs2root *gdq=new drs2root();
    drs2root* d2r = new drs2root("~/data/1509301458.drs");
    cout<<"Functions: ";
    cout<<"d2r->Find_event(i); ";
    cout<<"d2r->Next_event(); ";
    cout<<"d2r->Print_header(); ";
    cout<<endl;
}
