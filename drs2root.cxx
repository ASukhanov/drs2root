#include <TROOT.h>
#include <TSystem.h>
//#include <Bytes.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/stat.h> // for stat()

using namespace std;
using namespace  TMath;

#include "drs2root.h"

char gcobuf[256];
char* strnz(const char* source, size_t num)
{ strncpy(gcobuf, source, num); gcobuf[num]=0; return & gcobuf[0];}
#define CO(obj) strnz(obj,sizeof(obj))

drs2root::~drs2root()
{}

drs2root::drs2root(const Char_t *in, const Char_t *out)
{
    Int_t rc,ii;
  const Char_t *tname;
    struct stat statv;
    Char_t      oname[256];

  if((tname = gSystem->ExpandPathName(in))!=0)
  {
    strcpy(fname,tname);
  }
  printf("Opening %s\n",fname);
    fD = fopen(fname,"rb");
    if(fD==NULL)
    {
        perror("Could not open file ");
        perror(fname);
        return;
    }
    rc = stat(fname,&statv);
    if(rc!=0) {perror("Cannot fstat"); fsize = -1;}
    fsize = statv.st_size;
    printf("File opened %s[%d]\n",fname,fsize);
    strcpy(oname,fname);
    char *substr = strrchr(oname,'.');
    strcpy(substr+1,"root");
    printf("Output file %s\n",oname);
    if(ffile) {printf("deleting file\n");delete ffile;}
    ffile = new TFile(oname,"recreate");
    if(ffile == NULL) {printf("ERROR. Could not open %s\n",oname); return;}
    printf("File opened\n");

    if (fread(&fTime_header,sizeof(fTime_header),1,fD) != 1)
      {printf("ERROR in Time Header\n"); return;}
  cout<<"Got fTime_header["<<sizeof(fTime_header)<<"]"<<endl;
  //cout<<strncpy(fTime_header.stamp_time_Header<<",";
  cout<<CO(fTime_header.stamp_time_Header)<<",";
  cout<<CO(fTime_header.stamp_board_number)<<fTime_header.board_number;
  cout<<endl;
  for (ii=0; ii<kNCh; ii++)
  {
    cout<<CO(fTime_header.channel_time[ii].stamp_ch_header)<<":";
    if(fTime_header.channel_time[ii].stamp_ch_header[0] != 'C')
    {  
      cout<<endl<<"WARNING. Only "<<ii;
      cout<<" channels present, change kNCh in heder file and recompile"<<endl;
      return;
    }
    cout<<fTime_header.channel_time[ii].t[0];
    cout<<endl;
  }
  return;
}

Int_t drs2root::Find_event(Int_t ev)
{
  fpos = sizeof(fTime_header) + ev*sizeof(fEv);
  fseek(fD,fpos,SEEK_SET);
  return 0;
}

void drs2root::EventMinMax()
{
  Int_t ii,jj;
  struct Ch_Amplitude_t* ch;

  for(ii=0;ii<kNCh;ii++)
  {
    fMin[ii]=kMax;
    fMax[ii]=kMin;    
    ch = &fEv.ch[ii];
    for(jj=0;jj<kLCh;jj++)
    {
      if(ch->a[jj]<fMin[ii]) {fMin[ii]=ch->a[jj];fMinPos[ii]=jj;}
      if(ch->a[jj]>fMax[ii]) {fMax[ii]=ch->a[jj];fMaxPos[ii]=jj;}
    }
  }
}
Int_t drs2root::Next_event()
{
  if (fread(&fEv,sizeof(fEv),1,fD) != 1)
    {return 1;}
  EventMinMax();  
  return 0;
}
void drs2root::Print_header()
{
  cout<<CO(fEv.h.stamp_ev_header)<<":";
  cout<<fEv.h.evnum<<",";
  cout<<"date:"<<fEv.h.year;
  cout<<"-"<<fEv.h.month;
  cout<<"-"<<fEv.h.day;
  cout<<" "<<fEv.h.hour;
  cout<<":"<<fEv.h.minute;
  cout<<":"<<fEv.h.s<<".";
  cout<<CO(fEv.h.stamp_tno)<<" "<<fEv.h.tno<<",";
  cout<<endl;
  Int_t ii;
  for(ii=0;ii<kNCh;ii++)
  {
    cout<<"Max["<<ii<<"]="<<fMax[ii]<<"\t@"<<fMaxPos[ii]<<endl;
    cout<<"Min["<<ii<<"]="<<fMin[ii]<<"\t@"<<fMinPos[ii]<<endl;
  }
}
void drs2root::Print_event()
{
  Int_t ii,jj;
  struct Ch_Amplitude_t* ch;
  for(ii=0;ii<kNCh;ii++)
  {
    ch = &fEv.ch[ii];
    cout<<CO(fEv.ch[ii].stamp_ch_header)<<":"<<endl;
    for(jj=0;jj<kLCh;jj++)
    {
      cout<<ch->a[jj]<<",";
    }
    cout<<endl;
  }
}
