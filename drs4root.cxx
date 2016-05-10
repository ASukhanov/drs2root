// Version v2 2016-04-30. FILTERING disabled by default. Tested at FNAL.

//#include <stdlib.h> // for exit()
#include <TROOT.h>
#include <TSystem.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/stat.h> // for stat()
#include <math.h>

using namespace std;
using namespace  TMath;

#include "drs4root.h"
#define kCellWidth 200./double(kLCh)

int drs4root::baseline_npoints = 10;
int drs4root::bl_first = 5;
int drs4root::timing_calibration = 1;
int drs4root::regularize = 0;
double drs4root::threshold = 0.2;
int drs4root::gverb = 0;
double drs4root::invert[kNCh] = {1.,1.};
int drs4root::ghist = 0;
double drs4root::gthreshold_relative = 0.;

#ifdef FILTERING
#include "mfilter.h"
int drs4root::mf_shape=0;
int drs4root::mf_size=0;
int drs4root::gfilter_roi_length=10;
#endif

int gverb = drs4root::gverb;

char gcobuf[256];
char* strnz(const char* source, size_t num)
{ strncpy(gcobuf, source, num); gcobuf[num]=0; return & gcobuf[0];}
#define CO(obj) strnz(obj,sizeof(obj))
//'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// class implementation 
//
double sigma(int nn, double sum, double sum2)
{ return sqrt(1./double(nn-1)*(sum2-sum*sum/double(nn)));}

void printstat(const char *s, int nn, double sum, double sum2, double factor, const char *unit)
{  printf("%s%.1f+-%.1f %s", s, factor*sum/double(nn), factor*sigma(nn,sum,sum2), unit);}

void drs4root::Print_stat()
{
  // print statistics
  printf("ev#%5i,",fevcount);
  printstat(" dT=",ndt,sumdt,sumdt2,1000.,"ps"); 
  printstat(", A0=",fevcount,asum[0],asum2[0],1000.,"mV");
  printstat(", A1=",fevcount,asum[1],asum2[1],1000.,"mV");
  if(baseline_npoints > 0)
  {
    printstat(", BL0=",fevcount,bl_sum[0],bl_sum2[0],1000.,"mV");
    printstat(", BL1=",fevcount,bl_sum[1],bl_sum2[1],1000.,"mV");
  }
  printf("\n");
}

drs4root::~drs4root()
{}

drs4root::drs4root(const Char_t *in, const Char_t *out)
{
    Int_t rc,ii,ch;
  const Char_t *tname;
    struct stat statv;
    Char_t      oname[256];
    TFile *ffile = 0;

  if(strlen(in)==0) return;
  if((tname = gSystem->ExpandPathName(in))!=0) strcpy(fname,tname);
  printf("Opening %s\n",fname);
  fD = fopen(fname,"rb");
  if(fD==NULL)
  {
      perror("Could not open file ");
      perror(fname);
      return;
  }
  // input file
  rc = stat(fname,&statv);
  if(rc!=0) {perror("Cannot fstat"); fsize = -1;}
  fsize = statv.st_size;
  printf("File opened %s[%d]\n",fname,fsize);
#ifdef FILEFORMAT_PRE505
  printf("WARNING! unpacking assuming file format written with software version prior 5.0.5\n");
#endif

  // output file
  strcpy(oname,fname);
  char *substr = strrchr(oname,'.');
  strcpy(substr+1,"root");
  printf("Output file %s\n",oname);
  if(ffile) {printf("deleting file\n");delete ffile;}
  ffile = new TFile(oname,"recreate");
  if(ffile == NULL) {printf("ERROR. Could not open %s\n",oname); return;}
  printf("File opened\n");

  // print options
#ifdef FAST_CALIBRATION
  printf("Using fast calibration\n");
#endif 
  for(ch=0;ch<kNCh;ch++)
  {
    if(invert[ch]==1.) printf("Negative");
      else printf("Negative");
    printf(" pulse processing of ch %i\n",ch);
  }
  if(gthreshold_relative != 0.) printf("Threshold is relative = %.1f of the peak\n", gthreshold_relative);
  else printf("Threshold is absolute = %.1fV\n",threshold);
  printf("Verbosity=%i\n",gverb);
  if(baseline_npoints>1) 
    printf("Baseline subtraction for %i points starting at %i\n",baseline_npoints,bl_first);
#ifdef FILTERING
  if(mf_size > kMaxFilterLength) mf_size = kMaxFilterLength;
  if(mf_size && mf_shape) printf("Filtering of type %i[%i] enabled\n",mf_shape,mf_size);
  if(gfilter_roi_length) printf("Filter ROI is %i cells covering the peak\n",gfilter_roi_length);
  else printf("Filter ROI is full area\n");
#ifdef FILTER_NONUNIFORM
  printf("Non-Uniform-Sampled filtering\n");
#endif
#endif 
  if(regularize) printf("Cell Regularization is on!\n");
  
  // read time header
  if (fread(&fth,sizeof(fth),1,fD) != 1)
      {printf("ERROR in Time Header\n"); return;}
  cout<<"Got fth["<<sizeof(fth)<<"]"<<endl;
  //cout<<strncpy(fth.stamp_time_Header<<",";
  cout<<"Time header: "<<CO(fth.stamp_time_Header)<<",";
  if(fth.stamp_time_Header[0] != 'T')
  {
    cout<<"\nERROR. File format wrong. Time header should be 'TIME'. Change FILEFORMAT_PRE505 in drs4root.h accordingly and recompile\n";
    return;
  }
  cout<<CO(fth.stamp_board_number)<<fth.board_number;
  cout<<endl;
  for (ch=0; ch<kNCh; ch++)
  {
    cout<<CO(fth.ch[ch].stamp_ch_header)<<":";
    if(fth.ch[ch].stamp_ch_header[0] != 'C')
    {  
      cout<<endl<<"WARNING. Only "<<ch;
      cout<<" channels present, change kNCh in header file and recompile"<<endl;
      // event header found
      fseek(fD, -4, SEEK_CUR);
      return;
    }
    cout<<fth.ch[ch].t[0];
    cout<<endl;
    for(ii=0, sum_widths[ch]=0.; ii<kLCh; ii++) sum_widths[ch] += fth.ch[ch].t[ii];
    printf("Found timing calibration for channel #%d, sum_widths=%f\n", ch+1, sum_widths[ch]);
    if(timing_calibration==0) for (ii=0; ii<kLCh; ii++) fth.ch[ch].t[ii] = kCellWidth;
#ifdef FAST_CALIBRATION
    for (ii=1; ii<kLCh-1; ii++) fth.ch[ch].t[ii] += fth.ch[ch].t[ii-1];
    for (ii=kLCh-1; ii>0; ii--) fth.ch[ch].t[ii] = fth.ch[ch].t[ii-1];
    fth.ch[ch].t[0] = 0.;
#endif
    if(gverb&8) {for(ii=0;ii<kLCh;ii++) printf("%i:%03i ",ii,int(fth.ch[ch].t[ii]*1000.)); printf("\n");}
  }
  // create histograms\ and tree
  TString hname;
  if(ghist)
  {
    for(ch=0;ch<kNCh;ch++)
    {
      hname = "hch"; hname += ch;
      //cout<<"Creating hist d4r->fhch["<<ch<<"]"<<endl;
      fhch[ch] = new TH2S(hname,hname,kLCh*2,0,200.,4000,-0.1,1.);
      hname = "d4r_profile_"; hname += ch;
      fprofile[ch] = new TProfile(hname,hname,kLCh*10,-50.,50.,-0.1,1.);
    }
    ftree = new TTree("tree","drs2root tree");
    if(ftree==NULL){printf("failed to create tree\n");return;}
    ftree->Branch("peak",&fpeak,"fpeak[2]/D");
    ftree->Branch("peakpos",&fpeak_pos,"fpeak_pos[2]/D");
    ftree->Branch("LED_time",&fLED_time,"ftime[2]/D");
    ftree->Branch("baseline",&fbaseline,"fbaseline[2]/D");
    ftree->Branch("onecell",&fone_cell,"fone_cell[2]/D");
  #ifdef FILTERING
    ftree->Branch("peakMF",&fpkMF,"fpkMF[2]/D");
    ftree->Branch("peakposMF",&fpkpMF,"fpkpMF[2]/D");
    for(ch=0;ch<kNCh;ch++)
    {
      hname = "hchMF"; hname += ch;
      fhchMF[ch] = new TH2S(hname,hname,kLCh*2,0,200.,4000,-0.1,2.);
    }
  #endif    
  }
  Init();
  return;
}

Int_t drs4root::Skip_events(Int_t ev)
{
  fpos = sizeof(fth) + ev*sizeof(fEv);
  fseek(fD,fpos,SEEK_CUR);
  return 0;
}
void drs4root::Init()
{
  // initialize statistics
  int ch;
  fevcount = 0;
  ndt = 0;
  sumdt = sumdt2 = 0.;
  for(ch=0;ch<kNCh;ch++) 
  {
    fpeak[ch]=0.; 
    asum[ch]=0.; 
    asum2[ch]=0.; 
    bl_sum[ch]=0.; 
    bl_sum2[ch]=0.;
  }
#ifdef FILTERING
  mfilter_delete(); // to free the allocated memory
  mfilter_coeff = 0;
  //for(ch=0;ch<kNCh;ch++) {peak[ch]=0.; peakpos[ch]=0.;};
#endif
}

//void drs4root::First_event() { fseek(fD,firstev_pos,SEEK_SET); }

void drs4root::EventMinMax()
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
//,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#ifdef  FILTERING
//'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
void drs4root::Set_shape(Double_t *xx, Double_t *yy, Int_t size)
{
  double mf_max;
  Int_t ii;
  //double xf[kMaxFilterLength];
  //for(ii=0;ii<kMaxFilterLength;ii++) xf[ii] = double(ii)*kCellWidth;

  //'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  mf_size=mfilter_create(xx,yy,size,kCellWidth,mf_shape,mf_size,&mfilter_coeff,&mfilter_x);
  //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  if(mf_size > kMaxFilterLength) 
  {
    printf("ERROR, filter size %i > %i\n",mf_size,kMaxFilterLength);
    gSystem->Abort(1);
  }
  if(mf_size)
  {
    //for(ii=0;ii<mf_size+1;ii++) printf("->mfi_set(%.3f,%.3f)\n",mfilter_x[ii],mfilter_coeff[ii]);
    mfi_set(mfilter_x,mfilter_coeff,mf_size);
    printf("MFilter of type %i[%i] created:\n",mf_shape,mf_size);
    double l2=0.;
    for(ii=0;ii<mf_size;ii++) 
    {
      printf("%i:fx=%f,fy%f\n",ii,mfilter_x[ii],mfilter_coeff[ii]);
      l2 += mfilter_coeff[ii]*mfilter_coeff[ii];
    }
    mf_max=0.; 
    mf_idx_max=0;
    for(ii=0;ii<mf_size;ii++) if(mfilter_coeff[ii]>mf_max) {mf_max=mfilter_coeff[ii]; mf_idx_max=ii;}
    // The L2-normalized filter improves the signal/noise by factor of 1./mf_max. 
    printf("MFilter max=%f @ %i. Expected signal/noise improvement %.2f. L2=%f\n",mf_max,mf_idx_max,1./mf_max,l2);
    fgraph = new TGraph(mf_size,mfilter_x,mfilter_coeff);
  }
  else printf("Failed to create mfilter\n");
}
//,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#endif
//'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Int_t drs4root::Next_event()
{
  int ii,ch;
  double t1, t2, dt; 
  
  if (fread(&fEv,sizeof(fEv),1,fD) != 1)
    {return 1;}
  //EventMinMax();
  //return 0;  
  if (gverb&(1)) printf("sn=%i evn=%i range=%i\n",fevcount,fEv.h.evnum,fEv.h.range);
  fevcount++;
  for(ch=0;ch<kNCh;ch++)
  {
#ifdef FILEFORMAT_PRE505
    if (gverb&(2)) printf("ch%i:%4.4s trig cell %i timecal=%f",ch,fEv.ch[ch].stamp_ch_header,fEv.h.tno,fth.ch[ch].t[fEv.h.tno]);
#else
    if (gverb&(2)) printf("ch%i:%4.4s trig cell %i timecal=%f, scaler=%i\n",ch,fEv.ch[ch].stamp_ch_header,fEv.h.tno,fth.ch[ch].t[fEv.h.tno],fEv.ch[ch].scl);
#endif
    fpeak[ch] = 0.;
    for (ii=0 ; ii<1024 ; ii++) {
      // convert data to volts
      waveform[ch][ii] = (fEv.ch[ch].a[ii] / 65536. * (1.-2.*invert[ch]) + invert[ch] + fEv.h.range/1000.0 - 0.5);
      //waveform[ch][ii] = fEv.ch[ch].a[ii]/65536.;
      if(waveform[ch][ii]>fpeak[ch]) 
      {  fpeak[ch]=waveform[ch][ii]; fpeak_idx[ch]=ii;}
      if (gverb&8) printf("%i:%03i ",ii,int(waveform[ch][ii]*1000.));            
      // calculate time for this cell
  #ifdef FAST_CALIBRATION
      ftime[ch][ii] = fth.ch[ch].t[(ii+fEv.h.tno) % 1024] 
        - fth.ch[ch].t[fEv.h.tno];
      if (ftime[ch][ii] <0.) ftime[ch][ii] += sum_widths[ch];
  # else
      int jj;
      for (jj=0,ftime[ch][ii]=0 ; jj<ii ; jj++)
          ftime[ch][ii] += fth.ch[ch].t[(jj+fEv.h.tno) % 1024];
  #endif
      if (gverb&8) printf("%03i ",int(ftime[ch][ii]*1000.)); 
    }
    fpeak_pos[ch]=ftime[ch][fpeak_idx[ch]];
    if (gverb&(2|8|0x100)) 
      printf("ch%i, max %.3f @ %i/%.3f\n",ch,fpeak[ch],fpeak_idx[ch],fpeak_pos[ch]);
  }
  
  // align cell #0 of all channels
  t1 = ftime[0][(kLCh-fEv.h.tno) % kLCh];
  for (ch=1 ; ch<kNCh ; ch++) {
      t2 = ftime[ch][(kLCh-fEv.h.tno) % kLCh];
      dt = t1 - t2;
      if (gverb&8) printf("ch%i trig@%i, %f-%f, aligned for %f\n",ch,fEv.h.tno,t1,t2,dt);
      for (ii=0 ; ii<kLCh ; ii++) ftime[ch][ii] += dt;
      fpeak_pos[ch]=ftime[ch][fpeak_idx[ch]];
  }
  // baseline subtraction
  if(baseline_npoints > 1)
    for(ch=0;ch<kNCh;ch++)
    {
      for(fbaseline[ch]=0., ii=bl_first; ii<bl_first+baseline_npoints; ii++)
        fbaseline[ch] += waveform[ch][ii];
      fbaseline[ch] /= double(baseline_npoints);        
      for(ii=0;ii<kLCh;ii++) waveform[ch][ii] -= fbaseline[ch];
    }
  
  // histogram raw waveforms
  for(ch=0;ch<kNCh;ch++)
  if(fhch[ch]) fhch[ch]->FillN(kLCh,ftime[ch],waveform[ch],(double* )NULL);

  // calculate single cell noise
  if(baseline_npoints > 0)
    for (ch=0;ch<2;ch++)
    {
      fone_cell[ch] = waveform[ch][bl_first];
      bl_sum[ch] += fone_cell[ch]; bl_sum2[ch] += fone_cell[ch]*fone_cell[ch];
    }
    
  // recalculate amplitudes of first two channels by averaging around peak
  //int half_top_width=2;
  //double na = double(half_top_width*2+1);
  for (ch=0 ; ch<2 ; ch++)
  {
    //fpeak[ch]=0.;
    //for (ii=fpeak_idx[ch]-half_top_width; ii<=fpeak_idx[ch]+half_top_width; ii++) fpeak[ch] += waveform[ch][ii];
    //fpeak[ch] /=na;
    asum[ch] += fpeak[ch];
    asum2[ch] += fpeak[ch]*fpeak[ch];
    //printf("wma[%i]=%.3f @ %i/%.3f\n",ch,fpeak[ch],fpeak_idx[ch],fpeak_pos[ch]);
  }
  
  // DLED - Digital Leading Edge Discrimination
  //fLED_time[ch][0] = fLED_time[ch][1] = 0.;
  for (ch=0;ch<2;ch++)
  {
    if(gthreshold_relative != 0.) threshold = fpeak[ch]*gthreshold_relative;
    if (gverb&4) printf("Looking for threshold=%f crossing, trigcell %i\n",
      threshold,fEv.h.tno);
    fLED_time[ch] = 0.;
    for (ii=0 ; ii<1022 ; ii++)
        if (waveform[ch][ii] < threshold && waveform[ch][ii+1] >= threshold) 
        {
          fLED_time[ch] = (threshold-waveform[ch][ii])/(waveform[ch][ii+1]-waveform[ch][ii])*(ftime[ch][ii+1]-ftime[ch][ii])+ftime[ch][ii];
          if (gverb&4) printf("Cross[%i]: %f,%f @ %i,%f,%f t=%f\n", 
            ch, waveform[ch][ii], waveform[ch][ii+1], ii, ftime[ch][ii], ftime[ch][ii+1], ftime[ch][ch]);
          break;
        }
  }
  
  // histogram of the corrected timing
  for (ch=0;ch<2;ch++)
  {
    if(fprofile[ch])
      for(ii=0;ii<kLCh;ii++) fprofile[ch]->Fill(ftime[ch][ii]-fLED_time[ch],waveform[ch][ii]);
  }
  
  // calculate distance of peaks with statistics
  if (fLED_time[0] > 0. && fLED_time[1] > 0.) {
      ndt++;
      dt = fLED_time[1] - fLED_time[0];
      sumdt += dt;
      sumdt2 += dt*dt;
  }
    
#ifdef FILTERING
    // use the first channel of first event to build the matching filter
    if (fevcount==1 && mf_shape==kFShape_1st_event)                   
      Set_shape(ftime[0],waveform[0],kLCh);
    
    // FIR filtering
    //printf("mf_size=%i\n",mf_size);
    int roi_start; // beginning of the region of interest for filtering
    int roi_length;
    double *ob[kNCh];
    if (mf_size && mf_shape)
    for(ch=0;ch<2;ch++)
    {
      ob[ch] = waveform[ch];
      roi_length = gfilter_roi_length;
      roi_start = 0;
      #ifdef FILTER_NONUNIFORM
        if(roi_length == 0)  roi_length = kLCh - mf_size;
        else       roi_start = fpeak_idx[ch] - mf_idx_max - roi_length/2;
        //ob[ch] = NULL;// don't want it
        ob[ch] = &(waveform[ch][roi_start]);
      #endif
      if(gverb&0x100) printf("initial fpeak=%.3f @ %i/%.3f\n",fpeak[ch],fpeak_idx[ch],fpeak_pos[ch]);
      //in-place filter
      //printf("->filter#%i(%i[%i], %.3f,%.3f...)mf_idx_max=%i\n",ch,roi_start,roi_length,waveform[ch][roi_start],ftime[ch][roi_start],mf_idx_max);
      //'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      mfilter_filter(&(ftime[ch][roi_start]),&(waveform[ch][roi_start]),roi_length,ob[ch],&(fpkMF[ch]),&(fpkpMF[ch]),&fpkiMF[ch]);
      //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      fpkiMF[ch] += roi_start;// + mf_idx_max;
      
      if(gverb&0x40){for(ii=0;ii<kLCh;ii++) printf("%i,%.2f:%.4f ",ii,ftime[ch][ii],waveform[ch][ii]);}  
      if(gverb&0x100) printf("filtered peak=%.3f @ %i/%.3f\n",fpkMF[ch],fpkiMF[ch],fpkpMF[ch]);
      for (ch=0;ch<kNCh;ch++)
        if(fhchMF[ch]) fhchMF[ch]->FillN(kLCh,ftime[ch],waveform[ch],(double* )NULL);
    }
    if(gverb&0x100) 
      printf("dt=%f\n",fpkpMF[1] - fpkpMF[0]);
#endif
  
  // fill the tree
  if(ftree) ftree->Fill();
  if((fevcount%1000)==999)
    Print_stat();
  return 0;
}
void drs4root::Print_header()
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
  /*Int_t ii;
  for(ii=0;ii<kNCh;ii++)
  {
    cout<<"Max["<<ii<<"]="<<fMax[ii]<<"\t@"<<fMaxPos[ii]<<endl;
    cout<<"Min["<<ii<<"]="<<fMin[ii]<<"\t@"<<fMinPos[ii]<<endl;
  }*/
}
void drs4root::Print_event()
{
  Int_t ii,jj;
  struct Ch_Amplitude_t* ch;
  //Print_stat();
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
void drs4root::Loop(Int_t nn)
{
  fevcount = 0;
  Int_t ii;
  for(ii=0;ii<nn;ii++) if(Next_event()) break;
  return;
}
