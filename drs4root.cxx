#include <TROOT.h>
#include <TSystem.h>
//#include <Bytes.h>

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
#ifdef FILTERING
#include "mfilter.h"
int drs4root::mf_shape=0;
int drs4root::mf_size=0;
#endif

int gverb = drs4root::gverb;

char gcobuf[256];
char* strnz(const char* source, size_t num)
{ strncpy(gcobuf, source, num); gcobuf[num]=0; return & gcobuf[0];}
#define CO(obj) strnz(obj,sizeof(obj))

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
  printf("Threshold=%f\n",threshold);
  printf("Verbosity=%i\n",gverb);
  if(baseline_npoints>1) 
    printf("Baseline subtraction for %i points starting at %i\n",baseline_npoints,bl_first);
#ifdef FILTERING
  if(mf_size > kMaxFilterLength) mf_size = kMaxFilterLength;
  if(mf_size && mf_shape) printf("Filtering of type %i[%i] enabled\n",mf_shape,mf_size);
#endif 
  if(regularize) printf("Cell Regularization is on!\n");
  
  // read time header
  if (fread(&fth,sizeof(fth),1,fD) != 1)
      {printf("ERROR in Time Header\n"); return;}
  cout<<"Got fth["<<sizeof(fth)<<"]"<<endl;
  //cout<<strncpy(fth.stamp_time_Header<<",";
  cout<<CO(fth.stamp_time_Header)<<",";
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
  //firstev_pos = ftell(fD);
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
    wmax[ch]=0.; 
    asum[ch]=0.; 
    asum2[ch]=0.; 
    bl_sum[ch]=0.; 
    bl_sum2[ch]=0.;
  }
#ifdef FILTERING
  mfilter_delete(); // to free the allocated memory
  mfilter_coeff = 0;
  for(ch=0;ch<kNCh;ch++) {peak[kNCh]=0.; peakpos[kNCh]=0.;};
#endif
  if(ghist)
  {
    TString hname;
    //TString htitle;
    for(ch=0;ch<kNCh;ch++)
    {
      hname = "hch"; hname += ch;
      cout<<"Creating hist d4r->fhch["<<ch<<"]"<<endl;
      fhch[ch] = new TH2S(hname,hname,kLCh*2,0,200.,4000,-0.1,2.);
      //fhch[ch]->Print();
      hname = "hpos"; hname += ch;
      fhpos[ch] = new TH1S(hname,hname,1000,0.,100.);
    }
    fhdt = new TH1S("hdt","dT using FIR",2000,-10.,10.);
    fhdled = new TH1S("hdled","dT using Leading Edge Discriminator",2000,-10.,10.);
  }
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
Int_t drs4root::Next_event()
{
  int ii,ch;
  double t1, t2, tt[kNCh], dt, dv, v; 
  
  if (fread(&fEv,sizeof(fEv),1,fD) != 1)
    {return 1;}
  //EventMinMax();
  //return 0;  

  fevcount++;
  for(ch=0;ch<kNCh;ch++)
  {
    if (gverb&(2|8)) printf("ch%i trig cell %i timecal=%f\n",ch,fEv.h.tno,fth.ch[ch].t[fEv.h.tno]);
    for (ii=0 ; ii<1024 ; ii++) {
      // convert data to volts
      waveform[ch][ii] = (fEv.ch[ch].a[ii] / 65536. * (1.-2.*invert[ch]) + invert[ch] + fEv.h.range/1000.0 - 0.5);
      if(waveform[ch][ii]>wmax[ch]) {wmax[ch]=waveform[ch][ii]; cmax[ch]=ii;}
      if (gverb&2) printf("%i:%03i ",ii,int(waveform[ch][ii]*1000.));            
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
    if (gverb&(2|8)) printf("\nch%i, max %f @ %i\n",ch,wmax[ch],cmax[ch]);
  }
  
#ifdef FILTERING
  // use the first channel of first event to build the arbitr1.ary matching filter
  if (fevcount==1) 
  {
    mf_size=mfilter_create(waveform[0],kLCh,mf_shape,mf_size,&mfilter_coeff);
    if(mf_size)
    {
      printf("MFilter of type %i[%i] created:\n",mf_shape,mf_size);
      double l2=0.;
      for(ii=0;ii<mf_size;ii++) 
      {
        printf("%i:%f \n",ii,mfilter_coeff[ii]);
        l2 += mfilter_coeff[ii]*mfilter_coeff[ii];
      }
      double mf_max=0.; 
      int mf_pos=0;
      for(ii=0;ii<mf_size;ii++) if(mfilter_coeff[ii]>mf_max) {mf_max=mfilter_coeff[ii]; mf_pos=ii;}
      printf("MFilter max=%f @ %i. L2=%f\n",mf_max,mf_pos,l2);
    }
  }
#endif
  // align cell #0 of all channels
  t1 = ftime[0][(kLCh-fEv.h.tno) % kLCh];
  for (ch=1 ; ch<kNCh ; ch++) {
      t2 = ftime[ch][(kLCh-fEv.h.tno) % kLCh];
      dt = t1 - t2;
      if (gverb&8) printf("ch%i trig@%i, %f-%f, aligned for %f\n",ch,fEv.h.tno,t1,t2,dt);
      for (ii=0 ; ii<kLCh ; ii++) ftime[ch][ii] += dt;
  }
  // baseline subtraction
  if(baseline_npoints > 1)
    for(ch=0;ch<kNCh;ch++)
    {
      for(v=0., ii=bl_first; ii<bl_first+baseline_npoints; ii++)
        v += waveform[ch][ii];
      v /= double(baseline_npoints);        
      for(ii=0;ii<kLCh;ii++) waveform[ch][ii] -= v;
    }
  // regularization of the cell intervals to make them equidistant
  double time_shift, expected_time, prev;
  if(regularize)
  {
    for(ch=0;ch<2;ch++)
      for(prev=waveform[ch][0],ii=1;ii<kLCh;ii++)
      {
        dt=ftime[ch][ii]-ftime[ch][ii-1];
        expected_time = double(ii)*kCellWidth;
        time_shift = ftime[ch][ii] - expected_time;
        ftime[ch][ii] = expected_time;
        dv=waveform[ch][ii]-prev;
        prev = waveform[ch][ii];
        waveform[ch][ii] +=  dv/dt*time_shift;
      }
  }            
#ifdef FILTERING
  // FIR filtering
  //printf("mf_size=%i\n",mf_size);
  if (mf_size && mf_shape)
  for(ch=0;ch<2;ch++)
  {
    //in-place filter
    mfilter_filter(waveform[ch],kLCh,waveform[ch],&(peak[ch]),&(peakpos[ch]));
    //recalculate max
    for(ii=0;ii<kLCh-mf_size;ii++)
      if(waveform[ch][ii]>wmax[ch]) {wmax[ch]=waveform[ch][ii]; cmax[ch]=ii;}
    if(gverb&0x40){for(ii=0;ii<kLCh;ii++) printf("%i,%.2f:%.4f ",ii,ftime[ch][ii],waveform[ch][ii]);  printf("\nmax=%f @ %i\n",wmax[ch],cmax[ch]);}
    // calculate peak position using 3-point quadratic approximation
    double dymaxl, dymaxr, dx;
    dymaxl = waveform[ch][cmax[ch]] - waveform[ch][cmax[ch]-1];
    dymaxr = waveform[ch][cmax[ch]+1] - waveform[ch][cmax[ch]];
    dx = ftime[ch][cmax[ch]] - ftime[ch][ii-1];
    if(dymaxl-dymaxr == 0.) {cout<<"Computational error!\n";}
    fpeak_pos[ch] = dymaxl/(dymaxl-dymaxr)/dx + double(ftime[ch][cmax[ch]]); // TODO add the filter width
    if(gverb&0x200) 
      printf("ch%i: %.3f,%.3f-%.3f,%.3f-%.3f,%.3f, pos=%4f\n",ch,
             waveform[ch][cmax[ch]-1],ftime[ch][cmax[ch]-1],
             waveform[ch][cmax[ch]],ftime[ch][cmax[ch]],
             waveform[ch][cmax[ch]+1],ftime[ch][cmax[ch]+1],fpeak_pos[ch]);
  }
  if(gverb&0x200) printf("dt=%f\n",fpeak_pos[1] - fpeak_pos[0]);
#endif

  // calculate single cell noise
  if(baseline_npoints > 0)
    for (ch=0;ch<2;ch++)
    {
      v = waveform[ch][bl_first];
      bl_sum[ch] += v; bl_sum2[ch] += v*v;
    }
  // estimate amplitudes of first two channels by averaging around peak
  int half_top_width=2;
  double na = double(half_top_width*2+1);
  for (ch=0 ; ch<2 ; ch++)
  {
    wmax[ch]=0.;
    for (ii=cmax[ch]-half_top_width; ii<=cmax[ch]+half_top_width; ii++) wmax[ch] += waveform[ch][ii];
    wmax[ch] /=na;
    asum[ch] += wmax[ch];
    asum2[ch] += wmax[ch]*wmax[ch];
  }
  // DLED - Digital Leading Edge Discrimination
  if (gverb&4) printf("Looking for threshold=%f crossing, trigcell %i\n",
    threshold,fEv.h.tno);
  tt[0] = tt[1] = 0.;
  for (ch=0;ch<2;ch++)
  {
    for (ii=0 ; ii<1022 ; ii++)
        if (waveform[ch][ii] < threshold && waveform[ch][ii+1] >= threshold) 
        {
          tt[ch] = (threshold-waveform[ch][ii])/(waveform[ch][ii+1]-waveform[ch][ii])*(ftime[ch][ii+1]-ftime[ch][ii])+ftime[ch][ii];
          if (gverb&4) printf("Cross[%i]: %f,%f @ %i,%f,%f t=%f\n", 
            ch, waveform[ch][ii], waveform[ch][ii+1], ii, ftime[ch][ii], ftime[ch][ii+1], tt[ch]);
          break;
        }
  }
  // calculate distance of peaks with statistics
  if (tt[0] > 0 && tt[1] > 0) {
      ndt++;
      dt = tt[1] - tt[0];
      sumdt += dt;
      sumdt2 += dt*dt;
  }
  // histogram amplitudes
  if(ghist)
  {
    for (ch=0;ch<kNCh;ch++)
    {
      //printf("filling %i\n", ch);
      if(fhch[ch]) 
      {
        fhch[ch]->FillN(kLCh,ftime[ch],waveform[ch],(double* )NULL);
        //for(ii=0;ii<kLCh;ii++) fhch[ch]->Fill(ftime[ch][ii],waveform[ch][ii]);
      }
      fhpos[ch]->Fill(fpeak_pos[ch]);
    }
    fhdt->Fill(fpeak_pos[1]-fpeak_pos[0]);
    fhdled->Fill(dt);
  }
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
