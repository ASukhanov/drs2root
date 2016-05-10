// Version v2 2016-04-30. FILTERING disabled
// Version v3 2016-05-10. Reading files, written with drs-5.0.5

#ifndef drs4root_h
#define drs4root_h

#include <TFile.h>
#include <TH1S.h>
#include <TH2S.h>
#include <TProfile.h>
#include <TTree.h>
#include <TGraph.h>

//#define FILEFORMAT_PRE505 //v2/ uncomment this for reading files written using software prior to 505

//#define FILTERING

//#define FAST_CALIBRATION

enum
{
#ifdef FILTERING
  kMaxFilterLength = 100,
#endif
  kNCh = 2, //number of channels recorded in the file
  kLCh = 1024,
  kMax = 65536, // max amplitude
  kMin = 0
};

struct Event_Header_t
{
  char stamp_ev_header[4];
  UInt_t evnum;
  UShort_t year;
  UShort_t month;
  UShort_t day;
  UShort_t hour;
  UShort_t minute;
  UShort_t s;
  UShort_t ms;
  UShort_t range;
  char stamp_bno[2];
  UShort_t bno;
  char stamp_tno[2];
  UShort_t tno;
};
struct Ch_Amplitude_t
{
  char stamp_ch_header[4];
#ifndef FILEFORMAT_PRE505
  UInt_t scl; //v2/ Scaler in Hz
#endif
  UShort_t a[kLCh];
};
struct Event_t
{
  struct Event_Header_t h;
  struct Ch_Amplitude_t ch[kNCh];
};
struct Channel_Time_t
{
  char stamp_ch_header[4];
  Float_t t[kLCh];
};
struct Time_header_t
{
#ifndef FILEFORMAT_PRE505
  char file_header[4]; //v2
#endif
  char stamp_time_Header[4];
  char stamp_board_number[2];
  Short_t board_number;
  struct Channel_Time_t ch[kNCh];
};
  
class drs4root : public TFile
{
private:  
  int ndt;
  double sumdt, sumdt2;
  double sum_widths[kNCh];
  double asum[kNCh], asum2[kNCh];
  double bl_sum[kNCh], bl_sum2[kNCh];
  double waveform[kNCh][kLCh], ftime[kNCh][kLCh];
  long firstev_pos;
  //
  void EventMinMax();
  
public:
  static double threshold;
  static double invert[kNCh];
  static int baseline_npoints;
  static int bl_first;
  static int timing_calibration;
  static int regularize;
  static int gverb;
  static int ghist;
#ifdef FILTERING
  static int gfilter_roi_length;
#endif
  static double gthreshold_relative;

  struct Event_t fEv;
  struct Time_header_t fth;
  FILE    *fD;
  TFile   *ffile;
  Int_t   fsize;
  Char_t  fname[256];
  ULong_t fpos;
  Int_t   fevcount;
  UInt_t  fMin[kNCh], fMinPos[kNCh];
  UInt_t  fMax[kNCh], fMaxPos[kNCh];
  TTree   *ftree;
  Double_t fpeak[kNCh];
  Double_t fpeak_pos[kNCh];
  Double_t fLED_time[kNCh];
  Int_t   fpeak_idx[kNCh];
  Double_t fbaseline[kNCh];
  Double_t fone_cell[kNCh];
  TH2S  *fhch[kNCh];
  TProfile  *fprofile[kNCh]; // time-corrected profile of the waveform
  TGraph *fgraph;
   
#ifdef FILTERING
  Double_t fpkMF[kNCh], fpkpMF[kNCh];
  Int_t fpkiMF[kNCh];
  static int mf_shape, mf_size;
  double *mfilter_coeff,*mfilter_x;
  Int_t mf_idx_max;
  TH2S  *fhchMF[kNCh];
  //TH2S  *fhch_raw[kNCh];
#endif

  //drs4root();
  drs4root(const Char_t *in="", const Char_t *out="");
  ~drs4root();

  void Init();
  void Set_shape(Double_t *xx, Double_t *yy, Int_t size);
  TGraph *GetFilterGraph() {return fgraph;}; // graph of the filter coefficients
  Int_t Skip_events(Int_t ev);
  Int_t Next_event();
  void Print_stat();
  //void First_event(); //obsolete, just reuse init.C
  void Print_header();
  void Print_event();
  TProfile *GetProfile(Int_t ch) {return fprofile[ch];}
  TH2S *GetHistRaw(Int_t ch) {return fhch[ch];}
#ifdef  FILTERING
  TH2S *GetHistFiltered(Int_t ch) {return fhchMF[ch];}
  void SetFilter(Int_t shape=0, Int_t ntaps=10) {mf_shape = shape; mf_size = ntaps;}
#endif
  void Loop(Int_t nEvents=999999999);
  ClassDef(drs4root,0)
};
#endif
