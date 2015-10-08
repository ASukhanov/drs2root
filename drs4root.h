#ifndef drs4root_h
#define drs4root_h

#include <TFile.h>

#define FILTERING
#define FAST_CALIBRATION

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
  double wmax[kNCh];
  int cmax[kNCh];
  double asum[kNCh], asum2[kNCh];
  double bl_sum[kNCh], bl_sum2[kNCh];
  double waveform[kNCh][kLCh], ftime[kNCh][kLCh];
  long firstev_pos;
  //
  void EventMinMax();
  void Init();
  
public:
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

  //unsigned short voltage[kLCh];
  //double waveform[kNCh][kLCh], time[kNCh][kLCh];
  //float bin_width[kNCh][kLCh];
   
  static double threshold;
  static double invert[kNCh];
  static int baseline_npoints;
  static int bl_first;
  static int timing_calibration;
  static int regularize;
  static int nEvents;
  static int gverb;
#ifdef FILTERING
  static int mf_shape, mf_size;
  double *mfilter_coeff;
  double peak[kNCh],peakpos[kNCh];
#endif

  //drs4root();
  drs4root(const Char_t *in="", const Char_t *out="");
  ~drs4root();

  Int_t Skip_events(Int_t ev);
  Int_t Next_event();
  void Print_stat();
  //void First_event(); //obsolete, just reuse init.C
  void Print_header();
  void Print_event();
  ClassDef(drs4root,0)
};
#endif