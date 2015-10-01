#ifndef drs2root_h
#define drs2root_h

#include <TFile.h>



enum
{
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
  UShort_t reserved;
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
    struct Channel_Time_t channel_time[kNCh];
  };
  
class drs2root : public TFile
{
private:  
  void EventMinMax();
  
public:
  struct Event_t fEv;

  struct Time_header_t fTime_header;

    FILE    *fD;
    TFile   *ffile;
    Int_t   fsize;
    Char_t  fname[256];
    ULong_t fpos;
    Int_t   fevcount;
    UInt_t  fevnum;
    UInt_t  fMin[kNCh], fMinPos[kNCh];
    UInt_t  fMax[kNCh], fMaxPos[kNCh];

  drs2root(const Char_t *in="", const Char_t *out="");
  ~drs2root();

  Int_t Find_event(Int_t ev);
  Int_t Next_event();
  void Print_header();
  void Print_event();
  ClassDef(drs2root,0)
};
#endif