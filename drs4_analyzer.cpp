/*
   Name:           drs4_analyzer, based on read_binary.
   Created by:     Stefan Ritt <stefan.ritt@psi.ch>
   Date:           July 30th, 2014

   Purpose:        Example file to read binary data saved by DRSOsc.
 
   Compile and run it with:
 
      gcc drs4_analyzer.cpp mfilter.cpp -o drs4_analyzer -lm 
      ./drs4_analyzer <filename> <options>

   This program assumes that a pulse from a signal generator is split
   and fed into channels #1 and #2. It then calculates the time difference
   between these two pulses to show the performance of the DRS board
   for time measurements.

   $Id: read_binary.cpp 21495 2014-09-26 14:20:49Z ritt $
   
   2015-10-06 Version 4 by Andrei Sukhanov <sukhanov@bnl.gov>
   FAST_CALIBRATION: Makes it possible to apply the calibration ~100 times faster by storing not the cell widths but cell times.
*/
#include <stdio.h>
void usage()
{
  printf("Usage: drs4_analyzer file.drs options\n\n");
  printf("Process binary files saved by DRSOsc\n\n");
  printf("OPTIONS:\n");
  printf("  -nN  negative pulse processing of channel N\n");
  printf("  -c   disable timing calibration\n");
  printf("  -tV  set LED threshold to V {0.0:1.0]\n");
  printf("  -fSN enable FIR filtering, S is filter shape (a:first event, r:rectangular, t:triangle), N is filter length\n");
  printf("  -bN  N-point baseline subtraction (N>1), single cell noise is calculated when N>0\n");
  printf("  -r   regularization of cell intervals (make it equidistant), not yet commissioned\n");
  printf("  -eN  process V events\n");
  printf("  -vN  verbosity level\n");
}

#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>


#define FAST_CALIBRATION

#define FILTERING
#ifdef FILTERING
  #include "mfilter.h"
  enum{
    kMaxFilterLength = 100
  };
#endif

#define kCellWidth 200./1024.

typedef struct {
   char           time_header[4];
   char           bn[2];
   unsigned short board_serial_number;
} THEADER;

typedef struct {
    char           event_header[4];
   unsigned int   event_serial_number;
   unsigned short year;
   unsigned short month;
   unsigned short day;
   unsigned short hour;
   unsigned short minute;
   unsigned short second;
   unsigned short millisecond;
   unsigned short range;
   char           bs[2];
   unsigned short board_serial_number;
   char           tc[2];
   unsigned short trigger_cell;
} EHEADER;

int gverb = 0;

double sigma(int nn, double sum, double sum2)
{ return sqrt(1./double(nn-1)*(sum2-sum*sum/double(nn)));}

void printstat(const char *s, int nn, double sum, double sum2, double factor, const char *unit)
{  printf("%s%.1f+-%.1f %s", s, factor*sum/double(nn), factor*sigma(nn,sum,sum2), unit);}

/*-----------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
   THEADER th;
   EHEADER eh;
   char hdr[4];
   unsigned short voltage[1024];
   double waveform[4][1024], time[4][1024];
   float bin_width[4][1024];
   int i, j, ch, nEv, chn_index,ii;
   double t1, t2, tt[2], dt, dv;
   char filename[256];
   
   int ndt;
   double default_threshold=0.2;
   double threshold=default_threshold;
   double sumdt, sumdt2;
   
   int opt;
   double invert[4];
   for(ii=0;ii<4;ii++) invert[ii]=0.;
   int nEvents=1000000;
   double wmax[4];
   int cmax[4];
   double asum[4], asum2[4];
   double bl_sum[4], bl_sum2[4];
   int baseline_npoints=10;
   int bl_first=5;
   int timing_calibration = 1;
   int regularize=0;
   double sum_widths[4];
   double v;
   
#ifdef FILTERING
   int mf_shape=0, mf_size=0;
   double *mfilter_coeff=0;
   double peak[4],peakpos[4];
   for(ii=0;ii<4;ii++) {peak[ii]=0.; peakpos[ii]=0.;}
#endif
   
   while ((opt = getopt (argc, argv, "n:cv:t:e:f:rb:")) != -1)
     switch (opt)
     {
       case 'n':
         invert[atoi(optarg)] = 1.;
         printf("Negative pulse processing of ch %i\n",atoi(optarg));
         break;
       case 'v':
         gverb = atoi(optarg);
         printf("Verbosity=%i\n",gverb);
         break;
       case 't':
         threshold = atof(optarg);
         printf("Threshold=%f\n",threshold);
         break;
       case 'e':
         nEvents = atoi(optarg);
         printf("Number of events to process: %i\n",nEvents);
         break;
       case 'c':
         timing_calibration=0;
         break;
       case 'b':
         baseline_npoints = atoi(optarg); 
         break;
       case 'f':
#ifdef FILTERING
        switch (optarg[0]){
           case 'a': mf_shape = kFShape_arbitrary; break;
           case 'r': mf_shape = kFShape_rectangle; break;
           case 't': mf_shape = kFShape_triangle; break;
         }
         mf_size = kMaxFilterLength;
         if(strlen(optarg)>1) mf_size = atoi(optarg+1);
         if(mf_size > kMaxFilterLength) mf_size = kMaxFilterLength;
         printf("Filtering of type %i[%i] enabled\n",mf_shape,mf_size);
#endif
         break;
       case 'r':
         regularize=1;
         printf("Cell Regularization is on!\n");
         break;
       case '?':
         fprintf (stderr, "ERROR with option -%c\n", optopt);
         return 1;
       default:
         usage();
         return 1;
     }
   if(argc==optind) {usage();return 1;}

   if(threshold == default_threshold) printf("Using default threshold=%f\n",threshold);
#ifdef FAST_CALIBRATION
   printf("Using fast calibration\n");
#endif
    if(baseline_npoints>1) 
      printf("Baseline subtraction for %i points starting at %i\n",baseline_npoints,bl_first);
     
   // open the binary waveform file
   strcpy(filename, argv[optind]);
   FILE *f = fopen(filename, "r");
   if (f == NULL) {
      printf("Cannot find file \'%s\'\n", filename);
      return 0;
   }
   // read time header
   fread(&th, sizeof(th), 1, f);
   printf("Found data for board #%d\n", th.board_serial_number);

   // read time bin widths
   memset(bin_width, sizeof(bin_width), 0);
   for (ch=0 ; ch<5 ; ch++) {
      fread(hdr, sizeof(hdr), 1, f);
      if (hdr[0] != 'C') {
         // event header found
         fseek(f, -4, SEEK_CUR);
         break;
      }
      i = hdr[3] - '0' - 1;
      //&RA
      i=ch;
      fread(&bin_width[i][0], sizeof(float), 1024, f);
      for(ii=0, sum_widths[ch]=0.; ii<1024; ii++) sum_widths[ch] += bin_width[ch][ii];
      printf("Found timing calibration for channel #%d, sum_widths=%f\n", i+1, sum_widths[ch]);
      if(timing_calibration==0) for (ii=0; ii<1024; ii++) bin_width[ch][ii] = kCellWidth;
#ifdef FAST_CALIBRATION
      for (ii=1; ii<1024-1; ii++) bin_width[ch][ii] += bin_width[ch][ii-1];
      for (ii=1024-1; ii>0; ii--) bin_width[ch][ii] = bin_width[ch][ii-1];
      bin_width[ch][0] = 0.;
#endif
      if(gverb&8) {for(ii=0;ii<1024;ii++) printf("%i:%03i ",ii,int(bin_width[i][ii]*1000.)); printf("\n");}
   }
   // initialize statistics
   ndt = 0;
   sumdt = sumdt2 = 0;
   for(ch=0;ch<4;ch++) {wmax[ch]=0.; asum[ch]=0.; asum2[ch]=0.; bl_sum[ch]=0.; bl_sum2[ch]=0.;}
   
   // loop over all events in the data file
   for (nEv= 0 ; nEv<nEvents; nEv++) {
      // read event header
      i = (int)fread(&eh, sizeof(eh), 1, f);
      if (i < 1) break;
      // reach channel data
      for (ch=0 ; ch<5 ; ch++) {
         if (gverb&4) printf("Ch#%i, Trigger cell %i timecal=%f\n",
           ch,eh.trigger_cell,bin_width[ch][eh.trigger_cell]);
         i = (int)fread(hdr, sizeof(hdr), 1, f);
         if (i < 1) break;
         if (hdr[0] != 'C') {
            // event header found
            fseek(f, -4, SEEK_CUR);
            break;
         }
         //&RA/chn_index = hdr[3] - '0' - 1;
         chn_index = ch;
         fread(voltage, sizeof(short), 1024, f);

         for (i=0 ; i<1024 ; i++) {
            // convert data to volts
            waveform[chn_index][i] = (voltage[i] / 65536. * (1.-2.*invert[ch]) + invert[ch] + eh.range/1000.0 - 0.5);
            if(waveform[chn_index][i]>wmax[ch]) {wmax[ch]=waveform[chn_index][i]; cmax[ch]=i;}
            if (gverb&2) printf("%i:%03i ",i,int(waveform[chn_index][i]*1000.));            
            // calculate time for this cell
#ifdef FAST_CALIBRATION
            time[chn_index][i] = bin_width[chn_index][(i+eh.trigger_cell) % 1024] 
              - bin_width[chn_index][eh.trigger_cell];
            if (time[chn_index][i] <0.) time[chn_index][i] += sum_widths[ch];
# else
            for (j=0,time[chn_index][i]=0 ; j<i ; j++)
               time[chn_index][i] += bin_width[chn_index][(j+eh.trigger_cell) % 1024];
#endif
            if (gverb&8) printf("%03i ",int(time[chn_index][i]*1000.)); 
         }
         if (gverb&(2|8)) printf("\n");
      }
#ifdef FILTERING
      // use the first channel of first event to build the arbitrary matching filter
      if (mf_shape && nEv==0) 
      {
        mf_size=mfilter_create(waveform[0],1024,mf_shape,mf_size,&mfilter_coeff);
        //if(gverb&0x10) 
          printf("MFilter of type %i[%i] created:\n",mf_shape,mf_size);
        //if(gverb&0x20) 
        {
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
      t1 = time[0][(1024-eh.trigger_cell) % 1024];
      for (ch=1 ; ch<4 ; ch++) {
         t2 = time[ch][(1024-eh.trigger_cell) % 1024];
         dt = t1 - t2;
         if (gverb&8) printf("ch%i dt=%f\n",ch,dt);
         for (i=0 ; i<1024 ; i++) time[ch][i] += dt;
      }
      // baseline subtraction
      if(baseline_npoints > 1)
        for(ch=0;ch<2;ch++)
        {
          for(v=0., ii=bl_first; ii<bl_first+baseline_npoints; ii++)
            v += waveform[ch][ii];
          v /= double(baseline_npoints);        
          for(ii=0;ii<1024;ii++) waveform[ch][ii] -= v;
        }
      // regularization of the cell intervals to make them equidistant
      double time_shift, expected_time, prev;
      if(regularize)
      {
        for(ch=0;ch<2;ch++)
          for(prev=waveform[ch][0],ii=1;ii<1024;ii++)
          {
            dt=time[ch][ii]-time[ch][ii-1];
            expected_time = double(ii)*kCellWidth;
            time_shift = time[ch][ii] - expected_time;
            time[ch][ii] = expected_time;
            dv=waveform[ch][ii]-prev;
            prev = waveform[ch][ii];
            waveform[ch][ii] +=  dv/dt*time_shift;
          }
      }            
#ifdef FILTERING
      // FIR filtering
      if (mf_shape)
      for(ch=0;ch<2;ch++)
      {
        //in-place filter
        mfilter_filter(waveform[ch],1024,waveform[ch],&(peak[ch]),&(peakpos[ch]));
        //recalculate max
        for(ii=0;ii<1024;ii++)
          if(waveform[ch][ii]>wmax[ch]) {wmax[ch]=waveform[ch][ii]; cmax[ch]=ii;}
        if(gverb&0x40){for(ii=0;ii<1024;ii++) printf("%i,%.2f:%.4f ",ii,time[ch][ii],waveform[ch][ii]);  printf("\nmax=%f @ %i\n",wmax[ch],cmax[ch]);}
      }
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
        threshold,eh.trigger_cell);
      tt[0] = tt[1] = 0.;
      for (ch=0;ch<2;ch++)
      {
        for (i=0 ; i<1022 ; i++)
            if (waveform[ch][i] < threshold && waveform[ch][i+1] >= threshold) 
            {
              tt[ch] = (threshold-waveform[ch][i])/(waveform[ch][i+1]-waveform[ch][i])*(time[ch][i+1]-time[ch][i])+time[ch][i];
              if (gverb&4) printf("Cross[%i]: %f,%f @ %i,%f,%f t=%f\n", 
                ch, waveform[ch][i], waveform[ch][i+1], i, time[ch][i], time[ch][i+1], tt[ch]);
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
      // print statistics
      if(gverb&1 || (nEv%1000)==999)
      {
        printf("ev#%5i,",nEv+1);
        printstat(" dT=",nEv+1,sumdt,sumdt2,1000.,"ps"); 
        printstat(", A0=",nEv+1,asum[0],asum2[0],1000.,"mV");
        printstat(", A1=",nEv+1,asum[1],asum2[1],1000.,"mV");
        if(baseline_npoints > 0)
        {
          printstat(", BL0=",nEv+1,bl_sum[0],bl_sum2[0],1000.,"mV");
          printstat(", BL1=",nEv+1,bl_sum[1],bl_sum2[1],1000.,"mV");
        }
        printf("\n");
      }
  }   
        printf("ev#%5i,",nEv);
        printstat(" dT=",nEv,sumdt,sumdt2,1000.,"ps"); 
        printstat(", A0=",nEv,asum[0],asum2[0],1000.,"mV");
        printstat(", A1=",nEv,asum[1],asum2[1],1000.,"mV");
        if(baseline_npoints > 0)
        {
          printstat(", BL0=",nEv,bl_sum[0],bl_sum2[0],1000.,"mV");
          printstat(", BL1=",nEv,bl_sum[1],bl_sum2[1],1000.,"mV");
        }
        printf("\n");
}

