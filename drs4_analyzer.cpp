/*
   Name:           drs4_analyzer, based on read_binary.
   Created by:     Stefan Ritt <stefan.ritt@psi.ch>
   Date:           July 30th, 2014

   Purpose:        Example file to read binary data saved by DRSOsc.
 
   Compile and run it with:
 
      gcc drs4_analyzer.cpp mfilter.cpp -o drs4_analyzer -lm
 
      ./read_binary <filename>

   This program assumes that a pulse from a signal generator is split
   and fed into channels #1 and #2. It then calculates the time difference
   between these two pulses to show the performance of the DRS board
   for time measurements.

   $Id: read_binary.cpp 21495 2014-09-26 14:20:49Z ritt $
   
   2015-10-01 Version 2 by Andrei Sukhanov <sukhanov@bnl.gov>
   FAST_CALIBRATION: Makes it possible to apply the calibration ~100 times faster by storing not the cell widths but cell times.
*/
#include <stdio.h>
void usage()
{
  printf("proc_drs file.drs options\n\n");
  printf("Process binary files saved by DRSOsc\n\n");
  printf("OPTIONS:\n");
  printf("  -nN  negative pulse processing of channel N\n");
  printf("  -c   disable timing calibration\n");
  printf("  -tV  set LED threshold to V {0.0:1.0]\n");
  printf("  -eN  process V events\n");
  printf("  -vN  verbosity level\n");
}

#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define FAST_CALIBRATION

//#define MATCH_FILTERING
#ifdef MATCH_FILTERING
  int mfilter_create(const double *amplitude, int length);
  int mfilter_filter(double *amplitude, int length, double *peak_amplitude, double *peak_position);
  int mfilter_delete();
#endif

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
   double t1, t2, tt[2], dt;
   char filename[256];
   
   int ndt;
   double default_threshold=0.2;
   double threshold=default_threshold;
   double sumdt, sumdt2;
   
   int opt;
   double invert[4];
   for(ii=0;ii<4;ii++) invert[ii]=0.;
   int verb=0;
   int nEvents=1000000;
   double wmax[4];
   int cmax[4];
   double ampl[4], asigma[4];
   int timing_calibration = 1;
   double sum_widths[4];
   
   while ((opt = getopt (argc, argv, "n:cv:t:e:")) != -1)
     switch (opt)
     {
       case 'n':
         invert[atoi(optarg)] = 1.;
         printf("Negative pulse processing of ch %i\n",atoi(optarg));
         break;
       case 'v':
         verb = atoi(optarg);
         printf("Verbosity=%i\n",verb);
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
       case '?':
         fprintf (stderr, "ERROR with option -%c\n", optopt);
         return 1;
       default:
         usage();
         return 1;
     }

   if(threshold == default_threshold) printf("Using default threshold=%f\n",threshold);
#ifdef FAST_CALIBRATION
   printf("Using fast calibration\n");
#endif
     
   if(argc==optind) {usage();return 1;}

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
      if(timing_calibration==0) for (ii=0; ii<1024; ii++) bin_width[ch][ii] = 0.2;
#ifdef FAST_CALIBRATION
      for (ii=1; ii<1024-1; ii++) bin_width[ch][ii] += bin_width[ch][ii-1];
      for (ii=1024-1; ii>0; ii--) bin_width[ch][ii] = bin_width[ch][ii-1];
      bin_width[ch][0] = 0.;
#endif
      if(verb&8) {for(ii=0;ii<1024;ii++) printf("%i:%03i ",ii,int(bin_width[i][ii]*1000.)); printf("\n");}
   }
   
   // initialize statistics
   ndt = 0;
   sumdt = sumdt2 = 0;
   for(ch=0;ch<4;ch++) {wmax[ch]=0.; ampl[ch]=0.; asigma[ch]=0.;}
   
   // loop over all events in the data file
   for (nEv= 0 ; nEv<nEvents; nEv++) {
      // read event header
      i = (int)fread(&eh, sizeof(eh), 1, f);
      if (i < 1)
         break;
      // reach channel data
      for (ch=0 ; ch<5 ; ch++) {
         if (verb&2) printf("Ch#%i, Trigger cell %i timecal=%f\n",
           ch,eh.trigger_cell,bin_width[ch][eh.trigger_cell]);
         i = (int)fread(hdr, sizeof(hdr), 1, f);
         if (i < 1)
            break;
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
            if (verb&2) printf("%i:%03i ",i,int(waveform[chn_index][i]*1000.));            
            // calculate time for this cell
#ifdef FAST_CALIBRATION
            time[chn_index][i] = bin_width[chn_index][(i+eh.trigger_cell) % 1024] 
              - bin_width[chn_index][eh.trigger_cell];
            if (time[chn_index][i] <0.) time[chn_index][i] += sum_widths[ch];
# else
            for (j=0,time[chn_index][i]=0 ; j<i ; j++)
               time[chn_index][i] += bin_width[chn_index][(j+eh.trigger_cell) % 1024];
#endif
            if (verb&8) printf("%03i ",int(time[chn_index][i]*1000.)); 
         }
         if (verb&(2|8)) printf("\n");
      }
#ifdef MATCH_FILTERING
      // use first channel of first event to build the matching filter
      if (nEv==0) mfilter_create(waveform[0],1024);
#endif
      // align cell #0 of all channels
      t1 = time[0][(1024-eh.trigger_cell) % 1024];
      for (ch=1 ; ch<4 ; ch++) {
         t2 = time[ch][(1024-eh.trigger_cell) % 1024];
         dt = t1 - t2;
         if (verb&8) printf("ch%i dt=%f\n",ch,dt);
         for (i=0 ; i<1024 ; i++) time[ch][i] += dt;
      }
      
      tt[0] = tt[1] = 0.;
      //threshold = 0.3;
      
      // find peak in channels above threshold
      if (verb&4) printf("Looking for threshold=%f crossing, trigcell %i\n",
	threshold,eh.trigger_cell);
      for (ch=0;ch<2;ch++)
      {
        for (i=0 ; i<1022 ; i++)
            if (waveform[ch][i] < threshold && waveform[ch][i+1] >= threshold) 
            {
              tt[ch] = (threshold-waveform[ch][i])/(waveform[ch][i+1]-waveform[ch][i])*(time[ch][i+1]-time[ch][i])+time[ch][i];
              if (verb&4) printf("Cross[%i]: %f,%f @ %i,%f,%f t=%f\n", 
                ch, waveform[ch][i], waveform[ch][i+1], i, time[ch][i], time[ch][i+1], tt[ch]);
              break;
            }
      }

      // estimate amplitudes of first two channels by averaging around peak
      int half_top_width=2;
      double na = double(half_top_width*2+1);
      for (ch=0 ; ch<2 ; ch++)
      {
        wmax[ch]=0.;
        for (i=cmax[ch]-half_top_width; i<=cmax[ch]+half_top_width; i++) wmax[ch] += waveform[ch][i];
        wmax[ch] /=na;
        ampl[ch] += wmax[ch];
        asigma[ch] += wmax[ch]*wmax[ch];
      }
      
      if(verb&1 || (nEv%1000)==999) printf("Event #%d, t2-t1=%f, A0=%f, A1=%f\n", eh.event_serial_number, tt[1]-tt[0],wmax[0],wmax[1]);
            
      // calculate distance of peaks with statistics
      if (tt[0] > 0 && tt[1] > 0) {
         ndt++;
         dt = tt[1] - tt[0];
         sumdt += dt;
         sumdt2 += dt*dt;
      }
   }
   
   // print statistics
   printf("en#%i,dT = %1.3lfns +- %1.1lfps. ",nEv,
     sumdt/double(ndt), 1000*sqrt(1.0/(double(ndt)-1.)*(sumdt2-1.0/double(ndt)*sumdt*sumdt)));
   for (ch=0 ; ch<2 ; ch++)
     printf("A[%i]=%f +- %f. ",ch, ampl[ch]/double(nEv),
       sqrt(1.0/(double(nEv)-1.)*(asigma[ch]-1.0/double(nEv)*ampl[ch]*ampl[ch])));
   printf("\n");
   return 1;
}

