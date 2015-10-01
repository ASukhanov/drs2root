/*
   Name:           read_binary.cpp
   Created by:     Stefan Ritt <stefan.ritt@psi.ch>
   Date:           July 30th, 2014

   Purpose:        Example file to read binary data saved by DRSOsc.
 
   Compile and run it with:
 
      gcc -o read_binary read_binary.cpp -lm
 
      ./read_binary <filename>

   This program assumes that a pulse from a signal generator is split
   and fed into channels #1 and #2. It then calculates the time difference
   between these two pulses to show the performance of the DRS board
   for time measurements.

   $Id: read_binary.cpp 21495 2014-09-26 14:20:49Z ritt $
   
   2015-10-01 Version 2 by Andrei Sukhanov
*/
#include <stdio.h>
void usage()
{
  printf("proc_drs file.drs options\n\n");
  printf("Process binary files saved by DRSOsc\n\n");
  printf("OPTIONS:\n");
  printf("  -n   negative pulse processing\n");
  printf("  -c   disable timing calibration\n");
  printf("  -tV  set LED threshold to V {0.0:1.0]\n");
  printf("  -eV  process V events\n");
  printf("  -vV  verbosity level\n");
}

#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

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
   int i, j, ch, n, chn_index,ii;
   double t1, t2, dt;
   char filename[256];
   
   int ndt;
   double threshold, sumdt, sumdt2;
   
   int opt;
   double negative_pulse=0.;
   int verb=0;
   int nEvents=1000000;
   double wmax[4];
   int cmax[4];
   double ampl[4], asigma[4];
   int timing_calibration = 1;
   
   while ((opt = getopt (argc, argv, "ncv:t:e:")) != -1)
     switch (opt)
     {
       case 'n':
         negative_pulse = 1.;
         printf("Negative pulse processing\n");
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
      printf("Found timing calibration for channel #%d\n", i+1);
      //&RA
      i=ch;
      fread(&bin_width[i][0], sizeof(float), 1024, f);
      if(timing_calibration==0) for (ii=0; ii<1024; ii++) bin_width[ch][ii] = 0.2;
      if(verb&8) {for(ii=0;ii<1024;ii++) printf("%03i ",int(bin_width[i][ii]*1000.)); printf("\n");}
   }
   
   // initialize statistics
   ndt = 0;
   sumdt = sumdt2 = 0;
   for(ch=0;ch<4;ch++) {wmax[ch]=0.; ampl[ch]=0.; asigma[ch]=0.;}
   
   // loop over all events in the data file
   for (n= 0 ; n<nEvents; n++) {
      // read event header
      i = (int)fread(&eh, sizeof(eh), 1, f);
      if (i < 1)
         break;
      
      // reach channel data
      for (ch=0 ; ch<5 ; ch++) {
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
            waveform[chn_index][i] = (voltage[i] / 65536. * (1.-2.*negative_pulse) + negative_pulse + eh.range/1000.0 - 0.5);
            if(waveform[chn_index][i]>wmax[ch]) {wmax[ch]=waveform[chn_index][i]; cmax[ch]=i;}
            if (verb&2) printf("%03i ",int(waveform[chn_index][i]*1000.));            
            // calculate time for this cell
            for (j=0,time[chn_index][i]=0 ; j<i ; j++)
               time[chn_index][i] += bin_width[chn_index][(j+eh.trigger_cell) % 1024];
         }
         if (verb&2) printf("\n");
      }
      
      // align cell #0 of all channels
      t1 = time[0][(1024-eh.trigger_cell) % 1024];
      for (ch=1 ; ch<4 ; ch++) {
         t2 = time[ch][(1024-eh.trigger_cell) % 1024];
         dt = t1 - t2;
         for (i=0 ; i<1024 ; i++)
            time[ch][i] += dt;
      }
      
      t1 = t2 = 0;
      //threshold = 0.3;
      
      // find peak in channel 1 above threshold
      if (verb&4) printf("Looking for threshold=%f crossing\n",threshold);
      for (i=0 ; i<1022 ; i++)
          if (waveform[0][i] < threshold && waveform[0][i+1] >= threshold) {
            t1 = (threshold-waveform[0][i])/(waveform[0][i+1]-waveform[0][i])*(time[0][i+1]-time[0][i])+time[0][i];
            if (verb&4) printf("Cross0: %f @ %i, t1=%f\n", waveform[0][i], i, t1);
            break;
         }
      
      // find peak in channel 2 above threshold
      for (i=0 ; i<1022 ; i++)
         if (waveform[1][i] < threshold && waveform[1][i+1] >= threshold) {
            t2 = (threshold-waveform[1][i])/(waveform[1][i+1]-waveform[1][i])*(time[1][i+1]-time[1][i])+time[1][i];
            if (verb&4) printf("Cross1: %f @ %i, t2=%f, t2-t1=%f\n", waveform[1][i], i, t2, t2-t1);
            break;
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
      
      if(verb&1 || (n%1000)==999) printf("Event #%d, t2-t1=%f, A0=%f, A1=%f\n", eh.event_serial_number, t2-t1,wmax[0],wmax[1]);
            
      // calculate distance of peaks with statistics
      if (t1 > 0 && t2 > 0) {
         ndt++;
         dt = t2 - t1;
         sumdt += dt;
         sumdt2 += dt*dt;
      }
   }
   
   // print statistics
   printf("en#%i,dT = %1.3lfns +- %1.1lfps. ",n,
     sumdt/double(ndt), 1000*sqrt(1.0/(double(ndt)-1.)*(sumdt2-1.0/double(ndt)*sumdt*sumdt)));
   for (ch=0 ; ch<2 ; ch++)
     printf("A[%i]=%f +- %f. ",ch, ampl[ch]/double(n),
       sqrt(1.0/(double(n)-1.)*(asigma[ch]-1.0/double(n)*ampl[ch]*ampl[ch])));
   printf("\n");
   return 1;
}

