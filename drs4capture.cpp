/********************************************************************\

  Name:         drs4capture.cpp
  based on      drs_exam, created by:   Stefan Ritt
  Modified by   Andrei Sukhanov

  Contents:     Application to read out a DRS4
                evaluation board

  $Id: drs_exam.cpp 21308 2014-04-11 14:50:16Z ritt $
  Version v2 2016-05-12 argc,argv, non-calibrated data.
  The simplest way to compile: copy it to drsxxx/src/drs_exam; cd drsxxx; make

\********************************************************************/

#include <math.h>

#ifdef _MSC_VER

#include <windows.h>

#elif defined(OS_LINUX)

#define O_BINARY 0

#include <unistd.h>
#include <ctype.h>
#include <sys/ioctl.h>
#include <errno.h>

#define DIR_SEPARATOR '/'

#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>130.199.23.167
#include <getopt.h>
#include <time.h>

#include "strlcpy.h"
#include "DRS.h"

/*''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''*/
static void print_usage(const char *prog)
{
  printf("Usage: %s [-wcveh]\n", prog);
  puts("  -w --write    enable writing to disk\n"
  "  -c --nocalib    do not calibrate data\n"
  "  -vN --verbose    verbose\n"
  "  -eN --events    get N events\n"  
  "  -h --help     help\n"
  );
  exit(1);
}
int o_write=0;
int o_calibrate=1;
int o_verb=0;
int o_nev=10000;
//const char* o_filename="$DATA/data.txt";
const char* o_filename="/home/pi/data1/data/data.drs";
char gtxt[128];

static void parse_opts(int argc, char *argv[])
{
    static const struct option lopts[] = {
      { "write",  1, 0, 'w' },
      { "nocalib",   1, 0, 'c' },
      { "verb",   1, 0, 'v' },
      { "help",    0, 0, 'h' },
      { NULL, 0, 0, 0 }
    };
  int c;
  
  while (1) 
  {
    c = getopt_long(argc, argv, "wcv:n:h", lopts, NULL);
    
    if (c == -1)
      break;
    
    switch (c) {
      case 'w': o_write = 1; break;
      case 'c': o_calibrate = 0; break;
      case 'v': o_verb = atoi(optarg);break;
      case 'n': o_nev = atoi(optarg);break;
      case 'h':
      default:
        print_usage(argv[0]);
        break;
    }
  }
}
/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/

int main(int argc, char *argv[])
{
  int i, nev, nBoards;
  DRS *drs;
  DRSBoard *b;
  float time_array[8][1024];
  float wave_array[8][1024];
  unsigned short wave_ushort_array[8][1024];
  FILE  *ofile=NULL;
  
  /* do initial scan */
  drs = new DRS();
  
  /* show any found board(s) */
  for (i=0 ; i<drs->GetNumberOfBoards() ; i++) {
    b = drs->GetBoard(i);
    printf("Found DRS4 evaluation board, serial #%d, firmware revision %d\n", 
           b->GetBoardSerialNumber(), b->GetFirmwareVersion());
  }
  /* exit if no board found */
  nBoards = drs->GetNumberOfBoards();
  if (nBoards == 0) {
    printf("No DRS4 evaluation board found\n");
    return 0;
  }
  parse_opts(argc, argv);
  
  /* continue working with first board only */
  b = drs->GetBoard(0);
  
  /* initialize board */
  b->Init();
  
  /* set sampling frequency */
  b->SetFrequency(5, true);
  
  /* enable transparent mode needed for analog trigger */
  b->SetTranspMode(1);
  
  /* set input range to -0.5V ... +0.5V */
  b->SetInputRange(0);
  
  /* use following line to set range to 0..1V */
  //b->SetInputRange(0.5);
  
  /* use following line to turn on the internal 100 MHz clock connected to all channels  */
  //b->EnableTcal(1);
  
  /* use following lines to enable hardware trigger on CH1 at 50 mV positive edge */
  if (b->GetBoardType() >= 8) {        // Evaluaiton Board V4&5
    b->EnableTrigger(1, 0);           // enable hardware trigger
    b->SetTriggerSource(1<<0);        // set CH1 as source
  } else if (b->GetBoardType() == 7) { // Evaluation Board V3
    b->EnableTrigger(0, 1);           // lemo off, analog trigger on
    b->SetTriggerSource(0);           // use CH1 as source
  }
  b->SetTriggerLevel(0.05);            // 0.05 V
  b->SetTriggerPolarity(false);        // positive edge
  
  /* use following lines to set individual trigger elvels */
  //b->SetIndividualTriggerLevel(1, 0.1);
  //b->SetIndividualTriggerLevel(2, 0.2);
  //b->SetIndividualTriggerLevel(3, 0.3);
  //b->SetIndividualTriggerLevel(4, 0.4);
  //b->SetTriggerSource(15);
  
  b->SetTriggerDelayNs(0);             // zero ns trigger delay
  
  /* use following lines to enable the external trigger */
  //if (b->GetBoardType() == 8) {     // Evaluaiton Board V4
  //   b->EnableTrigger(1, 0);           // enable hardware trigger
  //   b->SetTriggerSource(1<<4);        // set external trigger as source
  //} else {                          // Evaluation Board V3
  //   b->EnableTrigger(1, 0);           // lemo on, analog trigger off
  // }
  
  if(o_write)
  {
    /* open file to save waveforms */
    ofile = fopen(o_filename, "w");
    if (ofile == NULL) {
      sprintf(gtxt,"ERROR: Cannot open file %s",o_filename);
      perror(gtxt);
      return 1;
    }
    printf("file %s opened\n",o_filename);
  }
  /* repeat o_nev times */
  time_t t_now,t_prev;
  t_prev = time(0);
  double t_diff;
  long fpos_now, fpos_prev=0;
  int evInterval=1000;
  for (nev=0 ; nev<o_nev ; nev++) {
    if((nev%evInterval)==evInterval-1) 
    { 
      t_now = time(0);
      t_diff = difftime(t_now,t_prev);
      fpos_now = ftell(ofile);
      printf("ev%06i %04.1f/s %04.2fMB/s\n",nev,(float)evInterval/t_diff,(float)(fpos_now-fpos_prev)/t_diff/1e6);
      memcpy(&t_prev,&t_now,sizeof(t_now));
      t_prev = time(0);
      fpos_prev = fpos_now;
    }
    /* start board (activate domino wave) */
    b->StartDomino();
    
    /* wait for trigger */
    if(o_verb&1) printf("Waiting for trigger...");

    fflush(stdout);
    while (b->IsBusy());
    
    /* read all waveforms */
    b->TransferWaves(0, 8);
    
    /* read time (X) array of first channel in ns */
    b->GetTime(0, 0, b->GetTriggerCell(0), time_array[0]);
    
    /* decode waveform (Y) array of first channel in mV */
    if(o_calibrate)
    {
      b->GetWave(0, 0, wave_array[0]);
    }else{
      b->GetRawWave(0, 0, wave_ushort_array[0], 0);
    }
    
    /* read time (X) array of second channel in ns
     *       Note: On the evaluation board input #1 is connected to channel 0 and 1 of
     *       the DRS chip, input #2 is connected to channel 2 and 3 and so on. So to
     *       get the input #2 we have to read DRS channel #2, not #1. */
    b->GetTime(0, 2, b->GetTriggerCell(0), time_array[1]);
    
    /* decode waveform (Y) array of second channel in mV */
    b->GetWave(0, 2, wave_array[1]);
    
    if(o_write)
    {
      /* Save waveform: X=time_array[i], Yn=wave_array[n][i] */
      if(o_calibrate)
      {
        fprintf(ofile, "Event #%d -calibrated-----------\n  t1[ns]  u1[mV]  t2[ns] u2[mV]\n", nev);
        for (i=0 ; i<1024 ; i++)
          fprintf(ofile, "%7.3f %7.1f %7.3f %7.1f\n", time_array[0][i], wave_array[0][i], time_array[1][i], wave_array[1][i]);
      }
      else
      {
        fprintf(ofile, "Event #%d -raw------------------\n  t1[ns]  u1      t2[ns] u2[mV]\n", nev);
        for (i=0 ; i<1024 ; i++)
          fprintf(ofile, "%7.3f %5i %7.3f %7.1f\n", time_array[0][i], wave_ushort_array[0][i], time_array[1][i], wave_array[1][i]);
      }
    }
    if(o_verb&1) printf("\rEvent #%d read successfully\n", nev);
  }
  if(o_write)  fclose(ofile);
  /* delete DRS object -> close USB connection */
  delete drs;
}
