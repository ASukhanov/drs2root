#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mfilter.h"

#define kFilterThreshold 0.1
//enum {double kFilterThreshold=0.1};

extern int gverb;

//global storage
int gmf_type=0;
int gmf_size=0;
double *gmf_out=0;
double *gmf_coeff=0;

int mfilter_create(const double *amplitude, const int size, const int type, const int fsize, double **coeff)
{
  int ii;
  double max=0;
  int imax,ifleft=0,ifright=0;
  double v,l2=0.;
  
  //init
  gmf_coeff = (double*) malloc(sizeof(double)*fsize);
  if(gmf_coeff==0) {printf("Could not allocate gmf_coeff\n"); exit(1);}
  
  gmf_type = type;
  //printf("filter of type %i[%i] allocated at %x\n", gmf_type, fsize, int(gmf_coeff));
  
  for(ii=0;ii<fsize;ii++) gmf_coeff[ii]=0.;
  
  switch (type)
  {
    case kFShape_arbitrary:
      for(ii=0;ii<size;ii++) if(amplitude[ii]>max) {max=amplitude[ii];imax=ii;}
      for(ii=imax-fsize/2;ii<imax;ii++) if(amplitude[ii]>kFilterThreshold) break;
      ifleft=ii;
      for(ii=imax+fsize/2;ii>imax;ii--) if(amplitude[ii]>kFilterThreshold) break;
      ifright=ii;
      gmf_size = ifright - ifleft;
      for(ii=0;ii<gmf_size;ii++) {gmf_coeff[ii]=amplitude[ii+ifleft];} 
      break;
    case kFShape_rectangle:
      gmf_size = fsize;
      for(ii=0;ii<gmf_size;ii++) gmf_coeff[ii]=1.;
      break;
    case kFShape_triangle:
      gmf_size = fsize;
      for(ii=0;ii<gmf_size/2;ii++) gmf_coeff[ii]=double(ii+1);
      for(ii=0;ii<=gmf_size/2;ii++) gmf_coeff[gmf_size/2+ii]=double(gmf_size/2-ii+1.);
      break;
  }
  
  // L2-norm
  for(l2=0., ii=0;ii<gmf_size;ii++) {l2 += gmf_coeff[ii]*gmf_coeff[ii];}
  for(ii=0;ii<gmf_size;ii++) {gmf_coeff[ii] /= sqrt(l2);}
  
  *coeff = gmf_coeff;
  //printf("returning %x into %x\n",int(gmf_coeff), int(coeff));

  //allocate working buffer;
  gmf_out = (double*) malloc(sizeof(double)*(size+gmf_size));
  if(gmf_out==0) {printf("Could not allocate gmf_out\n"); exit(1);}
  return gmf_size;
}

int mfilter_filter(const double *input, const int size, double *out, double *peak_amplitude, double *peak_position)
{
  int ii,jj;
  for(ii=0;ii<size;ii++)
  {
    for(gmf_out[ii]=0.,jj=0;jj<gmf_size;jj++)
      gmf_out[ii] += input[ii+jj]*gmf_coeff[jj];
  }
  if(gverb&0x80)
  {
    printf("mfilter_filter[%i]:\n",gmf_size);
    for(ii=0;ii<gmf_size;ii++) printf("%i:%f\n",ii,gmf_coeff[ii]);
    printf("input:\n"); for(ii=0;ii<size;ii++) printf("%i:%f ",ii,input[ii]); printf("\n");
    printf("output:\n"); for(ii=0;ii<size;ii++) printf("%i:%f ",ii,gmf_out[ii]); printf("\n");
  }
  for(ii=0;ii<size;ii++) out[ii] = gmf_out[ii];
}

int mfilter_delete()
{
  if(gmf_out) free(gmf_out);
  if(gmf_coeff) free(gmf_coeff);
  return 0;
}