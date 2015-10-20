#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mfilter.h"

// FILTER_NONUNIFORM is controlled in mfilter.h
#ifdef FILTER_NONUNIFORM
  #include <math.h>
  #include <gsl/gsl_errno.h>
  #include <gsl/gsl_spline.h>
#endif


#define kFilterThreshold 0.1
//enum {double kFilterThreshold=0.1};

//extern int gverb;

//global storage
int gmf_type=0;
int gmf_size=0;
int gmf_out_size=0;
double *gmf_out=0;
double *gmf_coeff=0;
double *gmf_x=0;

int mfilter_create(const double *xx, const double *yy, const int size, const double xstep, const int type, const int fsize, double **coeff, double **mfx)
{
  int ii;
  double max=0;
  int imax=0,ifleft=0,ifright=0;
  double l2=0.;
  
  //printf("Creating filter\n");
  //init
  gmf_coeff = (double*) malloc(sizeof(double)*fsize);
  gmf_x = (double*) malloc(sizeof(double)*fsize);
  if(gmf_coeff==0) {printf("Could not allocate gmf_coeff\n"); exit(1);}
  
  gmf_type = type;
  printf("filter of type %i[%i] allocated at %lx\n", gmf_type, fsize, long(gmf_coeff));
  //for(ii=0;ii<size;ii++) printf("xy=%.3f,%.3f\n",xx[ii],yy[ii]);
  
  for(ii=0;ii<fsize;ii++) gmf_coeff[ii]=0.;
  
  switch (type)
  {
    case kFShape_1st_event:
    case kFShape_external:
      for(ii=0;ii<size;ii++) if(yy[ii]>max) {max=yy[ii];imax=ii;}
      ii = imax-fsize/2;
      if(ii<1) ii = 0;
      for(;ii<imax;ii++) if(yy[ii]>kFilterThreshold) break;
      ifleft=ii;
      for(ii=imax; ii<imax+fsize/2; ii++) if(yy[ii]<kFilterThreshold) break;
      ifright=ii;
      gmf_size = ifright - ifleft;
      for(ii=0;ii<gmf_size;ii++) 
      {
        gmf_coeff[ii]=yy[ii+ifleft];
        gmf_x[ii]=xx[ii+ifleft] - xx[ifleft];
      } 
      break;
    case kFShape_rectangle:
      gmf_size = fsize+1;
      for(ii=0;ii<gmf_size;ii++) 
      {
        gmf_coeff[ii]=1.; 
        gmf_x[ii]=double(ii)*xstep;
      }
      gmf_size--;
      break;
    case kFShape_triangle:
      gmf_size = fsize+1;
      for(ii=0;ii<gmf_size/2;ii++) 
      {
        gmf_coeff[ii]=double(ii+1);
        gmf_x[ii]=double(ii)*xstep;
      }
      for(ii=0;ii<=gmf_size/2;ii++) 
      {
        gmf_coeff[gmf_size/2+ii]=double(gmf_size/2-ii+1.);
        gmf_x[ii]=double(ii)*xstep;
      }
      gmf_size--;
      break;
    default:
      gmf_size = 0;
      return gmf_size;
  }
  gmf_size--;// the rightmost point is used only for interpolation
  printf("gmf_size %i,ileft %i, max=%f @ %i\n",gmf_size,ifleft,max,imax);
  
  // L2-norm
  for(l2=0., ii=0;ii<gmf_size;ii++) {l2 += gmf_coeff[ii]*gmf_coeff[ii];}
  for(ii=0;ii<gmf_size;ii++) {gmf_coeff[ii] /= sqrt(l2);}
  
  //allocate working buffer;
  gmf_out_size = size+gmf_size;
  gmf_out = (double*) malloc(sizeof(double)*(gmf_out_size));
  if(gmf_out==0) {printf("Could not allocate gmf_out\n"); exit(1);}
  
  *coeff = gmf_coeff; *mfx = gmf_x;
  //printf("returning %x into %x\n",int(gmf_coeff), int(coeff));
  return gmf_size;
}

#ifndef FILTER_NONUNIFORM
  void mfi_set(double *xx, double *yy, int nn){} // xx not used for uniform samples
  
  int mfilter_filter(const double *xx, const double *yy, const int size, double *out, double *peak_yy, double *peak_position, int *peak_idx)
  {
    int ii,jj;
    *peak_idx = 0;
    *peak_yy = 0.;
    for(ii=0;ii<size;ii++)
    {
      for(gmf_out[ii]=0.,jj=0;jj<gmf_size;jj++)
        gmf_out[ii] += yy[ii+jj]*gmf_coeff[jj];
      if(ii+jj>size) continue;
      if(gmf_out[ii] > *peak_yy) {*peak_yy = gmf_out[ii]; *peak_idx = ii;}
    }
    if(out) for(ii=0;ii<size;ii++) out[ii] = gmf_out[ii];
    *peak_position = gmf_out[*peak_idx];
    /*if(gverb&0x80)
    {
      printf("mfilter_filter[%i]:\n",gmf_size);
      for(ii=0;ii<gmf_size;ii++) printf("%i:%f\n",ii,gmf_coeff[ii]);
      printf("yy:\n"); for(ii=0;ii<size;ii++) printf("%i:%f ",ii,yy[ii]); printf("\n");
      printf("output:\n"); for(ii=0;ii<size;ii++) printf("%i:%f ",ii,gmf_out[ii]); printf("\n");
    }*/
    return 0;
  }
#else
  //extern double mfilter_interpolated(double xx);
  ///double mfi_coeffy[kMaxFilterLength], mfi_coeffx[kMaxFilterLength];
  //int mfi_n;
  gsl_spline *gspline=NULL;
  gsl_interp_accel *gacc=NULL;

  double mfilter_interpolated(double xx)
  {
    double v;
    //printf("->mfi(%.3f) ",xx);
    if(xx < gmf_x[0] || gmf_x[gmf_size-1] < xx ) 
      printf("\nxx=%f out of range [%f..%f]\n",xx,gmf_x[0],gmf_x[gmf_size-1]);
    v = gsl_spline_eval (gspline, xx, gacc);
    //printf("=%f\n",v);
    return v;
  }
  void mfi_set(double *xx, double *yy, int nn)
  {
    int ii;
    nn++; // for interpolation need one additional point, the rightmost one
    //for(ii=0;ii<nn;ii++) {mfi_coeffy[ii] = yy[ii]; mfi_coeffx[ii] = double(ii);}
    //mfi_n = nn;

    gacc = gsl_interp_accel_alloc ();
    gspline = gsl_spline_alloc (gsl_interp_cspline, nn);
    gsl_spline_init (gspline, xx, yy, nn);
    printf("Filter interpolator [%i] is created:\n",nn);
    for(ii=0;ii<nn;ii++) {gmf_x[ii] = xx[ii];}
    //for(ii=0;ii<nn;ii++) printf("x=%.3f,y=%.3f\n",xx[ii],yy[ii]);
    //printf("Interpolator check:\n");
    //double x;
    //for(x=xx[0];x<xx[2];x += (xx[2]-xx[0])/10.)                  
    //  printf("x=%.3f,iy=%.3f\n",x,mfilter_interpolated(x));
  }

  int mfilter_filter(const double *xx, const double *yy, const int size, double *out, double *peak_yy, double *peak_position, int *peak_idx)
  {
    int ii,jj;
    *peak_idx = 0;
    *peak_yy = 0.;
    //printf("->filter(%.3f,%.3f,%i...\n",yy[0],xx[0],size);
    for(ii=0;ii<size;ii++)
    {
      for(gmf_out[ii]=0.,jj=0;jj<gmf_size-1;jj++)
      {
        //printf("ii,jj=%i,%i,x1=%.3f,x2=%.3f\n",ii,jj,xx[ii],xx[ii+jj]);
        gmf_out[ii] += yy[ii+jj]*mfilter_interpolated(xx[ii+jj]-xx[ii]);
      }
      //printf("gmf[%i/%.3f]=%.3f\n",ii,xx[ii],gmf_out[ii]);
      if(ii+gmf_size>gmf_out_size) continue;
      if(gmf_out[ii] > *peak_yy) {*peak_yy = gmf_out[ii]; *peak_idx = ii;}
    }
    *peak_position = xx[*peak_idx];
    if(out) for(ii=0;ii<size;ii++) out[ii] = gmf_out[ii];
    return 0;
  }
  /*
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
  */
#endif

int mfilter_delete()
{
  printf("Deleting filter\n");
  if(gmf_out != NULL) free(gmf_out); gmf_out = NULL;
  if(gmf_coeff != NULL) free(gmf_coeff); gmf_coeff = NULL;
  if(gmf_x != NULL) free(gmf_coeff); gmf_x = NULL;
#ifdef FILTER_NONUNIFORM
    if(gspline != NULL) gsl_spline_free (gspline);
    if(gacc != NULL) gsl_interp_accel_free (gacc);
#endif
  return 0;
}