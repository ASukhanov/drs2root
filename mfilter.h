// uncomment following if processing Non Uniform Samples
#define FILTER_NONUNIFORM

enum
{
  kFShape_1st_event = 1,
  kFShape_rectangle = 2,
  kFShape_triangle = 3,
  kFShape_external = 4
};

int mfilter_create(const double *x, const double *y, const int n, double xstep, const int type, const int ntaps, double **coefficients, double **coeff_x);
// Important, the size of x and y should be at least n+1


int mfilter_filter(const double *x, const double *y, const int size, double *out, double *peak_amplitude, double *peak_position, int *peak_idx);

//#ifdef FILTER_NONUNIFORM
//int mfilter_filter_nonuniform(const double *input, const double *xx, const int size, double *out, double *peak_amplitude, double *peak_position, int *peak_idx);
//#endif

void mfi_set(double *xx, double *yy, int nn);

int mfilter_delete();
