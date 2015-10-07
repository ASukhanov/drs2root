enum
{
  kFShape_arbitrary = 1,
  kFShape_rectangle = 2,
  kFShape_triangle = 3
};

int mfilter_create(const double *amplitude, const int size, const int type, const int fsize, double **mfilter_coeff);

int mfilter_filter(const double *input, const int size, double *out, double *peak_amplitude, 
                   double *peak_position);

int mfilter_delete();
