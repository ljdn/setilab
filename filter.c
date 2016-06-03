#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "filter.h"

int generate_low_pass(double Fs, double Fc,
		      int order, double coeffs[])
{
  double Ft;
  int n;

  assert(order>0 && !(order&0x1));
  assert(Fs>0 && Fc>0 && Fc<Fs/2);

  Ft = Fc/Fs;

  for (n=0;n<=order;n++) {
    if (n == order/2) {
      coeffs[n] = 2*Ft;
    } else {
      coeffs[n] = sin(2*M_PI*Ft*(n-order/2))/(M_PI*(n-order/2));
    }
  }

  return 0;
}

int generate_high_pass(double Fs, double Fc,
		       int order, double coeffs[])
{
  double Ft;
  int n;

  assert(order>0 && !(order&0x1));
  assert(Fs>0 && Fc>0 && Fc<Fs/2);

  Ft = Fc/Fs;

  for (n=0;n<=order;n++) {
    if (n == order/2) {
      coeffs[n] = 1 - 2*Ft;
    } else {
      coeffs[n] = -sin(2*M_PI*Ft*(n-order/2))/(M_PI*(n-order/2));
    }
  }
  return 0;
}

int generate_band_pass(double Fs, double Fcl, double Fch,
		       int order, double coeffs[])
{
  double Ftl, Fth;
  int n;

  assert(order>0 && !(order&0x1));
	printf("%f\n" "%f\n" "%f\n", Fs, Fcl, Fch);;
  assert(Fs>0 && Fcl>0 && Fcl<Fs/2 && Fch>0 && Fch<Fs/2);

  Ftl = Fcl/Fs;
  Fth = Fch/Fs;

  for (n=0;n<=order;n++) {
    if (n == order/2) {
      coeffs[n] = 2*(Fth-Ftl);
    } else {
      coeffs[n] = (sin(2*M_PI*Fth*(n-order/2))/(M_PI*(n-order/2))) -
	(sin(2*M_PI*Ftl*(n-order/2))/(M_PI*(n-order/2)));
    }
  }
  return 0;
}


int generate_band_stop(double Fs, double Fcl, double Fch,
		       int order, double coeffs[])
{
  double Ftl, Fth;
  int n;

  assert(order>0 && !(order&0x1));
  assert(Fs>0 && Fcl>0 && Fcl<Fs/2 && Fch>0 && Fch<Fs/2);

  Ftl = Fcl/Fs;
  Fth = Fch/Fs;

  for (n=0;n<=order;n++) {
    if (n == order/2) {
      coeffs[n] = 1-2*(Fth-Ftl);
    } else {
      coeffs[n] = (sin(2*M_PI*Ftl*(n-order/2))/(M_PI*(n-order/2))) -
	(sin(2*M_PI*Fth*(n-order/2))/(M_PI*(n-order/2)));
    }
  }
  return 0;
}

int hamming_window(int order, double coeffs[])
{
   int n;

  assert(order>0 && !(order&0x1));

  for (n=0;n<=order;n++) {
    coeffs[n] = coeffs[n] * (0.54 - 0.46 * cos(2*M_PI*n/order));
  }

  return 0;

}

// Simple (slow) convolution
// output must be same length as input.  coeffs assumed to be
int convolve(int length, double input_signal[],
	     int order, double coeffs[],
	     double output_signal[])
{
  int i,j;

  for (i=0;i<length;i++) {
    output_signal[i] = 0;
    for (j=order;j>=0;j--) {
      // Use coeff only if there is input signal
      // that matches, otherwise assume input signal
      // is zero (aperiodic model)
      if ((i-j)>=0 && (i-j)<length) {
	// Causal model, use inputs up to this point
	output_signal[i] += input_signal[i-j] * coeffs[j];
      }
    }
  }
  return 0;
}


// Simple (slow) convolution combined with power estimate for output
int convolve_and_compute_power(int length, double input_signal[],
			       int order, double coeffs[],
			       double *power)
{
  int i,j;
  double cur_sum;
  double pow_sum;

  pow_sum = 0;

  for (i=0;i<length;i++) {
    cur_sum = 0;
    for (j=order;j>=0;j--) {
      // Use coeff only if there is input signal
      // that matches, otherwise assume input signal
      // is zero (aperiodic model)
      if ((i-j)>=0 && (i-j)<length) {
	// Causal model, use inputs up to this point
	cur_sum += input_signal[i-j] * coeffs[j];
      }
    }
    pow_sum += cur_sum * cur_sum;
  }

  *power = pow_sum / length;

  return 0;
}

/* below taken from http://www.exstrom.com/journal/sigproc/liir.c */

/**********************************************************************
  binomial_mult - multiplies a series of binomials together and returns
  the coefficients of the resulting polynomial.

  The multiplication has the following form:

  (x+p[0])*(x+p[1])*...*(x+p[n-1])

  The p[i] coefficients are assumed to be complex and are passed to the
  function as a pointer to an array of doubles of length 2n.

  The resulting polynomial has the following form:

  x^n + a[0]*x^n-1 + a[1]*x^n-2 + ... +a[n-2]*x + a[n-1]

  The a[i] coefficients can in general be complex but should in most
  cases turn out to be real. The a[i] coefficients are returned by the
  function as a pointer to an array of doubles of length 2n. Storage
  for the array is allocated by the function and should be freed by the
  calling program when no longer needed.

  Function arguments:

  n  -  The number of binomials to multiply
  p  -  Pointer to an array of doubles where p[2i] (i=0...n-1) is
        assumed to be the real part of the coefficient of the ith binomial
        and p[2i+1] is assumed to be the imaginary part. The overall size
        of the array is then 2n.
*/

static double * binomial_mult (int n, double * p)
{
    int i, j;
    double *a;

    a = (double *)calloc( 2 * n, sizeof(double) );
    if( a == NULL ) return( NULL );

    for( i = 0; i < n; ++i )
    {
	for( j = i; j > 0; --j )
	{
	    a[2*j] += p[2*i] * a[2*(j-1)] - p[2*i+1] * a[2*(j-1)+1];
	    a[2*j+1] += p[2*i] * a[2*(j-1)+1] + p[2*i+1] * a[2*(j-1)];
	}
	a[0] += p[2*i];
	a[1] += p[2*i+1];
    }
    return( a );
}


/**********************************************************************
  ccof_bwlp - calculates the c coefficients for a butterworth lowpass
  filter. The coefficients are returned as an array of doubles

*/
static double * ccof_bwlp (int n)
{
    double * ccof;
    int m, i;

    ccof = (double*)calloc(n + 1, sizeof(double));
    if (!ccof) {
        return NULL;
    }

    ccof[0] = 1.0;
    ccof[1] = (double)n;
    m = n/2;
    for (i = 2; i <= m; ++i) {
        ccof[i] = (double)((n - i + 1) * ccof[i-1] / i);
        ccof[n-i]= ccof[i];
    }

    ccof[n-1] = (double)n;
    ccof[n] = 1.0;

    return ccof;
}



/**********************************************************************
  dcof_bwlp - calculates the d coefficients for a butterworth lowpass
  filter. The coefficients are returned as an array of doubles.

*/
static double * dcof_bwlp (int n, double fcf)
{
    int k;            // loop variables
    double theta;     // M_PI * fcf / 2.0
    double st;        // sine of theta
    double ct;        // cosine of theta
    double parg;      // pole angle
    double sparg;     // sine of the pole angle
    double cparg;     // cosine of the pole angle
    double a;         // workspace variable
    double *rcof;     // binomial coefficients
    double *dcof;     // dk coefficients

    rcof = (double *)calloc(2 * n, sizeof(double));
    if (!rcof) {
        return NULL;
    }

    theta = M_PI * fcf;
    st = sin(theta);
    ct = cos(theta);

    for (k = 0; k < n; ++k) {
        parg = M_PI * (double)(2*k+1)/(double)(2*n);
        sparg = sin(parg);
        cparg = cos(parg);
        a = 1.0 + st*sparg;
        rcof[2*k] = -ct/a;
        rcof[2*k+1] = -st*cparg/a;
    }

    dcof = binomial_mult(n, rcof);
    free(rcof);

    dcof[1] = dcof[0];
    dcof[0] = 1.0;
    for (k = 3; k <= n; ++k) {
        dcof[k] = dcof[2*k-2];
    }

    return dcof;
}

/**********************************************************************
  sf_bwlp - calculates the scaling factor for a butterworth lowpass filter.
  The scaling factor is what the c coefficients must be multiplied by so
  that the filter response has a maximum value of 1.

*/

static double sf_bwlp (int n, double fcf)
{
    int k;         // loop variables
    double omega;     // M_PI * fcf
    double fomega;    // function of omega
    double parg0;     // zeroth pole angle
    double sf;        // scaling factor

    omega = M_PI * fcf;
    fomega = sin(omega);
    parg0 = M_PI / (double)(2*n);

    sf = 1.0;
    for (k = 0; k < n/2; ++k) {
        sf *= 1.0 + fomega * sin((double)(2*k+1)*parg0);
    }

    fomega = sin(omega / 2.0);

    if (n % 2) {
        sf *= fomega + cos(omega / 2.0);
    }

    sf = pow(fomega, n) / sf;

    return sf;
}


/* butterworth lowpass filter. numerator coefficients
 * returned in b, denominator coefficients in a
 * n is order
 * fcf is cutoff freq
 *
 *
 * butter(n, fcf, b, a);
 *
 * NOTE: this has been validated against matlab's output
 */
void butter (int n, double fcf, double ** b, double ** a)
{
    int i = 0;
    double sf;

    /* get numerator coffs */
    *b = ccof_bwlp(n);
    if (!*b) {
        fprintf(stderr, "Error calculating numerator coffs\n");
        return;
    }

    /* get denominator coffs */
    *a = dcof_bwlp(n, fcf);
    if (!*a) {
        fprintf(stderr, "Error calculating denominator coffs\n");
        return;
    }


    sf = sf_bwlp(n, fcf);

    // scale the scaling factor for numerator coffs
    for (i = 0; i < n + 1; i++) {
        (*b)[i] *= sf;
    }

}


/* apply a butterworth filter (both sets of coefficients)
 * y = filter(b, a, x)
 */
void filter (int ord, double * a, double * b,
             int np, double * x, double * y)
{
    int i,j;

    y[0] = b[0] * x[0];

    for (i = 1; i < ord + 1; i++) {

        y[i] = 0.0;

        for (j = 0; j < i + 1; j++) {
            y[i] += b[j] * x[i - j];
        }

        for (j = 0; j < i; j++) {
            y[i] -=  a[j + 1] * y[i - j - 1];
        }
    }

    /* end of initial part */
    for (i = ord + 1; i < np + 1; i++) {

        y[i] = 0.0;

        for (j = 0; j < ord + 1; j++) {
            y[i] += b[j] * x[i - j];
        }

        for (j = 0; j < ord; j++) {
            y[i] -= a[j + 1] * y[i - j - 1];
        }
    }
}

/* y = filtfilt(b, a, x) */
void filtfilt (int ord, double * a, double * b,
               int np, double * x, double * y)
{

    int i;

    filter(ord, a, b, np, x, y);

    /* reverse the series */
    for (i = 0; i < np; i++) {
        x[i] = y[np - i - 1];
    }

    filter(ord, a, b, np, x, y);

    /* put it back */
    for (i = 0; i < np; i++) {
        x[i] = y[np - i - 1];
    }

    for (i = 0; i < np; i++) {
        y[i] = x[i];
    }
}
