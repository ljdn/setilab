#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"


void usage()
{
    printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands\n");
}

double avg_power(double *data, int num)
{
    int i;
    double ss;

    ss=0;
    for (i=0;i<num;i++) {
	ss += data[i]*data[i];
    }

    return ss/num;
}

double max_of(double *data, int num)
{
    double m=data[0];
    int i;

    for (i=1;i<num;i++) {
	if (data[i]>m) { m=data[i]; }
    }
    return m;
}

double avg_of(double *data, int num)
{
    double s=0;
    int i;

    for (i=0;i<num;i++) {
	s+=data[i];
    }
    return s/num;
}

void remove_dc(double *data, int num)
{
  int i;
  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (i=0;i<num;i++) {
    data[i] -= dc;
  }
}


int analyze_signal(signal *sig, int filter_order, int num_bands, double *lb, double *ub)
{
    double Fc, bandwidth;
    double filter_coeffs[filter_order+1];
    double signal_power;
    double band_power[num_bands];

    double start, end;

    unsigned long long tstart, tend;

    resources rstart, rend, rdiff;

    int band;

    Fc=(sig->Fs)/2;
    bandwidth = Fc / num_bands;

    remove_dc(sig->data,sig->num_samples);

    signal_power = avg_power(sig->data,sig->num_samples);

    printf("signal average power:     %lf\n", signal_power);

    get_resources(&rstart,THIS_PROCESS);
    start=get_seconds();
    tstart = get_cycle_count();

    for (band=0;band<num_bands;band++) {
	// Make the filter
	generate_band_pass(sig->Fs,
			   band*bandwidth+0.0001, // keep within limits
			   (band+1)*bandwidth-0.0001,
			   filter_order,
			   filter_coeffs);
	hamming_window(filter_order,filter_coeffs);

	// Convolve
	convolve_and_compute_power(sig->num_samples,
				   sig->data,
				   filter_order,
				   filter_coeffs,
				   &(band_power[band]));


    }
    tend = get_cycle_count();
    end = get_seconds();
    get_resources(&rend,THIS_PROCESS);

    get_resources_diff(&rstart, &rend, &rdiff);

    // Pretty print results
    double max_band_power = max_of(band_power,num_bands);
    double avg_band_power = avg_of(band_power,num_bands);
    int i;
    int wow=0;

#define MAXWIDTH 40

#define THRESHOLD 2.0

#define ALIENS_LOW   50000.0
#define ALIENS_HIGH  150000.0

    *lb=*ub=-1;

    for (band=0;band<num_bands;band++) {
      double band_low = band*bandwidth+0.0001;
      double band_high = (band+1)*bandwidth-0.0001;

      printf("%5d %20lf to %20lf Hz: %20lf ",
	     band, band_low, band_high, band_power[band]);

      for (i=0;i<MAXWIDTH*(band_power[band]/max_band_power);i++) {
	printf("*");
      }

      if ( (band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
	   (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {

	// band of interest

	if (band_power[band] > THRESHOLD * avg_band_power) {
	  printf("(WOW)");
	  wow=1;
	  if (*lb<0) { *lb=band*bandwidth+0.0001; }
	  *ub = (band+1)*bandwidth-0.0001;
	} else {
	  printf("(meh)");
	}
      } else {
	printf("(meh)");
      }

      printf("\n");
    }

    printf("Resource usages:\n"
	   "User time        %lf seconds\n"
	   "System time      %lf seconds\n"
	   "Page faults      %ld\n"
	   "Page swaps       %ld\n"
	   "Blocks of I/O    %ld\n"
	   "Signals caught   %ld\n"
	   "Context switches %ld\n",
	   rdiff.usertime,
	   rdiff.systime,
	   rdiff.pagefaults,
	   rdiff.pageswaps,
	   rdiff.ioblocks,
	   rdiff.sigs,
	   rdiff.contextswitches);


    printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\nNote that cycle count only makes sense if the thread stayed on one core\n", tend-tstart, cycles_to_seconds(tend-tstart), timing_overhead());
    printf("Analysis took %lf seconds by basic timing\n", end-start);

    return wow;

}

int main(int argc, char *argv[])
{
    signal *sig;
    double Fs;
    char sig_type;
    char *sig_file;
    int filter_order;
    int num_bands;
    double start, end;

    if (argc!=6) {
	usage();
	return -1;
    }

    sig_type = toupper(argv[1][0]);
    sig_file = argv[2];
    Fs = atof(argv[3]);
    filter_order = atoi(argv[4]);
    num_bands = atoi(argv[5]);

    assert(Fs>0.0);
    assert(filter_order>0 && !(filter_order & 0x1));
    assert(num_bands>0);

    printf("type:     %s\n"
	   "file:     %s\n"
	   "Fs:       %lf Hz\n"
	   "order:    %d\n"
	   "bands:    %d\n",
	   sig_type=='T' ? "Text" : sig_type=='B' ? "Binary" : sig_type=='M' ? "Mapped Binary" : "UNKNOWN TYPE",
	   sig_file,
	   Fs,
	   filter_order,
	   num_bands);

    printf("Load or map file\n");

    switch (sig_type) {
	case 'T':
	    sig = load_text_format_signal(sig_file);
	    break;

	case 'B':
	    sig = load_binary_format_signal(sig_file);
	    break;

	case 'M':
	    sig = map_binary_format_signal(sig_file);
	    break;

	default:
	    printf("Unknown signal type\n");
	    return -1;
    }

    if (!sig) {
	printf("Unable to load or map file\n");
	return -1;
    }

    sig->Fs=Fs;

    if (analyze_signal(sig,filter_order,num_bands,&start,&end)) {
	printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n",start,end,(end+start)/2.0);
    } else {
	printf("no aliens\n");
    }

    free_signal(sig);

    return 0;
}
