#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

//int num_bands;
int num_procs;
int num_threads;
pthread_t *tid;             // array of thread ids

typedef struct my_args {
  int thread;
  double bandwidth;
  int thread_size;
  int filter_order;
  double* filter_coeffs;
  signal* sig;
  signal* output;
  double* band_power;
} my_args;

void usage()
{
  printf("usage: p_band_scan text|bin|mmap signal_file Fs filter_order num_bands num_processsors num_threads\n");
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

// threads begin here
void *worker(void *input) {
  my_args *args = (my_args*) input;

  // make filter
  int i;
  int coeff = 0;
  printf("thread: %d\n" "size: %d\n", args->thread, args->thread_size);
  for(i=(args->thread)*(args->thread_size); i<(args->thread+1)*(args->thread_size); i++) {
    // int coeff = (i-1) * (args->filter_order+1);
    printf("band: %d\n", i);
    generate_band_pass(args->sig->Fs,
                      (i*(args->bandwidth))+0.0001,
                      (i+1)*(args->bandwidth)-0.0001,
                      (args->filter_order),
                      &(args->filter_coeffs[coeff]));
    hamming_window((args->filter_order), &(args->filter_coeffs[coeff]));

    // heavy lifting!
    convolve_and_compute_power(args->sig->num_samples,
                              args->sig->data,
                              args->filter_order,
                              &(args->filter_coeffs[coeff]),
                              &(args->band_power[i]));
    // increment to next batch of coeffs
    coeff++;
  }

  // done, exit with no return val
  pthread_exit(NULL);

}

int analyze_signal(signal *sig, int filter_order, int num_bands, int num_procs, int num_threads, double *lb, double *ub)
{
  double Fc, bandwidth, signal_power;
  // double filter_coeffs[(filter_order+1)*(num_bands+1)];
  double band_power[num_bands];
  signal *output[num_bands];
  long rc;

  double start, end;

  unsigned long long tstart, tend;

  resources rstart, rend, rdiff;

  tid = (pthread_t *) malloc(sizeof(pthread_t)*num_threads);

  if (!tid) {
    fprintf(stderr, "Cannot allocate memory\n");
    exit(-1);
  }

  Fc=(sig->Fs)/2;
  bandwidth = Fc / num_bands;

  int i;
  for(i=0; i<num_bands; i++) {
    output[i] = allocate_signal(sig->num_samples, sig->Fs, 0);
    if (!output[i]) {
      printf("Out of memory\n");
      return 0;
    }
  }

  remove_dc(sig->data,sig->num_samples);

  signal_power = avg_power(sig->data,sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  get_resources(&rstart,THIS_PROCESS);
  start=get_seconds();
  tstart = get_cycle_count();

  // if more bands than threads, separate into groups
  int howmany_threads;
  int thread_size;
  if (num_threads < num_bands) {
    howmany_threads = num_threads;
    thread_size = floor(num_bands/num_threads);
  } else {
    howmany_threads = num_bands;
    thread_size = 1;
  }


  for (int thread=0; thread<howmany_threads; thread++) {
    // printf("thread num: %d\n", thread);
    if (thread == howmany_threads-1) {
      thread_size = num_bands - (thread_size*thread);
    }
    // last thread_size gets the rest of the bands if it wasn't split evenly

    double filter_coeffs[(filter_order+1) * thread_size];
    // threads begin in worker function
    my_args my_data;
    my_data.thread = thread;
    my_data.bandwidth = bandwidth;
    my_data.thread_size = thread_size;
    my_data.filter_order = filter_order;
    my_data.filter_coeffs = filter_coeffs;
    my_data.sig = sig;
    my_data.output = *output;
    my_data.band_power = band_power;
    rc = pthread_create( &(tid[thread]),    // thread id
                         NULL,             // default attributes
                         worker,           // thread begins here
                         &my_data);        // argument
    if (rc != 0) {
      perror("Failed to start thread");
      exit(-1);
    }
    printf("thread (analyze signal): %d\n", thread);

  }

  // join threads
  for (int thread=0;thread<num_threads;thread++) {
    rc = pthread_join(tid[thread], NULL);
    if (rc != 0) {
      perror("Join failed");
      exit(-1);
    }
  }

  tend = get_cycle_count();
  end = get_seconds();
  get_resources(&rend,THIS_PROCESS);

  get_resources_diff(&rstart, &rend, &rdiff);

  // Pretty print results
  double max_band_power = max_of(band_power,num_bands);
  double avg_band_power = avg_of(band_power,num_bands);

  int wow=0;

  #define MAXWIDTH 40

  #define THRESHOLD 2.0

  #define ALIENS_LOW   50000.0
  #define ALIENS_HIGH  150000.0

  *lb=*ub=-1;

  for (int band=0;band<num_bands;band++) {
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

  if (argc!=8) {
    usage();
    return -1;
  }

  sig_type = toupper(argv[1][0]);
  sig_file = argv[2];
  Fs = atof(argv[3]);
  filter_order = atoi(argv[4]);
  num_bands = atoi(argv[5]);
  num_procs = atoi(argv[6]);
  num_threads = atoi(argv[7]);

  assert(Fs>0.0);
  assert(filter_order>0 && !(filter_order & 0x1));
  assert(num_bands>0);

  printf("type:     %s\n"
         "file:     %s\n"
         "Fs:       %lf Hz\n"
         "order:    %d\n"
         "bands:    %d\n"
         "processes: %d\n"
         "threads: %d\n",
         sig_type=='T' ? "Text" : sig_type=='B' ? "Binary" : sig_type=='M' ? "Mapped Binary" : "UNKNOWN TYPE",
         sig_file,
         Fs,
         filter_order,
         num_bands,
         num_procs,
         num_threads);

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

  if (analyze_signal(sig,filter_order,num_bands, num_procs, num_threads, &start,&end)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n",start,end,(end+start)/2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}
