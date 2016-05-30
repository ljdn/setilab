#define _GNU_SOURCE
#include <sched.h>    // for processor affinity
#include <unistd.h>   // unix standard apis
#include <pthread.h>  // pthread api

#include <stdlib.h>
#include <stdio.h>    


// The following code reads the intel/amd cycle counter
// This is the number of cycles since the machine was started
// Note that with billions of cycles per second, this is the
// most accurate timer on the machine
unsigned long long int rdtsc(void)
{
  unsigned a, d;

  // This is inline assembler - rdtsc is the specific instruction
  // we are assembling here
  __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));

  return ((unsigned long long)a) | (((unsigned long long)d) << 32);;
}

int     vector_len;   // length of vector we will sum
double *vector;       // the vector we will sum

int        num_threads;     // number of threads we will use
int        num_procs;       // number of processors we will use 
pthread_t *tid;             // array of thread ids
double    *partial_sum;     // partial sums, one for each processor

void *worker(void *arg)
{
  long myid = (long)arg;
  int i;
  int blocksize=vector_len / num_threads; // note: floor

  int mystart, myend;

  // This figures out the chunk of the vector I should
  // work on based on my id
  mystart = myid*blocksize;
  if (myid==(num_threads-1)) { // last processor
    // the last processor will take care of the leftover
    // elements of the vector, in case num_threads doesn't
    // divide vector_len
    myend = vector_len;
  } else {
    myend = (myid+1) * blocksize;
  }

  // Now I sum that chunk and put the result in partial_sum
  partial_sum[myid]=0.0;
  for (i=mystart;i<myend;i++) {
    partial_sum[myid] += vector[i];
  }

  // Done.  The master thread will sum up the partial sums
  pthread_exit(NULL);           // finish - no return value
  
}


int main(int argc, char *argv[]) 
{
  long i;
  long rc;
  double sequential_sum;
  double parallel_sum;
  unsigned long long start, sequential_time, parallel_time;
  
  if (argc!=4) {
    fprintf(stderr,"usage: parallel-sum-ex number-of-threads number-of-procs length-of-vector\n");
    exit(-1);
  }

  num_threads=atoi(argv[1]);      // number of threads
  num_procs=atoi(argv[2]);        // numer of processors to use
  vector_len=atoi(argv[3]);       // length of vector to sum

  vector = (double *) malloc(sizeof(double)*vector_len);
  tid = (pthread_t *) malloc(sizeof(pthread_t)*num_threads);
  partial_sum = (double *) malloc(sizeof(double)*num_threads);

  if (!vector || !tid || !partial_sum) { 
    fprintf(stderr, "cannot allocate memory\n");
    exit(-1);
  }
	
  // populate vector
  for (i=0;i<vector_len;i++) { 
    vector[i] = (double)i;
  }

  // sequential sum
  start=rdtsc();
  sequential_sum=0.0;
  for (i=0;i<vector_len;i++) { 
    sequential_sum+=vector[i];
  }
  sequential_time=rdtsc()-start;


  start=rdtsc();

  // launch threads
  // thread i will compute partial_sum[i], which will sum
  // vector[i*ceiling(vector_size/num_threads)  
  //          through
  //        (i+1)*floor(vector_size/num_threads) ]
  // the last thread will also handle the additional elements
  for (i=0;i<num_threads;i++) {
    rc=pthread_create( &(tid[i]), // thread id gets put here
		       NULL,      // use default attributes
		       worker,    // thread will begin in this function
		       (void*) i  // we'll give it i as the argument
	             );
    if (rc!=0) { 
      perror("Failed to start thread");
      exit(-1);
    }
  }

  // now we will join all the threads
  for (i=0;i<num_threads;i++) {
    rc=pthread_join(tid[i],NULL);   // 
    if (rc!=0) { 
	perror("join failed");
	exit(-1);
    }
  }

  // Now we know that the partial_sums are done, so we'll add them up
  parallel_sum=0.0;
  for (i=0;i<num_threads;i++) {
    parallel_sum+=partial_sum[i];
  }

  parallel_time=rdtsc()-start;

  printf("Sequential sum:   %lf (%llu cycles)\n"
	 "Parallel sum:     %lf (%llu cycles)\n", 
	 sequential_sum, sequential_time, parallel_sum, parallel_time);

  return 0;
}
