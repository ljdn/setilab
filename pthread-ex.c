#define _GNU_SOURCE
#include <sched.h>    // for processor affinity
#include <unistd.h>   // unix standard apis
#include <pthread.h>  // pthread api

#include <stdlib.h>
#include <stdio.h>    

long numproc;

pthread_t *tid;  // array of thread ids

void *worker(void *arg)
{
  long myid = (long)arg;
  long mytid = tid[myid];   // notice access to shared variable
  long mytime = 1+rand()%10; 

  printf("Hello from thread %ld (tid %ld)\n",myid,mytid);
  printf("Thread %ld is putting itself onto procesor %ld\n", myid, myid%numproc);

  // put ourselves on the desired processor
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(myid%numproc,&set);
  if (sched_setaffinity(0,sizeof(set),&set)<0) { // do it
    perror("Can't setaffinity");  // hopefully doesn't fail
    exit(-1);
  }

  printf("Thread %ld now sleeping for %ld seconds\n",myid,mytime);
  
  sleep(mytime);

 
  printf("Thread %ld done sleeping and now exiting\n",myid);
	 
  pthread_exit(NULL);           // finish - no return value
  
}

int main(int argc, char *argv[]) 
{
  long i;
  long rc;
  long num;
  long num_started;


  if (argc!=3 || atoi(argv[1])<1) { 
    fprintf(stderr,"usage: pthread-ex number-of-threads number-of-procs\n");
    exit(-1);
  }

  num=atoi(argv[1]);  // number of threads
  numproc=atoi(argv[2]); // numer of processors to use
  if (numproc>32) { 
    numproc=32;
  }

  tid = (pthread_t *) malloc(sizeof(pthread_t)*num);

  srand(time(NULL)); // seed random number generator


  printf("Starting %ld softare threads on %ld processors (hardware threads)\n",num, numproc);

  num_started=0;
  for (i=0;i<num;i++) {
    rc=pthread_create( &(tid[i]), // thread id gets put here
		       NULL,      // use default attributes
		       worker,    // thread will begin in this function
		       (void*) i  // we'll give it i as the argument
	             );
    if (rc==0) { 
      printf("Started thread %ld, tid %lu\n",i,tid[i]);
      num_started++;
    } else {
      printf("Failed to start thread %ld\n",i);
      perror("Failed to start thread");
      tid[i]=0xdeadbeef;
    }
  }

  printf("Finished starting threads (%ld started)\n",num_started);

  printf("Now joining\n");

  for (i=0;i<num;i++) {
    if (tid[i]!=0xdeadbeef) { 
      printf("Joining with %ld, tid %lu\n", i, tid[i]);
      rc=pthread_join(tid[i],NULL);   // 
      if (rc!=0) { 
	printf("Failed to join with %ld!\n",i);
	perror("join failed");
      } else { 
	printf("Done joining with %ld\n",i);
      }
    } else { 
      printf("Skipping %ld (wasn't started successfully)\n",i);
    }
  }

  printf("Done!\n");

  return 0;
}
