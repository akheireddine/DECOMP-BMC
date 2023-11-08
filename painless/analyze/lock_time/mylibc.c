#include "mylibc.h"

void logDiff(timespec start, timespec end)
{
   timespec delta;
   if ((end.tv_nsec-start.tv_nsec)<0) {
      delta.tv_sec = end.tv_sec-start.tv_sec-1;
      delta.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
   } else {
      delta.tv_sec = end.tv_sec-start.tv_sec;
      delta.tv_nsec = end.tv_nsec-start.tv_nsec;
   }
   printf("%d\n", delta.tv_sec*1000000000+delta.tv_nsec);
}


int pthread_mutex_lock(pthread_mutex_t * mutex) {
   timespec time1, time2;
   clock_gettime(CLOCK_MONOTONIC_RAW, &time1);
   int ret = (libc_pthread_mutex_lock())(mutex);
   clock_gettime(CLOCK_MONOTONIC_RAW, &time2);
   logDiff(time1, time2);
   return ret;
}
