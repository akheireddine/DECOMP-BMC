#ifndef MY_LIB_C
#define MY_LIB_C

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

#define LIBC ((void *)-1L) // Incorrect handle, so libc by default

int (*fn_lock_libc) (pthread_mutex_t * mutex) = NULL;

#define libc_pthread_mutex_lock() (fn_lock_libc != NULL ? fn_lock_libc : \
   (fn_lock_libc = (int (*)(pthread_mutex_t*))dlsym(LIBC, "pthread_mutex_lock")))

#endif //MY_LIB_C
