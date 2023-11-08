// Copyright (C) 2011 Microsoft Research
//RM Martins, 2012
/*
 * This file contains the old mpi commands
 * The simplest_worker has it own mpi commands
 * IN USE BUT IT IS BEING REFACTORED TO STATIC_WORKER.H
 */

#ifndef _MPI_COMMANDS_H_
#define _MPI_COMMANDS_H_

#include <vector>
#include <mpi.h>

#define MPI_MASTER 0

typedef enum 
{
  MPI_OK = 0,
  MPI_QUIT
} MPI_COMMAND;

typedef enum
{
  MPI_DATA_PROBLEM = 0,
  MPI_CMD,
  MPI_DATA_QUIT,
  MPI_DATA_CLAUSE,
  MPI_DATA_TRAIL,
  MPI_DATA_SHARED,
  MPI_DATA_INTERPOLANT,
  MPI_DATA_SHARED_VARIABLES,
  MPI_DATA_DISTRIBUTION,
  MPI_DATA_MASTER,
  MPI_DATA_DYNAMIC_SIZE,
  MPI_DATA_OCCURRENCES
} MPI_TAG;

class MPI_User {
public:
    MPI_User();
    ~MPI_User();

    void MPI_safe_send(const std::vector<int> &buf, int to, int tag);
    void MPI_safe_send(MPI_COMMAND cmd, int to, int tag=MPI_CMD);
      
    void MPI_wait_for(int * &buf, int from, MPI_Status &s);

protected:
    int* recv_buffer;
    int recv_buffer_size;
};

#endif