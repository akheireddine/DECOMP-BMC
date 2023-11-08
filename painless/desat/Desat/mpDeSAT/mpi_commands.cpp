// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011
/*
 * This file contains the old mpi commands
 * The simplest_worker has it own mpi commands
 * IN USE BUT IT IS BEING REFACTORED TO STATIC_WORKER.C
 */

#include <iostream>
#include <mpi.h>
#include "mpi_commands.h"


MPI_User::MPI_User() : 
    recv_buffer(0),
    recv_buffer_size(0)
{
}

MPI_User::~MPI_User() {
    if (recv_buffer)
        delete [] recv_buffer;
}

void MPI_User::MPI_safe_send(const std::vector<int> &buf, int to, int tag)
{  
	while (MPI_Send((void*) buf.data(), buf.size(), MPI_INT, to, tag, MPI_COMM_WORLD)!=MPI_SUCCESS)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		std::cout << "Worker " << rank << " send failure." << std::endl;
	}
}

void MPI_User::MPI_safe_send(MPI_COMMAND cmd, int to, int tag)
{  
	std::vector<int> tmp(1, cmd);
	MPI_safe_send(tmp, to, tag);
}

void MPI_User::MPI_wait_for(int * & buf, int from, MPI_Status &s)
{
	/* New version that breaks the old static version: mpdesat_worker.cpp, mpdesat_node.cpp and mpdesat_master.cpp

	while (MPI_Probe(from, MPI_ANY_TAG, MPI_COMM_WORLD, &s)!=MPI_SUCCESS)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		std::cout << "Node " << rank << ": Probe failure while waiting for command." << std::endl;
	}

	if (recv_buffer_size <= s.count) {
		if (recv_buffer != 0)
			delete [] recv_buffer;
		recv_buffer = new int[s.count+1];
		recv_buffer_size = s.count;
	}	

	while (MPI_Recv(recv_buffer + 1, s.count, MPI_INT, s.MPI_SOURCE, s.MPI_TAG, MPI_COMM_WORLD, &s)!=MPI_SUCCESS)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		std::cout << "Node " << rank << ": Receive failure while waiting for command." << std::endl;
	}

    recv_buffer[0] = s.count;
    buf = recv_buffer;
	
	*/

	while (MPI_Probe(from, MPI_ANY_TAG, MPI_COMM_WORLD, &s)!=MPI_SUCCESS)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		std::cout << "Node " << rank << ": Probe failure while waiting for command." << std::endl;
	}

	buf = new int[s.count];

	while (MPI_Recv(buf, s.count, MPI_INT, s.MPI_SOURCE, s.MPI_TAG, MPI_COMM_WORLD, &s)!=MPI_SUCCESS)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		std::cout << "Node " << rank << ": Receive failure while waiting for command." << std::endl;
	}
}