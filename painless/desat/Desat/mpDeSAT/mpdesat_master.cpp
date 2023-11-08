// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011
// RM Martins, 2012
/*
 * This file contains the master from the static topologies: flat, tree and tree-flat
 * WARNING: due to some changes in the code it may have some memory problems
 * IN USE BUT IT IS BEING REFACTORED TO STATIC_WORKER.C		
 */

#include <mpi.h>

#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <minisat1p.h>
#include <decomposition_batch.h>

#include "mpi_commands.h"
#include "mpdesat_master.h"
#include "returnvalue.h"

#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib") //added

MPDeSATMaster::MPDeSATMaster(ExpressionManager &m, int size, int rank, const char *filename, int c, InterpolationMode i, int master, std::vector<int> worker) :
MPDeSATNode(m, rank, master, worker),
	size(size),
	cores(c),
	filename(filename),
	n_workers(size-1),
	d(NULL),
	globalSolver(NULL),
	sendingTo(-1),
	bufferSize(10485760)
{
	assert(rank==MPI_MASTER);

	globalSolver = new MiniSAT_1p(m, false, MiniSAT_1p::LIFTING);    
	d = new BatchDecomposition();
	d->setPartitions(n_workers);
	assumptions.resize(n_workers);
	interpolants.resize(n_workers, m.mkNil());
	for (unsigned i=0; i< n_workers; i++)
		occurrences.push_back(new VariableOccurrence());

}

MPDeSATMaster::~MPDeSATMaster(void)
{
	if (globalSolver) delete globalSolver;
	if (d) delete d;
}

//send methods
void MPDeSATMaster::sendStatus(bool status)
{
	std::vector<signed> buf;
	status ? buf.push_back(1) : buf.push_back(0);

	for(unsigned id = 0; id < workersId.size(); id++)
		MPI_safe_send(buf,workersId[id],MPI_DATA_QUIT);
}

void MPDeSATMaster::sendSharedVariables()
{
	std::vector<signed> buf;
	buf.push_back(0);
	for(unsigned j = 1; j <= maxVar; j++)
	{
		if(sharedVariables.isShared(j)) buf.push_back(1);
		else buf.push_back(0);
	}

	buf[0] = buf.size();

	for(unsigned id = 0; id < workersId.size(); id++)
		MPI_safe_send(buf,workersId[id],MPI_DATA_SHARED_VARIABLES);
}

void MPDeSATMaster::sendTrail(void)
{

	sendBuffer.clear();
	sendBuffer.push_back(0);

	for (int v=1; v<=(signed)maxVar; v++)
	{
		if (sharedVariables.isShared(v))
		{        
			ModelValue mv = globalSolver->get(v);
			if (mv != M_UNDEF){
				if(mv==M_TRUE) sendBuffer.push_back(v);
				else sendBuffer.push_back(-v);
			} 
		}

	}
	sendBuffer[0] = sendBuffer.size();

	for(unsigned id=0; id < workersId.size(); id++)
		MPI_safe_send(sendBuffer, workersId[id], MPI_DATA_TRAIL);

	sendBuffer.clear();

}

//initialization methods
void MPDeSATMaster::setClauseMax(unsigned n)
{

	d->setClauseMax(n);
	maxClauses = n;

	std::vector<int> problem_statment(2);
	problem_statment[0] = maxVar;
	problem_statment[1] = n;

	//sending to all workers the information regarding the number of variables and number of clauses of the problem
	for(unsigned id = 1; id <= n_workers; id++)
		MPI_safe_send(problem_statment, id, MPI_DATA_PROBLEM);

}

bool MPDeSATMaster::addClause(const std::vector<signed> &literals)
{ 
	d->setPartitions(workersId.size());
	unsigned w = d->where(literals);

	//printClause(literals);
	if(sendBuffer.empty()){
		sendBuffer.push_back(0); //reserve position for the size of the array
		sendBuffer.push_back(0); //reserve position for the number of clauses in the array
		num_clauses = 0;
	}

	//TO DO: improve the way the splitting of the clause buffer is done
	if((sendingTo != w && sendingTo !=-1) || sendBuffer.size() > bufferSize )
	{
		//must send my buffer since I will change worker or I reached the predefined limit for the buffer
		//NOTE: the buffer size can be larger than the predefined limit, where the surplus is the size-1 of the clause that was last read
		sendBuffer[0] = sendBuffer.size();
		sendBuffer[1] = num_clauses;
		MPI_safe_send(sendBuffer, workersId[sendingTo], MPI_DATA_CLAUSE);
		sendBuffer.clear();

	} 

	if(sendBuffer.empty()){
		sendBuffer.push_back(0); //reserve position for the size of the array
		sendBuffer.push_back(0); //reserve position for the number of clauses in the array
		num_clauses = 0;
	}
	num_clauses++;
	sendBuffer.push_back(literals.size());
	sendingTo = w;
	for(unsigned i=0; i < literals.size(); i++){
		sendBuffer.push_back(literals[i]);

		std::vector<unsigned>::iterator it = std::find(workersId.begin(), workersId.end(), workersId[sendingTo]);
		unsigned workerPosition = std::distance( workersId.begin(), it);

		occurrences[workerPosition]->setOccurs(literals[i]);	
	}

	return true;
}

void MPDeSATMaster::setVariableMax(unsigned n)
{
	if (n<maxVar) return;

	maxVar = n;
	globalSolver->setVariableMax(n);  
	d->setVariableMax(n);

}

void MPDeSATMaster::stopReadDimacsFile()
{
	//If sendBuffer is not empty then there are still some clauses in cache that need to be sent to the last worker
	if(!sendBuffer.empty()){
		sendBuffer[0] = sendBuffer.size();
		sendBuffer[1] = num_clauses;
		MPI_safe_send(sendBuffer, workersId[sendingTo], MPI_DATA_CLAUSE);
		sendBuffer.clear();
	}

	//Master flaging workers that all clauses have been read and distributed
	std::vector<signed> lit(1);
	lit[0] = 0; //this is done by sending an empty array (only one member that describes the size 0 of the array)
	for(unsigned id=0; id < workersId.size(); id++)
		MPI_safe_send(lit, workersId[id], MPI_DATA_CLAUSE);
}

//solving methods
bool MPDeSATMaster::solve(void)
{
	clock_t start = clock();

	ModelValue res = M_UNDEF;

	clock_t before = clock();
	readDimacsFile(filename);

	stopReadDimacsFile(); //flags the workers that there is no more clauses to be read

	mpi_loadTime += clock() - before;

	for(unsigned i=0; i < workersId.size(); i++)
	{
		sharedVariables.add(occurrences[i]);
		//for(unsigned j=1; j <= maxVar; j++)
		//	if(sharedVariables.occurs(j,i)) std::cout << "contains Node " << workersId[i] << " has variable " << j << std::endl;
	}

	sharedVariables.update();

	//for(unsigned j=1; j <= maxVar; j++)
	//	if(sharedVariables.isShared(j)) std::cout << "variable " << j << " is shared " << std::endl;

	sendSharedVariables();

	statsSharedVariables();

	rounds = 0;

	before = clock();
	MPI_Barrier(MPI_COMM_WORLD); //all processes are ready to start solving the formula
	mpi_idleTime += clock() - before;

	std::cout << "Solving with " << maxVar << " variables, " << maxClauses << " clauses and " << size << " nodes." << std::endl;

	while (res==M_UNDEF)
	{

		rounds++;
		
		//PROCESS_MEMORY_COUNTERS pmc;
		//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc));
		//saveMemoryUsed((double)pmc.WorkingSetSize);

		if(mpi_verbose > 0)
			std::cout << std::endl << "Node " << rank << " Round " <<  rounds << std::endl;

		if (!solveGlobals())
		{
			res = M_FALSE;
			sendStatus(true);

		} else {
			sendStatus(false);	

			//sends the trail to each worker
			MPI_Status s;
			sendTrail();

			//receives the interpolant or solution from each worker

			interpolants.clear();
			bool all_sat = true;
			models.clear();

			for(unsigned id = 0; id < workersId.size(); id++)
			{
				//buf_interpolant contains a solution or an interpolant
				//if the first position is a 0 then it contains a solution, otherwise it represents the size of the array and contains an interpolant

				int * buf_interpolant;
				before = clock();
				MPI_wait_for(buf_interpolant,workersId[id],s);
				mpi_idleTime += clock() - before;

				assert(s.MPI_TAG == MPI_DATA_INTERPOLANT);
				if(buf_interpolant[0] == 0) 
				{

					interpolants.push_back(m.mkNil());

					std::vector<signed> worker_model;
					//worker_model.push_back(0); //variables number only start at position 1

					if(mpi_verbose > 0) std::cout << "Node " << rank << " received from Node " << workersId[id] << "(Model) : ";

					unsigned size = buf_interpolant[1];

					for(unsigned i = 2; i < size; i++)
					{

						if(mpi_verbose > 0)
							std::cout << buf_interpolant[i] << " ";

						worker_model.push_back(buf_interpolant[i]);

					}

					if(mpi_verbose > 0) std::cout << std::endl;
					models.push_back(worker_model);
					solutions_imported++;

				} else {

					interpolants.push_back(receiveInterpolant(buf_interpolant));
					if(mpi_verbose > 0) 
					{
						std::string str = m.toString(interpolants.back());
						std::cout << "Node " << rank << " received from Node " << workersId[id] << "(ITP) : " << str.c_str() << std::endl;
					}
					all_sat = false;	
				}				
			}



			if(all_sat)
			{
				all_sat_iterations++;
				assert(models.size() == workersId.size());
				if (!findDisagreement())
				{
					res = M_TRUE;
					sendStatus(true); 
				}
				interpolants.clear();
			}
			else if (!importInterpolants()){
				res = M_TRUE;     
				sendStatus(true);
			}

		}

	}

	mpi_totalTime += clock() - start;
	if(mpi_stats > 0) printStats();
	else showTime(mpi_totalTime,"total",mpi_totalTime);

	//synchronization at the end of the solving process
	MPI_Barrier(MPI_COMM_WORLD);

	return res==M_TRUE;
}

bool MPDeSATMaster::solveGlobals(void)
{
	verbose = 0;
	clock_t before = clock();
	if(verbose > 1)
	{
		std::cout << "Finding global assignment..." << std::endl;
		std::cout << "Trail size: " << trail.size() << std::endl;

		for (unsigned i=0; i<trail.size(); i++)
		{
			std::cout << trail[i];
			if (reasons[i]) std::cout << "!";
			std::cout << " ";
		}
		std::cout << std::endl;
	}

	bool r=globalSolver->solve(trail);  

	while (!r && trail.size()>0)
	{ 
		std::vector<signed> temp;
		((MiniSAT_1p*)globalSolver)->getConflict(temp);

		signed last = 0;
		bool reason = true;
		bool inConflict= false;

		if(verbose > 1)
		{
			std::cout << "GLOBAL CONFLICT:";
			if (temp.size()==0)
				std::cout << " EMPTY";
			else
				for (unsigned i=0; i<temp.size(); i++)
					std::cout << temp[i] << " ";
			std::cout << std::endl;

			std::cout << "TRAIL:";
			for (unsigned i=0; i<trail.size(); i++)
			{
				std::cout << trail[i];
				if (reasons[i]) std::cout << "!";
				std::cout << " ";
			}
			std::cout << std::endl;
		}

		inConflict = std::find(temp.begin(), temp.end(), -last)!=temp.end();    
		while (trail.size()>0 && !inConflict)
		{
			last = trail.back();
			trail.pop_back();
			reason = reasons.back();
			reasons.pop_back();      
			inConflict = std::find(temp.begin(), temp.end(), -last)!=temp.end();
		}

		if (inConflict && !reason)
		{      
			trail.push_back(-last);
			reasons.push_back(true);
		}

		if (trail.size()==0 && reason)
			r = false;
		else
		{      
			if(verbose > 1)
			{
				std::cout << "TRAIL:";
				for (unsigned i=0; i<trail.size(); i++)
				{
					std::cout << trail[i];
					if (reasons[i]) std::cout << "!";
					std::cout << " ";
				}
				std::cout << std::endl;
			}

			r=globalSolver->solve(trail);
		}
	}
	mpi_solvingTime += clock() - before;

	if(verbose > 1)
	{
		if (!r)
		{  
			std::cout << "No global assignment." << std::endl;
		}
		else 
		{

			std::cout << "Node " << rank << " Round " << rounds << " Global Assignment:";
			for (int v=1; v<=(signed)maxVar; v++)
				if (sharedVariables.isShared(v))
				{        
					ModelValue mv = globalSolver->get(v);
					if (mv != M_UNDEF){
						if(mv==M_TRUE) std::cout << v;
						else std::cout << -v;
						std::cout << " ";
					} 
				} 
				std::cout << std::endl;
		}
	}

	return r;
}

bool MPDeSATMaster::importInterpolants(void)
{
	clock_t before = clock();
	bool did_something = false;

	for (unsigned i=0; i< workersId.size(); i++)
	{
		CExpression itp = interpolants[i];
		if (m.isNil(itp)){
			continue;
		}

		if (!m.isTrue(itp))
		{
			if(mpi_verbose > 0)
			{
				std::string t = m.toString(itp);
				std::cout << "Global interpolant import: " <<  t.c_str() << std::endl;
			}

			//dot file is useful for debugging purposes

			//std::string filename = "filename";
			//std::stringstream convert; // stringstream used for the conversion
			//convert << rounds;
			//convert << i;
			//filename += convert.str();
			//m.toDot(itp,filename);

			((MiniSAT_1p*)globalSolver)->clearNewClauses();
			//m.sizeRBC(itp,rank);
			globalSolver->addConstraint(itp);      

			unsigned before_clauses = ((MiniSAT_1p*)globalSolver)->numClauses();
			unsigned before_variables = ((MiniSAT_1p*)globalSolver)->nVars();
			((MiniSAT_1p*)globalSolver)->addNewClauses();

			unsigned extra_variables = ((MiniSAT_1p*)globalSolver)->nVars() - before_variables;
			unsigned extra_clauses = ((MiniSAT_1p*)globalSolver)->numClauses() - before_clauses;

			if(extra_variables > variables_interpolants) variables_interpolants = extra_variables;
			if(extra_clauses > clauses_interpolants) clauses_interpolants = extra_clauses;

			interpolants_imported++;
			did_something=true;
		}

		interpolants[i] = m.mkNil();
	}

	return did_something;
}