// Copyright (C) 2011 Microsoft Research
// RM Martins, 2012 
/*
 * This file contains experiments with a 4 stage version of mpDeSAT:
 * (i)   sequential solving; 
 * (ii)  problem splitting;
 * (iii) exchanging the master to an empty node when a limit is reached;
 * (iv)  exchanging the master between nodes when all nodes have formulas.
 * WARNING: Deprecated file!
 */


#include <mpi.h>

#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <decomposition_batch.h>

#include "mpi_commands.h"
#include "mpdesat_dynworker.h"
#include "returnvalue.h"

//#include <windows.h>
//#include <psapi.h>
//#pragma comment(lib, "psapi.lib") //added


MPDeSATDynWorker::MPDeSATDynWorker(ExpressionManager &m, int size, int rank, const char *filename, int c, InterpolationMode i, int master, std::vector<int> worker) :
    MPDeSATNode(m, rank, master, worker),
	size(size),
	cores(c),
	filename(filename),
	workerSolver(NULL),
    sequentialSolver(NULL),
    globalSolver(NULL),
	quit(false),
	solving_stage(false),
	dynamic_timeout(false),
	mpi_startTime(0),
	overflow_limit(100000),
	control_clauses(0),
	control_factor(1.3),
	changeToMaster(false),
	sequential(false),
	phase_transition(true),
	dynamic_rounds(0),
	dynamic_limit(100),
	dynamic_increase(10),
	dynamic_total(0)
{
	workersId.clear();

	if (rank == 0) 
    {
		for (int j=1; j < cores; j++)
			idleWorkers.push_back(j);

        sequentialSolver = new MiniSAT_1p(m, false, MiniSAT_1p::LIFTING);    
	}

	d = new BatchDecomposition();

	if (i != NONE)
    {
		workerSolver = new Partition(m, masterSharedVariables, rank, verbosity);
		workerSolver->setInterpolationMode(i, sharedVariables);
	} 
    else 
        workerSolver = new Partition(m, masterSharedVariables, rank , verbosity, false);

	workerOccurrences = new VariableOccurrence();
	myVariables = new VariableOccurrence();
}

MPDeSATDynWorker::~MPDeSATDynWorker(void)
{
	if (workerSolver) delete workerSolver;
    if (sequentialSolver) delete sequentialSolver;
    if (globalSolver) delete globalSolver;
}

//receive methods
void MPDeSATDynWorker::receiveSharedVariables(void)
{

	MPI_Status s;
	int * shared_variables = NULL;

	clock_t	before = clock();
	MPI_wait_for(shared_variables, masterId, s);
	//std::cout << "Node " << rank << " received msg from " << s.MPI_SOURCE << std::endl;
	mpi_idleTime += clock()-before;

	assert(s.MPI_TAG == MPI_DATA_SHARED_VARIABLES);

	int size = shared_variables[0];
	//std::cout << "Node " << rank << " size: " << size << std::endl;
	for(int i=1; i < size; i++)
	{
		//std::cout << "Node " << rank << " setting variable to shared " << shared_variables[i] << std::endl;
		masterSharedVariables.setShared(shared_variables[i], true);
		//sharedVariables.setShared(shared_variables[i], true);
	}
}

void MPDeSATDynWorker::receiveClauses(void)
{
	MPI_Status s;
	int * clauses = NULL;

	clock_t before = clock();
	MPI_wait_for(clauses, masterId, s);
	mpi_idleTime += clock()-before;

	//clauses contain an array with clauses: 1st element is the size of the array, 2nd element is the number of clauses in the array
	//3rd element is the size of the clause, then the literals that compose that clause (repeat until end of array)
	//std::cout << "Node " << rank << " got clauses from " << s.MPI_SOURCE << std::endl;
	assert(s.MPI_TAG == MPI_DATA_CLAUSE);

	int size = clauses[0]; //first position contains the size of the array
	int num_clauses = clauses[1]; //second position contains the number of clauses in the array

	//std::cout << "Node " << rank << "size: " << size << std::endl;
	//std::cout << "Node " << rank << "num_clauses: " << num_clauses << std::endl;


	for (int i=2; i < size; i++)
	{
		int clauseSize = clauses[i];
		//std::cout << "clauseSize: " <<  clauses[i] << std::endl;
		std::vector<int> literals(clauseSize);

		//std::cout << "clause: ";
		for(int j=0; j < clauseSize; j++)
		{

			literals[j] = clauses[++i];
			//greenOcurrences->setOccurs(literals[j]);
			//std::cout << literals[j] << " ";
		}
		//std::cout << std::endl;

		//std::cout << "Node " << rank << " adding a clause " << std::endl;
		//printClause(literals);
		workerSolver->addClause(literals);	
	}
	//std::cout << "Node " << rank << " end of receive clauses" << std::endl;

    std::cout << "[" << rank << "] received " << num_clauses << "clauses." << std::endl;
}

void MPDeSATDynWorker::receiveProblem(void)
{
	MPI_Status s;

	int * problem_statment = NULL;
	//root worker sends the information about the formula to all workers
	clock_t before = clock();
	MPI_wait_for(problem_statment, masterId, s);
	mpi_idleTime += clock()-before;
	//std::cout << "Node " << rank << " got problem statment from " << s.MPI_SOURCE << std::endl;

	assert(s.MPI_TAG == MPI_DATA_PROBLEM);

	maxVar = problem_statment[0];
	maxClauses = problem_statment[1];

    std::cout << "[" << rank << "] received problem statement." << std::endl;
}


void MPDeSATDynWorker::receiveTrail(void)
{
	MPI_Status s;
	int * assump = NULL;

	clock_t before = clock();
	MPI_wait_for(assump, masterId,s);
	mpi_idleTime = clock()-before;

	assert(s.MPI_TAG == MPI_DATA_TRAIL);

	for(unsigned t=0; t < trail_shared.size(); t++)
		sharedVariables.setShared(trail_shared[t],false);

	//clearing the current trail and reasons
	trail.clear();
	reasons.clear();
	assumptions.clear();

	trail_master.clear();
	trail_shared.clear(); //keeps track of the shared variables that are enforce by the master

	if (mpi_verbose > 0) 
        std::cout << "Node " << rank  << " received trail: ";

	for (int i = 1; i < assump[0]; i++)
	{
		if(mpi_verbose > 0) std::cout << assump[i] << " ";
		//if(occurrences->occurs(assump[i])) //only considers as assumptions the variables that occur in this partition
		assumptions.push_back(assump[i]);
		trail.push_back(assump[i]);
		reasons.push_back(true);
		trail_master.push_back(assump[i]);

		if (!sharedVariables.isShared(assump[i]))
        {
			sharedVariables.setShared(assump[i],true);
			trail_shared.push_back(assump[i]);
		}
	}

	if(mpi_verbose > 0) std::cout << std::endl;
}

bool MPDeSATDynWorker::receiveTermination(void)
{
	
//resume:

	MPI_Status s;
	bool termination = false;
	
	int * cmd = NULL;

	clock_t before = clock();
	//std::cout << "Node " << rank << " is waiting for a termination signal " << std::endl;
	MPI_wait_for(cmd, masterId, s);
	mpi_idleTime += clock()-before;

	assert(s.MPI_TAG == MPI_DATA_QUIT || 
           s.MPI_TAG == MPI_DATA_MASTER || 
           s.MPI_TAG == MPI_DATA_DYNAMIC_SIZE || 
           s.MPI_TAG == MPI_DATA_OCCURRENCES);

	//if(s.MPI_TAG == MPI_DATA_QUIT) std::cout << "Node " << rank << " is receiving a quit msg from " << s.MPI_SOURCE << std::endl;
	//if(s.MPI_TAG == MPI_DATA_MASTER) std::cout << "Node " << rank << " is receiving a master msg from "  << s.MPI_SOURCE  << std::endl;
	//if(s.MPI_TAG == MPI_DATA_DYNAMIC_SIZE) std::cout << "Node " << rank << " is receiving a dynamic size msg from "  << s.MPI_SOURCE << std::endl;
	//if(s.MPI_TAG == MPI_DATA_OCCURRENCES) std::cout << "Node " << rank << " is receiving a data occurrences msg from "  << s.MPI_SOURCE  << std::endl;

	if(s.MPI_TAG == MPI_DATA_QUIT)
	{
	    cmd[0] == 0 ? termination = false : termination = true;	
	} 
    else if (s.MPI_TAG == MPI_DATA_MASTER) 
    {	
	    masterId = cmd[0];
	    //std::cout << "Node " << rank << " changing the masterId to " << masterId << std::endl;	
	    //goto resume;

	    //std::cout << "Node " << rank << " phase: " << phase_transition << std::endl;
	
		sendBuffer.clear();
		sendBuffer.push_back(0);

		//std::cout << "Node " << rank << "variableOcurrs size: " << variableOccurs.size() << std::endl;

		//std::set<int>::iterator it;
		//for(it = variableOccurs.begin(); it != variableOccurs.end(); it++)
		//{
		//	int v = *it;
		//	sendBuffer.push_back(v);
		//}

		for(unsigned t=1 ; t <= maxVar; t++)
		{
			if(workerSolver->occurs(t))
				sendBuffer.push_back(t);
		}

		sendBuffer[0] = sendBuffer.size();
		//std::cout << "Termination Node " << rank << " sending shared variables to " << masterId << std::endl;
		//std::cout << "Node " << rank << " sending shared variables " << sendBuffer[0] << std::endl;
		MPI_safe_send(sendBuffer,masterId,MPI_DATA_SHARED_VARIABLES);
		sendBuffer.clear();	

	    receiveTermination();	
	} 
    else if (s.MPI_TAG == MPI_DATA_DYNAMIC_SIZE) 
    {
		//sendBuffer.clear();
		//sendBuffer.push_back(((MiniSAT_1p*)workerSolver->getSolver())->numClauses());

		//MPI_safe_send(sendBuffer, masterId, MPI_DATA_DYNAMIC_SIZE);

		//int * dynamic_size;

		//MPI_wait_for(dynamic_size, masterId, s);
		//assert(s.MPI_TAG == MPI_DATA_MASTER);
		
		//std::cout << "Node " << rank << " receiving the new master to be " << cmd[0] << std::endl;
		masterId = cmd[0];
			
		if (masterId == rank)
        {		
			//std::cout << "Node " << rank << " is going to be the new master " << std::endl;
			//bool changeToMaster = true;

			trail.clear();
			reasons.clear();
					
			unsigned size = cmd[1];
			if(size > 3){
				bool goto_trail = true;

				//for(unsigned t=1; t < size; t++)
				//	std::cout << dynamic_size[t] << " ";
				//std::cout << std::endl;

				for(unsigned t=1; t < size; t++)
				{
					signed v = cmd[t];
					if(v == 0)
					{
						goto_trail = false;
						continue;
					}

					if(goto_trail){
						trail.push_back(v);
						reasons.push_back(false);
					}
					else reasons[v-1] = true;

				}
			}
			
			dynamic_rounds = 0;
			dynamic_total++;
			return termination;

		} 
        else 
        {
			//std::cout << "Node " << rank << " going to waiting for termination status " << std::endl;
			dynamic_rounds = 0;
			dynamic_total++;
			receiveTermination();
		}	
	} 
    else 
    {	
		//int * dynamic_occurrences;
		//MPI_wait_for(dynamic_occurrences, masterId, s);
		//std::cout << "Node " << rank << " receiving a request regarding the occurrences" << std::endl;

		if (cmd[0] == 0)
        {
			//requesting this worker to send its own occurrences

			//std::cout << "Node " << rank << " must send its own occurrences" << std::endl;
			
			variableOccurs.clear();
			unitClauses.clear();

			std::set<int> occurs;

			//workerSolver->clearVariableOccurrence();
			
			//std::cout << "Node " << rank << " variables " << ((MiniSAT_1p*)workerSolver->getSolver())->numVars() << " clauses " << ((MiniSAT_1p*)workerSolver->getSolver())->numClauses();

			((MiniSAT_1p*)workerSolver->getSolver())->getUsedVariables(variableOccurs,maxVar);
			((MiniSAT_1p*)workerSolver->getSolver())->getUnitClauses(unitClauses,maxVar);

			for(std::set<int>::iterator it = variableOccurs.begin(); it != variableOccurs.end(); it++)
			{
				int v = *it;
				occurs.insert(v);
			}

			for(std::set<int>::iterator it = unitClauses.begin(); it != unitClauses.end(); it++)
			{
				int v = *it;
				unsigned x = v < 0 ? x = -v : x = v;
				if(variableOccurs.find(x)==variableOccurs.end())
				{
					occurs.insert(v);
				} 
			}
						
			sendBuffer.clear();
			sendBuffer.push_back(0);

			for(unsigned v=1; v <= maxVar; v++)
			{
				if(std::find(occurs.begin(),occurs.end(),v)!=occurs.end())
					sendBuffer.push_back(v);
			}
			sendBuffer[0] = sendBuffer.size();			

			//std::cout << "Node " << rank << " sending the occurrences to " << masterId << std::endl;
			MPI_safe_send(sendBuffer, masterId, MPI_DATA_OCCURRENCES);

			receiveTermination();
		
		} 
        else 
        {
		
			//receiving the variable occurrences of all nodes
			//std::cout << "Node " << rank << " receiving the variable occurrences of all nodes" << std::endl;

			unsigned size = cmd[0];
			unsigned current_worker = 0;

			for (unsigned i=0; i< (unsigned)cores; i++)
				changeOccurrences.push_back(new VariableOccurrence());

			//for(unsigned t = 0 ; t < size; t++)
			//	std::cout << cmd[t] << " ";
			//std::cout << std::endl;

			for(unsigned t = 2; t < size; t++)
            {				
				if (cmd[t] == 0)
				{
					current_worker++;
					continue;
				}

				changeOccurrences[current_worker]->setOccurs(cmd[t]);			
			}

			//for(unsigned t=0 ; t < (unsigned)cores; t++)
			//{
			//	for(unsigned v=1; v <= maxVar; v++)
			//	{
			//		if(changeOccurrences[t]->occurs(v)){
			//			std::cout << "Node " << rank << " Partition " << t << " contains variable " << v << std::endl;
			//		}
			//	}
			//}

			phase_transition = false;

			receiveTermination();
		}
	}
	
	return termination;
}

//send methods
void MPDeSATDynWorker::sendTrail(void)
{
	sendBuffer.clear();
	sendBuffer.push_back(0);

	//sends the trail of the master downwards the tree
	for(unsigned t=0; t < trail_master.size(); t++)
		sendBuffer.push_back(trail_master[t]);

	for (int v=1; v<=(signed)maxVar; v++)
	{
		//TO DO: check for both phases in the same iteration
		if (sharedVariables.isShared(v)
			&& std::find(trail_master.begin(),trail_master.end(),v)==trail_master.end()
			&& std::find(trail_master.begin(),trail_master.end(),-v)==trail_master.end() 
			)
		{        
			//ModelValue mv = workerSolver->get(v);
			ModelValue mv = globalSolver->get(v);
			if (mv != M_UNDEF){
				if(mv==M_TRUE) sendBuffer.push_back(v);
				else sendBuffer.push_back(-v);
			} 
		}
	}


	sendBuffer[0] = sendBuffer.size();

	//std::cout << "Node " << rank << " sent trail: ";
	//for(unsigned t=1; t < sendBuffer.size(); t++)
	//	std::cout << sendBuffer[t] << " ";
	//std::cout << std::endl;
	
	for(unsigned id=0; id < workersId.size(); id++)
		MPI_safe_send(sendBuffer, workersId[id], MPI_DATA_TRAIL);

	sendBuffer.clear();
}

void MPDeSATDynWorker::sendTermination(bool status)
{
	sendBuffer.clear();

	status ? sendBuffer.push_back(1) : sendBuffer.push_back(0);

	for(unsigned id = 0; id < workersId.size(); id++)
		MPI_safe_send(sendBuffer,workersId[id],MPI_DATA_QUIT);

	sendBuffer.clear();
}

void MPDeSATDynWorker::sendSolution(void)
{
	//solution has been found
	sendBuffer.clear();

	sendBuffer.push_back(0);
	sendBuffer.push_back(0);

	if(mpi_verbose > 0) std::cout << "Node " << rank << " sending model to Node " << masterId << "(Model) : ";

	for(signed var = 1; var <= (signed)maxVar; var++)
	{
		int mvi = workerSolver->get(var);

		switch(mvi)
		{
		case M_TRUE:
			if(mpi_verbose > 0) std::cout << var << " ";
			if(myVariables->occurs(var) && masterSharedVariables.isShared(var)) sendBuffer.push_back(var);
			break;
		case M_FALSE:
			if(mpi_verbose > 0) std::cout << -var << " ";
			if(myVariables->occurs(var) && masterSharedVariables.isShared(var)) sendBuffer.push_back(-var);
			break;
		case M_UNDEF:
			break;
		default: 
			throw std::exception("Variable value not defined");
		}
	}

	if(mpi_verbose > 0) std::cout << std::endl << std::endl;

	sendBuffer[1] = sendBuffer.size();

	solutions_exported++;
	MPI_safe_send(sendBuffer,masterId,MPI_DATA_INTERPOLANT);

	sendBuffer.clear();
}


//solving methods
bool MPDeSATDynWorker::solveGlobals(void)
{
	//PROCESS_MEMORY_COUNTERS pmc;
	//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
	//std::cout << "Node " << rank << " [memory - before globals] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;

	unsigned num_calls = 0;

	if (verbose > 1)
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
	num_calls++;

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

			//((MiniSAT_1p*)globalSolver)->cleanLearnts();
			r=globalSolver->solve(trail);
			num_calls++;
		}
	}

	//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
	//std::cout << "Node " << rank << " [memory - after globals] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;

	//std::cout << "Node " << rank << " [memory - calls] : " << num_calls << std::endl;

	//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
	//std::cout << "Node " << rank << " [memory - after clear globals] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;

	dynamic_rounds++;

	return r;

}

void MPDeSATDynWorker::solveSolution(void)
{
	//solution has been found
	sendBuffer.clear();

	sendBuffer.push_back(0);
	sendBuffer.push_back(0);

	for (unsigned i=0; i< workersId.size(); i++)
	{
		for(unsigned t=0; t < models[i].size(); t++)
		{      

			int mvi = models[i][t];
			if(std::find(sendBuffer.begin(), sendBuffer.end(), mvi)==sendBuffer.end()) 
				sendBuffer.push_back(mvi);				
		}
	}

	sendBuffer[1] = sendBuffer.size();

	if(mpi_verbose > 0)
	{
		std::cout << "Node " << rank << " sending model to Node " << masterId << "(Model) : ";
		for(signed var=2; var < (signed)sendBuffer.size(); var++)
		{
			std::cout << sendBuffer[var] << " ";

		}
		std::cout << std::endl;
	}

	solutions_exported++;
	MPI_safe_send(sendBuffer,masterId,MPI_DATA_INTERPOLANT);

	sendBuffer.clear();
}

void MPDeSATDynWorker::solvePartitions()
{
	//PROCESS_MEMORY_COUNTERS pmc;
	//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
	//std::cout << "Node " << rank << " [memory - after partitions] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;

	workerInterpolant = m.mkNil();

	clock_t before = clock();
	Expression t = workerSolver->getInterpolant(assumptions);
	mpi_solvingTime += clock()-before;

	sendBuffer.clear();

	if (!workerSolver->em().isTrue(t))
	{
		workerInterpolant = m.duplicate(t, workerSolver->em());

		//std::string filename = "filename";
		//std::stringstream convert; // stringstream used for the conversion
		//convert << rank;
		//convert << interpolants_exported;
		//filename += convert.str();
		//m.toDot(workerInterpolant,filename);

		sendBuffer.push_back(0); //reserve space for the size of the buffer
		//m.sizeRBC(workerInterpolant, rank);

		bool overflow_problem = m.sizeRBC(workerInterpolant, overflow_limit);

		if (overflow_problem)
        {
			workerInterpolant = m.mkNil();

			vec<Lit> &cfs = ((MiniSAT_1p*)workerSolver)->conflict;
			std::vector<signed> cf_tmp;
			for (int i=0; i<cfs.size(); i++)        
				cf_tmp.push_back(litToLiteral(cfs[i]));


			if (cf_tmp.size()==0)
				workerInterpolant = m.mkFalse();
			else if (cf_tmp.size()==1)
				workerInterpolant = m.mkLiteral(cf_tmp[0]);
			else
			{
				workerInterpolant = m.mkOr(m.mkLiteral(cf_tmp[0]), m.mkLiteral(cf_tmp[1]));
				for (unsigned i=2; i<cf_tmp.size(); i++)
				{
					workerInterpolant = m.mkOr(workerInterpolant, m.mkLiteral(cf_tmp[i]));
				}

			}
		}

		m.toMPI(sendBuffer,workerInterpolant);
		sendBuffer[0] = sendBuffer.size();

		if(mpi_verbose > 0)
		{
			std::string str = m.toString(workerInterpolant);
			std::cout << "Node " << rank << " interpolant export (ITP): " << str.c_str() << std::endl;
		}

		interpolants_exported++;
		interpolantStats(workerInterpolant);
		sizeInterpolant += sendBuffer.size();
		if(sendBuffer.size() > maxSizeInterpolant){
			maxSizeInterpolant = sendBuffer.size();
			//std::cout << "Node " << rank << " generates an interpolant of size: " << maxSizeInterpolant << std::endl;
		}
		MPI_safe_send(sendBuffer,masterId,MPI_DATA_INTERPOLANT);

	} else {

		sendSolution();

	}

	dynamic_rounds++;

	sendBuffer.clear();
	workerInterpolant = m.mkNil();
	t = m.mkNil();
	//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
	//std::cout << "Node " << rank << " [memory - after partitions] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;
}

bool MPDeSATDynWorker::importInterpolants(void)
{
	clock_t before = clock();
	bool did_something = false;

	for (unsigned i=0; i< workersId.size(); i++)
	{
		CExpression itp = interpolants[i];
		if (m.isNil(itp))
        {
			continue;
		}

		if (!m.isTrue(itp))
		{
			if (verbosity>0)
			{
				std::string t = m.toString(itp);
				std::cout << "Node " << rank << " Global interpolant import: " <<  t.c_str() << std::endl;
			}

			//dot file is useful for debugging purposes

			//std::string filename = "filename";
			//std::stringstream convert; // stringstream used for the conversion
			//convert << rounds;
			//convert << i;
			//filename += convert.str();
			//m.toDot(itp,filename);

			//PROCESS_MEMORY_COUNTERS pmc;
			//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
		

			//std::cout << "Node " << rank << " [memory - before ITP->clauses] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;
			
			//((MiniSAT_1p*)workerSolver->getSolver())->clearNewClauses();
			//(workerSolver->getSolver())->addConstraint(itp);
			//((MiniSAT_1p*)workerSolver->getSolver())->addNewClauses();

			((MiniSAT_1p*)globalSolver)->clearNewClauses();
			globalSolver->addConstraint(itp);
			((MiniSAT_1p*)globalSolver)->addNewClauses();
			
			//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
			//std::cout << "Node " << rank << " [memory - ITP->clauses] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;
			//std::cout << "Node " << rank << " [memory - variables, clauses, learnts] : ( " << ((MiniSAT_1p*)workerSolver->getSolver())->numClauses() << 
			//	" , " << ((MiniSAT_1p*)workerSolver->getSolver())->numVars() << " , " << ((MiniSAT_1p*)workerSolver->getSolver())->nLearnts() << " )" << std::endl;  

			did_something=true;

		}

		interpolants[i] = m.mkNil();
	}

	mpi_solvingTime = clock() - before;
	interpolants.clear();
	return did_something;
}

void MPDeSATDynWorker::sendSharedVariablesDyn(void)
{
	sendBuffer.clear();

	sendBuffer.push_back(0);


	for(unsigned j = 1; j <= maxVar; j++)
	{
		if(sharedVariables.isShared(j))
			sendBuffer.push_back(j);
	}

	sendBuffer[0] = sendBuffer.size();
	
	for(unsigned id = 0; id < workersId.size(); id++)
		MPI_safe_send(sendBuffer,workersId[id],MPI_DATA_SHARED_VARIABLES);

	sendBuffer.clear();
}

bool MPDeSATDynWorker::sequential_solving(void)
{
	//STAGE 1:: sequential solving
	//TO DO: parallel reading for large files

	bool res = false;

	if (rank != MPI_MASTER)
    {

		while(waitsForWork());

	} 
    else 
    {
		sequential = true;
		masterId = 0;
		readDimacsFile(filename);
		
		clock_t before = clock();
		bool disable_sequential = false;
		if (limit_conflicts <= 0 || limit_clauses <= 0 || limit_time <= 0)
			disable_sequential = true;
		
		if (!disable_sequential)
        {
			((MiniSAT_1p*)sequentialSolver)->setLimitResources(limit_conflicts,limit_clauses,limit_time); //run minisat with limited resources
				res = sequentialSolver->solve();
		}

		double solving_time = (clock() - before)/(double)CLOCKS_PER_SEC;
		sequential = false;

		mpi_nodeTime = clock();

		bool resources = ((MiniSAT_1p*)sequentialSolver)->getResources();

		if (disable_sequential || resources)
		{
			////OLD VERSION
			////sequential solving did not solve the formula
			////split the formula between the different problem nodes
			//unsigned next_worker = idleWorkers.back();
			//idleWorkers.pop_back();

			//workingWorkers.push_back(next_worker);
			//workersId.push_back(next_worker);

			//sendBuffer.clear();
			////send information regarding the node that will send half of clauses
			//sendBuffer.push_back(2); //problem splitting stage
			//sendBuffer.push_back(rank);
			//MPI_safe_send(sendBuffer,next_worker,MPI_DATA_MASTER);

			//splitHalf(next_worker);
			std::cout << "Sequential Solving [stage 1] : " << (clock() - mpi_startTime)/(double)CLOCKS_PER_SEC << std::endl;
			splitN(cores-1); //TO DO: change the method so that we can choose how many formula workers we want
			std::cout << "Formula Decomposition [stage 2] : " << (clock() - mpi_startTime)/(double)CLOCKS_PER_SEC << std::endl;
			std::cout << "Formula Decomposition [stage 2] (active workers) : " << cores - idleWorkers.size() << std::endl;
			std::cout << "Formula Decomposition [stage 2] (idle workers) : " << idleWorkers.size() << std::endl;
		} 
        else 
        {
			quit = true;
			//formula is "easy" and was solved by sequential solving
			sendBuffer.clear();
			sendBuffer.push_back(1); //sequential stage
			//sending termination signal to all the idle workers
			for(unsigned t=0; t < idleWorkers.size(); t++)
				MPI_safe_send(sendBuffer,idleWorkers[t],MPI_DATA_MASTER);
		}		
	}
	
	return res;
}

void MPDeSATDynWorker::workerToMaster(void)
{
	if (globalSolver) delete globalSolver;
	globalSolver = new MiniSAT_1p(m, false, MiniSAT_1p::LIFTING);
	while (globalSolver->numVars() < workerSolver->numVars()+1) globalSolver->addVar();

	workersId.clear();
	for(unsigned t=0; t < (unsigned)cores; t++)
	{
		if(rank != t) workersId.push_back(t);
	}

	unsigned global_clauses = ((MiniSAT_1p*)workerSolver->getSolver())->numClauses();

	std::vector<signed> clause;
	for(unsigned t=0; t < global_clauses; t++)
	{
		((MiniSAT_1p*)workerSolver->getSolver())->removeLastToClause(clause);
		globalSolver->addClause(clause);
	}

	((MiniSAT_1p*)workerSolver->getSolver())->cleanLearnts();
	//workerSolver->newSolver();
	//while(workerSolver->numVars() < maxVar+1) workerSolver->addVar();

	trail.clear();
	trail_master.clear();
	assumptions.clear();
	reasons.clear();
	trail_shared.clear();
}

void MPDeSATDynWorker::masterToWorker(void)
{
	//masterSharedVariables.cleanShared();

	for(unsigned v = 1 ; v <= maxVar ; v++)
	{
		if(sharedVariables.isShared(v)){
			masterSharedVariables.setShared(v, true);
		}
	}

	unsigned global_vars = ((MiniSAT_1p*)globalSolver)->numVars();
	unsigned global_clauses = ((MiniSAT_1p*)globalSolver)->numClauses();
	for(unsigned t=0; t < global_vars; t++) workerSolver->addVar();

	std::vector<signed> clause;
	for(unsigned t=0; t < global_clauses; t++)
	{
		((MiniSAT_1p*)globalSolver)->removeLastToClause(clause);
		workerSolver->addClause(clause);
	}

	((MiniSAT_1p*)globalSolver)->cleanLearnts();
	if (globalSolver) delete globalSolver;
	globalSolver = new MiniSAT_1p(m, false, MiniSAT_1p::LIFTING);
	while(globalSolver->numVars() < maxVar+1) globalSolver->addVar();

	//sharedVariables.cleanShared();	
}

bool MPDeSATDynWorker::solve(void)
{
	clock_t mpi_startTime = clock();

	ModelValue res = M_UNDEF;
	double cutoff_timeout_sequential = ceil(sec_timeout*((double)pct_timeout/100));
	bool r = false;

	r = sequential_solving();

	//STAGE 2:: split problem
	if (quit)
    {
		std::cout << "Node " << rank << " [solved] : variables " <<  ((MiniSAT_1p*)workerSolver->getSolver())->numVars() << " clauses " << ((MiniSAT_1p*)workerSolver->getSolver())->numClauses() << std::endl;
		std::cout << "Node " << rank << " [solved] : idle Time " << mpi_idleTime /(double)CLOCKS_PER_SEC << std::endl;
		std::cout << "Node " << rank << " [solved] : CPU Time " << (clock() - mpi_nodeTime - mpi_idleTime) /(double)CLOCKS_PER_SEC << std::endl;
		if (rank!=0) std::cout << "Node " << rank << " [solved] : efficiency 1 [not used]" << std::endl;
		else std::cout << "Node " << rank << " [solved] : efficiency 1 [sequential]" << std::endl;

		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0)
        {
			std::cout << "Node " << rank << " formula has been solved by the sequential solver." << std::endl;
			std::cout << "Total Time: " <<  (clock() - mpi_startTime)/(double)CLOCKS_PER_SEC << std::endl;
		}

		return r;
	} 
	
	//else {
	//
	//	//get a new master to all the nodes
	//	if(rank == 0)
	//	{
	//	workersId.push_back(rank); //include itself in the worker list of the next master

	//	unsigned next_master = idleWorkers.back();
	//	idleWorkers.pop_back();

	//	//std::cout << "Node " << rank << " updating master idle " << std::endl;
	//	updateMasterIdle(next_master);

	//	//std::cout << "Node " << rank << " changing idle " << std::endl;
	//	changeMaster(next_master, 5);

	//	//std::cout << "Node " << rank << " updating master" << std::endl;
	//	updateMaster(next_master);

	//	workersId.clear();
	//	
	//	} else if (rank != masterId){
	//	
	//		int * new_master;
	//		MPI_Status s;
	//		//std::cout << "Node " << rank << " waiting for the new master " << std::endl;
	//		MPI_wait_for(new_master,masterId,s);
	//		assert(s.MPI_TAG == MPI_DATA_MASTER);
	//		
	//		masterId = new_master[0];
	//		//std::cout << "Node " << rank << " updated the new master to " << masterId << std::endl;
	//		delete new_master;
	//	}
	//
	//
	//}

	//else if(!solving_stage) {
	//	split_problem((int)ceil(sec_timeout*((double)pct_timeout2/100)));
	//}

	
	//if(rank !=0 && masterId!=rank && !quit)
	//{
	//	//sends request to node 0
	//	sendBuffer.clear();
	//	sendBuffer.push_back(0);
	//	std::cout << "Node " << rank << " sending stop request to process 0" << std::endl;
	//	MPI_safe_send(sendBuffer,0,MPI_DATA_MASTER);

	//}

	//if(rank == 0)
	//{

	//	workersId.push_back(rank); //include itself in the worker list of the next master

	//	assert(!idleWorkers.empty());
	//	unsigned next_master = idleWorkers.back();
	//	idleWorkers.pop_back();

	//	updateMasterIdle(next_master);

	//	changeMaster(next_master);
	//	updateMaster(next_master);
	//	workersId.clear();
	//}

	//if(rank != 0 && masterId!=rank && !quit)
	//{
	//	int * new_master;
	//	MPI_Status s;
	//	MPI_wait_for(new_master,MPI_MASTER,s);
	//	masterId = new_master[1];	

	//}

	while(workerSolver->numVars() < maxVar+1) workerSolver->addVar(); //creates dummy variables if necessary
	
	globalSolver = new MiniSAT_1p(m, false, MiniSAT_1p::LIFTING);    
	while(globalSolver->numVars() < maxVar+1) globalSolver->addVar();

	if (rank == masterId && !solving_stage)
	{
		//((MiniSAT_1p*)sequentialSolver)->printDB();

		((MiniSAT_1p*)sequentialSolver)->getUnitClauses(variableOccurs, maxVar); //get literals that are assigned at decision level 0
		//((MiniSAT_1p*)sequentialSolver)->getUsedVariables(variableOccurs, maxVar); //get literals that are assigned at decision level 0
		
		//std::cout << "Node " << rank << " sent " << variableOccurs.size() << " unit clauses." << std::endl;
		//sendBuffer.clear();
		//sendBuffer.push_back(0);

		//std::set<int>::iterator it;
		//for(it = variableOccurs.begin(); it != variableOccurs.end(); it++)
		//{
		//	int v = *it;
		//	unsigned x = v < 0 ? x = -v : x = v;
		//	sendBuffer.push_back(v);
		//	//workerSolver->insertVariableOccurrence(x);
		//}

		//sendBuffer[0] = sendBuffer.size();
		//for(unsigned t=0; t < workersId.size(); t++)
		//{
		//	MPI_safe_send(sendBuffer, workersId[t], MPI_DATA_SHARED_VARIABLES);
		//}

		//sendBuffer.clear();
		
		((MiniSAT_1p*)sequentialSolver)->cleanLearnts();
		if (sequentialSolver) delete sequentialSolver;
        sequentialSolver = 0;

		std::cout << "MASTER: " << rank << " has " << workersId.size() << " workers" << std::endl;
		
		//for(unsigned i=0; i < workersId.size(); i++)
			//std::cout << "Node " << rank << " has worker " << workersId[i] << std::endl;

		for (unsigned i=0; i< workersId.size(); i++)
			masterOccurrences.push_back(new VariableOccurrence());

		assumptions.resize(workersId.size());
		interpolants.resize(workersId.size(), m.mkNil());
		models.resize(workersId.size());
	
		MPI_Status s;

		for(unsigned t=0; t < workersId.size(); t++)
		{
			int * variableOccurrences = NULL;
			std::cout << "Node " << rank << " is waiting for variable occurrences from " << workersId[t] << std::endl;
			clock_t before = clock();
			MPI_wait_for(variableOccurrences, workersId[t], s);
			mpi_idleTime += clock()-before;
            std::cout << "Node " << rank << " received variable occurrences from " << workersId[t] << std::endl;

			assert(s.MPI_TAG == MPI_DATA_SHARED_VARIABLES);
			
			int size = variableOccurrences[0];
			for(int i = 1; i < size; i++)
				masterOccurrences[t]->setOccurs(variableOccurrences[i]);
		}

		for(unsigned i=0; i < workersId.size(); i++)
			sharedVariables.add(masterOccurrences[i]);
		
		sharedVariables.update();
		
		for(std::set<int>::iterator it = variableOccurs.begin(); it != variableOccurs.end(); it++)
		{
			int v = *it;
			unsigned x = v < 0 ? x = -v : x = v;

			//check if the variable occurs in at least one more partition
			bool forceShared = false;
			for(unsigned i=0; i < workersId.size(); i++)
			{
				if(masterOccurrences[i]->occurs(v))
				{
					forceShared = true;
					break;
				}
			}

			if(forceShared && !sharedVariables.isShared(x)){
				sharedVariables.setShared(x, true);
			}
			
			std::vector<signed> tmp(1);
			tmp[0] = v;
			globalSolver->addClause(tmp);
			//workerSolver->insertVariableOccurrence(x);

		}
		
		//for(unsigned v=1; v <= maxVar; v++)
		//	if(!sharedVariables.isShared(v)) variableOccurs.insert(v);
		
		sendSharedVariablesDyn();
		statsSharedVariables();


		//std::cout << "Node " << rank << " clearing the variable occurrence " << std::endl;
		//workerSolver->clearVariableOccurrence();
	} 
    else if (!quit && !solving_stage)
    {
		//get unit clauses from the sequential solver
		//int * unit_clauses;
		//MPI_Status s;

		//MPI_wait_for(unit_clauses, masterId, s);
		//assert(s.MPI_TAG == MPI_DATA_SHARED_VARIABLES);

		//unsigned size = unit_clauses[0];
		//for(unsigned t = 1 ; t < size ; t++){
		//	std::vector<signed> unit(1);
		//	unit[0] = unit_clauses[t];
		//	workerSolver->addClause(unit);
		//	variableOccurs.insert(unit_clauses[t]);
		//}
		
		 //keep learnt clauses from the previous iterations? (may change the set of shared variables)

		//std::cout << "WORKER " << rank << std::endl;
		
		//if(rank == 16) std::cout << "get used variables" << std::endl;
		//if(rank == 16) std::cout << "variables: " << ((MiniSAT_1p*)workerSolver->getSolver())->numVars() << " clauses " << 
		//	((MiniSAT_1p*)workerSolver->getSolver())->numClauses() << std::endl;
		
		((MiniSAT_1p*)workerSolver->getSolver())->getUsedVariables(variableOccurs, maxVar);
		sendBuffer.clear();
		sendBuffer.push_back(0);

		//if(rank == 16) std::cout << "Node " << rank << "variableOcurrs size: " << variableOccurs.size() << std::endl;

		workerSolver->clearVariableOccurrence();
		std::set<int>::iterator it;
		for(it = variableOccurs.begin(); it != variableOccurs.end(); it++)
		{
			int v = *it;
			//if(rank == 16) std::cout << "Node " << rank << " has variable " << v << std::endl;
			unsigned x = v < 0 ? x = -v : x = v;
			sendBuffer.push_back(x);
			workerSolver->insertVariableOccurrence(x);
			
			//if(rank == 16) std::cout << "inserted Node " << rank << " has variable " << v << std::endl;
		}

		sendBuffer[0] = sendBuffer.size();
		//if(rank == 16) std::cout << "Node " << rank << " sending shared variables to " << masterId << std::endl;
		MPI_safe_send(sendBuffer,masterId,MPI_DATA_SHARED_VARIABLES);
		sendBuffer.clear();

		//if(rank == 16) std::cout << "Node " << rank << " receiving shared variables" << std::endl;
		receiveSharedVariables();

		control_clauses = ((MiniSAT_1p*)workerSolver->getSolver())->numClauses();
	}


	((MiniSAT_1p*)workerSolver->getSolver())->cleanLearnts();

	//while(workerSolver->numVars() < maxVar+1) workerSolver->addVar(); //creates dummy variables if necessary
	//if(rank == 16 || rank == 17) std::cout << "Node " << rank << " cleaning " << std::endl;
	//((MiniSAT_1p*)workerSolver->getSolver())->cleanLearnts(); //keep learnt clauses from the previous iterations? (may change the set of shared variables)
	 std::cout << "Node " << rank << " Solving with " << workerSolver->numVars() << " variables and " << workerSolver->numClauses() << " clauses." << std::endl;
	 //std::cout << "Node " << rank << " has master " << masterId << std::endl;

	//((MiniSAT_1p*)workerSolver)->resetProof();
	//if(rank == 0){
	//	((MiniSAT_1p*)workerSolver->getSolver())->printDB();
	//}

	 //for(unsigned t=1; t <= maxVar; t++){
		// if(workerSolver->variablesOccurrence().occurs(t))
		//	 std::cout << "Node " << rank << " has variable " << t << std::endl;
	 //}

	//std::cout << "Node " << rank << " workers: ";
	//for(unsigned t=0; t < workersId.size(); t++)
	//	std::cout << workersId[t] << " ";
	//std::cout << std::endl;

	//((MiniSAT_1p*)workerSolver->getSolver())->setLimitResources(0,0);
	//if(rank != masterId) ((MiniSAT_1p*)workerSolver->getSolver())->setLimitResources(0,0,0);
	//else {
	//	std::cout << "Node " << rank << " setting limit resources for the search " << std::endl;
	//	((MiniSAT_1p*)workerSolver->getSolver())->setLimitResources(limit_conflicts,0,0);
	//}

	((MiniSAT_1p*)workerSolver->getSolver())->setLimitResources(0,0,0); //reset limit resources for all nodes

	if (solving_stage)
    {
		masterId = rank;
		control_clauses = maxClauses;
	
		masterOccurrences.clear();

		for (unsigned i=0; i< workersId.size(); i++)
			masterOccurrences.push_back(new VariableOccurrence());

		assumptions.resize(workersId.size());
		interpolants.resize(workersId.size(), m.mkNil());
		models.resize(workersId.size());

		MPI_Status s;

		//std::cout << "Node " << rank  << " workers: ";
		//for(unsigned t=0; t < workersId.size(); t++)
		//	std::cout << workersId[t] << " " ;
		//std::cout << std::endl;

		for (unsigned t=0; t < workersId.size(); t++)
		{
			int * variableOccurrences = NULL;
			//std::cout << "Node " << rank << " is waiting for variable occurrences from " << workersId[t] << std::endl;
			clock_t before = clock();
			MPI_wait_for(variableOccurrences, workersId[t], s);
			mpi_idleTime += clock() - before;

			assert(s.MPI_TAG == MPI_DATA_SHARED_VARIABLES);

			int size = variableOccurrences[0];
			for(int i = 1; i < size; i++)
				masterOccurrences[t]->setOccurs(variableOccurrences[i]);				
		}

		for(unsigned i=0; i < workersId.size(); i++)
			sharedVariables.add(masterOccurrences[i]);

		sharedVariables.update();
		
		//for(unsigned v=1; v <= maxVar; v++)
		//	if(sharedVariables.isShared(v)) variableOccurs.insert(v);	
	
		workerSolver->clearVariableOccurrence();

		statsSharedVariables();
	}

	//mpi_nodeTime = clock();
	//PROCESS_MEMORY_COUNTERS pmc;

	if (rank == masterId) std::cout << "Node " << rank << " [master - time] : " << (clock() - mpi_startTime)/(double)CLOCKS_PER_SEC << std::endl;
	if (rank == masterId) std::cout << "Node " << rank << " [master - control] : " << control_clauses << std::endl;
	
	//if(rank == masterId) {
	//	std::cout << "Node " << rank << " variables " << ((MiniSAT_1p*)globalSolver)->numVars() << " clauses " << ((MiniSAT_1p*)globalSolver)->numClauses() << std::endl;
	//	((MiniSAT_1p*)globalSolver)->printDB();
	//}


	variableOccurs.clear();
	unitClauses.clear();
	//((MiniSAT_1p*)workerSolver->getSolver())->getUsedVariables(variableOccurs,maxVar);
	//((MiniSAT_1p*)workerSolver->getSolver())->getUnitClauses(unitClauses,maxVar);

	if (rank == 0)
	{
	    ((MiniSAT_1p*)globalSolver)->getUsedVariables(variableOccurs,maxVar);
	    ((MiniSAT_1p*)globalSolver)->getUnitClauses(unitClauses,maxVar);

		for (std::set<int>::iterator it = variableOccurs.begin(); it != variableOccurs.end(); it++)
	    {
		    int v = *it;
		    //unsigned x = v < 0 ? x = -v : x = v;
		    assert( v > 0);
		    myVariables->setOccurs(v);
	    }

	    for (std::set<int>::iterator it = unitClauses.begin(); it != unitClauses.end(); it++)
	    {
		    int v = *it;
		    unsigned x = v < 0 ? x = -v : x = v;
		    if(variableOccurs.find(x)==variableOccurs.end())
		    {
			    myVariables->setOccurs(v);
		    } 

	    }
	} 
    else 
    {	
		for(unsigned v=1 ; v <= maxVar; v++)
		{
			if(workerSolver->occurs(v)){
				myVariables->setOccurs(v);
			}
		}
	}

	unsigned num_myVariables = 0;
	for(unsigned v=1; v <= maxVar; v++)
	{
		if(myVariables->occurs(v))
			num_myVariables++;
	}

	std::cout << "Node " << rank << " has " << num_myVariables << std::endl;


	while (!quit)
	{
		rounds++;

		if(rank == masterId)
		{
			
			//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
			//std::cout << "Node " << rank << " [memory - MB] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;

			//std::cout << "Node " << rank << " rounds " << rounds << std::endl;
			if (!solveGlobals())
			{
				////NYI: old idea
				//if(dynamic_timeout)
				//{
				//	std::cout << "Node " << rank << " dynamic timeout" << std::endl;
				//	if(idleWorkers.size()!=0)
				//	{
				//		//dynamic_split stage
				//		
				//		workersId.push_back(rank); //include itself in the worker list of the next master

				//		assert(!idleWorkers.empty());
				//		unsigned next_master = idleWorkers.back();
				//		idleWorkers.pop_back();

				//		std::cout << "Node " << rank << " updating master idle " << std::endl;
				//		updateMasterIdle(next_master);

				//		std::cout << "Node " << rank << " changing master " << std::endl;
				//		changeMaster(next_master, 3);

				//		std::cout << "Node " << rank << " update master " << std::endl;
				//		updateMaster(next_master);

				//		((MiniSAT_1p*)workerSolver->getSolver())->setLimitResources(0,0,0);

				//	} else {
				//		//dynamic_balancing stage
				//		std::cout << "Node " << rank << " dynamic balancing stage " << std::endl;
				//		((MiniSAT_1p*)workerSolver->getSolver())->setLimitResources(0,0,0);

				//	}

				//	dynamic_timeout = false;

				//} else {

					std::cout << "Node " << rank << " [solution - rounds] : " << dynamic_total << std::endl;
					std::cout << "Node " << rank << " [solution - answer] : UNSAT " << std::endl;
					std::cout << "Node " << rank << " [solution - time] : " << (clock() - mpi_startTime)/(double)CLOCKS_PER_SEC << std::endl;
					

					//((MiniSAT_1p*)workerSolver->getSolver())->printDB();

					res = M_FALSE;
					quit = true;
					sendStatus(true);

				//}

			} else {


				sendStatus(false);	

				//sends the trail to each worker
				MPI_Status s;
				sendTrail();

				//receives the interpolant or solution from each worker
				interpolants.clear();
				bool all_sat = true;
				models.clear();

				//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
				//std::cout << "Node " << rank << " [memory - import - begin] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;

				//std::cout << "Node " << rank << " waiting for interpolant or solution " << std::endl;
				for(unsigned id = 0; id < workersId.size(); id++)
				{
					//buf_interpolant contains a solution or an interpolant
					//if the first position is a 0 then it contains a solution, otherwise it represents the size of the array and contains an interpolant
					int * buf_interpolant;

					//std::cout << "Node " << rank << " waiting for interpolant or solution from " << workersId[id] << std::endl;
					clock_t before = clock();
					MPI_wait_for(buf_interpolant,workersId[id],s);
					mpi_idleTime += clock() - before;

					assert(s.MPI_TAG == MPI_DATA_INTERPOLANT);
					if(buf_interpolant[0] == 0) 
					{

						interpolants.push_back(m.mkNil());

						std::vector<signed> worker_model;
						//worker_model.push_back(0); //variables number only start at position 1
						//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
						//std::cout << "Node " << rank << " [from " << workersId[id] << " ] - [memory - MODEL] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;

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
						//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
						//std::cout << "Node " << rank <<  " [before receive ITP] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;
						interpolants.push_back(receiveInterpolant(buf_interpolant));
						//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
						//std::cout << "Node " << rank <<  " [after receive ITP] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;
						
						if(m.isFalse(interpolants.back()))
							std::cout << "Node " << rank << " [ false ITP ] : from " << workersId[id] << std::endl;
						
						if(mpi_verbose > 0) 
						{
							std::string str = m.toString(interpolants.back());
							std::cout << "Node " << rank << " received from Node " << workersId[id] << "(ITP) : " << str.c_str() << std::endl;
						}
						all_sat = false;
						//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
						//std::cout << "Node " << rank <<  " [from " << workersId[id] << " ] - [memory - ITP] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;
					}					
				}

				if(all_sat)
				{
					all_sat_iterations++;
					assert(models.size() == workersId.size());
					//if (phase_transition && !findDisagreement())
					//{
					//	std::cout << "Node " << rank << " [solution - rounds] : " << dynamic_total << std::endl;
					//	std::cout << "Node " << rank << " [solution - answer] : SAT " << std::endl;
					//	std::cout << "Node " << rank << " [solution - time] : " << (clock() - mpi_startTime)/(double)CLOCKS_PER_SEC << std::endl;

					//	res = M_TRUE;
					//	quit = true;
					//	sendStatus(true); 
					//}
					
					//if (!phase_transition && !findDisagreementDyn())
					if (!findDisagreementDyn())
					{
						std::cout << "Node " << rank << " [solution - rounds] : " << dynamic_total << std::endl;
						std::cout << "Node " << rank << " [solution - answer] : SAT " << std::endl;
						std::cout << "Node " << rank << " [solution - time] : " << (clock() - mpi_startTime)/(double)CLOCKS_PER_SEC << std::endl;

						res = M_TRUE;
						quit = true;
						sendStatus(true); 
					}

					interpolants.clear();
				}
				else if (!importInterpolants()){
					res = M_TRUE;
					quit = true;
					sendStatus(true);
				}


				//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
				//std::cout << "Node " << rank << " [memory - import - end] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;

				//Current worker checks if the formula has become too large
				//if(control_factor !=0 && ((MiniSAT_1p*)globalSolver)->numClauses() > control_factor*control_clauses)
				//std::cout << "Node " << rank << " dynamic rounds " << dynamic_rounds << std::endl;
				if(dynamic_rounds == dynamic_limit)
				{
					dynamic_limit += dynamic_increase;
					dynamic_total++;

					//checks if there are idle nodes
					if(!idleWorkers.empty()){
						//get an empty node and switch roles with it
						workersId.push_back(rank); //include itself in the worker list of the next master

						unsigned next_master = idleWorkers.back();
						idleWorkers.pop_back();

						//std::cout << "Node " << rank << " updating master idle " << std::endl;
						updateMasterIdle(next_master);
						
						//std::cout << "Node " << rank << " changing idle " << std::endl;
						changeMaster(next_master, 3);

						//std::cout << "Node " << rank << " updating master" << std::endl;
						updateMaster(next_master);
						
						workersId.clear();

						sendBuffer.clear();
						sendBuffer.push_back(0);

						variableOccurs.clear();
						unitClauses.clear();
						//((MiniSAT_1p*)workerSolver->getSolver())->getUsedVariables(variableOccurs,maxVar);
						//((MiniSAT_1p*)workerSolver->getSolver())->getUnitClauses(unitClauses,maxVar);
						
						((MiniSAT_1p*)globalSolver)->getUsedVariables(variableOccurs,maxVar);
						((MiniSAT_1p*)globalSolver)->getUnitClauses(unitClauses,maxVar);

						workerSolver->clearVariableOccurrence();
						
						for(std::set<int>::iterator it = variableOccurs.begin(); it != variableOccurs.end(); it++)
						{
							int v = *it;
							//unsigned x = v < 0 ? x = -v : x = v;
							assert( v > 0);
							workerSolver->insertVariableOccurrence(v);
							sendBuffer.push_back(v);
						}

						for(std::set<int>::iterator it = unitClauses.begin(); it != unitClauses.end(); it++)
						{
							int v = *it;
							unsigned x = v < 0 ? x = -v : x = v;
							if(variableOccurs.find(x)==variableOccurs.end())
							{
								workerSolver->insertVariableOccurrence(x);
								sendBuffer.push_back(x);
							} 

							std::vector<signed> tmp(1);
							tmp[0] = v;
							workerSolver->addClause(tmp);

						}

						sendBuffer[0] = sendBuffer.size();
						MPI_safe_send(sendBuffer,masterId,MPI_DATA_SHARED_VARIABLES);
						sendBuffer.clear();

						//((MiniSAT_1p*)workerSolver->getSolver())->deleteProof();
						//((MiniSAT_1p*)workerSolver->getSolver())->resetProof();

						//((MiniSAT_1p*)workerSolver->getSolver())->cleanLearnts(); // should we clean learnts?
						masterToWorker();
						dynamic_rounds = 0;

					} else {
						//disable this stage for now
						//control_factor = 0;

						//std::cout << "Node " << rank << " DYNAMIC BALANCING" << std::endl;

						//dynamic_balancing

						if(phase_transition){
							//at this point each node sends their occurrences

							//std::cout << "Node " << rank << " sending request to all workers " << std::endl;

							sendBuffer.clear();
							sendBuffer.push_back(0);
							for(unsigned t=0; t < workersId.size(); t++)
								MPI_safe_send(sendBuffer, workersId[t], MPI_DATA_OCCURRENCES); //send request to all workers
							sendBuffer.clear();


							//save the current worker occurrences
							variableOccurs.clear();
							unitClauses.clear();

							((MiniSAT_1p*)globalSolver)->getUsedVariables(variableOccurs,maxVar); //TO DO: this is too expensive, keep the variable occurrences from each solver
							((MiniSAT_1p*)globalSolver)->getUnitClauses(unitClauses,maxVar);

							//std::cout << "Node " << rank << " creating changeOccurrences!" << std::endl;
							for (unsigned i=0; i< (unsigned)cores; i++)
								changeOccurrences.push_back(new VariableOccurrence());

							for(std::set<int>::iterator it = variableOccurs.begin(); it != variableOccurs.end(); it++)
							{
								int v = *it;
								changeOccurrences[rank]->setOccurs(v);
							}

							for(std::set<int>::iterator it = unitClauses.begin(); it != unitClauses.end(); it++)
							{
								int v = *it;
								unsigned x = v < 0 ? x = -v : x = v;
								if(variableOccurs.find(x)==variableOccurs.end())
								{
									changeOccurrences[rank]->setOccurs(v);
								} 

							}

							//get all the occurrences from all workers
							//std::cout << "Node " << rank << " getting all the occurrences from all workers " << std::endl;

							for(unsigned t=0; t < (unsigned)cores; t++)
							{
								if(rank != t){
									int * worker_occurrences;
									MPI_Status s;
									//std::cout << "Node " << rank << " getting the occurrences from " << t << std::endl;
									MPI_wait_for(worker_occurrences,t, s);

									unsigned size = worker_occurrences[0];
									for(unsigned i=1; i < size; i++)
										changeOccurrences[t]->setOccurs(worker_occurrences[i]);									
								}
							}


							sendBuffer.clear();
							sendBuffer.push_back(0);
							for(unsigned t=0; t < (unsigned)cores; t++)
							{
								sendBuffer.push_back(0);
								for(unsigned v=1; v <= maxVar; v++)
								{
									if(changeOccurrences[t]->occurs(v))
										sendBuffer.push_back(v);
								}
							}
							sendBuffer[0] = sendBuffer.size();

							//sharedVariables.cleanShared();
							

							//send the occurrences to all workers - this information will be used if they become masters

							std::cout << "Node " << rank << " sending the occurrences to all workers! " << std::endl;

							for(unsigned t=0; t < workersId.size(); t++)
								MPI_safe_send(sendBuffer,workersId[t],MPI_DATA_OCCURRENCES);

							
							phase_transition = false;
						
						}

						variableOccurs.clear();
						unitClauses.clear();

						((MiniSAT_1p*)globalSolver)->getUsedVariables(variableOccurs,maxVar); //TO DO: this is too expensive, keep the variable occurrences from each solver
						((MiniSAT_1p*)globalSolver)->getUnitClauses(unitClauses,maxVar);

						//std::cout << "Node " << rank << " creating changeOccurrences!" << std::endl;
						for(std::set<int>::iterator it = variableOccurs.begin(); it != variableOccurs.end(); it++)
						{
							int v = *it;
							myVariables->setOccurs(v);
						}

						for(std::set<int>::iterator it = unitClauses.begin(); it != unitClauses.end(); it++)
						{
							int v = *it;
							unsigned x = v < 0 ? x = -v : x = v;
							if(variableOccurs.find(x)==variableOccurs.end())
							{
									myVariables->setOccurs(v);
							} 

						}

						unsigned num_myVariables = 0;
						for(unsigned v=1; v <= maxVar; v++)
						{
							if(myVariables->occurs(v))
								num_myVariables++;
						}
						
						std::cout << "Node " << rank << " has " << num_myVariables << std::endl;

						//asking for the next master

						//sendBuffer.clear();
						//sendBuffer.push_back(0);
						//for(unsigned t=0; t < workersId.size(); t++)
						//	MPI_safe_send(sendBuffer,workersId[t],MPI_DATA_DYNAMIC_SIZE);

						//sendBuffer.clear();
						//
						//unsigned next_master = 0;
						//unsigned min = INT_MAX;

						//for(unsigned t=0; t < workersId.size(); t++)
						//{

						//	int * dynamic_size;
						//	MPI_Status s;
						//	MPI_wait_for(dynamic_size, workersId[t], s);
						//	assert(s.MPI_TAG == MPI_DATA_DYNAMIC_SIZE);

						//	if((unsigned)dynamic_size[0] < min){
						//		min = dynamic_size[0];
						//		next_master = workersId[t];
						//	}
						//}

						//assert( min != INT_MAX );

						//std::cout << "Node " << rank << " found that the new master is " << next_master << " and has a clause size of " << min << std::endl;

						
						unsigned next_master = (rank + 1) % cores;
						//std::cout << "Node " << rank <<  " next master is " << next_master << std::endl;

						sendBuffer.clear();
						sendBuffer.push_back(next_master);

						//std::cout << "Node " << rank << " updating the new master on all workers " << std::endl;

						//std::cout << "Node " << rank << " trail : " << trail.size() << " reasons : " << reasons.size() << " masterTrail: " << trail_master.size() << 
						//	"sharedTrail: " << trail_shared.size() << std::endl;

						std::vector<signed> masterTrailReason;
						masterTrailReason.push_back(next_master);
						masterTrailReason.push_back(0);

						masterTrailReason.push_back(0);

						//std::vector<signed> reasonsSend;

						//for(unsigned t=0; t < trail.size(); t++)
						//{
						//	if(changeOccurrences[next_master]->occurs(trail[t]))
						//	{
						//	
						//	masterTrailReason.push_back(trail[t]);
						//	if(reasons[t]) reasonsSend.push_back(masterTrailReason.size()-2);
						//	}
						//}
						//masterTrailReason.push_back(0);
						//for(unsigned t=0; t < reasonsSend.size(); t++)
						//{
						//	masterTrailReason.push_back(reasonsSend[t]);
						//}
						//
						//masterTrailReason[1] = masterTrailReason.size();

						for(unsigned t=0; t < workersId.size(); t++)
						{
							//if(workersId[t] == next_master) MPI_safe_send(masterTrailReason,next_master,MPI_DATA_MASTER); 
							//else MPI_safe_send(sendBuffer,workersId[t],MPI_DATA_MASTER);
							//std::cout << "Node " << rank << " sending the new master to " << workersId[t] << std::endl;
							if(workersId[t] == next_master) MPI_safe_send(masterTrailReason,next_master,MPI_DATA_DYNAMIC_SIZE); 
							else MPI_safe_send(sendBuffer,workersId[t],MPI_DATA_DYNAMIC_SIZE);
						}

						masterTrailReason.clear();

						masterId = next_master;
						workersId.clear();
						sendBuffer.clear();

						//masterToWorker();

						
						for(unsigned v = 1 ; v <= maxVar ; v++)
						{
							if(sharedVariables.isShared(v)){
								masterSharedVariables.setShared(v, true);
							}
						}

						unsigned global_vars = ((MiniSAT_1p*)globalSolver)->numVars();
						unsigned global_clauses = ((MiniSAT_1p*)globalSolver)->numClauses();
						for(unsigned t=0; t < global_vars; t++) workerSolver->addVar();

						std::vector<signed> clause;
						for(unsigned t=0; t < global_clauses; t++)
						{
							((MiniSAT_1p*)globalSolver)->removeLastToClause(clause);
							workerSolver->addClause(clause);
						}

						((MiniSAT_1p*)globalSolver)->cleanLearnts();
						if (globalSolver) delete globalSolver;
						globalSolver = new MiniSAT_1p(m, false, MiniSAT_1p::LIFTING);
						while(globalSolver->numVars() < maxVar+1) globalSolver->addVar();

						//std::cout << "sending solving request to new master" << std::endl;
						
						//std::cout << "Node " << rank << " updating master" << std::endl;
						//updateMaster(next_master);

						dynamic_rounds = 0;

					}

				}

			}

		} 
        else 
        {
			quit = receiveTermination();

			if (!quit && rank!=masterId)
            {
				receiveTrail();
				solvePartitions();
			}

			if (rank == masterId)
            {
				//std::cout << "Node " << rank << " is changing from worker to master" << std::endl;
				
				changeToMaster = false;
				workersId.clear();
				//masterOccurrences.clear();
				//for (unsigned i=0; i< (unsigned)cores-1; i++)
				//	masterOccurrences.push_back(new VariableOccurrence());

				//for(unsigned t=0; t < (unsigned)cores; t++)
				//{
				//	if(rank != t){
				//		workersId.push_back(t);
				//		masterOccurrences.push_back(changeOccurrences[t]);
				//	}
				//}

				workerToMaster();
				control_clauses = globalSolver->numClauses();

				masterSharedVariables.cleanShared();

				sharedVariables.cleanShared();
				sharedVariables.cleanOccurrences();

				for(unsigned i=0; i < (unsigned)cores; i++)
					sharedVariables.add(changeOccurrences[i]);

				sharedVariables.update();
	
				statsSharedVariables();

				idleWorkers.clear();

				//std::cout << "Node " << rank << " [master - time] : " << (clock() - mpi_startTime)/(double)CLOCKS_PER_SEC << std::endl;
				//std::cout << "Node " << rank << " [master - control] : " << control_clauses << std::endl;

				//for(unsigned t=0 ; t < (unsigned)cores; t++)
				//{
				//	for(unsigned v=1; v <= maxVar; v++)
				//	{
				//		if(sharedVariables.occurs(v,t)){
				//			std::cout << "Partition " << t << " contains variable " << v << std::endl;
				//		}
				//	}
				//}

			}
		}

	}

	if (rank == masterId)
	{
		std::cout << "Total Workers: " << workersId.size()+1 << std::endl;
		std::cout << "Total Idle: " << cores - workersId.size()-1 << std::endl;

		//for(unsigned t=0 ; t < (unsigned)cores; t++)
		//{
		//	for(unsigned v=1; v <= maxVar; v++)
		//	{
		//		if(sharedVariables.occurs(v,t)){
		//			std::cout << "Partition " << t << " contains variable " << v << std::endl;
		//		}
		//	}
		//}

		sendBuffer.clear();
		sendBuffer.push_back(1); //sequential stage
		//sending termination signal to all the idle workers
		for(unsigned t=0; t < idleWorkers.size(); t++)
			MPI_safe_send(sendBuffer,idleWorkers[t],MPI_DATA_MASTER);

		//sendBuffer.clear();
		//sendBuffer.push_back(idleWorkers.size());
		//for(unsigned t=0; t < idleWorkers.size(); t++)
		//	sendBuffer.push_back(idleWorkers[t]);

		////std::cout << "Node " << rank << " is sending the idle Workers to node 0" << std::endl;
		//MPI_safe_send(sendBuffer,0,MPI_DATA_MASTER);
		//sendBuffer.clear();
		std::cout << "Total Time: " <<  (clock() - mpi_startTime)/(double)CLOCKS_PER_SEC << std::endl;
	}


	////OLD VERSION
	//if(rank == 0){

	//	idleWorkers.clear();
	//	int * receivIdleWorkers = NULL;
	//	MPI_Status s;

	//	//std::cout << "Node " << rank << " is receiving the idle Workers" << std::endl;
	//	MPI_wait_for(receivIdleWorkers,masterId,s);

	//	for(int t=0; t < receivIdleWorkers[0]; t++)
	//		idleWorkers.push_back(receivIdleWorkers[t+1]);

	//	//std::cout << "Node " << rank << " is sending the termination signal" << std::endl;

	//	sendBuffer.clear();
	//	sendBuffer.push_back(1);
	//	//sending termination signal to all the idle workers
	//	for(unsigned t=0; t < idleWorkers.size(); t++)
	//	{
	//		//std::cout << "Node " << rank << "is sending the termination signal to " << idleWorkers[t] << std::endl;
	//		MPI_safe_send(sendBuffer,idleWorkers[t],MPI_DATA_MASTER);
	//	}
	//}

	//synchronization at the end of the solving process
	if(rank != masterId) std::cout << "Node " << rank << " [solved] : variables " <<  ((MiniSAT_1p*)workerSolver->getSolver())->numVars() << " clauses " << ((MiniSAT_1p*)workerSolver->getSolver())->numClauses() << std::endl;
	else std::cout << "Node " << rank << " [solved] : variables " <<  ((MiniSAT_1p*)globalSolver)->numVars() << " clauses " << ((MiniSAT_1p*)globalSolver)->numClauses() << std::endl;
	std::cout << "Node " << rank << " [solved] : idle Time " << mpi_idleTime /(double)CLOCKS_PER_SEC << std::endl;
	std::cout << "Node " << rank << " [solved] : CPU Time " << (clock() - mpi_nodeTime - mpi_idleTime) /(double)CLOCKS_PER_SEC << std::endl;
	std::cout << "Node " << rank << " [solved] : efficiency " << std::setprecision(5) <<  
	((clock() - mpi_nodeTime - mpi_idleTime) /(double)CLOCKS_PER_SEC) / ((clock() - mpi_startTime)/(double)CLOCKS_PER_SEC) << std::endl;
	
	MPI_Barrier(MPI_COMM_WORLD);

	

	return res==M_TRUE;

}

bool MPDeSATDynWorker::findDisagreementDyn(void)
{

	bool res = false;

	std::vector<signed> new_trail;
	std::vector<unsigned> global_model(maxVar+1,M_UNDEF);
	std::vector<unsigned> preferencesTrue(maxVar+1,0);
	std::vector<unsigned> preferencesFalse(maxVar+1,0);

	for (unsigned i=0; i< workersId.size(); i++)
	{
		unsigned w = workersId[i];

		for (unsigned t=0; t < models[i].size(); t++)
		{
			signed mvi = models[i][t];

			//if (sharedVariables.occurs(mvi, w)) //it is required to check if it is shared?
			//{
				//std::cout << "variable " << mvi << " partition " << i << std::endl;

				if (mvi < 0)
                {
					preferencesFalse[-mvi]++;
					if(global_model[-mvi] == M_UNDEF || global_model[-mvi] == M_FALSE){
						global_model[-mvi] = M_FALSE;
					} else {
						//conflict in the set of models
						global_model[-mvi] = 4; //conflict

					}
				} 
                else 
                {
					preferencesTrue[mvi]++;
					if(global_model[mvi] == M_UNDEF || global_model[mvi] == M_TRUE){
						global_model[mvi] = M_TRUE;
					} else {
						//conflict in the set of models
						global_model[mvi] = 4; //conflict

					}

				}
			//}
		}
	}

	for (signed t=1; t < (signed)global_model.size(); t++)
    {
		if (global_model[t] == 4)
        {

			signed v = t;

			if(preferencesTrue[t] > preferencesFalse[t])
				v = -t;

			if (std::find(trail.begin(), trail.end(), v) == trail.end())
				new_trail.push_back(v);
		}
	}


	for (unsigned i=0; i<new_trail.size(); i++)
	{
		trail.push_back(new_trail[i]);
		reasons.push_back(false);
	}

	//std::cout << "Trail: ";
	//for(unsigned t=0; t < trail.size(); t++)
	//{
	//	std::cout << trail[t] << " ";
	//}
	//std::cout << std::endl;

	global_model.clear();
	preferencesFalse.clear();
	preferencesTrue.clear();

	return new_trail.size()>0;

}

bool MPDeSATDynWorker::waitsForWork(void)
{
	bool change_master = false;
	MPI_Status s;
	int * masterInformation = NULL;

	//worker waits until information regarding the new master
	//std::cout << "Node " << rank << " is waiting for work" << std::endl;
	MPI_wait_for(masterInformation, masterId, s);
	assert (s.MPI_TAG == MPI_DATA_MASTER);
	//std::cout << "Node " << rank << " is getting work from " << s.MPI_SOURCE << std::endl;

	//updating master information
	if (masterInformation[0] == 4)
    {
		std::cout << "Node " << rank << " is changing master to " << masterInformation[1] << std::endl;
		change_master = true;
		masterId = masterInformation[1];		
		return change_master;
	}

	//solving stage
	if (masterInformation[0] == 3 || masterInformation[0] == 5)
	{
		//this is a new master and should go directly to the solving stage

		if (masterInformation[0] == 3) solving_stage = true;
		std::cout << "Node " << rank << " is the new master" << std::endl;

		int * nodeInformation = NULL;
		MPI_wait_for(nodeInformation, masterId, s);
		std::cout << "Node " << rank << " got nodeInformation from " << s.MPI_SOURCE << std::endl;

		workersId.clear();
		idleWorkers.clear();

		//receive the information regarding the working/idle nodes
		for (int t=0; t < cores; t++)
		{
			if (t != rank)
			{
				if(nodeInformation[t] == 0) 
					idleWorkers.push_back(t);
				if(nodeInformation[t] == 1)
					workersId.push_back(t);
			}
		}		
	}

	//sequential stage
	//the formula has been solved without the participation of this node
	if (masterInformation[0] == 1)
		quit = true;
	else
    {
		//std::cout << "Node " << rank << " setting the new master Id to " << masterInformation[1] << std::endl;
		masterId = masterInformation[1];

		//worker receives the problem statment from the master
		std::cout << "Node " << rank << " receiving the problem" << std::endl;
		receiveProblem();

		//worker receives the clauses from the master
		std::cout << "Node " << rank << " receiving clauses" << std::endl;
		receiveClauses();

		//problem splitting stage
		if(masterInformation[0] == 3 || masterInformation[0] == 5) masterId = rank;
	}

	std::cout << "Node " << rank << " [active] : " << (clock() - mpi_startTime)/(double)CLOCKS_PER_SEC << std::endl;
	mpi_nodeTime = clock();	

	return change_master;
}

void MPDeSATDynWorker::splitHalf(unsigned to_worker)
{
	//send information regarding the problem statment
	std::vector<int> problem_statment(2);
	problem_statment[0] = maxVar;
	problem_statment[1] = maxClauses;
	MPI_safe_send(problem_statment, to_worker, MPI_DATA_PROBLEM);

	//send half of the clauses
	sendBuffer.clear();
	((MiniSAT_1p*)workerSolver->getSolver())->splitDB(sendBuffer);
	MPI_safe_send(sendBuffer, to_worker, MPI_DATA_CLAUSE);

	//std::cout << "Node " << rank << " splitting my formula into half and giving it to " << to_worker << " (size: " << sendBuffer[0] << " , clauses: " << sendBuffer[1] << ") " << std::endl;
	sendBuffer.clear();
}

void MPDeSATDynWorker::splitN(unsigned workers)
{
	//assumes that the formula fits the memory of the sequential solver
		
	unsigned problem_workers = workers;
	//if(idleWorkers.size() % 2 == 0) problem_workers = idleWorkers.size()/2;
	//else problem_workers = (idleWorkers.size()-1)/2;

	//unsigned clauses = ((MiniSAT_1p*)workerSolver->getSolver())->numClauses();
	unsigned clauses = ((MiniSAT_1p*)sequentialSolver)->numClauses();
	std::cout << "Node " << rank << " has " << clauses << " clauses." << std::endl;

	d->setClauseMax(clauses);
	d->setPartitions(problem_workers);

	control_clauses = (unsigned)ceil((double)clauses/problem_workers);

	std::vector<signed> literals;
	
	unsigned next_worker = idleWorkers.back();
	idleWorkers.pop_back();

	workingWorkers.push_back(next_worker);
	workersId.push_back(next_worker);

	sendBuffer.clear();
	
	//send information regarding the node that will send half of clauses
	sendBuffer.push_back(2); //problem splitting stage
	sendBuffer.push_back(rank);
	MPI_safe_send(sendBuffer,next_worker,MPI_DATA_MASTER);
	sendBuffer.clear();
	
	std::vector<int> problem_statment(2);
	unsigned w = d->where(literals);
	
	//while (((MiniSAT_1p*)workerSolver->getSolver())->numClauses() > 0)
	while (((MiniSAT_1p*)sequentialSolver)->numClauses() > 0)	
	{
		unsigned next = d->where(literals);
		//((MiniSAT_1p*)workerSolver->getSolver())->removeLast(sendBuffer);
		((MiniSAT_1p*)sequentialSolver)->removeLast(sendBuffer);
		
        if (next != w)
        {
            std::cout << "[" << rank << "] sending problem statement to " << next_worker << std::endl;
			problem_statment[0] = maxVar;
			problem_statment[1] = maxClauses;
			MPI_safe_send(problem_statment, next_worker, MPI_DATA_PROBLEM);
			            
			MPI_safe_send(sendBuffer, next_worker, MPI_DATA_CLAUSE);
			sendBuffer.clear();
		
			//if(next == problem_workers-1){
			//	//std::cout << "Node " << rank << " keeping the clauses from the last partition " << std::endl;
			//	break;
			//}

			next_worker = idleWorkers.back();
			idleWorkers.pop_back();

			workingWorkers.push_back(next_worker);
			workersId.push_back(next_worker);
		
			//send information regarding the node that will send half of clauses
			sendBuffer.push_back(2); //problem splitting stage
			sendBuffer.push_back(rank);
			MPI_safe_send(sendBuffer,next_worker,MPI_DATA_MASTER);
			sendBuffer.clear();
		} 
		
		w = next;
	}

	//sending clauses to the last worker
	problem_statment[0] = maxVar;
	problem_statment[1] = maxClauses;
	MPI_safe_send(problem_statment, next_worker, MPI_DATA_PROBLEM);
	problem_statment.clear();

	MPI_safe_send(sendBuffer, next_worker, MPI_DATA_CLAUSE);
	sendBuffer.clear();
}

void MPDeSATDynWorker::changeMaster(unsigned next_master, unsigned stage)
{
	//std::cout << "Node " << rank << " changing master to " << next_master << std::endl;
	sendBuffer.clear();
	sendBuffer.push_back(stage); //request the node to go to the solving stage
	sendBuffer.push_back(rank);
	
	std::cout << "sending solving request to new master" << std::endl;
	MPI_safe_send(sendBuffer,next_master,MPI_DATA_MASTER);
	sendBuffer.clear();

	for(int t=0; t < cores; t++)
		sendBuffer.push_back(0);

	for(unsigned t=0; t < workersId.size(); t++)
		sendBuffer[workersId[t]] = 1;

	sendBuffer[rank] = 1;

	//if(rank!=0) workersId.clear();

	std::cout << "sending node information to " << next_master << std::endl;
	MPI_safe_send(sendBuffer,next_master,MPI_DATA_MASTER);
	sendBuffer.clear();

	//send information regarding the problem statment
	std::vector<int> problem_statment(2);
	problem_statment[0] = maxVar;
	problem_statment[1] = ((MiniSAT_1p*)globalSolver)->numClauses();
		//((MiniSAT_1p*)workerSolver->getSolver())->numClauses();
	//problem_statment[1] = maxClauses;
	//std::cout << "sending problem statment to " << next_master << std::endl;
	MPI_safe_send(problem_statment, next_master, MPI_DATA_PROBLEM);

	//send empty set of clauses		
	sendBuffer.push_back(0);
	sendBuffer.push_back(0);
	//std::cout << "sending set of clauses to " << next_master << std::endl;
	MPI_safe_send(sendBuffer, next_master, MPI_DATA_CLAUSE);
	sendBuffer.clear();
}

void MPDeSATDynWorker::updateMasterIdle(unsigned next_master)
{
	sendBuffer.clear();
	sendBuffer.push_back(4); //Nodes that are idle just update their reference to the current master and return to an idle position
	sendBuffer.push_back(next_master);
	for(unsigned t=0; t < idleWorkers.size(); t++)
	{
		MPI_safe_send(sendBuffer,idleWorkers[t],MPI_DATA_MASTER);
	}
	sendBuffer.clear();
}

void MPDeSATDynWorker::updateMaster(unsigned next_master)
{
	sendBuffer.clear();
	sendBuffer.push_back(next_master);
	for(unsigned t=0; t < workersId.size(); t++)
	{
		std::cout << "update Master: " << workersId[t] << " to " << next_master << std::endl;
		if(rank != workersId[t]) MPI_safe_send(sendBuffer,workersId[t],MPI_DATA_MASTER);
		else masterId = next_master;

	}
	sendBuffer.clear();
}

void MPDeSATDynWorker::split_problem(double cutoff_timeout_sequential)
{
	//stage 2: split the formula into a set of problem nodes

	bool partition_formula = false;
	//unsigned active_nodes = 0;
	while (!partition_formula)
	{
		//current node tries to solve the formula 
		clock_t start = clock();
		//std::cout << "Node " << rank << " is searching in a subset of the formula" << std::endl;
		//((MiniSAT_1p*)workerSolver->getSolver())->setTimeout((int)cutoff_timeout_sequential);
		((MiniSAT_1p*)workerSolver->getSolver())->setLimitResources(limit_conflicts,limit_clauses,0);
		workerSolver->solve();
		double solving_time = (clock() - start)/(double)CLOCKS_PER_SEC;
		//std::cout << "Node " << rank << " finished searching on a subset of the formula: " << solving_time << std::endl;

		//if(solving_time > cutoff_timeout_sequential)
		if(((MiniSAT_1p*)workerSolver->getSolver())->getResources())
		{
			//std::cout << "Node " << rank << " splitting my formula into half" << std::endl;
			//timeout on the current formula -> splits the formula

			if(rank != 0)
			{
				//sends request to node 0
				sendBuffer.clear();
				sendBuffer.push_back(1);
				//std::cout << "Node " << rank << " sending request to process 0" << std::endl;
				MPI_safe_send(sendBuffer,0,MPI_DATA_MASTER);

				//receive next_worker from node 0
				int * next_worker_rec;
				MPI_Status s;
				//std::cout << "Node " << rank << " waiting for the next worker from node 0" << std::endl;
				MPI_wait_for(next_worker_rec, 0, s);

				//std::cout << "Node " << rank << " next_worker_rec:1 " << next_worker_rec[0] << " next_worker_rec:2 " << next_worker_rec[1] << std::endl;
				assert(next_worker_rec[0] == 2); //problem splitting stage
				unsigned next_worker = next_worker_rec[1];				

				if(next_worker == 0) partition_formula = true;
				else splitHalf(next_worker);
			} 
            else 
            {
				std::vector<unsigned> new_working;
				std::vector<unsigned> extra_working;

				//checks for pending requests so that the master information is updated
				for(unsigned t=0; t < workingWorkers.size(); t++){
					int * request;
					MPI_Status s;
					MPI_wait_for(request,workingWorkers[t],s);

					if(request[0] == 1){

						unsigned next_worker = 0;
						if(!idleWorkers.empty())
						{
							next_worker = idleWorkers.back();
							idleWorkers.pop_back();
							extra_working.push_back(next_worker);
						} else {
							//next_worker = 0; //there are no more available nodes (we must keep one node to be the master in the next stage)
							new_working.push_back(workingWorkers[t]);
						}

						sendBuffer.clear();
						sendBuffer.push_back(2);
						sendBuffer.push_back(workingWorkers[t]);
						//std::cout << "Node " << rank << " sending master " << workingWorkers[t] << " to " << next_worker << std::endl; 
						MPI_safe_send(sendBuffer,next_worker,MPI_DATA_MASTER); //sends the master information to the idle worker

						sendBuffer.clear();
						sendBuffer.push_back(2);
						sendBuffer.push_back(next_worker);
						//std::cout << "Node " << rank << " sending worker " << next_worker << " to " << workingWorkers[t] << std::endl;
						MPI_safe_send(sendBuffer,workingWorkers[t],MPI_DATA_MASTER); //send the worker information to the new master

						if(next_worker !=0)
							workersId.push_back(next_worker);
					} 
                    else 
                    {
						//active_nodes++;
						new_working.push_back(workingWorkers[t]);
					}					
				}

				for(unsigned t=0; t < extra_working.size(); t++)
					workingWorkers.push_back(extra_working[t]);

				for(unsigned t=0; t < new_working.size(); t++)
				{
					std::vector<unsigned>::iterator it = std::find(workingWorkers.begin(), workingWorkers.end(), new_working[t]);
					unsigned removePosition = std::distance( workingWorkers.begin(), it);
					workingWorkers.erase(workingWorkers.begin()+removePosition);
				}

				//splits the formula
				unsigned next_worker;
				if (idleWorkers.empty())
                {
					partition_formula = true;
				} 
                else 
                {
					next_worker = idleWorkers.back();
					idleWorkers.pop_back();
					workersId.push_back(next_worker);

					workingWorkers.push_back(next_worker);

					sendBuffer.clear();
					//send information regarding the node that will send half of clauses
					sendBuffer.push_back(2); //problem splitting stage
					sendBuffer.push_back(rank);
					//std::cout << "Node " << rank << " sending master " << rank << " to " << next_worker << std::endl;
					MPI_safe_send(sendBuffer,next_worker,MPI_DATA_MASTER);

					splitHalf(next_worker);		
				}
			}
		} 
        else 
        { 
			//std::cout << "Node " << rank << " last problem spitting iteration: " << solving_time << std::endl;
			partition_formula = true;
		}

	}

	if (rank == 0) 
    {
		while(!workingWorkers.empty())
		{
			std::vector<unsigned> new_working;

			//checks for pending requests to that the master information is updated
			for (unsigned t=0; t < workingWorkers.size(); t++)
            {
				int * request;
				MPI_Status s;
				MPI_wait_for(request,workingWorkers[t],s);

				if (request[0] == 1)
                {
					unsigned next_worker;
					next_worker = idleWorkers.back();
					idleWorkers.pop_back();
					workingWorkers.push_back(next_worker);

					sendBuffer.clear();
					sendBuffer.push_back(2); //problem splitting stage
					sendBuffer.push_back(next_worker);
					//std::cout << "Node " << rank << " sending master " << next_worker << " to " << workingWorkers[t] << std::endl;
					MPI_safe_send(sendBuffer,workingWorkers[t],MPI_DATA_MASTER);
				} 
                else 
                {
					//active_nodes++;
					new_working.push_back(workingWorkers[t]);
				}				
			}

			for(unsigned t=0; t < new_working.size(); t++)
			{
				std::vector<unsigned>::iterator it = std::find(workingWorkers.begin(), workingWorkers.end(), new_working[t]);
				unsigned removePosition = std::distance( workingWorkers.begin(), it);
				workingWorkers.erase(workingWorkers.begin()+removePosition);
			}
		}
	}
}

void MPDeSATDynWorker::sendStatus(bool status)
{
	std::vector<signed> buf;
	status ? buf.push_back(1) : buf.push_back(0);

	for(unsigned id = 0; id < workersId.size(); id++)
		MPI_safe_send(buf,workersId[id],MPI_DATA_QUIT);
}