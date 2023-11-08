// Copyright (C) 2011 Microsoft Research
// RM Martins, 2012
/*
 * This file contains the worker from the static topologies: flat, tree and tree-flat
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
#include "mpdesat_worker.h"
#include "returnvalue.h"

#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib") //added


MPDeSATWorker::MPDeSATWorker(ExpressionManager &m, int size, int rank, const char *filename, int c, InterpolationMode i, int master, std::vector<int> worker) :
MPDeSATNode(m, rank, master, worker),
	size(size),
	cores(c),
	filename(filename),
	n_workers(size-1),
	d(NULL),
	workerSolver(NULL),
	sendingTo(-1),
	bufferSize(10485760)
{

	workerSolver = new Partition(m, masterSharedVariables, rank, verbosity);
	if(i!= NONE) workerSolver->setInterpolationMode(i, sharedVariables);
	
	d = new BatchDecomposition();
	d->setPartitions(n_workers);
	assumptions.resize(n_workers);
	interpolants.resize(n_workers, m.mkNil());

	workerOcurrences = new VariableOccurrence();

	for (unsigned t=0; t< n_workers; t++)
		masterOccurrences.push_back(new VariableOccurrence());

}

MPDeSATWorker::~MPDeSATWorker(void)
{
	if (workerSolver) delete workerSolver;
	if (d) delete d;
}

//receive methods
void MPDeSATWorker::receiveSharedVariables(MPI_Status &s)
{

	int * shared_variables = NULL;

	clock_t	before = clock();
	MPI_wait_for(shared_variables, masterId, s);
	mpi_idleTime += clock()-before;

	assert(s.MPI_TAG == MPI_DATA_SHARED_VARIABLES);

	masterSharedVariables.update(shared_variables);	
}

void MPDeSATWorker::receiveClauses(MPI_Status &s)
{

	int * clauses = NULL;
	MPI_wait_for(clauses, masterId, s);

	unsigned num_current_clauses = 0;

	//clauses contain an array with clauses: 1st element is the size of the array, 2nd element is the number of clauses in the array
	//3rd element is the size of the clause, then the literals that compose that clause (repeat until end of array)

	assert(s.MPI_TAG == MPI_DATA_CLAUSE);

	int size = clauses[0]; //first position contains the size of the array
	int num_clauses = clauses[1]; //second position contains the number of clauses in the array

	if (size == 0)
    {	
		s.MPI_TAG = MPI_DATA_QUIT; //there are no more clauses to be read
		sendBuffer.clear();
		sendBuffer.push_back(0);

		for(unsigned i=0; i < workersId.size(); i++)
			MPI_safe_send(sendBuffer,workersId[i], MPI_DATA_CLAUSE);
	} else {
		d->setClauseMax(num_clauses);
		d->setPartitions(workersId.size());

		for(int i=2; i < size; i++)
		{

			int clauseSize = clauses[i];
			std::vector<int> literals(clauseSize);
			sendBuffer.clear();
			sendBuffer.push_back(0);
			sendBuffer.push_back(0);
			sendBuffer.push_back(clauseSize);

			unsigned w = d->where(literals);

			for(int j=0; j < clauseSize; j++)
			{
				literals[j] = clauses[++i];
				sendBuffer.push_back(literals[j]);

				if(n_workers != 0){

					//preserves the order of the workers occurrences so that each of them will get associated with the correct worker
					std::vector<unsigned>::iterator it = std::find(workersId.begin(), workersId.end(), workersId[w]);
					unsigned workerPosition = std::distance(workersId.begin(), it);

					masterOccurrences[workerPosition]->setOccurs(literals[j]);

				}
			}

			num_current_clauses++;

			if(n_workers == 0){
				//printClause(literals, rank);
				workerSolver->addClause(literals);	
			}
			else {

				//TO DO: use an array instead of sending the clause one by one
				sendBuffer[0] = sendBuffer.size();
				sendBuffer[1] = num_current_clauses;
				num_current_clauses = 0;
				//printClause(sendBuffer, rank);
				MPI_safe_send(sendBuffer, workersId[w], MPI_DATA_CLAUSE);

			}

		}


	}

	sendBuffer.clear();	
}

void MPDeSATWorker::receiveProblem(MPI_Status &s)
{
	clock_t before = clock();

	int * problem_statment = NULL;
	//root worker sends the information about the formula to all workers
	MPI_wait_for(problem_statment, MPI_MASTER, s);

	assert(s.MPI_TAG == MPI_DATA_PROBLEM);

	maxVar = problem_statment[0];
	maxClauses = problem_statment[1];
	
	while(s.MPI_TAG != MPI_DATA_QUIT)
		receiveClauses(s);

	mpi_idleTime += clock()-before;

}

void MPDeSATWorker::receiveTrail(void)
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

	if(mpi_verbose > 0) std::cout << "Node " << rank  << " received trail: ";

	for(int i = 1; i < assump[0]; i++)
	{
		if(mpi_verbose > 0) std::cout << assump[i] << " ";
		//if(occurrences->occurs(assump[i])) //only considers as assumptions the variables that occur in this partition
		assumptions.push_back(assump[i]);
		trail.push_back(assump[i]);
		reasons.push_back(true);
		trail_master.push_back(assump[i]);

		if(!sharedVariables.isShared(assump[i])){
			sharedVariables.setShared(assump[i],true);
			trail_shared.push_back(assump[i]);
		}

	}

	if(mpi_verbose > 0) std::cout << std::endl;
}

bool MPDeSATWorker::receiveTermination()
{
	MPI_Status s;
	bool termination = false;

	int * cmd = NULL;

	clock_t before = clock();
	MPI_wait_for(cmd,masterId,s);
	mpi_idleTime += clock()-before;

	assert(s.MPI_TAG = MPI_DATA_QUIT);

	cmd[0] == 0 ? termination = false : termination = true;
	
	return termination;
}

//send methods
void MPDeSATWorker::sendTrail(void)
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
			ModelValue mv = workerSolver->get(v);
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

void MPDeSATWorker::sendTermination(bool status)
{
	sendBuffer.clear();

	status ? sendBuffer.push_back(1) : sendBuffer.push_back(0);

	for(unsigned id = 0; id < workersId.size(); id++)
		MPI_safe_send(sendBuffer,workersId[id],MPI_DATA_QUIT);

	sendBuffer.clear();
}

void MPDeSATWorker::sendSolution(void)
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
			sendBuffer.push_back(var);
			break;
		case M_FALSE:
			if(mpi_verbose > 0) std::cout << -var << " ";
			sendBuffer.push_back(-var);
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
Expression MPDeSATWorker::solveWorker(void)
{

	clock_t before = clock();
	workerInterpolant = m.mkNil();
	Expression t = workerSolver->getInterpolant(trail);

	if(verbose > 1)
	{
		std::cout << "Node " << rank << "Finding global assignment..." << std::endl;
		std::cout << "Node " << rank << " Trail size: " << trail.size() << std::endl;

		for (unsigned i=0; i<trail.size(); i++)
		{
			std::cout << trail[i];
			if (reasons[i]) std::cout << "!";
			std::cout << " ";
		}
		std::cout << std::endl;
	}

	bool r=workerSolver->em().isTrue(t);

	while (!r && trail.size()>0)
	{ 

		std::vector<signed> temp;
		((MiniSAT_1p*)workerSolver->getSolver())->getConflict(temp);

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

			r=workerSolver->solve(trail);

		}
	}
	mpi_solvingTime += clock()-before;
	
	//t is already false when the unsatisfiability does not depend on the assumptions!
	//there is no need to check if the unsatisfiability depends on the assumptions

		
	return t;

}

void MPDeSATWorker::solveSolution(void)
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

void MPDeSATWorker::solvePartitions()
{
	//PROCESS_MEMORY_COUNTERS pmc;
	//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc));
	//saveMemoryUsed((double)pmc.WorkingSetSize);

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

	sendBuffer.clear();
	workerInterpolant = m.mkNil();
	t = m.mkNil();

}

bool MPDeSATWorker::importInterpolants(void)
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


			((MiniSAT_1p*)workerSolver->getSolver())->clearNewClauses();
			(workerSolver->getSolver())->addConstraint(itp);
			((MiniSAT_1p*)workerSolver->getSolver())->addNewClauses();

			did_something=true;

		}

		interpolants[i] = m.mkNil();
	}

	mpi_solvingTime = clock() - before;
	return did_something;
}

bool MPDeSATWorker::waitInterpolantSolution(bool quit)
{
	bool sat = true;

	//PROCESS_MEMORY_COUNTERS pmc;
	//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc));
	//saveMemoryUsed((double)pmc.WorkingSetSize);

	clock_t before;

	MPI_Status s;

	interpolants.clear();
	bool all_sat = true;
	models.clear();

	//PROCESS_MEMORY_COUNTERS pmc;
	//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc));
	//std::cout << "Node " << rank << " WorkingSetSize: " << pmc.WorkingSetSize / 1048576 << " MB" << std::endl;

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

			if(mpi_verbose > 0) std::cout << "Node " << rank << " received from Node " << workersId[id] << "(Model) : ";;

			unsigned size = buf_interpolant[1];

			for(unsigned i = 2; i < size; i++)
			{

				if(mpi_verbose > 0)
					std::cout << buf_interpolant[i] << " ";

				worker_model.push_back(buf_interpolant[i]);

			}

			if(mpi_verbose > 0) std::cout << std::endl;
			models.push_back(worker_model);

		} else {

			interpolants.push_back(receiveInterpolant(buf_interpolant));
			if(mpi_verbose > 0) 
			{
				std::string str = m.toString(interpolants.back());
				std::cout << "Node " << rank << " received from Node " << workersId[id] << " (ITP) : " << str.c_str() << std::endl;
			}
			all_sat = false;	
		}		
	}


	if(all_sat)
	{

		//Node received only models from its children

		all_sat_iterations++;
		assert(models.size() == workersId.size());

		clock_t before = clock();

		if (!findDisagreement()){
			//build a solution from the models of the children
			solveSolution();
		}
		else {

			//The models that the node received do not agree with each other
			//worker sends status, and a new trail to the children

			//sendInterpolantSolution(quit);
			sat = false;

		}

		clock_t mpi_solvingTime = clock() - before;

		interpolants.clear();


	}
	else if (!importInterpolants())
	{
		sendSolution();
	} else {

		//formula is UNSAT : send a new trail to children or an interpolant to master
		//sendInterpolantSolution(quit);
		sat = false;

	}

	return sat;

}

bool MPDeSATWorker::sendInterpolantSolution(bool quit)
{
	bool sat = false;
	//PROCESS_MEMORY_COUNTERS pmc;
	//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc));
	//saveMemoryUsed((double)pmc.WorkingSetSize);

	//std::cout << "Node " << rank << " WorkingSetSize: " << pmc.WorkingSetSize / 1048576 << " MB" << std::endl;

	sendBuffer.clear();
	Expression t = solveWorker();

	if(!workerSolver->em().isTrue(t)){

		//formula is UNSAT 

		workerInterpolant = m.mkNil();
		workerInterpolant = m.duplicate(t, workerSolver->em()); //why duplicate in the MPI version?
		t = m.mkNil();

		//std::string filename = "filename;"
		//std::stringstream convert; // stringstream used for the conversion
		//convert << rank;
		//convert << interpolants_exported;
		//filename += convert.str();
		//m.toDot(workerInterpolant,filename);

		//std::cout << "file has been written" << std::endl;

		if(mpi_verbose) std::cout << "Node " << rank << " interpolant export (ITP) " << m.toString(workerInterpolant) << std::endl;

		//std::vector<signed> buf;
		sendBuffer.push_back(0); //reserve space for the size of the buffer
		//m.toBuffer(sendBuffer,workerInterpolant);

		//m.sizeRBC(workerInterpolant, rank);
		m.toMPI(sendBuffer, workerInterpolant);
		sendBuffer[0] = sendBuffer.size();

		interpolants_exported++;
		sizeInterpolant += sendBuffer.size();
		if(sendBuffer.size() > maxSizeInterpolant){
			maxSizeInterpolant = sendBuffer.size();
			//std::cout << "Node " << rank << " generates an interpolant of size: " << maxSizeInterpolant << std::endl;
		}
		MPI_safe_send(sendBuffer,masterId,MPI_DATA_INTERPOLANT);
		sendBuffer.clear();
		

	} else {
		//send the status and trail to the children

		t = m.mkNil();

		rounds++;
		if(mpi_verbose > 0)
		{
			std::cout << "Node " << rank << " Round " <<  rounds << std::endl << std::endl;
		}
		sendTermination(quit);
		sendTrail();

		interpolants.clear();
		//waitInterpolantSolution(quit);
		sat = true;

	}

	return sat;
}

bool MPDeSATWorker::solve(void)
{
	ModelValue res = M_UNDEF;

	MPI_Status s;
	clock_t start = clock();

	clock_t before = clock();
	receiveProblem(s); //receive information regarding the number of variables and clauses of the formula
	mpi_loadTime += clock()-before;

	for(unsigned i=0; i < workersId.size(); i++)
		sharedVariables.add(masterOccurrences[i]); //stores the occurrences of the literals for each partition

	sharedVariables.update(); //updates the shared variables from the bottom level

	receiveSharedVariables(s); //from the top level

	if(n_workers !=0) sendSharedVariables(); //to the bottom level

	statsSharedVariables();

	before = clock();
	MPI_Barrier(MPI_COMM_WORLD); //synchronization of the topology to start solving the problem
	mpi_idleTime += clock()-before;

	rounds = 0;

	//TO DO: for trivial instances that are solved during reading we can send a termination flag across the network

	bool quit = false;

	//TO DO: only send assumptions that refer to variables of that partition
	while(workerSolver->numVars() < maxVar+1) workerSolver->addVar(); //creates dummy variables if necessary

	//	clock_t before;

	//for(unsigned t=1; t <= maxVar; t++)
	//{
	//	if(sharedVariables.isShared(t))
	//		std::cout << "Node " << rank << " variable " << t << " is shared " << std::endl;

	//	if(masterSharedVariables.isShared(t))
	//		std::cout << "Node " << rank << " variable " << t << " is master shared " << std::endl;

	//}

	
	 //for(unsigned t=1; t <= maxVar; t++)
		//if(workerSolver->variablesOccurrence().occurs(t)) std::cout << "Node " <<  rank << " has variable " << t << std::endl;

	while(!quit)
	{

		//PROCESS_MEMORY_COUNTERS pmc;
		//GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc));
		//saveMemoryUsed((double)pmc.WorkingSetSize);

		quit = receiveTermination();

		if(!quit){

			receiveTrail();

			if(n_workers !=0){

				rounds++;
				if(mpi_verbose > 0){
					std::cout << "Node " << rank << " Round " <<  rounds << std::endl << std::endl;
				}

				bool res = true;
				while(res){
					res = sendInterpolantSolution(quit);
					if(res) res = !waitInterpolantSolution(quit);
				}

			} else solvePartitions();

		}

	}

	if(n_workers != 0) sendTermination(quit);

	mpi_totalTime += clock() - start;
	if(mpi_stats > 0) printStats();

	MPI_Barrier(MPI_COMM_WORLD); //synchronization of the topology for a clean exit

	return true;
}