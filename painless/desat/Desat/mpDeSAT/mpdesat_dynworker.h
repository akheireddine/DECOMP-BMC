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

#ifndef _MPDESAT_DYNWORKER_H_
#define _MPDESAT_DYNWORKER_H_

#include <time.h>

#include <vector>

#include <minisat1p.h>
#include "mpdesat_node.h"
#include <interpolation_mode.h>
#include <decomposition.h>
#include <partition.h>
#include <shared_variables.h>

#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib") //added

class MPDeSATDynWorker : public MPDeSATNode
{
public:
	MPDeSATDynWorker(ExpressionManager &m, int size, int rank, const char *filename, int c, InterpolationMode i, int master, std::vector<int> worker);
	virtual ~MPDeSATDynWorker(void);

	virtual bool addClause(const std::vector<signed> &literals) { 

		if (sequential)
        {
		    //    if(literals.size() == 1)
		    //{
		    //	std::cout << "UNIT CLAUSE: " << literals[0] << std::endl;
		    //	//variableOccurs.insert(literals[0]);
		    //}
		
		    for (unsigned i=0; i<literals.size(); i++)
		    {
			        unsigned var = (literals[i]<0) ? -literals[i] : literals[i];
			        while (var>=numVars()) sequentialSolver->addVar();
		    }
		    //printClause(literals);

	        return sequentialSolver->addClause(literals);
	    }
		else 
            return workerSolver->addClause(literals);
	}

	virtual bool addUnit(signed l) { throw std::exception("NYI: addUnit"); }
	virtual void setVariableMax(unsigned n) { 	
		maxVar = n;
		while(workerSolver->numVars() < maxVar+1) workerSolver->addVar(); //creates dummy variables if necessary
	} 
	virtual void setClauseMax(unsigned n) { maxClauses = n; }

	virtual unsigned numClauses(void) const { return sequentialSolver->numClauses(); }
	virtual unsigned numVars(void) const { return sequentialSolver->numVars(); }
	virtual signed addVar(void) { return sequentialSolver->addVar(); }

	virtual bool solve(void);
	virtual bool solve(const std::vector<signed> &assumptions) { throw std::exception("NYI: solve with assumptions"); }

	virtual ModelValue get(signed l) const { throw std::exception("NYI: get ModelValue"); }  

	virtual Expression getInterpolant(const std::vector<signed> &A) { throw std::exception("NYI: getInterpolant"); }
	virtual bool addConstraint(const Expression &e) { throw std::exception("NYI: addConstraint @ worker.h"); }
	virtual signed addExtension(const Expression &e) { throw std::exception("NYI: addExtension"); }

	virtual Expression getModel(void) const { throw std::exception("NYI: getModel"); }

	virtual void setVerbose(int v) { mpi_verbose = v; }

protected:
	int size;
	int cores;
	const char *filename;
	VariableOccurrence* myVariables;
	VariableOccurrence* workerOccurrences;
	std::vector<VariableOccurrence*> masterOccurrences;
	Expression workerInterpolant;
	unsigned n_workers;
	Decomposition *d;
	Partition *workerSolver;
	std::vector<Partition*> partitions;
	std::vector<signed> trail_master;
	std::vector<signed> trail_shared;
	std::vector<signed> trail_root;
	int sendingTo;
	unsigned bufferSize;

	//dynamic 
	int pct_timeout;
	int pct_timeout1;
	int pct_timeout2;
	int pct_timeout3;
	int sec_timeout;
	bool quit;
	std::vector<unsigned> idleWorkers;
	bool solving_stage;
	bool sequential;

	std::set<int> variableOccurs;
	std::set<int> unitClauses;
	
	SharedVariables workerSharedVariables;
	bool dynamic_timeout;
	clock_t mpi_startTime;
	clock_t mpi_nodeTime;
	std::vector<unsigned> workingWorkers;
	int limit_conflicts;
	int limit_clauses;
	int limit_time;
	unsigned overflow_limit;
	unsigned control_clauses;
	double control_factor;
	unsigned dynamic_rounds;
	unsigned dynamic_limit;
	unsigned dynamic_increase;
	unsigned dynamic_total;

	bool changeToMaster;
	bool phase_transition;

	std::vector<VariableOccurrence*> changeOccurrences;
	
	SATSolver *sequentialSolver;
	SATSolver *globalSolver;

	//receive methods
	void receiveProblem(void);
	void receiveClauses(void);
	void receiveSharedVariables(void);
	void receiveTrail(void);
	bool receiveTermination(void);

	//send methods
	void sendTrail(void);
	void sendStatus(bool status);
	void sendTermination(bool status);

	//solving methods
	bool importInterpolants(void);
	Expression solveWorker(void);
	void sendSolution(void);
	void solvePartitions(void);
	void solveSolution(void);
	bool waitInterpolantSolution(bool quit);
	bool sendInterpolantSolution(bool quit);

	//dynamic methods
	bool waitsForWork(void);
	void splitHalf(unsigned next_master);
	void splitN(unsigned workers);
	void split_problem(double cutoff_timeout_sequential);
	void updateMasterIdle(unsigned next_master);
	void changeMaster(unsigned next_master, unsigned stage);
	void updateMaster(unsigned next_master);
	bool solveGlobals(void);
	void sendSharedVariablesDyn(void);
	void masterToWorker(void);
	void workerToMaster(void);
	bool findDisagreementDyn(void);

	//dynamic stages
	bool sequential_solving(void);

public:

	void setResources(int conflicts, int clauses, int time) 
	{
		limit_conflicts = conflicts;
		limit_clauses = clauses;
		limit_time = time;
	}

	unsigned getMaster(void) { return masterId; }

	
	void printMemoryStats(void)
	{

		PROCESS_MEMORY_COUNTERS pmc;
		GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) );
		std::cout << "Node " << rank << " [memory - MB] : " << (double)pmc.WorkingSetSize/1048576 << std::endl;
		std::cout << "Node " << rank << " [memory - formula] : " << "( " << ((MiniSAT_1p*)workerSolver->getSolver())->numVars() << " , " << 
			((MiniSAT_1p*)workerSolver->getSolver())->numClauses() << " , " <<  ((MiniSAT_1p*)workerSolver->getSolver())->nLearnts() << " )" << std::endl;
		std::cout << "Node " << rank << " [memory - time] : " << (clock() - mpi_startTime)/(double)CLOCKS_PER_SEC << std::endl;

	}
	
};

#endif