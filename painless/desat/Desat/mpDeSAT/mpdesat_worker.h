// Copyright (C) 2011 Microsoft Research
// RM Martins, 2012
/*
 * This file contains the worker from the static topologies: flat, tree and tree-flat
 * WARNING: due to some changes in the code it may have some memory problems
 * IN USE BUT IT IS BEING REFACTORED TO STATIC_WORKER.H		
 */

#ifndef _MPDESAT_WORKER_H_
#define _MPDESAT_WORKER_H_

#include <time.h>

#include <vector>

#include  "mpdesat_node.h"
//#include <satsolver.h>
#include <interpolation_mode.h>
#include <decomposition.h>
#include <partition.h>
#include <shared_variables.h>

class MPDeSATWorker : public MPDeSATNode
{
public:
	MPDeSATWorker(ExpressionManager &m, int size, int rank, const char *filename, int c, InterpolationMode i, int master, std::vector<int> worker);
	virtual ~MPDeSATWorker(void);

	virtual bool addClause(const std::vector<signed> &literals) { throw std::exception("NYI: addClause"); }
	virtual bool addUnit(signed l) { throw std::exception("NYI: addUnit"); }
	virtual void setVariableMax(unsigned n) { throw std::exception("NYI: setVariableMax"); }
	virtual void setClauseMax(unsigned n) { throw std::exception("NYI: setClauseMax"); }

	virtual unsigned numClauses(void) const { throw std::exception("NYI: numClauses"); }
	virtual unsigned numVars(void) const { throw std::exception("NYI: numVars"); }
	virtual signed addVar(void) { throw std::exception("NYI: addVar"); }

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
	VariableOccurrence* workerOcurrences;
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

	//receive methods
	void receiveProblem(MPI_Status &s);
	void receiveClauses(MPI_Status &s);
	void receiveSharedVariables(MPI_Status &s);
	void receiveTrail(void);
	bool receiveTermination();

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

};

#endif