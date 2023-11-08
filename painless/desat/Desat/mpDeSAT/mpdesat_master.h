// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011
// RM Martins, 2012
/*
 * This file contains the master from the static topologies: flat, tree and tree-flat
 * WARNING: due to some changes in the code it may have some memory problems
 * IN USE BUT IT IS BEING REFACTORED TO STATIC_WORKER.H		
 */

#ifndef _MPDESAT_MASTER_H_
#define _MPDESAT_MASTER_H_

#include <time.h>

#include <vector>

#include  "mpdesat_node.h"
#include <interpolation_mode.h>
#include <decomposition.h>
#include <partition.h>
#include <shared_variables.h>

class MPDeSATMaster : public MPDeSATNode
{
public:
	MPDeSATMaster(ExpressionManager &m, int size, int rank, const char *filename, int c, InterpolationMode i, int master, std::vector<int> worker);
	virtual ~MPDeSATMaster(void);

	virtual bool addClause(const std::vector<signed> &literals);
	virtual bool addUnit(signed l) { throw std::exception("NYI: addUnit"); }
	virtual void setVariableMax(unsigned n);
	virtual void setClauseMax(unsigned n);

	virtual unsigned numClauses(void) const { throw std::exception("NYI: numClauses"); }
	virtual unsigned numVars(void) const { throw std::exception("NYI: numVars"); }
	virtual signed addVar(void) { throw std::exception("NYI: addVar"); }

	virtual bool solve(void);
	virtual bool solve(const std::vector<signed> &assumptions) { throw std::exception("NYI: solve with assumptions"); }

	virtual ModelValue get(signed l) const { throw std::exception("NYI: get ModelValue"); }  

	virtual Expression getInterpolant(const std::vector<signed> &A) { throw std::exception("NYI: getInterpolant"); }
	virtual bool addConstraint(const Expression &e) { throw std::exception("NYI: addConstraint"); }
	virtual signed addExtension(const Expression &e) { throw std::exception("NYI: addExtension"); }

	virtual Expression getModel(void) const { throw std::exception("NYI: getModel"); }
	virtual void setVerbose(int v) { mpi_verbose = v; }


protected:

	int size;
	int cores;
	const char *filename;
	unsigned n_workers;
	Decomposition *d;
	SATSolver *globalSolver;
	std::vector<Partition*> partitions;
	int sendingTo;
	unsigned bufferSize;
	unsigned num_clauses;
	std::vector<VariableOccurrence*> occurrences;

	//initialization methods
	void stopReadDimacsFile(void);

	//send methods
	void sendTrail(void);
	void sendStatus(bool status);
	void sendSharedVariables();

	//solving methods
	bool importInterpolants(void);
	bool solveGlobals(void);

};

#endif