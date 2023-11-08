// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011
// RM Martins, 2012
/*
 * This file contains common code to master and the worker of the static topologies (mpdesat_master and mpdesat_worker)
 * WARNING: due to some changes in the code it may have some memory problems
 * IN USE BUT IT IS BEING REFACTORED TO STATIC_WORKER.H		
 */

#ifndef _MPDESAT_NODE_H_
#define _MPDESAT_NODE_H_

#include <time.h>

#include <vector>

#include <satsolver.h>
#include <interpolation_mode.h>
#include <decomposition.h>
#include <partition.h>
#include <shared_variables.h>

class MPDeSATNode : public SATSolver, public MPI_User
{
public:
	MPDeSATNode(ExpressionManager &m, unsigned rank, int master,  std::vector<int> worker);
	virtual ~MPDeSATNode(void);

	virtual bool addClause(const std::vector<signed> &literals){ throw std::exception("NYI: addClause"); };
	virtual bool addUnit(signed l) { throw std::exception("NYI: addUnit"); }
	virtual void setVariableMax(unsigned n) { throw std::exception("NYI: setVariableMax"); };
	virtual void setClauseMax(unsigned n) { throw std::exception("NYI: setClauseMax"); };

	virtual unsigned numClauses(void) const { throw std::exception("NYI: numClauses"); }
	virtual unsigned numVars(void) const { throw std::exception("NYI: numVars"); }
	virtual signed addVar(void) { throw std::exception("NYI: addVar"); }

	virtual bool solve(void) {throw std::exception ("NYI: solve"); };
	virtual bool solve(const std::vector<signed> &assumptions) { throw std::exception("NYI: solve with assumptions"); }

	virtual ModelValue get(signed l) const { throw std::exception("NYI: get ModelValue"); }  

	virtual Expression getInterpolant(const std::vector<signed> &A) { throw std::exception("NYI: getInterpolant"); }
	virtual bool addConstraint(const Expression &e) { throw std::exception("NYI: addConstraint @ master"); }
	virtual signed addExtension(const Expression &e) { throw std::exception("NYI: addExtension"); }

	virtual Expression getModel(void) const { throw std::exception("NYI: getModel"); }
	virtual void setVerbose(int v) { mpi_verbose = v; }

public:

	//time statistics
	clock_t mpi_loadTime;
	clock_t mpi_solvingTime;
	clock_t mpi_idleTime;
	clock_t mpi_totalTime;

	//search statistics 
	unsigned interpolants_imported;
	unsigned rounds;
	unsigned solutions_imported;
	unsigned all_sat_iterations;
	unsigned clauses_interpolants;
	unsigned variables_interpolants;
	unsigned print_rounds;
	unsigned maxSizeInterpolant;
	unsigned sizeInterpolant;
	unsigned interpolants_exported;
	unsigned solutions_exported;

	//memory warning 
	double memory_used;
	bool memory_flag;

	//interpolant statistics
	unsigned interpolant_unit_clause;
	unsigned interpolant_binary_clause;
	unsigned interpolant_ternary_clause;
	unsigned interpolant_other_clause;
	unsigned interpolant_clauses;
	unsigned interpolant_not_clause;
	unsigned interpolant_size_clauses;


protected:
	unsigned rank;

	int mpi_verbose;
	int mpi_stats;
	int verbose;

	unsigned masterId;
	std::vector<unsigned> workersId;

	unsigned maxVar;
	unsigned maxClauses;

	unsigned n_workers;
	unsigned n_cores;

	SharedVariables sharedVariables;  
	SharedVariables masterSharedVariables;

	std::vector<signed> assumptions;
	std::vector<Expression> interpolants;
	std::vector<signed> trail;
	std::vector<bool> reasons;  
	std::vector< std::vector<signed> > models;

	bool have_error, have_bad_alloc;
	std::exception exception;
	std::bad_alloc ba_exception;

	std::vector<int> recvBuffer;
	std::vector<int> sendBuffer;

	std::vector<double> speedup_time;

	//stats methods
	void printClause(const std::vector<signed> &literals);
	void saveMemoryUsed(double mem);
	void showTime(clock_t totalTime, const char *name, clock_t time);
	void statsSharedVariables(void);
	void interpolantStats(Expression t);

	//send methods
	void sendSharedVariables(void);

	//receive methods
	Expression receiveInterpolant(int * buf);

	//solving methods
	bool findDisagreement(void);

public: 

	void setStats(int v) { mpi_stats = v; }
	void printStats(void);

};

#endif