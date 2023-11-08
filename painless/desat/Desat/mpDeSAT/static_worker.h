// Copyright (C) 2011 Microsoft Research
// RM Martins, 2012
/*
 * This file contains the new refactored version of mpDeSAT-flat, mpDeSAT-tree and mpDeSAT-tree-flat using a state machine
 * WARNING: It is still under development and should not be used
 */

#ifndef _MPDESAT_STATIC_WORKER_H_
#define _MPDESAT_STATIC_WORKER_H_

#include <time.h>
#include <sstream>

#include <vector>
#include <minisat1p.h>
#include <decomposition.h>
#include <mpi.h>
#include <interpolator.h>
#include <partition.h>

#define MPI_DEFAULT_MASTER 0;

class StaticWorker : public SATSolver
{
public:
    StaticWorker(ExpressionManager &em, const char *fn, int initial_nodes, InterpolationMode interpolation, int masterID, std::vector<int> workersID);
	virtual ~StaticWorker(void);

	virtual bool addClause(const std::vector<signed> &literals);
    virtual bool addUnit(signed l) { throw std::exception("NYI: addUnit"); }
	virtual void setVariableMax(unsigned n);
	virtual void setClauseMax(unsigned n);

	virtual unsigned numClauses(void) const;
	virtual unsigned numVars(void) const;
	virtual signed addVar(void);

	virtual bool solve(void);
	virtual bool solve(const std::vector<signed> &assumptions) { throw std::exception("NYI: solve with assumptions"); }

	virtual ModelValue get(signed l) const { throw std::exception("NYI: get ModelValue"); }  

	virtual Expression getInterpolant(const std::vector<signed> &A) { throw std::exception("NYI: getInterpolant"); }
	virtual bool addConstraint(const Expression &e) { throw std::exception("NYI: addConstraint"); }
	virtual signed addExtension(const Expression &e) { throw std::exception("NYI: addExtension"); }

	virtual Expression getModel(void) const { throw std::exception("NYI: getModel"); }    

    inline bool isMaster(void) const { return rank == MPI_DEFAULT_MASTER; }
    inline bool isLeaf(void) const { return workers.empty(); }

protected:

	inline unsigned position( std::vector<int> vec, int value) const
	{
		int pos_worker = -1;
		int from_worker = value;
		std::vector<int>::iterator it = std::find(vec.begin(), vec.end(), from_worker);
		pos_worker = std::distance(vec.begin(), it);
		assert(pos_worker!= -1);

		return (unsigned)pos_worker;
	
	}

	inline void say(int lvl, const char* format, ...) const
	{
		if (verbosity >= lvl) 
		{
			printf("[%d] ", rank);
			va_list args;
			va_start(args, format);
			vprintf(format, args);
			va_end(args);
			printf("\n");
			fflush(stdout);
		}
	}

	inline void saybuffer(int lvl, std::vector<int> buf, const char* format, ...){
	
		if (verbosity >= lvl)
		{
			printf("[%d] ", rank);
			va_list args;
			va_start(args, format);
			vprintf(format, args);
			va_end(args);

			for (unsigned i = 0; i < buf.size(); i++)
				printf("%d ",buf[i]);
			printf("\n");
			fflush(stdout);
		}
	
	}

    const char * filename;

    typedef enum { INITIALIZING=0, IDLE, SHARED_VARIABLES, SOLVING_MASTER, SOLVING_WORKER, DEAD, DEAD_LEAF } State;
    typedef enum { INIT_DATA=0, CLAUSE_DATA, TRAIL_DATA, MODEL_DATA, CONFLICT_DATA, SHARED_VARIABLES_DATA, TERMINATE_DATA } Tag;
    int rank;
    int size;
    int active_nodes;
    State state;
    bool result;
	
	unsigned overflow_limit;
	unsigned interpolants_exported;

	InterpolationMode imode;

	int master;
	std::vector<int> workers;

	unsigned maxVar, maxClauses;

	SATSolver * globalSolver;
	Partition * solver;

    Decomposition * d;

    SharedVariables sharedVariables; 

	Expression workerInterpolant;

    MPI_Status status;
    std::vector<int> sendbuf, recvbuf;
    std::vector<int> assumptions;

    unsigned numModels, numConflicts, numWorkers;

    void getMessage(void);

    void safe_send(int to, int tag);
    Tag receive();

    void solveMaster(void);
	void resolveMaster(void);

	bool findDisagreement(void);

	void resolveWorker(void);
    void solveWorker(void);
    void conflictAnalysis(void);
    void resolve(void);
	void buildSharedVariables(void);
	Expression fromBufferToITP(void);
	void importSharedVariables(void);


    std::vector<ModelValue> tmp_model;
    std::vector<std::vector<signed> > models;
    //std::vector<std::vector<signed> > conflicts;
	std::vector<Expression> interpolants;

	

    std::vector<signed> trail;    
    std::vector<bool> reasons;
    std::vector<unsigned> level;

	std::vector<signed> tmp_clause;

    std::stringstream output;

	std::vector<VariableOccurrence*> occurrences;

};

#endif
