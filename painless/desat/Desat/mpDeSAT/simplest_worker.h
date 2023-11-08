// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger & RM Martins, 2012
/*
 * This file contains the latest version of our solver (mpDeSAT-memory):
 * (i)   it does not use interpolation;
 * (ii)  it uses a flat structure;
 * (iii) it can start the solving process with "w" nodes and will split the formula automatically if a cutoff is reached
 * (iv)  further testing to the memory limit of MiniSAT should be done (including the way we split the database)
 * (v) ... still a lot of work needs to be done so it is competitive with MiniSAT1.14p
 * WARNING: Do not use HALFLEARNT clause since it makes the SAT solver performs much worse
 */

#ifndef _MPDESAT_SIMPLEST_WORKER_H_
#define _MPDESAT_SIMPLEST_WORKER_H_

#include <time.h>
#include <sstream>

#include <vector>
#include <minisat1p.h>
#include <decomposition.h>

class SimplestWorker : public SATSolver
{
public:
    SimplestWorker(const char *fn, int initial_nodes, int working_nodes, int mlimit);
	virtual ~SimplestWorker(void);

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

    inline bool isMaster(void) const { return rank == master; }
	void getStats(void);

protected:
    ExpressionManager em;
    
    const char * filename;

    typedef enum { INITIALIZING=0, IDLE, SOLVING_MASTER, SOLVING_WORKER, DEAD } State;
    typedef enum { INIT_DATA=0, CLAUSE_DATA, LEARNT_DATA, HALFLEARNT_DATA, TRAIL_DATA, MODEL_DATA, CONFLICT_DATA, MEMORY_DATA, SPLIT_DATA, TERMINATE_DATA} Tag;
	typedef enum { ORIGINAL_CLAUSE=0, LEARNT_CLAUSE, HALFLEARNT_CLAUSE } ClauseType;
	int rank;
    int size;
    int master;
    int active_nodes;
	int total_nodes;
    State state;
    bool result;
	unsigned rounds;
	unsigned master_changed;
	clock_t startTime;
	int memory_limit;
	bool memory_exhausted;
	bool active_status;


    unsigned maxVar, maxClauses;
    MiniSAT_1p solver;

    Decomposition * d;

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

    MPI_Status status;
    std::vector<int> sendbuf, recvbuf;
    std::vector<int> assumptions;

    void getMessage(void);

    void safe_send(int to, int tag);
    Tag receive();

    void solveMaster(void);
    void solveWorker(void);
    void conflictAnalysis(void);
    void resolve(void);
	void splitWorker(void);
	void splitMaster(void);

	bool safe_addClause(const std::vector<signed> &literals, int type);
	bool safe_solve(const std::vector<signed> &trail);

    std::vector<int> tmp_from;
    std::vector<ModelValue> tmp_model;
    std::vector<std::vector<signed> > models;
    std::vector<std::vector<signed> > conflicts;
    unsigned numConflicts, numModels;

	std::vector<signed> tmp_clause;
	std::vector<signed> tmp_trail;
	std::vector<signed> worker_clause;

	std::vector<bool> variableOccurs;

    std::stringstream output;



};

#endif
