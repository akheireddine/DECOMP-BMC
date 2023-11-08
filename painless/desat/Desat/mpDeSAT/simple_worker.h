// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger & RM Martins, 2012
/*
 * This file contains the simple worker:
 * (i) it does not use interpolation and it was the start of the new framework
 * WARNING: Deprecated file! Should be used only if we want to extend it from scratch into a different direction.
 */

#ifndef _MPDESAT_SIMPLE_WORKER_H_
#define _MPDESAT_SIMPLE_WORKER_H_

#include <time.h>
#include <sstream>

#include <vector>
#include <minisat1p.h>
#include <decomposition.h>

class SimpleWorker : public SATSolver
{
public:
    SimpleWorker(const char *fn, int initial_nodes);
	virtual ~SimpleWorker(void);

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

protected:
    ExpressionManager em;
    
    const char * filename;

    typedef enum { INITIALIZING=0, IDLE, SOLVING_MASTER, SOLVING_WORKER, DEAD } State;
    typedef enum { INIT_DATA=0, CLAUSE_DATA, TRAIL_DATA, MODEL_DATA, CONFLICT_DATA, TERMINATE_DATA } Tag;
    int rank;
    int size;
    int master;
    int active_nodes;
    State state;
    bool result;


    unsigned maxVar, maxClauses;
    MiniSAT_1p solver;

    Decomposition * d;

    inline void say(const char* format, ...) const
    {
      if (verbosity<=0) return;
      printf("[%d] ", rank);
      va_list args;
      va_start(args, format);
      vprintf(format, args);
      va_end(args);
      printf("\n");
      fflush(stdout);
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

    std::vector<ModelValue> tmp_model;
    std::vector<std::vector<signed> > models;
    std::vector<std::vector<signed> > conflicts;

    std::vector<signed> trail;    
    std::vector<std::vector<signed> > reasons;
    std::vector<unsigned> level;

    std::stringstream output;
};

#endif
