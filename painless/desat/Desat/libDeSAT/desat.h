// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _DESAT_H_
#define _DESAT_H_

#include <time.h>
#include <omp.h>

#include <vector>
#include <unordered_set>

#include "satsolver.h"
#include "decomposition.h"
#include "translation.h"
#include "interpolation_mode.h"
#include "decomposition_mode.h"
#include "interpolator.h"
#include "partition.h"

namespace Desat
{
  class DeSAT : public SATSolver
  {
  public:
    DeSAT(ExpressionManager &m, unsigned partitions, DecompositionMode decomposition, unsigned cores, std::string filename = "");
    virtual ~DeSAT(void);

    clock_t globalTime;
    clock_t partitionsTime;
    clock_t importTime;
    clock_t lastIterationTime;

    virtual bool addClause(const std::vector<signed> &literals);
    virtual bool addClause(const std::vector<signed> &literals, int w);
    virtual bool addUnit(signed l);
    virtual void setVariableMax(unsigned n);
    virtual void setClauseMax(unsigned n);
    virtual void setCallbackExportRoot(void (* exportFunction)(void *, int, std::vector<int> &)){ throw std::runtime_error("NYI: exportFunction"); };
    virtual void setCallbackImportClauses(bool (*impt)(void *, std::vector<int> &, int &)){ throw std::runtime_error("NYI: importFunction"); };

    virtual void setIssuer(void * issuer){ throw std::runtime_error("NYI: setIssuer"); };

    virtual void setConstraintMode(ConstraintMode constr);

    inline void interrupt() { early_stop = true; }
    inline void clearInterrupt() { early_stop = false; }

    virtual signed addVar(void);
    virtual unsigned numClauses(void) const;
    virtual unsigned numVars(void);

    virtual bool solve(void);
    virtual bool solve(const std::vector<signed> &assumptions);

    inline virtual ModelValue get(signed l);

    virtual Expression getInterpolant(const std::vector<signed> &A);
    virtual bool addConstraint(CExpression &e);
    virtual signed addExtension(CExpression &e);

    virtual int computeLBD(const std::vector<signed> &cls){ throw std::runtime_error("NYI: computeLBD"); };

    void printInterpolant(int partition, std::ofstream &op);

    virtual Expression getModel(void);
    virtual std::vector<int> getModelVector(void);
    virtual std::vector<int> getFinalModel(void);

    virtual void clearNewClauses(void);
    virtual void addNewClauses(void);
    virtual void getConflict(std::vector<signed> &out);
    virtual std::vector<std::vector<signed>> getNewClauses(void);

    inline SATSolver *getGlobalSolver() { return globalSolver; };

    virtual void setPhase(const int var, const bool phase);


    virtual void setVerbose(int v);

    void setInterpolator(InterpolationMode i);
    void setDecompositor(DecompositionMode i);

    virtual void setInterrupt();
    virtual void unsetInterrupt();

    std::unordered_set<signed> getInterpolantsVariables(int i);
    void extractVariablesITP(CExpression itp, std::unordered_set<signed> &vars);

  public:
    unsigned interpolants_imported;
    unsigned rounds;
    unsigned solutions_imported;
    unsigned all_sat_found;

    InterpolationMode interpolationMode;
    DecompositionMode decompositionMode;

    SharedVariables sharedVariables;

    void showDistribution(void) const;

    bool solveGlobals(void);
    bool solvePartitions(void);
    bool solvePartition(int pid);
    bool importInterpolants(std::vector<std::vector<int>> &itp_to_export);
    bool importInterpolants(std::vector<std::vector<int>> &itp_cls_export, int i);

    bool findDisagreement(void);
    void setAssumption(void);

  protected:
    unsigned n_partitions;
    unsigned n_cores;
    unsigned maxVar;
    std::string filename;

    bool early_stop;

    Decomposition *d;
    SATSolver *globalSolver;
    std::vector<Partition *> partitions;

    void init(unsigned partitions, unsigned cs);

    std::vector<signed> assumptions;
    std::vector<Expression> interpolants;
    std::vector<signed> trail;
    std::vector<bool> reasons;
    std::vector<signed> temp;

    bool have_error, have_bad_alloc;
    std::exception exception;
    std::bad_alloc ba_exception;
  };

} // namespace Desat

#endif
