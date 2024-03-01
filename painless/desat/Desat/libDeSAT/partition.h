// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _PARTITION_H_
#define _PARTITION_H_

#include <vector>

#include <expression.h>
#include <ExpressionManager.h>

#include "cadical.h"
#include "satsolver.h"
#include "interpolation_mode.h"
#include "interpolator.h"

namespace Desat
{
  class Partition : public SATSolver
  {
  public:
    Partition(ExpressionManager &em, SharedVariables &sharedVariables, unsigned id, unsigned verbosity = 0, bool proof = true);
    ~Partition(void);

    void setInterpolationMode(InterpolationMode interpolationMode);
    void setInterpolationMode(InterpolationMode interpolationMode, SharedVariables &masterSharedVariables);
    void setInterpolationMode(ExpressionManager &em, InterpolationMode interpolationMode, SharedVariables &masterSharedVariables);

    Expression getInterpolant(const std::vector<signed> &assumptions);

    virtual bool addClause(const std::vector<signed> &literals);
    virtual bool addUnit(signed l) { throw std::runtime_error("NYI: addUnit"); };
    virtual void setVariableMax(unsigned n)
    {
      v.reserve(n + 1);
      solver->setVariableMax(n);
    };
    virtual void setClauseMax(unsigned n) { solver->setClauseMax(n); };

    virtual void setCallbackExportRoot(void (*expt)(void *, int, std::vector<int> &)) { throw std::runtime_error("NYI: exportFunction"); };
    virtual void setCallbackImportClauses(bool (*impt)(void *, std::vector<int> &, int &)){ throw std::runtime_error("NYI: importFunction"); };

    virtual void setIssuer(void * issuer){ throw std::runtime_error("NYI: setIssuer"); };

    virtual unsigned numClauses(void) const { return solver->numClauses(); };
    virtual unsigned numVars(void) { return solver->numVars(); };
    virtual signed addVar(void) { return solver->addVar(); };

    virtual bool solve(void) { return solver->solve(); }
    virtual bool solve(const std::vector<signed> &assumptions) { return solver->solve(assumptions); }

    virtual ModelValue get(signed l) { return solver->get(l); }

    virtual void setConstraintMode(ConstraintMode constr){ solver->setConstraintMode(constr); }

    virtual bool addConstraint(CExpression &e) { return solver->addConstraint(e); }
    virtual signed addExtension(CExpression &e) { throw std::runtime_error("NYI: addExtension"); }

    virtual int computeLBD(const std::vector<signed> &cls) { throw std::runtime_error("NYI: computeLBD");}

    virtual Expression getModel(void) { throw std::runtime_error("NYI: getModel"); }
    virtual std::vector<int> getModelVector(void) { throw std::runtime_error("NYI: getModelVector"); }
    virtual std::vector<int> getFinalModel(void) { throw std::runtime_error("NYI: getFinalModel"); }

    virtual void clearNewClauses(void) { throw std::runtime_error("NYI: clearNewClauses"); };
    virtual void addNewClauses(void) { throw std::runtime_error("NYI: addNewClauses"); };
    virtual void clearAssumed(void) { throw std::runtime_error("NYI: clearAssumed"); };
    virtual std::vector<std::vector<signed>> getNewClauses(void) { throw std::runtime_error("NYI: getNewClauses"); };
    virtual void getConflict(std::vector<signed> &out) { throw std::runtime_error("NYI: getConflict"); };

    virtual void setInterrupt();
    virtual void unsetInterrupt();

    virtual void setPhase(const int var, const bool phase);

    const VariableOccurrence &variablesOccurrence(void) const { return v; }

    ExpressionManager &em(void) { return m; }

    SATSolver *getSolver(void) { return solver; }

    void insertVariableOccurrence(unsigned t) { v.setOccurs(t); }

    bool occurs(unsigned x) { return v.occurs(x); }

    void clearVariableOccurrence(void) { v.clear(); }

    void newSolver(void)
    {
      delete solver;
      solver = new Cadical(m, true);
    }

  protected:
    unsigned id;
    ExpressionManager &m;
    SATSolver *solver;
    Interpolator *interpolator;
    VariableOccurrence v;
    SharedVariables &sharedVariables;

    inline void print(const char *format, ...) const
    {
      if (verbosity <= 0)
        return;
      va_list args;
      va_start(args, format);
      vprintf(format, args);
      va_end(args);
    }
  };

} // namespace Desat

#endif
