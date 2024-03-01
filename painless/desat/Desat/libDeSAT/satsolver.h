// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _SATSOLVER_H_
#define _SATSOLVER_H_

#include <stdarg.h>

#include <cassert>

#include <expression.h>

#include "dimacs_parser.h"
#include "interpolator.h"
#include "model.h"
#include "conversion_mode.h"

namespace Desat {

typedef enum { ORIGINAL_CLAUSE = 0, LEARNT_CLAUSE, HALFLEARNT_CLAUSE } ClauseType;

class SATSolver : public DimacsParser
{
public:  
  SATSolver(ExpressionManager &m) : m(m), verbosity(0), interpolator(NULL) {}
  virtual ~SATSolver(void) {}

  virtual bool addClause(const std::vector<signed> &literals) = 0;
  virtual void setVariableMax(unsigned n) = 0;
  virtual void setClauseMax(unsigned n) = 0;

  virtual void setCallbackExportRoot(void (*expt)(void *, int, std::vector<int> &)) = 0;
  virtual void setCallbackImportClauses(bool (*impt)(void *, std::vector<int> &, int &)) = 0;

  virtual void setIssuer(void * issuer) = 0;

  virtual unsigned numClauses(void) const = 0;
  virtual unsigned numVars(void) = 0;
  virtual signed addVar(void) = 0;

  virtual void clearNewClauses(void) = 0;
  virtual void addNewClauses(void) = 0;
  virtual std::vector<std::vector<signed> > getNewClauses(void) = 0;

  virtual bool solve(void) = 0;
  virtual bool solve(const std::vector<signed> &assumptions) = 0;

  virtual ModelValue get(signed l) = 0;  

  virtual void setInterpolator(const Interpolator *i) { interpolator = i; }
  virtual Expression getInterpolant(const std::vector<signed> &A) = 0;
  virtual bool addConstraint(CExpression &e) = 0;
  virtual signed addExtension(CExpression &e) = 0;

  virtual int computeLBD(const std::vector<signed> &cls) = 0;

  virtual Expression getModel(void) = 0;
  virtual std::vector<int> getModelVector(void) = 0;
  virtual std::vector<int> getFinalModel(void) = 0;

  virtual void setInterrupt() = 0;
  virtual void unsetInterrupt() = 0;

  virtual void getConflict(std::vector<signed> &out) = 0;

  virtual void setVerbose(int v) { verbosity=v; }
  virtual void setConstraintMode(ConstraintMode constr)=0;

  virtual void getStats(unsigned long &confl, unsigned long &propag, unsigned long &restart, unsigned long &decision){};

  virtual void setPhase(const int var, const bool phase)=0;

  virtual int returnVerbose(void) { return verbosity; }

  SATSolver(const SATSolver &other) : m(other.m) { assert(false); }


protected:
  int verbosity;
  ExpressionManager &m;
  const Interpolator *interpolator;  

  void print(const char* format, ...) const;
};


// Inlines

inline void SATSolver::print(const char* format, ...) const
{
  if (verbosity<=0) return;
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}



}

#endif