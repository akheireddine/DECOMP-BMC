// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _MINISAT1P_H_
#define _MINISAT1P_H_

#include <vector>
#include <map>

#include "../MiniSat/Solver.h"
#include "../MiniSat/Proof.h"

#include "satsolver.h"
#include "interpolator.h"

#include "expression.h"
#include "ExpressionManager.h"

namespace Desat
{

  class MiniSAT_1p : public SATSolver, public Solver
  {
  public:
    MiniSAT_1p(ExpressionManager &m, int nb_clauses=0, bool proof = true, ModelMode mm = NORMAL);
    MiniSAT_1p(ExpressionManager &m, SharedVariables &sharedVariables, int nb_clauses=0, bool proof = true, ModelMode mm = NORMAL);

    virtual ~MiniSAT_1p(void);

    virtual bool addClause(const std::vector<signed> &literals);
    virtual bool addUnit(signed l);
    virtual void setVariableMax(unsigned n);
    virtual void setClauseMax(unsigned n);
    virtual void setConstraintMode(ConstraintMode constr);

    virtual void setCallbackExportRoot(void (*expt)(void *, int, std::vector<int> &));
    virtual void setCallbackImportClauses(bool (*impt)(void *, std::vector<int> &, int &));

    virtual void setIssuer(void *issuer);

    virtual signed addVar(void);
    virtual unsigned numClauses(void) const;
    virtual unsigned numVars(void) const;

    virtual bool solve(void);
    virtual bool solve(const std::vector<signed> &assumptions);

    inline virtual ModelValue get(signed l) const;
    virtual Expression getInterpolant(const std::vector<signed> &beta);
    virtual Expression getModel(void) const;
    virtual std::vector<int> getModelVector(void) const;
    virtual std::vector<int> getFinalModel(void) const;

    virtual bool addConstraint(CExpression &e);
    virtual signed addExtension(CExpression &e);

    int computeLBD(const std::vector<signed> &cls){ return Solver::computeLBDInt(cls); }

    virtual void getStats(unsigned long &confl, unsigned long &propag, unsigned long &restart, unsigned long &decision);

    virtual void setVerbose(int v);

    virtual void setInterrupt();
    virtual void unsetInterrupt();

    virtual void setPhase(const int var, const bool phase);

    void getConflict(std::vector<signed> &out)
    {
      out.clear();
      for (int i = 0; i < conflict.size(); i++)
      {
        const Lit &l = conflict[i];
        signed x = sign(l) ? -var(l) : var(l);
        out.push_back(x);
      }
    }

    bool addLearntClause(const std::vector<signed> &literals);
    bool addHalfLearntClause(const std::vector<signed> &literals);

    int n_clauses;
  private:
    vec<Lit> temp_lits;

  protected:
    vec<Lit> assumptions;
    std::map<unsigned, Expression> interpolants;
    std::vector<unsigned> interpolate_stack;
    ConstraintMode constraintMode;

    VariableOccurrence v;
    SharedVariables &sharedVariables;
    std::vector<bool> local_variables;

    class InterpolationTraverser : public ProofTraverser
    {
    public:
      InterpolationTraverser(void);
      virtual void root(const vec<Lit> &c);
      virtual void chain(const vec<ClauseId> &cs, const vec<Lit> &xs);
      virtual void deleted(ClauseId c);
      virtual void done(ClauseId e);
      virtual ~InterpolationTraverser();

      class ResolutionChainInfo
      {
      public:
        vec<ClauseId> cs;
        std::vector<signed> xs;
#ifdef _DEBUG
        vec<Lit> resolvent;
#endif
      };
      vec<ResolutionChainInfo> resChainInfo;

    protected:
      void resolve(vec<Lit> &main, vec<Lit> &other, Lit x) const;
    };

    InterpolationTraverser iTraverser;

    Expression interpolate(void);
    Expression interpolate(int cid, const vec<Lit> &assumptions);

    typedef std::map<CExpression, Lit> ClausifyCache;
    typedef std::map<CExpression, vec<vec<Lit>>> ExplosionCache;
    ClausifyCache clausifyCache;
    ExplosionCache explosionCache;

    // FIXME these 3 vectors were marked CExpresssion but this failed
    std::vector<Expression> clausifyCNFStack;
    std::vector<Expression> clausifyExtendStack;
    std::vector<Expression> clausifyExplodeStack;

    std::vector<bool> assumed;
    std::vector<signed> tmp_unit;

    Lit clausify_clause(CExpression &e);
    Lit clausify_cube(CExpression &e);
    Lit clausify_cnf(CExpression &e);
    Lit extend(CExpression &e);
    Lit explode(CExpression &e);

    Lit pre_extend(CExpression &e);
    bool checkCnf(CExpression &e);
    bool checkConj(CExpression &e);
    void directCNF(CExpression &e);
    void getClause(CExpression &e, vec<Lit> &clause);

    virtual bool resolve_conflict(Clause *confl) { return Solver::resolve_conflict(confl); }

    bool addClause(const vec<Lit> &literals, int type)
    {

      switch (type)
      {
      case ORIGINAL_CLAUSE: //ORIGINAL
        Solver::addClause(literals);
        break;
      case LEARNT_CLAUSE: //LEARNT
        Solver::addLearntClause(literals);
        break;
      case HALFLEARNT_CLAUSE: //HALFLEARNT
        Solver::addHalfLearntClause(literals);
        break;
      default:
        throw std::runtime_error("unknown clause type");
      }

      return Solver::okay();
    }

    bool addUnit(Lit l)
    {
      Solver::addUnit(l);
      return Solver::okay();
    }

    void getLits(CExpression &e, vec<Lit> &lits);

    vec<vec<Lit>> new_clauses;

  public:
    void removeLastToClause(std::vector<signed> &clause);

    void clearNewClauses(void) { new_clauses.clear(); }
    void addNewClauses(void);
    void clearAssumed(void) { assumed.clear(); }
    void addNewClauses(vec<vec<Lit>> interpolant_clauses);
    std::vector<std::vector<signed>> getNewClauses(void);
    //void setTimeout(int seconds) { Solver::setTimeout(seconds); }
    void splitDB(std::vector<signed> &sendBuffer);
    void removeLast(std::vector<signed> &sendBuffer);
    void resetProof(void)
    {
      //delete proof;
      if (proof != NULL)
        delete proof;
      InterpolationTraverser iTraverser;
      Solver::proof = new Proof(iTraverser);
    }
    void deleteProof(void)
    {
      if (proof != NULL)
        delete proof;
      Solver::proof = NULL;
    }
    //void getUsedVariables(std::set<int> &vars);
  };

}

#endif
