// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC_NEW
#include <stdlib.h>
#include <crtdbg.h>
#endif

#include <iostream>

#include "../MiniSat/Solver.h"
#include "../MiniSat/Sort.h"

#include "minisat1p.h"
#include "interpolator.h"

using namespace Desat;

#define MINI_LIT(lit) lit > 0 ? Lit(lit, false) : Lit((-lit), true)

#define INT_LIT(lit) sign(lit) ? -(var(lit)) : (var(lit))

SharedVariables s_minisat;

MiniSAT_1p::MiniSAT_1p(ExpressionManager &m, int nb_clauses, bool proof, ModelMode mm) : SATSolver(m),
                                                                                         sharedVariables(s_minisat),
                                                                                         iTraverser(),
                                                                                         tmp_unit(1, 0),
                                                                                         n_clauses(nb_clauses),
                                                                                         // constraintMode(PLAISTED_GREENBAUM_EXTENSION)
                                                                                         constraintMode(EXPLOSION)
// constraintMode(TSEITIN_EXTENSION)
{
  if (proof)
    Solver::proof = new Proof(iTraverser);
  Solver::newVar(); // Don't use variable 0
  if (mm == LIFTING)
    Solver::model_lifting = true;
  else if (mm == CHECKING)
    Solver::partial_models = true;

  sharedVariables.add(&v);
}

MiniSAT_1p::MiniSAT_1p(ExpressionManager &m, SharedVariables &sharedVariables, int nb_clauses, bool proof, ModelMode mm) : SATSolver(m),
                                                                                                                           sharedVariables(sharedVariables),
                                                                                                                           iTraverser(),
                                                                                                                           tmp_unit(1, 0),
                                                                                                                           n_clauses(nb_clauses),
                                                                                                                           // constraintMode(PLAISTED_GREENBAUM_EXTENSION)
                                                                                                                           constraintMode(EXPLOSION)
{
  if (proof)
    Solver::proof = new Proof(iTraverser);
  Solver::newVar(); // Don't use variable 0
  if (mm == LIFTING)
    Solver::model_lifting = true;
  else if (mm == CHECKING)
    Solver::partial_models = true;

  sharedVariables.add(&v);
}

MiniSAT_1p::~MiniSAT_1p(void)
{
  if (proof)
    delete proof;
}

void MiniSAT_1p::setVerbose(int v)
{
  Solver::verbosity = (v > 1) ? 1 : 0;
  SATSolver::setVerbose(v);
}

void MiniSAT_1p::setConstraintMode(ConstraintMode constr)
{
  constraintMode = constr;
}

void MiniSAT_1p::setInterrupt()
{
  Solver::interrupt();
}

void MiniSAT_1p::unsetInterrupt()
{
  Solver::clearInterrupt();
}

void MiniSAT_1p::setCallbackExportRoot(void (*expt)(void *, int, std::vector<int> &))
{
  Solver::exportClauseCallback = expt;
}

void MiniSAT_1p::setCallbackImportClauses(bool (*impt)(void *, std::vector<int> &, int &))
{
  Solver::importClauseCallback = impt;
}

void MiniSAT_1p::setIssuer(void *issuer)
{
  Solver::issuer = issuer;
}

bool MiniSAT_1p::addLearntClause(const std::vector<signed> &literals)
{
  if (!Solver::okay())
    return false;

  temp_lits.clear();

  for (unsigned i = 0; i < literals.size(); i++)
  {
    const signed &l = literals[i];
    bool sgn = (l < 0);
    signed v = (sgn) ? -l : l;
    assert(v < Solver::nVars());
    temp_lits.push(Lit(v, sgn));
  }

  return addClause(temp_lits, LEARNT_CLAUSE);
}

bool MiniSAT_1p::addHalfLearntClause(const std::vector<signed> &literals)
{
  if (!Solver::okay())
    return false;

  temp_lits.clear();

  for (unsigned i = 0; i < literals.size(); i++)
  {
    const signed &l = literals[i];
    bool sgn = (l < 0);
    signed v = (sgn) ? -l : l;
    assert(v < Solver::nVars());
    temp_lits.push(Lit(v, sgn));
  }

  return addClause(temp_lits, HALFLEARNT_CLAUSE);
}

bool MiniSAT_1p::addClause(const std::vector<signed> &literals)
{
  if (!Solver::okay())
    return false;

  temp_lits.clear();

  for (unsigned i = 0; i < literals.size(); i++)
  {
    const signed &l = literals[i];
    bool sgn = (l < 0);
    signed var = (sgn) ? -l : l;
    v.setOccurs(var);
    assert(var < Solver::nVars());
    temp_lits.push(Lit(var, sgn));
  }

  return addClause(temp_lits, ORIGINAL_CLAUSE);
}

bool MiniSAT_1p::addUnit(signed l)
{
  if (!Solver::okay())
    return false;

  print("New Unit: (%d)\n", l);
  bool sgn = (l < 0);
  signed v = (sgn) ? -l : l;
  assert(v < Solver::nVars());
  Lit ml(v, sgn);
  return addUnit(ml);
}

void MiniSAT_1p::setVariableMax(unsigned n)
{
  while (n >= (unsigned)Solver::nVars())
    Solver::newVar();
}

void MiniSAT_1p::setClauseMax(unsigned n)
{
}

signed MiniSAT_1p::addVar(void)
{
  Solver::newVar();
  return Solver::nVars() - 1;
}

unsigned MiniSAT_1p::numClauses(void) const
{
  return Solver::nClauses();
}

unsigned MiniSAT_1p::numVars(void) const
{
  return Solver::nVars();
}

bool MiniSAT_1p::solve(void)
{
  if (!Solver::okay())
    return false;

  bool r = Solver::solve();

  return r;
}

bool MiniSAT_1p::solve(const std::vector<signed> &assumptions)
{
  if (!Solver::okay())
  {
    //std::cout << "TRIVIAL UNSATISFIABILITY" << std::endl;
    return false;
  }

  Solver::cancelUntil(0);

  vec<Lit> a;
  for (unsigned i = 0; i < assumptions.size(); i++)
  {
    const signed &l = assumptions[i];
    bool sgn = (l < 0);
    Lit x(sgn ? -l : l, sgn);
    a.push(x);
  }

  Solver::simplifyDB(); //simplifyDB is already included in the solving method of MiniSAT
  bool r = Solver::solve(a);

  if (!r)
  {
    //std::cout << "Unsatisfiability depends on:";
    //for (int i=0; i<Solver::conflict.size(); i++)
    //  std::cout << " " << (sign(Solver::conflict[i]) ? "-" : "") << var(Solver::conflict[i]);
    //std::cout << std::endl;
  }

  return r;
}

void MiniSAT_1p::getStats(unsigned long &confl, unsigned long &propag, unsigned long &restart, unsigned long &decision){
  confl = Solver::stats.conflicts;
  propag = Solver::stats.propagations;
  restart = Solver::stats.starts;
  decision = Solver::stats.decisions;
}

void MiniSAT_1p::setPhase(const int var, const bool phase){
  Solver::setPolarity(var, phase);
}

Expression MiniSAT_1p::getInterpolant(const std::vector<signed> &beta)
{
  // (*this) is A
  // We want an interpolant I, i.e.,
  // A -> I and -(I and beta)
  // and I is over the shared variables.

  // TODO: Save the search state?

  if (!Solver::okay())
  {
    return interpolate();
  }

  //std::cout << "Assumptions: ";
  assumptions.clear();
  for (unsigned i = 0; i < beta.size(); i++)
  {
    const signed &l = beta[i];
    bool sgn = (l < 0);
    int v = (sgn) ? -l : l;
    assumptions.push(Lit(v, sgn));
    //if (sgn) std::cout << "-";
    //std::cout << v ;
  }
  //std::cout << std::endl;

  Solver::simplifyDB();
  bool r = Solver::solve(assumptions);

  if (r)
  {
    // The assignment in beta actually works.
#ifdef _DEBUG
    for (int i = 0; i < assumptions.size(); i++)
    {
      const Lit &l = assumptions[i];
      const lbool &v = Solver::model[var(l)];
      assert(v == D_Undef ||
             (sign(l) && v == D_False) ||
             (!sign(l) && v == D_True));
    }
#endif

    return m.mkTrue();
  }
  else
  {

    if (SATSolver::verbosity > 0)
    {
      std::cout << "Unsatisfiability depends on:";
      for (int i = 0; i < Solver::conflict.size(); i++)
        std::cout << " " << (sign(Solver::conflict[i]) ? "-" : "") << var(Solver::conflict[i]);
      std::cout << std::endl;
    }

    return interpolate();
  }

  // TODO: Restore the search state?
}

Expression MiniSAT_1p::interpolate(void)
{
  if (Solver::conflict_id == ClauseId_NULL) // Conflict under assumption
  {
    return m.mkFalse();
  }

  return interpolate(Solver::conflict_id, assumptions);
}

Expression MiniSAT_1p::interpolate(int cid, const vec<Lit> &assumptions)
{
  iTraverser.done(cid);

  // #ifdef _DEBUG
  // vec<Lit> &c = iTraverser.resChainInfo[cid].resolvent;
  // for (int i=0; i<c.size(); i++)
  // {
  //   const Lit &l = c[i];
  //   bool ok=false;
  //   for (int j=0; j<assumptions.size(); j++)
  //     if (assumptions[j] == ~l) { ok = true; break; }
  //   assert(ok);
  // }
  // #endif

  if (assumed.size() < (unsigned)(2 * nVars() + 1))
    assumed.resize(2 * nVars() + 1, false);

  for (int i = 0; i < assumptions.size(); i++)
  {
    assumed[index(assumptions[i])] = true;
  }

  interpolants.clear();
  interpolate_stack.clear();
  interpolate_stack.push_back(cid);

  while (!interpolate_stack.empty())
  {
    unsigned i = interpolate_stack.back();

    if (interpolants.find(i) != interpolants.end())
      interpolate_stack.pop_back(); // Cached
    else
    {
      const InterpolationTraverser::ResolutionChainInfo &r = iTraverser.resChainInfo[i];
      bool children_done = true;

      for (int j = 0; j < r.cs.size(); j++)
      {
        if (interpolants.find(r.cs[j]) == interpolants.end())
        {
          children_done = false;
          interpolate_stack.push_back(r.cs[j]);
        }
      }

      if (children_done)
      {
        if (r.cs.size() == 0)
        {
          //print("ROOT (%d):", i);
          //for (unsigned j=0; j<r.xs.size(); j++)
          //  print(" %d", r.xs[j]);
          //print("\n");

          interpolants[i] = interpolator->root(r.xs, false);
          //std::cout << "ROOT B: " << m.toString(interpolants[i]) << std::endl;

          // For every assumption in the root, get rid of it!
          for (unsigned j = 0; j < r.xs.size(); j++)
          {
            const signed &x = r.xs[j];

            if (assumed[index(MINI_LIT(-x))])
            {
//std::cout << "Assumed: " << -x << std::endl;
#ifdef _DEBUG
              vec<Lit> &r = iTraverser.resChainInfo[i].resolvent;
              bool check = false;
              for (int ri = 0; ri < r.size(); ri++)
                if (r[ri] == MINI_LIT(x))
                {
                  check = true;
                  break;
                }
              assert(check);
#endif

              tmp_unit[0] = -x;
              Expression right = interpolator->root(tmp_unit, true);
              interpolants[i] = interpolator->resolve(interpolants[i], right, x);
            }
          }

          //std::cout << "ROOT(" << i << ") ITP: " << m.toString(interpolants[i]) << std::endl;
        }
        else
        {
          //print("RES from: %d\n", r.cs[0]);

          Expression t = interpolants[r.cs[0]];
          assert(!m.isNil(t));

          for (unsigned j = 0; j < r.xs.size(); j++)
          {
            //print("RES to: %d via %d\n", r.cs[j+1], r.xs[j]);
            t = interpolator->resolve(t, interpolants[r.cs[j + 1]], r.xs[j]);
          }

          //std::cout << "Partial(" << i << "," << t << "): " << m.toString(t) << std::endl;
          interpolants[i] = t;
        }

        interpolate_stack.pop_back();
      }
    }
  }

  // Clean up
  for (int i = 0; i < assumptions.size(); i++)
    assumed[index(assumptions[i])] = false;

  // std::cout << "Final(" << cid << "," << interpolants[cid] << "): " << m.toString(interpolants[cid]) << std::endl;

  return interpolants[cid];
}

bool MiniSAT_1p::addConstraint(CExpression &e)
{
  Lit top = lit_Undef;

  ClausifyCache::const_iterator it = clausifyCache.find(e);
  if (it != clausifyCache.end())
    top = it->second;
  else
  {
    if (m.isLiteral(e))
    {
      signed l = m.getLiteral(e);
      Lit q = MINI_LIT(l);
      top = m.isNegative(e) ? ~q : q;
    }
    //else if (m.isClause(e))
    //  top=clausify_clause(e);
    //else if (m.isCube(e))
    //  top=clausify_cube(e);
    //else if (m.isCNF(e))
    //  top=clausify_cnf(e);
    else
    {
      switch (constraintMode)
      {
      case TSEITIN_EXTENSION:
      case PLAISTED_GREENBAUM_EXTENSION:
        top = pre_extend(e);
        break;
      case EXPLOSION:
        top = explode(e);
        break;
      default:
        throw std::runtime_error("Unknown constraint translation mode.");
      }
    }
  }

  if (top != lit_Undef)
  {
    addUnit_tmp[0] = top;
    new_clauses.push(addUnit_tmp);
  }

  return true;
}

signed MiniSAT_1p::addExtension(CExpression &e)
{
  Lit q = extend(e);
  return (sign(q)) ? -var(q) : var(q);
}

void MiniSAT_1p::getLits(CExpression &e, vec<Lit> &lits)
{
  lits.clear();

  std::vector<signed> temp;
  m.getLiterals(e, temp);
  for (unsigned i = 0; i < temp.size(); i++)
    lits.push(MINI_LIT(temp[i]));
}

std::vector<std::vector<signed>> MiniSAT_1p::getNewClauses()
{
  std::vector<std::vector<signed>> new_clauses_literal;
  for (size_t i = 0; i < new_clauses.size(); i++)
  {
    //    vec<Lit> sorted_cls;
    std::vector<signed> cls;
    //    new_clauses[i].copyTo(sorted_cls);
    //    sortUnique(sorted_cls);

    for (int j = 0; j < new_clauses[i].size(); j++)
    {
      cls.emplace_back(INT_LIT(new_clauses[i][j]));
    }
    new_clauses_literal.emplace_back(cls);
  }
  return new_clauses_literal;
}

Lit MiniSAT_1p::pre_extend(CExpression &e)
{
  // Give it to the solver if already in CNF
  if (checkCnf(e))
  {
    directCNF(e);
    return lit_Undef;
  }
  // Otherwise perform cnfization
  else
  {
    // map< enodeid_t, int > enodeid_to_incoming_edges;
    // computeIncomingEdges( f, enodeid_to_incoming_edges ); // Compute incoming edges for f and children
    // f = rewriteMaxArity( f, enodeid_to_incoming_edges );  // Rewrite f with maximum arity for operators
    Lit res = extend(e);
    return res;
  }
}

//
// Check whether a formula is in cnf
//
bool MiniSAT_1p::checkCnf(CExpression &e)
{
  bool res = checkConj(e) || m.isClause(e);
  return res;
}

//
// Check if a formula is a conjunction of clauses
//
bool MiniSAT_1p::checkConj(CExpression &e)
{
  if (!m.isAnd(e))
    return false;

  size_t sz = m.nChildren(e);
  for (unsigned i = 0; i < sz; i++)
  {
    CExpression c = m.getChild(e, i);
    if (!checkConj(c) && !m.isClause(c))
      return false;
  }
  return true;
}

void MiniSAT_1p::getClause(CExpression &e, vec<Lit> &clause)
{
  if (m.isClause(e))
  {
    std::vector<signed> lits;
    m.getLiterals(e, lits);
    for (unsigned i = 0; i < lits.size(); i++)
    {
      Lit l = MINI_LIT(lits[i]);
      clause.push(l);
    }
  }
  else
    assert(false);
}

//
// Give the formula to the solver
//
void MiniSAT_1p::directCNF(CExpression &e)
{
  //
  // A unit clause
  //
  if (m.isLiteral(e))
  {
    Lit res = MINI_LIT(m.getLiteral(e));
  }
  //
  // A clause
  //
  else if (m.isOr(e))
  {
    vec<Lit> big_clause;
    getClause(e, big_clause);
    new_clauses.push(big_clause);
    return;
  }
  //
  // Conjunction
  //
  else if (m.isAnd(e))
  {
    size_t sz = m.nChildren(e);
    for (unsigned i = 0; i < sz; i++)
    {
      CExpression c = m.getChild(e, i);
      directCNF(c);
    }
  }
}

Lit MiniSAT_1p::extend(CExpression &e)
{
  //std::cout << "Extend" << std::endl;
  clausifyExtendStack.clear();
  clausifyExtendStack.push_back(e);

  while (!clausifyExtendStack.empty())
  {
    CExpression cur = clausifyExtendStack.back();

    if (clausifyCache.find(cur) != clausifyCache.end())
    {
      clausifyExtendStack.pop_back();
    }
    else
    {
      if (m.isTrue(cur))
      {
        Solver::newVar();
        Lit nv = Lit(nVars() - 1, false);
        // addUnit(nv);
        addUnit_tmp[0] = nv;
        new_clauses.push(addUnit_tmp);
        clausifyCache[cur] = nv;
        clausifyExtendStack.pop_back();
      }
      else if (m.isFalse(cur))
      {
        Solver::newVar();
        Lit nv = Lit(nVars() - 1, false);
        //addUnit(~nv);
        addUnit_tmp[0] = ~nv;
        new_clauses.push(addUnit_tmp);
        clausifyCache[cur] = nv;
        clausifyExtendStack.pop_back();
      }
      else if (m.isLiteral(cur))
      {
        clausifyCache[cur] = MINI_LIT(m.getLiteral(cur));
        clausifyExtendStack.pop_back();
      }
      else if (m.isAnd(cur))
      {
        bool all_children_done = true;
        size_t sz = m.nChildren(cur);
        for (unsigned i = 0; i < sz; i++)
        {
          CExpression c = m.getChild(cur, i);
          if (clausifyCache.find(c) == clausifyCache.end())
          {
            all_children_done = false;
            clausifyExtendStack.push_back(c);
          }
        }

        if (all_children_done)
        {
          Solver::newVar();
          Lit nv = Lit(nVars() - 1, false);
          vec<Lit> big(1, nv);
          vec<Lit> lits(2, ~nv);

          for (unsigned i = 0; i < sz; i++)
            big.push(~clausifyCache[m.getChild(cur, i)]);

          for (int i = 1; i < big.size(); i++)
          {
            lits[1] = ~big[i];
            new_clauses.push(lits);
          }

          if (constraintMode == TSEITIN_EXTENSION)
            new_clauses.push(big);

          clausifyCache[cur] = nv;
          clausifyExtendStack.pop_back();
        }
      }
      else if (m.isOr(cur))
      {
        bool all_children_done = true;
        size_t sz = m.nChildren(cur);

        //std::cout << "EXTENDING OR: " << m.toString(cur) << std::endl;
        for (unsigned i = 0; i < sz; i++)
        {
          CExpression c = m.getChild(cur, i);
          //std::cout << "CHILD: " << m.toString(c) << std::endl;
          if (clausifyCache.find(c) == clausifyCache.end())
          {
            all_children_done = false;
            clausifyExtendStack.push_back(c);
          }
        }

        if (all_children_done)
        {
          Solver::newVar();
          Lit nv = Lit(nVars() - 1, false);
          vec<Lit> big(1, ~nv);

          for (unsigned i = 0; i < sz; i++)
            big.push(clausifyCache[m.getChild(cur, i)]);

          if (constraintMode == TSEITIN_EXTENSION)
          {
            vec<Lit> lits(2, nv);
            for (int i = 1; i < big.size(); i++)
            {
              lits[1] = ~big[i];
              new_clauses.push(lits);
            }
          }

          new_clauses.push(big);
          clausifyCache[cur] = nv;
          clausifyExtendStack.pop_back();
        }
      }
      else if (m.isNil(cur))
        throw std::runtime_error("Cannot clausify NIL.");
      else
        throw std::runtime_error("Unexpected expression type.");
    }
  }

  assert(clausifyCache.find(e) != clausifyCache.end());
  Lit res = clausifyCache[e];
  // clausifyCache.clear();
  return res;
}

Lit MiniSAT_1p::explode(CExpression &e)
{
  clausifyExplodeStack.clear();
  clausifyExplodeStack.push_back(e);

  int limit_clauses = 100; //0.20 * n_clauses;
  std::vector<CExpression *> exp_tmp;

  while (!clausifyExplodeStack.empty())
  {
    if (explosionCache.find(e) != explosionCache.end() && explosionCache[e].size() > limit_clauses)
    {
      for (auto exp : exp_tmp)
      {
        explosionCache.erase(*exp);
      }
      Lit top = extend(e);
      return top;
    }
    CExpression cur = clausifyExplodeStack.back();

    if (explosionCache.find(cur) != explosionCache.end())
    {
      clausifyExplodeStack.pop_back();
    }
    else
    {
      if (m.isTrue(cur))
      {
        explosionCache[cur].push(vec<Lit>());
        clausifyExplodeStack.pop_back();
      }
      else if (m.isFalse(cur))
      {
        explosionCache[cur].clear();
        clausifyExplodeStack.pop_back();
      }
      else if (m.isLiteral(cur) || m.isClause(cur))
      {
        explosionCache[cur].push(vec<Lit>());
        getLits(cur, explosionCache[cur].last());
        clausifyExplodeStack.pop_back();
      }
      else if (m.isAnd(cur))
      {
        bool all_children_done = true;
        size_t sz = m.nChildren(cur);
        for (unsigned i = 0; i < sz; i++)
        {
          CExpression c = m.getChild(cur, i);
          if (explosionCache.find(c) == explosionCache.end())
          {
            all_children_done = false;
            clausifyExplodeStack.push_back(c);
            exp_tmp.emplace_back(&c);
          }
        }

        if (all_children_done)
        {
          explosionCache[cur].clear();
          for (unsigned i = 0; i < sz; i++)
          {
            CExpression c = m.getChild(cur, i);
            for (int j = 0; j < explosionCache[c].size(); j++)
              explosionCache[cur].push(explosionCache[c][j]);
          }
          clausifyExplodeStack.pop_back();
        }
      }
      else if (m.isOr(cur))
      {
        bool all_children_done = true;
        size_t sz = m.nChildren(cur);
        for (unsigned i = 0; i < sz; i++)
        {
          CExpression c = m.getChild(cur, i);
          if (explosionCache.find(c) == explosionCache.end())
          {
            all_children_done = false;
            clausifyExplodeStack.push_back(c);
            exp_tmp.emplace_back(&c);
          }
        }

        if (all_children_done)
        {
          assert(m.nChildren(cur) == 2);
          vec<vec<Lit>> &result = explosionCache[cur];
          result.clear();

          const vec<vec<Lit>> &a = explosionCache[m.getChild(cur, 0)];
          const vec<vec<Lit>> &b = explosionCache[m.getChild(cur, 1)];
          unsigned a_sz = a.size();
          unsigned b_sz = b.size();

          if (a_sz * b_sz > limit_clauses)
          {
            for (auto exp : exp_tmp)
            {
              explosionCache.erase(*exp);
            }
            Lit top = extend(e);
            return top;
          }

          // Copy a b_sz times
          for (unsigned j = 0; j < a_sz; j++)
            for (unsigned k = 0; k < b_sz; k++)
              result.push(a[j]);

          for (unsigned k = 0; k < b_sz; k++)
            for (unsigned j = 0; j < a_sz; j++)
            {
              vec<Lit> &q = result[(k * a_sz) + j];

              for (int m = 0; m < b[k].size(); m++)
              {
                bool found = false;
                for (int n = 0; n < q.size(); n++)
                  if (q[n] == b[k][m])
                  {
                    found = true;
                    break;
                  }

                if (!found)
                  q.push(b[k][m]);
              }
            }

          // std::cout << "RESULT: " << std::endl;
          // for (int i = 0; i < result.size(); i++)
          // {
          //   std::cout << "C: ";
          //   for (int j = 0; j < result[i].size(); j++)
          //     std::cout << " " << (sign(result[i][j]) ? "-" : "") << var(result[i][j]);
          //   std::cout << std::endl;
          // }

          clausifyExplodeStack.pop_back();
        }
      }
      else if (m.isNil(cur))
        throw std::runtime_error("Cannot clausify NIL.");
      else
        throw std::runtime_error("Unexpected expression type.");
    }
  }

  assert(explosionCache.find(e) != explosionCache.end());
  for (int i = 0; i < explosionCache[e].size(); i++)
    new_clauses.push(explosionCache[e][i]);

  return lit_Undef;
}

Lit MiniSAT_1p::clausify_clause(CExpression &e)
{
  //std::cout << "CCLause" << std::endl;

  assert(m.isClause(e));
  temp_lits.clear();
  getLits(e, temp_lits);

  //if (temp_lits.size()==1)
  //  addUnit(temp_lits[0]);
  //else
  //  addClause(temp_lits);
  new_clauses.push(temp_lits);

  return lit_Undef;
}

Lit MiniSAT_1p::clausify_cube(CExpression &e)
{
  //std::cout << "CCube" << std::endl;

  assert(m.isCube(e));
  temp_lits.clear();
  getLits(e, temp_lits);

  for (int i = 0; i < temp_lits.size(); i++)
  {
    //addUnit(temp_lits[i]);
    addUnit_tmp[0] = temp_lits[i];
    new_clauses.push(addUnit_tmp);
  }

  return lit_Undef;
}

Lit MiniSAT_1p::clausify_cnf(CExpression &e)
{
  //std::cout << "CCNF" << std::endl;

  assert(m.isCNF(e));

  clausifyCNFStack.clear();
  clausifyCNFStack.push_back(e);

  while (!clausifyCNFStack.empty())
  {
    CExpression q = clausifyCNFStack.back();
    clausifyCNFStack.pop_back();

    if (m.isAnd(q))
    {
      for (unsigned i = 0; i < m.nChildren(q); i++)
      {
        CExpression child = m.getChild(q, i);
        clausifyCNFStack.push_back(child);
      }
    }
    else
    {
      assert(m.isClause(q));
      clausify_clause(q);
    }
  }

  return lit_Undef;
}

ModelValue MiniSAT_1p::get(signed l) const
{
  bool sgn = (l < 0);
  signed v = (sgn) ? -l : l;
  lbool mv = Solver::model[v];
  // std::cout << "MODEL=" << (mv==D_True) << std::endl;
  if (!sgn && mv == D_True || sgn && mv == D_False)
  {
    return M_TRUE;
  }
  else if (!sgn && mv == D_False || sgn && mv == D_True)
  {
    return M_FALSE;
  }
  else
  {
    return M_UNDEF;
  }
}

Expression MiniSAT_1p::getModel(void) const
{
  std::vector<Expression> children;

  for (int i = 1; i < Solver::nVars(); i++)
  {
    ModelValue v = get(i);
    switch (v)
    {
    case M_TRUE:
      children.push_back(m.mkLiteral(i));
      break;
    case M_FALSE:
      children.push_back(m.mkLiteral(-i));
      break;
    default:
        /* ignore */;
    }
  }

  if (children.size() == 0)
    return m.mkTrue();
  else if (children.size() == 1)
    return children[0];
  else
  {
    Expression r = m.mkAnd(children[0], children[1]);
    for (unsigned i = 2; i < children.size(); i++)
      r = m.mkAnd(r, children[i]);
    return r;
  }
}

std::vector<int> MiniSAT_1p::getModelVector(void) const
{
  std::vector<int> model;
  for (int i = 1; i < Solver::nVars(); i++)
  {
    ModelValue vm = get(i);
    if (v.occurs(i) && !sharedVariables.isShared(i))
    {
      if (vm == M_TRUE)
        model.push_back(i);
      else if (vm == M_FALSE)
        model.push_back(-i);
    }
  }
  return model;
}

std::vector<int> MiniSAT_1p::getFinalModel(void) const
{
  throw std::runtime_error("NYI: getFinalModel");
  return {};
}

MiniSAT_1p::InterpolationTraverser::InterpolationTraverser(void)
{
}

MiniSAT_1p::InterpolationTraverser::~InterpolationTraverser()
{
}

void MiniSAT_1p::InterpolationTraverser::root(const vec<Lit> &c)
{
  resChainInfo.push();
  ResolutionChainInfo &r = resChainInfo.last();
  // We abuse r.xs to save the clause literals
  r.cs.clear();
  r.xs.clear();

  //c.copyTo(r.xs);

  for (int i = 0; i < c.size(); i++)
    r.xs.push_back(INT_LIT(c[i]));

#ifdef _DEBUG
  c.copyTo(r.resolvent);
#endif
}

void MiniSAT_1p::InterpolationTraverser::chain(const vec<ClauseId> &cs, const vec<Lit> &xs)
{
  resChainInfo.push();
  ResolutionChainInfo &r = resChainInfo.last();
  cs.copyTo(r.cs);
  // xs.copyTo(r.xs);

  for (int i = 0; i < xs.size(); i++)
    r.xs.push_back(INT_LIT(xs[i]));

#ifdef _DEBUG
  vec<Lit> &c = resChainInfo.last().resolvent;
  resChainInfo[cs[0]].resolvent.copyTo(c);
  for (int i = 0; i < xs.size(); i++)
  {
    // resolve(c, clauses[cs[i+1]], xs[i]);
    ClauseId next_id = cs[i + 1];
    resolve(c, resChainInfo[next_id].resolvent, xs[i]);
  }
#endif
}

void MiniSAT_1p::InterpolationTraverser::deleted(ClauseId c)
{
  //std::cout << "DELETE" << std::endl;
}

void MiniSAT_1p::InterpolationTraverser::done(ClauseId e)
{
#ifdef _DEBUG
  vec<Lit> &c = resChainInfo[e].resolvent;
  //std::cout << "FINAL:";
  //for (int i = 0; i < c.size(); i++)
  //  std::cout << " " << (sign(c[i]) ? "-":"") << var(c[i]);
  //std::cout << std::endl;
#endif
}

void MiniSAT_1p::InterpolationTraverser::resolve(vec<Lit> &main, vec<Lit> &other, Lit x) const
{
  //std::cout << "RESOLUTION: " << std::endl;

  //std::cout << " MAIN:";
  //for (int i = 0; i < main.size(); i++)
  //  std::cout << " " << (sign(main[i]) ? "-":"") << var(main[i]);
  //std::cout << std::endl;

  //std::cout << " OTHER:";
  //for (int i = 0; i < other.size(); i++)
  //  std::cout << " " << (sign(other[i]) ? "-":"") << var(other[i]);
  //std::cout << std::endl;

  //std::cout << " OVER: " << (sign(x) ? "-":"") << var(x) << std::endl;

  Lit p;
  bool ok1 = false, ok2 = false;
  for (int i = 0; i < main.size(); i++)
  {
    if (main[i] == x)
    {
      ok1 = true, p = main[i];
      main[i] = main.last();
      main.pop();
      break;
    }
  }

  assert(ok1);

  for (int i = 0; i < other.size(); i++)
  {
    if (other[i] != ~x)
      main.push(other[i]);
    else
    {
      if (p != ~other[i])
        throw std::runtime_error("PROOF ERROR! Resolved on variable with SAME polarity in both clauses."); // ... %d\n", x+1);
      ok2 = true;
    }
  }

  if (!ok1 || !ok2)
    throw std::runtime_error("PROOF ERROR! Resolved on missing variable."); // ... %d\n", x+1);

  sortUnique(main);
}

void MiniSAT_1p::addNewClauses(void)
{
  int sz = new_clauses.size();
  for (int i = 0; i < sz; i++)
  {
    //std::cout << "Adding:";
    //for (int j=0;j<new_clauses[i].size();j++)
    //  std::cout << " " << (sign(new_clauses[i][j])?"-":"") << var(new_clauses[i][j]);
    //std::cout << std::endl;

    Solver::newClause(new_clauses[i], false, ClauseId_NULL);
  }
  new_clauses.clear();
}

void MiniSAT_1p::addNewClauses(vec<vec<Lit>> interpolant_clauses)
{
  int sz = interpolant_clauses.size();
  for (int i = 0; i < sz; i++)
  {
    // std::cout << "Adding:";
    //for (int j=0;j<interpolant_clauses[i].size();j++)
    //   std::cout << " " << (sign(interpolant_clauses[i][j])?"-":"") << var(interpolant_clauses[i][j]);
    // std::cout << std::endl;

    Solver::newClause(interpolant_clauses[i], false, ClauseId_NULL);
  }
  new_clauses.clear();
  interpolant_clauses.clear();
}

void MiniSAT_1p::splitDB(std::vector<signed> &sendBuffer)
{
  //TO DO: break the buffer into pieces
  //TO DO: splitDB sends sets of clauses to the different workers

  sendBuffer.clear();
  sendBuffer.push_back(0);
  sendBuffer.push_back(0);

  vec<vec<Lit>> split_clauses;
  Solver::splitDB(split_clauses);

  int sz = split_clauses.size();
  sendBuffer[1] = sz;

  for (int i = 0; i < sz; i++)
  {

    sendBuffer.push_back(split_clauses[i].size());
    for (int j = 0; j < split_clauses[i].size(); j++)
      if (sign(split_clauses[i][j]))
        sendBuffer.push_back(-var(split_clauses[i][j]));
      else
        sendBuffer.push_back(var(split_clauses[i][j]));
  }

  sendBuffer[0] = sendBuffer.size();
}

void MiniSAT_1p::removeLastToClause(std::vector<signed> &clause)
{
  clause.clear();
  vec<vec<Lit>> split_clauses;
  Solver::removeLast(split_clauses);

  int sz = split_clauses.size();

  for (int i = 0; i < sz; i++)
  {

    for (int j = 0; j < split_clauses[i].size(); j++)
      if (sign(split_clauses[i][j]))
        clause.push_back(-var(split_clauses[i][j]));
      else
        clause.push_back(var(split_clauses[i][j]));
  }
}

void MiniSAT_1p::removeLast(std::vector<signed> &sendBuffer)
{
  if (sendBuffer.size() < 2)
  {
    sendBuffer.push_back(0);
    sendBuffer.push_back(0);
  }

  vec<vec<Lit>> split_clauses;
  Solver::removeLast(split_clauses);

  int sz = split_clauses.size();
  sendBuffer[1] = sz;

  for (int i = 0; i < sz; i++)
  {

    sendBuffer.push_back(split_clauses[i].size());
    for (int j = 0; j < split_clauses[i].size(); j++)
      if (sign(split_clauses[i][j]))
        sendBuffer.push_back(-var(split_clauses[i][j]));
      else
        sendBuffer.push_back(var(split_clauses[i][j]));
  }

  sendBuffer[0] = sendBuffer.size();
  sendBuffer[1]++;
}
