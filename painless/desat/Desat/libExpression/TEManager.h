// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _TE_MANAGER_H_
#define _TE_MANAGER_H_

#include <iostream>

#include <vector>
#include <unordered_map>

#include "TE.h"

class TrivialExpressionManager
{
public:
  TrivialExpressionManager(void);
  ~TrivialExpressionManager(void);

  typedef std::unordered_map<const TrivialExpression*, TrivialExpression*> Cache;
  std::vector<const TrivialExpression *> tStack;  

  NilTE te_nil;
  TrueTE te_true;
  FalseTE te_false;

  class Index1Set : public std::vector<LiteralTE*> {};  
  Index1Set literals, nliterals;
  std::vector<TrivialExpression*> ands, ors;

  TrivialExpression* mkNil(void)
  {
    return &te_nil;
  }

  TrivialExpression* mkTrue(void)
  {
    return &te_true;
  }

  TrivialExpression* mkFalse(void)
  {
    return &te_false;
  }

  TrivialExpression* mkLiteral(signed l)
  {
    bool sgn = (l<0);
    unsigned v = (sgn) ? -l : l;

    if (sgn)
    {
      if (nliterals.size() <= v)
      {
        size_t old_size = nliterals.size();
        nliterals.resize(v+1);

        for (unsigned i=old_size; i<nliterals.size(); i++)
          nliterals[i] = new LiteralTE(-(signed)i);
      }
    }
    else
    {
      if (literals.size() <= v)
      {
        size_t old_size = literals.size();
        literals.resize(v+1);

        for (unsigned i=old_size; i<literals.size(); i++)
          literals[i] = new LiteralTE(i);
      }
    }  

    return (sgn) ? nliterals[v] : literals[v];
  }

  TrivialExpression* mkAnd(const TrivialExpression *a, const TrivialExpression *b)
  {
    ands.push_back(new AndTE(a, b));
    return ands.back();
  }

  TrivialExpression* mkAnd(std::vector<const TrivialExpression*> &children)
  {
    assert(children.size()>2);
    TrivialExpression *last = mkAnd(children[0], children[1]);
    for (unsigned i=2; i<children.size(); i++)
      last = mkAnd(last, children[1]);
    return last;
  }

  TrivialExpression* mkOr(const TrivialExpression *a, const TrivialExpression *b)
  {
    ors.push_back(new OrTE(a, b));
    return ors.back();
  }

  TrivialExpression* mkOr(std::vector<const TrivialExpression*> &children)
  {
    assert(children.size()>2);
    TrivialExpression *last = mkOr(children[0], children[1]);
    for (unsigned i=2; i<children.size(); i++)
      last = mkOr(last, children[1]);
    return last;
  }
  
  TrivialExpression* mkEq(const TrivialExpression *a, const TrivialExpression *b)
  {
    assert(false);
    return NULL;
  }

  std::string toString(const TrivialExpression *a) const
  {
    return a->toString();
  }  

  bool isTrue(const TrivialExpression *a) const
  {
    return a->isTrue();
  }

  bool isFalse(const TrivialExpression *a) const
  {
    return a->isFalse();
  }

  bool isNil(const TrivialExpression *a) const
  {
    return a->isNil();
  }

  bool isAnd(const TrivialExpression *a) const
  {
    return a->isAnd();
  }

  bool isOr(const TrivialExpression *a) const
  {
    return a->isOr();
  }

  bool isLiteral(const TrivialExpression *a) const
  {
    return a->isLiteral();
  }

  signed getLiteral(const TrivialExpression *a) const
  {
    return a->getLiteral();    
  }

  unsigned nChildren(const TrivialExpression *a) const
  {
    return a->nChildren();
  }

  const TrivialExpression *getChild(const TrivialExpression *a, unsigned i) const
  {
    return a->getChild(i);
  }

  bool isNegative(const TrivialExpression *a) const
  {
    return false;
  }

  bool isClause(const TrivialExpression *a) const
  {
    return a->isClause();
  }

  bool isCube(const TrivialExpression *a) const
  {
    return a->isCube();
  }

  bool isCNF(const TrivialExpression *a) const
  {
    return a->isCNF();
  }

  std::vector<const TrivialExpression*> getLiteralsStack;

  void getLiterals(const TrivialExpression *a, std::vector<signed> &lits)
  {
    getLiteralsStack.clear();
    getLiteralsStack.push_back(a);
    lits.clear();

    while(!getLiteralsStack.empty())
    {
      const TrivialExpression *q = getLiteralsStack.back();
      getLiteralsStack.pop_back();
  
      if(isAnd(q) || isOr(q))
      {        
        for (unsigned i=0; i<nChildren(q); i++)
          getLiteralsStack.push_back(getChild(q, i));
      }
      else if (isLiteral(q))
        lits.push_back(getLiteral(q));
    }
  }

  TrivialExpression *mkOr(const std::vector<TrivialExpression*> &literals)
  {
    TrivialExpression *res;

    if (literals.size()==0)
      res = mkTrue();
    else if (literals.size()==1)
      res = literals[0];
    else
    {
      res = mkOr(literals[0], literals[1]);
      for (unsigned i=2; i<literals.size(); i++)
        res = mkOr(res, literals[i]);
    }

    return res;
  }  

  TrivialExpression* duplicate(const TrivialExpression *a, const TrivialExpressionManager &other);
  TrivialExpression* fromBuffer(std::vector<signed> &buffer, unsigned inx);
  void toBuffer(const TrivialExpression *t, std::vector<signed> &buffer) const;
};

#endif
