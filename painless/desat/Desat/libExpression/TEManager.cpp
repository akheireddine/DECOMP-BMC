// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#include "TEManager.h"

TrivialExpressionManager::TrivialExpressionManager(void)
{
}

TrivialExpressionManager::~TrivialExpressionManager(void)
{
  for (unsigned i=0; i<literals.size(); i++)
    delete literals[i];
  for (unsigned i=0; i<nliterals.size(); i++)
    delete nliterals[i];
  for (unsigned i=0; i<ands.size(); i++)
    delete ands[i];
  for (unsigned i=0; i<ors.size(); i++)
    delete ors[i];
}

void TrivialExpressionManager::toBuffer(const TrivialExpression *t, std::vector<signed> &buffer) const
{
  buffer.push_back(t->type);

  if (isAnd(t) || isOr(t))
  {
    buffer.push_back(nChildren(t));
    for (unsigned i=0; i<nChildren(t); i++)
      toBuffer(getChild(t, i), buffer);
  }
}

TrivialExpression* TrivialExpressionManager::fromBuffer(std::vector<signed> &buffer, unsigned inx)
{
  assert(buffer.size()>inx);
  TrivialExpression::E_TYPE type = (TrivialExpression::E_TYPE) buffer[inx];
  assert (type==TrivialExpression::AND || type==TrivialExpression::OR);

  std::vector<const TrivialExpression*> children;

  unsigned s = buffer[++inx];
  children.clear();
  inx++;
  for (unsigned i=0; i<s; i++)
    children.push_back(fromBuffer(buffer, inx+i));

  switch(type)
  {
  case TrivialExpression::AND: return mkAnd(children);
  case TrivialExpression::OR: return mkOr(children);
  default:
    throw std::runtime_error("UNSUPPORTED");
  }
}

TrivialExpression* TrivialExpressionManager::duplicate(const TrivialExpression *a, const TrivialExpressionManager &other)
{
  tStack.clear();
  tStack.push_back(a);

  Cache cache;

  while(!tStack.empty())
  {
    const TrivialExpression *q = tStack.back();

    Cache::const_iterator it = cache.find(q);
    if (it!=cache.end())
      tStack.pop_back();
    else
    {
      if (other.isNil(q))
      {
        cache[q] = mkNil();
        tStack.pop_back();
      }
      else if (other.isTrue(q))
      {
        cache[q] = mkTrue();
        tStack.pop_back();
      }
      else if (other.isFalse(q))
      {
        cache[q] = mkFalse();
        tStack.pop_back();
      }
      else if (other.isLiteral(q))
      {
        cache[q] = mkLiteral(other.getLiteral(q));
        tStack.pop_back();
      }
      else
      {
        bool children_done = true;

        for (unsigned i=0; i<other.nChildren(q); i++)
        {
          const TrivialExpression *c = other.getChild(q, i);
          if (cache.find(c)==cache.end())
          {
            tStack.push_back(c);
            children_done = false;
          }
        }


        if (children_done)
        {
          cache[q] = cache[other.getChild(q, 0)];
          for (unsigned i=1; i<other.nChildren(q); i++)
          {
            if (other.isAnd(q))
              cache[q] = mkAnd(cache[q], cache[other.getChild(q, i)]);
            else if (other.isOr(q))
              cache[q] = mkOr(cache[q], cache[other.getChild(q, i)]);
          }

          tStack.pop_back();
        }
      }

    }
  }

  assert(cache.find(a)!=cache.end());
  return cache[a];
}
