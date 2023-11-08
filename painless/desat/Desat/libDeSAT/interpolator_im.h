// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _INTERPOLATOR_IM_H_
#define _INTERPOLATOR_IM_H_

#include "interpolator.h"

// Inverse McMillan-Interpolator
class InterpolatorIM : public Interpolator {
public:
  InterpolatorIM(ExpressionManager &m, const SharedVariables &v) : Interpolator(m,v) {}

  virtual Expression root(const std::vector<signed> &clause, bool from_B=false) const;
  virtual Expression resolve(const Expression a, const Expression b, signed pivot) const;
};


// Inlines

inline Expression InterpolatorIM::root(const std::vector<signed> &clause, bool from_B) const
{
  if (!from_B)
    return m.mkFalse();  

  std::vector<Expression> children;

  for (unsigned i=0; i<clause.size(); i++)
  {
    const signed &l=clause[i];
    
    if (v.isShared(l))
    {
      //std::cout << "SHARED: " << l << std::endl;      
      children.push_back(m.mkLiteral( - l));
    }
    else
    {
      //std::cout << "UNSHARED: " << l << std::endl;
    }
  }

  switch (children.size())
  {
  case 0: return m.mkTrue(); break;
  case 1: return children[0]; break;
  default: 
    {
      Expression r = m.mkAnd(children[0], children[1]);
      for (unsigned i=2; i<children.size(); i++)
        r = m.mkAnd(r, children[i]);
      return r;
    }
  }
}

inline Expression InterpolatorIM::resolve(Expression a, Expression b, signed pivot) const
{  
  //std::cout << "GET A: " << m.toString(a) << std::endl;
  //std::cout << "GET B: " << m.toString(b) << std::endl;

  // Note: pivot cannot be from B since it is asummed to be an assignment.
  return m.mkOr(a, b);
}

#endif