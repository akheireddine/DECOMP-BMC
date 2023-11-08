// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _INTERPOLATOR_M_H_
#define _INTERPOLATOR_M_H_

#include "interpolator.h"

// McMillan-Interpolator
class InterpolatorM : public Interpolator {
public:
  InterpolatorM(ExpressionManager &m, const SharedVariables &v) : Interpolator(m,v) {}  

  virtual Expression root(const std::vector<signed> &clause, bool from_B=false) const;
  virtual Expression resolve(const Expression a, const Expression b, signed pivot) const;
};


// Inlines

inline Expression InterpolatorM::root(const std::vector<signed> &clause, bool from_B) const
{

  if (from_B)
    return m.mkTrue();  

  std::vector<Expression> children;

  for (unsigned i=0; i<clause.size(); i++)
  {
    const signed &l=clause[i];

	if (v.isShared(l))
    {
      //std::cout << "SHARED: " << l << std::endl;      
      children.push_back(m.mkLiteral(l));
    }
    else
    {
      //std::cout << "UNSHARED: " << l << std::endl;
    }
	
  }  
  
  switch (children.size())
  {
  case 0: return m.mkFalse(); break;
  case 1: return children[0]; break;
  default: 
    {
      Expression r = m.mkOr(children[0], children[1]);
      for (unsigned i=2; i<children.size(); i++)
        r = m.mkOr(r, children[i]);
      return r;
    }
  }
}

inline Expression InterpolatorM::resolve(Expression a, Expression b, signed pivot) const
{  
  //std::cout << "GET A: " << m.toString(a) << std::endl;
  //std::cout << "GET B: " << m.toString(b) << std::endl;
  //std::cout << "PIVOT: " << over << " (" << t.isShared(over) << ")" << std::endl;
	if (v.isShared(pivot))
    return m.mkAnd(a, b);
  else
    return m.mkOr(a, b);
}

#endif