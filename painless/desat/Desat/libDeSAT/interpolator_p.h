// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _INTERPOLATOR_P_H_
#define _INTERPOLATOR_P_H_

#include "interpolator.h"

// Pudlak-Interpolator
class InterpolatorP : public Interpolator {
public:
  InterpolatorP(ExpressionManager &m, const SharedVariables &v) : Interpolator(m,v) {}

  virtual Expression root(const std::vector<signed> &clause, bool from_B=false) const;
  virtual Expression resolve(const Expression a, const Expression b, signed pivot) const;
};


// Inlines

inline Expression InterpolatorP::root(const std::vector<signed> &clause, bool from_B) const
{
  if (from_B)
    return m.mkTrue();
  else
    return m.mkFalse();
}

inline Expression InterpolatorP::resolve(Expression a, Expression b, signed pivot) const
{
  // Precondition: The corresponding resolvent for `a' contains `over,
  // while the resolvent for `b' contains `-over'. The phase of the pivot
  // makes a crucial difference here.

  if (v.isShared(pivot))
  {   
    Expression or1 = m.mkOr(m.mkLiteral(pivot), a);
    Expression or2 = m.mkOr(m.mkLiteral(-pivot), b);
    return m.mkAnd(or1, or2);
  }
  else
    return m.mkOr(a, b);
}

#endif