// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include <set>

#include <ExpressionManager.h>

#include "shared_variables.h"

class Interpolator
{
public:
  Interpolator(ExpressionManager &m, const SharedVariables &v) : m(m),v(v) {}
  virtual ~Interpolator() = default;

  virtual Expression root(const std::vector<signed> &clause, bool from_B=false) const = 0;
  virtual Expression resolve(const Expression a, const Expression b, signed pivot) const = 0;

protected:
  ExpressionManager &m;
  const SharedVariables &v;
};

#endif
