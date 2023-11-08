// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _DECOMPOSITION_CYCLE_H_
#define _DECOMPOSITION_CYCLE_H_

#include "decomposition.h"

class CycleDecomposition : public Decomposition
{
public:
  CycleDecomposition(void);
  virtual ~CycleDecomposition(void);

  virtual int where(const std::vector<signed> &clause);

protected:
  unsigned last;
};

#endif