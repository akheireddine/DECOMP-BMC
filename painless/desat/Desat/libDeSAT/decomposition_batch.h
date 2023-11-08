// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _DECOMPOSITION_BATCH_H_
#define _DECOMPOSITION_BATCH_H_

#include "decomposition.h"

class BatchDecomposition : public Decomposition
{
public:
  BatchDecomposition(void);
  virtual ~BatchDecomposition(void);

  virtual int where(const std::vector<signed> &clause);

protected:
  unsigned count;
};

#endif